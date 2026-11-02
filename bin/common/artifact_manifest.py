"""Track files created by AIDAmri processing stages.

The tracker takes a recursive snapshot of one modality directory when a stage
starts and updates a dataset-local JSON manifest when the Python process exits.
Only newly created paths are considered owned outputs. Existing files that were
modified are recorded for provenance, but are deliberately not safe to delete.
"""

from __future__ import annotations

import atexit
import json
import os
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

try:
    import fcntl
except ImportError:  # pragma: no cover - AIDAmri currently targets Linux
    fcntl = None


MANIFEST_SUFFIX = "_aidamri_manifest.json"
LOCK_FILENAME = ".aidamri_manifest.lock"
MANIFEST_VERSION = 1
MANIFEST_FILE_MODE = 0o644
MODES = ("anat", "dwi", "func", "t2map")
STAGES = ("preprocessing", "registration", "processing")
MANIFEST_TIMEZONE = ZoneInfo("Europe/Berlin")


def _now() -> str:
    return datetime.now(MANIFEST_TIMEZONE).isoformat()


def manifest_filename(folder: str | os.PathLike[str], mode: str | None = None) -> str:
    """Build the dataset-local manifest name from the BIDS-style path."""
    root = Path(folder).resolve()
    # Use path components instead of image names so every file in one dataset
    # contributes to the same subject/session/modality manifest.
    subject = next((part for part in reversed(root.parts) if part.startswith("sub-")), None)
    if subject is None:
        raise ValueError(f"Cannot derive subject ID from modality folder: {root}")
    session = next(
        (part for part in reversed(root.parts) if part.startswith("ses-")),
        "ses-none",
    )
    mode = mode or root.name
    if mode not in MODES:
        raise ValueError(f"Unsupported modality: {mode!r}")
    return f".{subject}_{session}_{mode}{MANIFEST_SUFFIX}"


def manifest_path(
    folder: str | os.PathLike[str], mode: str | None = None
) -> Path:
    root = Path(folder).resolve()
    return root / manifest_filename(root, mode)


def _is_manifest_metadata(file_name: str) -> bool:
    return file_name == LOCK_FILENAME or (
        file_name.startswith(".")
        and (
            file_name.endswith(MANIFEST_SUFFIX)
            or f"{MANIFEST_SUFFIX}." in file_name
        )
    )


def _canonical_stage(stage: str) -> str:
    aliases = {
        "preprocess": "preprocessing",
        "preprocessing": "preprocessing",
        "register": "registration",
        "registration": "registration",
        "process": "processing",
        "processing": "processing",
    }
    try:
        return aliases[stage.strip().lower()]
    except (AttributeError, KeyError) as exc:
        raise ValueError(f"Unsupported processing stage: {stage!r}") from exc


def _path_signature(path: Path) -> tuple[int, int]:
    stat = path.stat(follow_symlinks=False)
    return stat.st_size, stat.st_mtime_ns


def snapshot(folder: str | os.PathLike[str]) -> dict[str, dict[str, tuple[int, int] | None]]:
    """Return files and directories below *folder* using relative POSIX paths."""
    root = Path(folder).resolve()
    files: dict[str, tuple[int, int]] = {}
    directories: dict[str, None] = {}

    for current_root, dir_names, file_names in os.walk(root, followlinks=False):
        current = Path(current_root)

        for dir_name in list(dir_names):
            path = current / dir_name
            relative = path.relative_to(root).as_posix()
            if path.is_symlink():
                files[relative] = _path_signature(path)
                dir_names.remove(dir_name)
            else:
                directories[relative] = None

        for file_name in file_names:
            if _is_manifest_metadata(file_name):
                continue
            path = current / file_name
            relative = path.relative_to(root).as_posix()
            try:
                files[relative] = _path_signature(path)
            except FileNotFoundError:
                # A concurrently running tool may move a temporary output while
                # the snapshot is being collected. The next snapshot sees the
                # final path.
                continue

    return {"files": files, "directories": directories}


def _empty_manifest(mode: str, initial_stage: str) -> dict:
    return {
        "version": MANIFEST_VERSION,
        "mode": mode,
        "tracking_started_at": _now(),
        "tracking_started_before_stage": initial_stage,
        "updated_at": _now(),
        "stages": {},
    }


def load_manifest(
    folder: str | os.PathLike[str], mode: str | None = None
) -> dict | None:
    path = manifest_path(folder, mode)
    if not path.is_file():
        return None
    with path.open("r", encoding="utf-8") as stream:
        manifest = json.load(stream)
    if not isinstance(manifest, dict):
        raise ValueError(f"Invalid manifest structure in '{path}'")
    if manifest.get("version") != MANIFEST_VERSION:
        raise ValueError(
            f"Unsupported manifest version in '{path}': "
            f"{manifest.get('version')!r}"
        )
    if manifest.get("mode") not in MODES:
        raise ValueError(f"Invalid manifest mode in '{path}'")
    stages = manifest.get("stages")
    if not isinstance(stages, dict):
        raise ValueError(f"Invalid stages section in '{path}'")
    for stage, stage_data in stages.items():
        if stage not in STAGES or not isinstance(stage_data, dict):
            raise ValueError(f"Invalid stage {stage!r} in '{path}'")
        for field in ("created_files", "created_directories", "modified_files"):
            paths = stage_data.get(field, [])
            if not isinstance(paths, list) or not all(isinstance(path, str) for path in paths):
                raise ValueError(
                    f"Invalid {field!r} list for stage {stage!r} in '{path}'"
                )
    return manifest


def _write_manifest_unlocked(folder: Path, manifest: dict) -> None:
    manifest["updated_at"] = _now()
    path = manifest_path(folder, manifest.get("mode"))
    file_descriptor, temp_name = tempfile.mkstemp(
        prefix=f"{path.name}.",
        suffix=".tmp",
        dir=folder,
        text=True,
    )
    try:
        # mkstemp defaults to 0600. Make the final manifest readable outside
        # root-owned containers before atomically moving it into place.
        os.fchmod(file_descriptor, MANIFEST_FILE_MODE)
        stream = os.fdopen(file_descriptor, "w", encoding="utf-8")
        file_descriptor = None
        with stream:
            json.dump(manifest, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_name, path)
    finally:
        if file_descriptor is not None:
            os.close(file_descriptor)
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass


def write_manifest(folder: str | os.PathLike[str], manifest: dict) -> None:
    """Atomically write a manifest while serializing concurrent updates."""
    root = Path(folder).resolve()
    # The stable lock file coordinates the parent batch process and child script.
    lock_path = root / LOCK_FILENAME
    with lock_path.open("a", encoding="utf-8") as lock_stream:
        if fcntl is not None:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX)
        try:
            _write_manifest_unlocked(root, manifest)
        finally:
            if fcntl is not None:
                fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)


def _merge_changes(
    folder: Path,
    mode: str,
    stage: str,
    created_files: set[str],
    created_directories: set[str],
    modified_files: set[str],
) -> None:
    lock_path = folder / LOCK_FILENAME
    with lock_path.open("a", encoding="utf-8") as lock_stream:
        if fcntl is not None:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX)
        try:
            manifest = load_manifest(folder, mode) or _empty_manifest(mode, stage)
            manifest_mode = manifest.get("mode")
            if manifest_mode != mode:
                raise ValueError(
                    f"Manifest mode mismatch in '{folder}': "
                    f"expected {mode!r}, found {manifest_mode!r}"
                )

            stage_data = manifest.setdefault("stages", {}).setdefault(
                stage,
                {"created_files": [], "created_directories": [], "modified_files": []},
            )
            stage_data["created_files"] = sorted(
                set(stage_data.get("created_files", [])) | created_files
            )
            stage_data["created_directories"] = sorted(
                set(stage_data.get("created_directories", [])) | created_directories
            )
            stage_data["modified_files"] = sorted(
                set(stage_data.get("modified_files", [])) | modified_files
            )
            stage_data["last_observed_at"] = _now()
            _write_manifest_unlocked(folder, manifest)
        finally:
            if fcntl is not None:
                fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)


@dataclass
class OutputTracker:
    folder: Path
    mode: str
    stage: str
    before: dict
    _finalized: bool = False

    @classmethod
    def start(cls, folder: str | os.PathLike[str], mode: str, stage: str) -> "OutputTracker":
        root = Path(folder).resolve()
        if not root.is_dir():
            raise NotADirectoryError(f"Output folder does not exist: {root}")
        if mode not in MODES:
            raise ValueError(f"Unsupported modality: {mode!r}")
        tracker = cls(root, mode, _canonical_stage(stage), snapshot(root))
        atexit.register(tracker.finalize)
        return tracker

    def finalize(self) -> None:
        if self._finalized:
            return
        self._finalized = True
        try:
            atexit.unregister(self.finalize)
        except Exception:
            pass

        after = snapshot(self.folder)
        before_files = self.before["files"]
        after_files = after["files"]
        # Only paths absent from the initial snapshot are safe reset candidates.
        created_files = set(after_files) - set(before_files)
        modified_files = {
            path
            for path in set(after_files) & set(before_files)
            if after_files[path] != before_files[path]
        }
        created_directories = set(after["directories"]) - set(self.before["directories"])

        _merge_changes(
            self.folder,
            self.mode,
            self.stage,
            created_files,
            created_directories,
            modified_files,
        )


def start_output_tracking(
    folder: str | os.PathLike[str], mode: str, stage: str
) -> OutputTracker:
    """Start tracking outputs and finalize automatically at process exit."""
    return OutputTracker.start(folder, mode, stage)
