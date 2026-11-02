"""Reset AIDAmri modality folders using runtime-generated output manifests."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path, PurePosixPath


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.artifact_manifest import (
    LOCK_FILENAME,
    MANIFEST_SUFFIX,
    MODES,
    load_manifest,
    manifest_filename,
    manifest_path,
    snapshot,
    write_manifest,
)


STAGE_ORDER = {
    "base": 0,
    "preprocessing": 1,
    "registration": 2,
    "processing": 3,
}
NIFTI_SUFFIXES = (".nii", ".nii.gz")
MUTED_SUMMARY_SUFFIXES = (
    "_T2w.json",
    "_EPI.json",
    "_dwi.json",
    ".bvec",
    ".bval",
)
STROKE_MASK_SUFFIX = "Stroke_mask.nii.gz"
PROTECTED_UNMANAGED_SUFFIXES = MUTED_SUMMARY_SUFFIXES + (STROKE_MASK_SUFFIX,)


class ResetError(RuntimeError):
    """Raised when a safe manifest-based reset cannot be planned."""


def normalize_phase(phase):
    aliases = {
        "base": "base",
        "grundzustand": "base",
        "convert2nifti": "base",
        "preprocess": "preprocessing",
        "preprocessing": "preprocessing",
        "register": "registration",
        "registration": "registration",
        "process": "processing",
        "processing": "processing",
    }
    try:
        return aliases[phase.strip().lower()]
    except (AttributeError, KeyError) as exc:
        raise ValueError(f"Unsupported phase: {phase!r}") from exc


def is_bids_like(project_root):
    """Return whether *project_root* contains a BIDS-like subject structure."""
    project_root = Path(project_root)
    if not project_root.is_dir():
        print(f"Error: project_root '{project_root}' is not a directory.")
        return False

    subject_dirs = [
        path for path in project_root.iterdir()
        if path.name.startswith("sub-") and path.is_dir()
    ]
    if not subject_dirs:
        print("Error: No 'sub-*' directories found below project_root.")
        return False

    modalities = {"anat", "dwi", "fmap", "func", "t2map"}
    for subject_dir in subject_dirs:
        subject_children = [path for path in subject_dir.iterdir() if path.is_dir()]
        if any(path.name in modalities for path in subject_children):
            return True
        for session_dir in subject_children:
            if not session_dir.name.startswith("ses-"):
                continue
            if any(
                path.is_dir() and path.name in modalities
                for path in session_dir.iterdir()
            ):
                return True

    print(
        "Error: No typical BIDS modality directories were found below "
        "the 'sub-*' directories."
    )
    return False


def has_any_brkraw(project_root):
    for _, dir_names, _ in os.walk(project_root):
        if "brkraw" in dir_names:
            return True
    return False


def find_modality_dirs(project_root, mode):
    modality_dirs = []
    for root, dir_names, _ in os.walk(project_root):
        if os.path.basename(root) == mode:
            modality_dirs.append(Path(root).resolve())
            dir_names.clear()
    return sorted(modality_dirs)


def validate_project_manifests(project_root, mode):
    """Validate manifests in existing folders of the selected modality only."""
    if mode not in MODES:
        raise ValueError(f"Unsupported modality: {mode!r}")
    errors = []
    for folder in find_modality_dirs(project_root, mode):
        try:
            expected_path = manifest_path(folder, mode)
        except ValueError as exc:
            errors.append(str(exc))
            continue
        if not expected_path.is_file():
            errors.append(f"Manifest is missing: {expected_path}")
            continue
        try:
            manifest = load_manifest(folder, mode)
        except (OSError, ValueError) as exc:
            errors.append(f"Manifest is invalid: {expected_path}: {exc}")
            continue
        if manifest.get("mode") != mode:
            errors.append(
                f"Manifest mode mismatch: {expected_path} "
                f"(found {manifest.get('mode')!r}, expected {mode!r})"
            )
    return errors


def _manifest_path(folder, relative_path):
    """Resolve a manifest path lexically without following dataset symlinks."""
    if not isinstance(relative_path, str):
        raise ResetError(f"Manifest path is not a string: {relative_path!r}")
    path = PurePosixPath(relative_path)
    if not path.parts or path.is_absolute() or ".." in path.parts:
        raise ResetError(f"Unsafe path in manifest: {relative_path!r}")
    return folder.joinpath(*path.parts)


def _stage_is_after(stage, target_phase):
    if stage not in STAGE_ORDER:
        raise ResetError(f"Unknown stage in manifest: {stage!r}")
    return STAGE_ORDER[stage] > STAGE_ORDER[target_phase]


def _is_nifti_path(path):
    return str(path).endswith(NIFTI_SUFFIXES)


def _build_base_restore_plan(folder, manifest):
    """Find and validate preprocessing backups required for a base reset."""
    restore_files = []
    errors = []
    preprocessing = manifest.get("stages", {}).get("preprocessing", {})
    created_files = set(preprocessing.get("created_files", []))
    modified_niftis = {
        path
        for path in preprocessing.get("modified_files", [])
        if _is_nifti_path(path)
    }

    backup_paths = {
        path
        for path in created_files
        if len(PurePosixPath(path).parts) == 2
        and PurePosixPath(path).parts[0] == "brkraw"
        and _is_nifti_path(path)
    }

    # Every modified preprocessing input needs its matching original backup.
    for relative_destination in modified_niftis:
        backup_relative = f"brkraw/{PurePosixPath(relative_destination).name}"
        if backup_relative not in backup_paths:
            errors.append(
                f"No registered brkraw backup for modified input: "
                f"{folder / relative_destination}"
            )

    for relative_source in sorted(backup_paths):
        source = _manifest_path(folder, relative_source)
        destination = folder / PurePosixPath(relative_source).name
        if not source.is_file() or source.is_symlink():
            errors.append(f"Backup is missing or is not a regular file: {source}")
            continue
        if destination.is_symlink() or destination.is_dir():
            errors.append(
                f"Restore destination is not a regular file path: {destination}"
            )
            continue
        restore_files.append((source, destination))

    if modified_niftis and not restore_files:
        errors.append(f"No usable brkraw NIfTI backup found in: {folder / 'brkraw'}")

    return restore_files, errors


def _base_restore_errors(plans):
    return [
        f"{plan['folder']}: {error}"
        for plan in plans
        for error in plan.get("restore_errors", [])
    ]


def _known_manifest_directories(manifest):
    """Return explicit and implicit directory paths referenced by a manifest."""
    known_directories = set()
    referenced_paths = set()

    for stage_data in manifest.get("stages", {}).values():
        created_directories = set(stage_data.get("created_directories", []))
        known_directories.update(created_directories)
        referenced_paths.update(created_directories)
        referenced_paths.update(stage_data.get("created_files", []))
        referenced_paths.update(stage_data.get("modified_files", []))

    # A directory containing a tracked path is implicitly known even if the
    # directory itself existed before tracking began.
    for relative_path in referenced_paths:
        for parent in PurePosixPath(relative_path).parents:
            if parent != PurePosixPath("."):
                known_directories.add(parent.as_posix())

    return known_directories


def _is_protected_unmanaged_file(path):
    return path.name.endswith(PROTECTED_UNMANAGED_SUFFIXES)


def _is_stroke_mask(path):
    return path.name.endswith(STROKE_MASK_SUFFIX)


def build_reset_plan(project_root, mode, phase, delete_unmanaged=False):
    """Build a non-mutating deletion plan from all matching manifests."""
    project_root = Path(project_root).resolve()
    phase = normalize_phase(phase)
    plans = []

    for folder in find_modality_dirs(project_root, mode):
        warnings = []
        if not (folder / "brkraw").is_dir():
            warnings.append("No direct brkraw directory; modality folder is skipped.")
            plans.append(
                {
                    "folder": folder,
                    "skip": True,
                    "warnings": warnings,
                    "restore_errors": (
                        ["A base reset requires a direct brkraw directory."]
                        if phase == "base"
                        else []
                    ),
                }
            )
            continue

        try:
            manifest = load_manifest(folder, mode)
        except (OSError, ValueError) as exc:
            warnings.append(f"Manifest cannot be read: {exc}")
            plans.append({"folder": folder, "skip": True, "warnings": warnings})
            continue

        if manifest is None:
            warnings.append(
                f"No {manifest_filename(folder, mode)} found; "
                "nothing is deleted for safety."
            )
            plans.append({"folder": folder, "skip": True, "warnings": warnings})
            continue
        if manifest.get("mode") != mode:
            warnings.append(
                f"Manifest belongs to {manifest.get('mode')!r}, expected {mode!r}; "
                "folder is skipped."
            )
            plans.append({"folder": folder, "skip": True, "warnings": warnings})
            continue

        files = set()
        directories = set()
        selected_stages = []
        missing_files = set()
        missing_directories = set()
        modified_files = []
        restore_files = []
        restore_errors = []

        # Classify current files against all recorded stage ownership.
        all_created_files = {
            path
            for stage_data in manifest.get("stages", {}).values()
            for path in stage_data.get("created_files", [])
        }
        all_modified_files = {
            path
            for stage_data in manifest.get("stages", {}).values()
            for path in stage_data.get("modified_files", [])
        }
        current_snapshot = snapshot(folder)
        current_files = set(current_snapshot["files"])
        unknown_files = sorted(
            _manifest_path(folder, relative_path)
            for relative_path in current_files - all_created_files - all_modified_files
        )
        known_directories = _known_manifest_directories(manifest)
        unknown_directories = sorted(
            (
                _manifest_path(folder, relative_path)
                for relative_path in (
                    set(current_snapshot["directories"]) - known_directories
                )
            ),
            key=lambda path: len(path.parts),
            reverse=True,
        )
        protected_unmanaged_files = [
            path for path in unknown_files if _is_protected_unmanaged_file(path)
        ]
        unmanaged_files_to_delete = [
            path for path in unknown_files if path not in protected_unmanaged_files
        ]
        protected_directory_paths = {
            parent.as_posix()
            for path in protected_unmanaged_files
            for parent in PurePosixPath(_relative(folder, path)).parents
            if parent != PurePosixPath(".")
        }
        unmanaged_directories_to_delete = [
            path
            for path in unknown_directories
            if _relative(folder, path) not in protected_directory_paths
        ]

        initial_stage = manifest.get("tracking_started_before_stage")
        if phase == "base":
            if initial_stage != "preprocessing":
                restore_errors.append(
                    "Manifest tracking did not start before preprocessing; "
                    "a complete base restore cannot be verified."
                )
            planned_restores, validation_errors = _build_base_restore_plan(
                folder, manifest
            )
            restore_files.extend(planned_restores)
            restore_errors.extend(validation_errors)

        restored_destinations = {
            _relative(folder, destination) for _, destination in restore_files
        }
        if initial_stage in STAGE_ORDER:
            earliest_safe_phase = STAGE_ORDER[initial_stage] - 1
            if STAGE_ORDER[phase] < earliest_safe_phase:
                warnings.append(
                    "Manifest tracking started only before stage "
                    f"{initial_stage!r}. Older untracked artifacts are preserved; "
                    "the requested target state may be incomplete."
                )

        for stage, stage_data in manifest.get("stages", {}).items():
            if not _stage_is_after(stage, phase):
                continue
            selected_stages.append(stage)

            modified = stage_data.get("modified_files", [])
            modified_files.extend(
                (stage, _manifest_path(folder, relative_path))
                for relative_path in modified
                if relative_path not in restored_destinations
            )

            for relative_path in stage_data.get("created_files", []):
                path = _manifest_path(folder, relative_path)
                if os.path.lexists(path):
                    if path.is_dir() and not path.is_symlink():
                        warnings.append(
                            f"Path recorded as a file is now a directory and is preserved: {path}"
                        )
                    else:
                        files.add(path)
                else:
                    missing_files.add(relative_path)

            for relative_path in stage_data.get("created_directories", []):
                path = _manifest_path(folder, relative_path)
                if os.path.lexists(path):
                    if path.is_dir() and not path.is_symlink():
                        directories.add(path)
                    else:
                        warnings.append(
                            f"Path recorded as a directory is no longer a directory and is preserved: {path}"
                        )
                else:
                    missing_directories.add(relative_path)

        # Do not warn about modified paths that this reset already handles by
        # restoring their original or deleting them as owned stage artifacts.
        handled_modified_paths = restored_destinations | {
            _relative(folder, path) for path in files
        }
        modified_files = [
            (stage, path)
            for stage, path in modified_files
            if _relative(folder, path) not in handled_modified_paths
            and os.path.lexists(path)
        ]

        plans.append(
            {
                "folder": folder,
                "skip": False,
                "manifest": manifest,
                "phase": phase,
                "selected_stages": selected_stages,
                "files": sorted(files),
                "directories": sorted(
                    directories,
                    key=lambda path: len(path.parts),
                    reverse=True,
                ),
                "missing_files": missing_files,
                "missing_directories": missing_directories,
                "restore_files": restore_files,
                "restore_errors": restore_errors,
                "restored_destinations": restored_destinations,
                "metadata_files": (
                    [manifest_path(folder, mode), folder / LOCK_FILENAME]
                    if phase == "base"
                    else []
                ),
                "modified_files": sorted(
                    modified_files,
                    key=lambda item: (item[0], item[1]),
                ),
                "unknown_files": unknown_files,
                "unknown_directories": unknown_directories,
                "protected_unmanaged_files": protected_unmanaged_files,
                "unmanaged_files_to_delete": unmanaged_files_to_delete,
                "unmanaged_directories_to_delete": unmanaged_directories_to_delete,
                "delete_unmanaged": delete_unmanaged,
                "deleted_unmanaged_files": [],
                "removed_unmanaged_directories": [],
                "runtime_warnings": [],
                "warnings": warnings,
            }
        )

    return plans


def _normal_stroke_mask_folders(plans):
    return {
        plan["folder"]
        for plan in plans
        if not plan.get("delete_unmanaged")
        and any(_is_stroke_mask(path) for path in plan.get("unknown_files", []))
    }


def print_reset_plan(plans, mode, phase):
    print(f"\nReset plan for mode={mode}, target phase={phase}")
    print("Only later-stage artifacts registered in the manifest will be deleted.")
    stroke_mask_folders = _normal_stroke_mask_folders(plans)
    if stroke_mask_folders:
        folder_word = "folder" if len(stroke_mask_folders) == 1 else "folders"
        print(
            f"\n[INFO] Stroke mask files are present in "
            f"{len(stroke_mask_folders)} selected {folder_word}. "
            "They are preserved; individual paths are shown below but omitted "
            "from the final summary."
        )
    if any(plan.get("delete_unmanaged") for plan in plans):
        print(
            "\n[DANGER] --delete-unmanaged is active. Files not known to the "
            "manifest are selected for permanent deletion."
        )
        print(
            "Protected BIDS sidecars, bvec/bval files, and *Stroke_mask.nii.gz "
            "files are preserved; other unmanaged user data can be deleted."
        )

    for plan in plans:
        print(f"\n{plan['folder']}")
        for warning in plan["warnings"]:
            print(f"  [WARNING] {warning}")
        for error in plan.get("restore_errors", []):
            print(f"  [ERROR] Base restore blocked: {error}")
        if plan["skip"]:
            continue
        for source, destination in plan["restore_files"]:
            print(f"  Restore original NIfTI: {source} -> {destination}")
        restore_sources = {source for source, _ in plan["restore_files"]}
        for stage, path in plan["modified_files"]:
            print(
                f"  [WARNING] Existing file modified by {stage} cannot be "
                f"restored: {path}"
            )
        if plan.get("delete_unmanaged"):
            for path in plan.get("protected_unmanaged_files", []):
                print(
                    "  [INFO] Protected from --delete-unmanaged; preserved: "
                    f"{path}"
                )
            for path in plan.get("unmanaged_files_to_delete", []):
                print(f"  [DANGER] Delete unmanaged file: {path}")
            for path in plan.get("unmanaged_directories_to_delete", []):
                print(f"  [DANGER] Remove unmanaged directory if empty: {path}")
        else:
            for path in plan["unknown_files"]:
                print(f"  [INFO] Not managed by the manifest; preserved: {path}")
        for path in plan["files"]:
            if path not in restore_sources:
                print(f"  Delete file: {path}")
        for path in plan["directories"]:
            print(f"  Remove directory if empty: {path}")
        for path in plan.get("metadata_files", []):
            print(f"  Delete reset metadata: {path}")
        if (
            not plan["restore_files"]
            and not plan["files"]
            and not plan["directories"]
            and not plan.get("metadata_files")
            and not (
                plan.get("delete_unmanaged")
                and (
                    plan.get("unmanaged_files_to_delete")
                    or plan.get("unmanaged_directories_to_delete")
                )
            )
        ):
            print("  No artifacts to delete.")


def _folder_identity(folder):
    subject = next(
        (part for part in reversed(folder.parts) if part.startswith("sub-")),
        "sub-UNKNOWN",
    )
    session = next(
        (part for part in reversed(folder.parts) if part.startswith("ses-")),
        "ses-none",
    )
    return subject, session, folder.name


def _summary_unknown_files(plan):
    """Return unmanaged files that should be repeated in the final summary."""
    return [
        path
        for path in plan.get("unknown_files", [])
        if not path.name.endswith(MUTED_SUMMARY_SUFFIXES)
        and not _is_stroke_mask(path)
    ]


def _summary_protected_unmanaged_files(plan):
    return [
        path
        for path in plan.get("protected_unmanaged_files", [])
        if not path.name.endswith(MUTED_SUMMARY_SUFFIXES)
    ]


def print_issue_summary(plans):
    """Print all actionable warnings and information grouped by dataset."""
    stroke_mask_folders = _normal_stroke_mask_folders(plans)
    affected_plans = [
        plan
        for plan in plans
        if (
            plan.get("warnings")
            or plan.get("restore_errors")
            or plan.get("modified_files")
            or (
                plan.get("delete_unmanaged")
                and (
                    plan.get("unmanaged_files_to_delete")
                    or plan.get("unmanaged_directories_to_delete")
                    or _summary_protected_unmanaged_files(plan)
                )
            )
            or (
                not plan.get("delete_unmanaged")
                and _summary_unknown_files(plan)
            )
            or plan.get("runtime_warnings")
        )
    ]

    print("\nWarning and information summary")
    if stroke_mask_folders:
        folder_word = "folder" if len(stroke_mask_folders) == 1 else "folders"
        print(
            f"[INFO] Stroke mask files are present in "
            f"{len(stroke_mask_folders)} selected {folder_word} and were preserved."
        )
    if not affected_plans:
        if not stroke_mask_folders:
            print("No warnings or reportable unmanaged files found.")
        return

    for plan in affected_plans:
        subject, session, mode = _folder_identity(plan["folder"])
        print(f"\n{subject} | {session} | {mode}")
        for warning in plan.get("warnings", []):
            print(f"  [WARNING] {warning}")
        for error in plan.get("restore_errors", []):
            print(f"  [ERROR] Base restore blocked: {error}")
        for stage, path in plan.get("modified_files", []):
            print(
                f"  [WARNING] Modified by {stage}; original content cannot be "
                f"restored: {path}"
            )
        if plan.get("delete_unmanaged"):
            deleted_unmanaged = set(plan.get("deleted_unmanaged_files", []))
            removed_unmanaged = set(
                plan.get("removed_unmanaged_directories", [])
            )
            for path in plan.get("unmanaged_files_to_delete", []):
                action = (
                    "deleted"
                    if path in deleted_unmanaged
                    else "selected for deletion"
                )
                print(f"  [WARNING] Unmanaged file {action}: {path}")
            for path in _summary_protected_unmanaged_files(plan):
                print(
                    "  [INFO] Protected from --delete-unmanaged; preserved: "
                    f"{path}"
                )
            for path in plan.get("unmanaged_directories_to_delete", []):
                action = (
                    "removed"
                    if path in removed_unmanaged
                    else "selected for removal if empty"
                )
                print(f"  [WARNING] Unmanaged directory {action}: {path}")
        else:
            for path in _summary_unknown_files(plan):
                print(
                    f"  [INFO] Not managed by the manifest; review manually: {path}"
                )
        for warning in plan.get("runtime_warnings", []):
            print(f"  [WARNING] {warning}")

    if any(
        plan.get("unmanaged_files_to_delete")
        or plan.get("unmanaged_directories_to_delete")
        for plan in affected_plans
    ):
        print(
            "\n[DANGER] Review every unmanaged path above. These files are "
            "not recoverable through the manifest after deletion."
        )
    else:
        print(
            "\nReview the listed files and directories and intervene manually "
            "if necessary."
        )


def _relative(folder, path):
    return path.relative_to(folder).as_posix()


def _update_manifest_after_reset(plan, deleted_files, removed_directories):
    # Remove only entries that were successfully deleted or were already absent.
    manifest = plan["manifest"]
    folder = plan["folder"]
    deleted_relative = {_relative(folder, path) for path in deleted_files}
    removed_relative = {_relative(folder, path) for path in removed_directories}
    deleted_relative |= plan["missing_files"]
    removed_relative |= plan["missing_directories"]
    restored_relative = plan.get("restored_destinations", set())

    for stage in plan["selected_stages"]:
        stage_data = manifest.get("stages", {}).get(stage, {})
        stage_data["created_files"] = [
            path for path in stage_data.get("created_files", [])
            if path not in deleted_relative
        ]
        stage_data["created_directories"] = [
            path for path in stage_data.get("created_directories", [])
            if path not in removed_relative
        ]
        stage_data["modified_files"] = [
            path for path in stage_data.get("modified_files", [])
            if path not in restored_relative and path not in deleted_relative
        ]
        if not any(
            stage_data.get(key)
            for key in ("created_files", "created_directories", "modified_files")
        ):
            manifest["stages"].pop(stage, None)

    write_manifest(folder, manifest)


def _restore_base_files(plan):
    """Atomically replace working NIfTIs while retaining backups until success."""
    restored_files = set()

    for source, destination in plan.get("restore_files", []):
        temporary_path = None
        try:
            # Recheck immediately before copying in case the dataset changed
            # after the reset plan was displayed.
            if not source.is_file() or source.is_symlink():
                raise OSError("backup is missing or is not a regular file")
            if destination.is_symlink() or destination.is_dir():
                raise OSError("restore destination is not a regular file path")

            file_descriptor, temporary_name = tempfile.mkstemp(
                prefix=f".{destination.name}.restore.",
                suffix=".tmp",
                dir=destination.parent,
            )
            os.close(file_descriptor)
            temporary_path = Path(temporary_name)

            # The temporary file lives beside the destination, so os.replace()
            # performs one atomic replacement even when the destination exists.
            shutil.copy2(source, temporary_path)
            with temporary_path.open("rb") as stream:
                os.fsync(stream.fileno())
            os.replace(temporary_path, destination)
            restored_files.add(destination)
            print(f"Restored original NIfTI: {source} -> {destination}")
        except OSError as exc:
            warning = (
                f"Original NIfTI could not be restored: {source} -> "
                f"{destination}: {exc}"
            )
            plan["runtime_warnings"].append(warning)
            print(f"[ERROR] {warning}")
            return False, restored_files
        finally:
            if temporary_path is not None:
                try:
                    temporary_path.unlink()
                except FileNotFoundError:
                    pass

    return True, restored_files


def _delete_base_metadata(plan):
    """Delete the manifest and its lock after a completed base reset."""
    deleted_count = 0
    succeeded = True

    for path in plan.get("metadata_files", []):
        try:
            path.unlink()
            deleted_count += 1
            print(f"Deleted reset metadata: {path}")
        except FileNotFoundError:
            continue
        except OSError as exc:
            succeeded = False
            warning = f"Reset metadata could not be deleted: {path}: {exc}"
            plan["runtime_warnings"].append(warning)
            print(f"[WARNING] {warning}")

    return succeeded, deleted_count


def _delete_unmanaged_paths(plan):
    """Delete explicitly requested unmanaged files and then empty directories."""
    deleted_files = 0
    removed_directories = 0
    succeeded = True

    for path in plan.get("unmanaged_files_to_delete", []):
        try:
            path.unlink()
            deleted_files += 1
            plan["deleted_unmanaged_files"].append(path)
            print(f"Deleted unmanaged file: {path}")
        except FileNotFoundError:
            continue
        except OSError as exc:
            succeeded = False
            warning = f"Unmanaged file could not be deleted: {path}: {exc}"
            plan["runtime_warnings"].append(warning)
            print(f"[WARNING] {warning}")

    # Directories are never removed recursively. Deepest paths are attempted
    # first, after their unmanaged files have been deleted individually.
    for path in plan.get("unmanaged_directories_to_delete", []):
        try:
            path.rmdir()
            removed_directories += 1
            plan["removed_unmanaged_directories"].append(path)
            print(f"Removed empty unmanaged directory: {path}")
        except FileNotFoundError:
            continue
        except OSError as exc:
            succeeded = False
            warning = (
                f"Unmanaged directory could not be removed and is preserved: "
                f"{path}: {exc}"
            )
            plan["runtime_warnings"].append(warning)
            print(f"[WARNING] {warning}")

    return succeeded, deleted_files, removed_directories


def apply_reset_plan(plans):
    restore_errors = _base_restore_errors(plans)
    if restore_errors:
        print("\nBase reset aborted: Required originals cannot be restored.")
        print_issue_summary(plans)
        return False

    deleted_files_count = 0
    deleted_dirs_count = 0
    restored_files_count = 0
    deleted_metadata_count = 0
    deleted_unmanaged_files_count = 0
    removed_unmanaged_dirs_count = 0
    reset_succeeded = True

    # Restore all modality folders before deleting any backup or derived file.
    for plan in plans:
        if plan["skip"] or not plan.get("restore_files"):
            continue
        restore_succeeded, restored_files = _restore_base_files(plan)
        restored_files_count += len(restored_files)
        if not restore_succeeded:
            print(
                "\nBase reset stopped before deleting registered artifacts. "
                "All brkraw backups were preserved."
            )
            print_issue_summary(plans)
            return False

    for plan in plans:
        if plan["skip"]:
            continue
        deleted_files = set()
        removed_directories = set()

        for path in plan["files"]:
            try:
                path.unlink()
                deleted_files.add(path)
                deleted_files_count += 1
                print(f"Deleted: {path}")
            except FileNotFoundError:
                deleted_files.add(path)
            except OSError as exc:
                warning = f"File could not be deleted: {path}: {exc}"
                plan["runtime_warnings"].append(warning)
                print(f"[WARNING] {warning}")

        if plan.get("delete_unmanaged"):
            (
                unmanaged_succeeded,
                unmanaged_files_count,
                unmanaged_dirs_count,
            ) = _delete_unmanaged_paths(plan)
            deleted_unmanaged_files_count += unmanaged_files_count
            removed_unmanaged_dirs_count += unmanaged_dirs_count
            reset_succeeded = reset_succeeded and unmanaged_succeeded

        # Tracked directories are also removed only when empty.
        for path in plan["directories"]:
            try:
                path.rmdir()
                removed_directories.add(path)
                deleted_dirs_count += 1
                print(f"Removed empty directory: {path}")
            except FileNotFoundError:
                removed_directories.add(path)
            except OSError:
                warning = (
                    f"Directory is not empty and is preserved: {path}. "
                    "Unmanaged content was not deleted."
                )
                plan["runtime_warnings"].append(warning)
                print(f"[WARNING] {warning}")

        _update_manifest_after_reset(plan, deleted_files, removed_directories)
        if plan.get("phase") == "base":
            metadata_succeeded, metadata_count = _delete_base_metadata(plan)
            deleted_metadata_count += metadata_count
            reset_succeeded = reset_succeeded and metadata_succeeded

    print("\nDone!")
    if any(plan.get("phase") == "base" for plan in plans):
        print(f"Restored original NIfTI files: {restored_files_count}")
        print(f"Deleted reset metadata files: {deleted_metadata_count}")
    print(f"Deleted files: {deleted_files_count}")
    print(f"Removed empty directories: {deleted_dirs_count}")
    if any(plan.get("delete_unmanaged") for plan in plans):
        print(f"Deleted unmanaged files: {deleted_unmanaged_files_count}")
        print(f"Removed unmanaged directories: {removed_unmanaged_dirs_count}")
    print_issue_summary(plans)
    return reset_succeeded


def reset_folder(
    project_root,
    mode="anat",
    phase="base",
    dry_run=False,
    delete_unmanaged=False,
):
    """Plan and optionally apply a safe manifest-based project reset."""
    if mode not in MODES:
        raise ValueError("mode must be 'anat', 'dwi', 'func', or 't2map'")
    manifest_errors = validate_project_manifests(project_root, mode)
    if manifest_errors:
        raise ResetError(
            "Not all existing folders of the selected modality have a valid manifest:\n"
            + "\n".join(manifest_errors)
        )
    phase = normalize_phase(phase)
    plans = build_reset_plan(
        project_root,
        mode,
        phase,
        delete_unmanaged=delete_unmanaged,
    )
    print_reset_plan(plans, mode, phase)
    restore_errors = _base_restore_errors(plans)
    if restore_errors:
        raise ResetError(
            "Base reset cannot safely restore all original NIfTI files:\n"
            + "\n".join(restore_errors)
        )
    if dry_run:
        print_issue_summary(plans)
    else:
        apply_reset_plan(plans)
    return plans


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Reset anat, dwi, func, or t2map folders using manifests generated "
            "by the processing scripts."
        ),
        epilog=(
            "Example: %(prog)s --input /path/to/project --mode anat "
            "--phase preprocessing --dry-run"
        ),
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the BIDS-like project directory",
    )
    parser.add_argument(
        "-m",
        "--mode",
        required=True,
        choices=MODES,
        help="Modality to reset",
    )
    parser.add_argument(
        "-p",
        "--phase",
        required=True,
        choices=("base", "preprocessing", "registration", "processing"),
        help="Processing phase to reset to",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be deleted without changing anything",
    )
    parser.add_argument(
        "--delete-unmanaged",
        action="store_true",
        help=(
            "Permanently delete files and empty directories not known to the "
            "manifest (dangerous)"
        ),
    )
    parser.add_argument(
        "-y",
        "--yes",
        action="store_true",
        help="Skip the confirmation prompt",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    project_root = Path(args.input).expanduser().resolve()

    print(f"Checking BIDS structure below: {project_root}")
    if not is_bids_like(project_root):
        print("Aborting because the BIDS structure is missing or invalid.")
        return 1
    manifest_errors = validate_project_manifests(project_root, args.mode)
    if manifest_errors:
        print(
            f"Error: Not all existing {args.mode} folders contain a valid "
            f"*{MANIFEST_SUFFIX}:"
        )
        for error in manifest_errors:
            print(f"  - {error}")
        print("Aborted: No changes were made.")
        return 1
    if not has_any_brkraw(project_root):
        print("Error: No 'brkraw' directory was found below project_root.")
        return 1

    try:
        plans = build_reset_plan(
            project_root,
            args.mode,
            args.phase,
            delete_unmanaged=args.delete_unmanaged,
        )
    except ResetError as exc:
        print(f"Error: {exc}")
        return 1

    print_reset_plan(plans, args.mode, args.phase)
    restore_errors = _base_restore_errors(plans)
    if restore_errors:
        print("\nBase reset aborted: Required originals cannot be restored.")
        print("No changes were made.")
        print_issue_summary(plans)
        return 1

    if args.dry_run:
        print("\nDry-run: No changes were made.")
        print_issue_summary(plans)
        return 0

    has_actions = any(
        not plan["skip"]
        and (
            plan.get("restore_files")
            or plan["files"]
            or plan["directories"]
            or plan.get("metadata_files")
            or (
                plan.get("delete_unmanaged")
                and (
                    plan.get("unmanaged_files_to_delete")
                    or plan.get("unmanaged_directories_to_delete")
                )
            )
        )
        for plan in plans
    )
    if not has_actions:
        print("\nNo deletion operations are required.")
        print_issue_summary(plans)
        return 0

    if not args.yes:
        expected_confirmation = (
            "DELETE UNMANAGED" if args.delete_unmanaged else "Yes"
        )
        prompt = (
            "[DANGER] Type 'DELETE UNMANAGED' to delete all displayed "
            "unmanaged paths: "
            if args.delete_unmanaged
            else "Type 'Yes' to apply this reset plan: "
        )
        try:
            confirmation = input(prompt)
        except (EOFError, KeyboardInterrupt):
            print("\nAborted: No changes were made.")
            print_issue_summary(plans)
            return 1
        if confirmation != expected_confirmation:
            print("Aborted: No changes were made.")
            print_issue_summary(plans)
            return 1

    return 0 if apply_reset_plan(plans) else 1


if __name__ == "__main__":
    sys.exit(main())
