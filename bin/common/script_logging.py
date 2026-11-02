"""Mirror script output to a log file during standalone execution."""

import atexit
import os
import shlex
import sys
from pathlib import Path


DISABLE_LOG_ENV = "AIDAMRI_DISABLE_SCRIPT_LOG"
GIT_INFO_FILENAME = "AIDAmri_git_information.txt"
BUILD_GIT_INFO_PATH = Path("/aida/build") / GIT_INFO_FILENAME
_LOG_FILES = []


def _git_information_paths(output_dir):
    paths = [BUILD_GIT_INFO_PATH]
    if output_dir is not None:
        output_path = Path(output_dir)
        # The proc folder contains sub-*/ses-*/<modality>, regardless of its name.
        for folder in (output_path, *output_path.parents):
            if folder.name.startswith("ses-") and folder.parent.name.startswith("sub-"):
                paths.append(folder.parent.parent / GIT_INFO_FILENAME)
                break
        else:
            paths.append(output_path / GIT_INFO_FILENAME)
    return paths


def append_git_information(log_file=None, output_dir=None):
    """Finish a step log with the build's Git information, even after errors."""
    if log_file is None:
        # Batch runs redirect both streams into the same log file.
        sys.stdout.flush()
        sys.stderr.flush()
        log_file = sys.stdout
    if log_file.closed:
        return
    for path in _git_information_paths(output_dir):
        try:
            contents = path.read_text(encoding="utf-8")
            break
        except (OSError, UnicodeError) as exc:
            contents = f"WARNING: Could not read AIDAmri Git information: {exc}\n"
    log_file.write("\n" + contents)
    if not contents.endswith("\n"):
        log_file.write("\n")
    log_file.flush()


def build_script_log_path(output_dir, log_name):
    """Prefix step logs with the enclosing sub/session/modality, if available."""
    output_path = Path(os.path.abspath(output_dir))
    for folder in (output_path, *output_path.parents):
        if folder.name not in {"anat", "dwi", "func", "t2map"}:
            continue
        session = folder.parent.name
        subject = folder.parent.parent.name
        if subject.startswith("sub-") and session.startswith("ses-"):
            log_name = f"{subject}_{session}_{folder.name}_{log_name}"
            break
    return os.path.join(output_dir, log_name)


def script_logging_disabled():
    """Return whether batchProc.py is responsible for capturing script output."""
    return any(
        os.environ.get(name, "").lower() in {"1", "true", "yes"}
        for name in (DISABLE_LOG_ENV, "AIDAMRI_DISABLE_PROCESS_LOG")
    )


def log_explicit_cli_options(logger=None):
    """Record only the actual CLI arguments, and only for standalone runs."""
    if script_logging_disabled():
        return

    message = (
        f"Script: {os.path.basename(sys.argv[0])}\n"
        f"Explicit command line options: {shlex.join(sys.argv[1:])}"
    )
    if logger is None:
        print(message, flush=True)
    else:
        logger.info(message)


class _TeeStream:
    def __init__(self, stream, log_file):
        self.stream = stream
        self.log_file = log_file

    def write(self, text):
        self.stream.write(text)
        self.log_file.write(text)
        return len(text)

    def flush(self):
        self.stream.flush()
        self.log_file.flush()

    def isatty(self):
        return self.stream.isatty()

    def fileno(self):
        return self.stream.fileno()

    @property
    def encoding(self):
        return self.stream.encoding

    def __getattr__(self, name):
        return getattr(self.stream, name)


def get_terminal_stream():
    """Return stdout without our log mirrors, or None for redirected output."""
    stream = sys.stdout
    while isinstance(stream, _TeeStream):
        stream = stream.stream
    return stream if stream.isatty() else None


def setup_script_logging(output_dir, log_name):
    """Enable file logging unless batchProc.py already captures the output."""
    output_dir = os.path.abspath(output_dir)
    if script_logging_disabled():
        atexit.register(append_git_information, output_dir=output_dir)
        return

    log_file = open(build_script_log_path(output_dir, log_name), "w", encoding="utf-8")
    _LOG_FILES.append(log_file)
    atexit.register(append_git_information, log_file, output_dir=output_dir)
    sys.stdout = _TeeStream(sys.stdout, log_file)
    sys.stderr = _TeeStream(sys.stderr, log_file)
    log_explicit_cli_options()
