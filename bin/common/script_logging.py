"""Mirror script output to a log file during standalone execution."""

import os
import shlex
import sys


DISABLE_LOG_ENV = "AIDAMRI_DISABLE_SCRIPT_LOG"
_LOG_FILES = []


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


def setup_script_logging(output_dir, log_name):
    """Enable file logging unless batchProc.py already captures the output."""
    if script_logging_disabled():
        return

    log_file = open(os.path.join(output_dir, log_name), "w", encoding="utf-8")
    _LOG_FILES.append(log_file)
    sys.stdout = _TeeStream(sys.stdout, log_file)
    sys.stderr = _TeeStream(sys.stderr, log_file)
    log_explicit_cli_options()
