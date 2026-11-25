#!/bin/bash
set -e
cd -- "$(dirname -- "$0")"
PYTHON_BIN=""
for candidate in python3.12 python3.11 python3.10 python3; do
  if command -v "$candidate" >/dev/null 2>&1 && "$candidate" -c 'import sys;sys.exit(not ((3,10) <= sys.version_info[:2] <= (3,12)))'; then
    PYTHON_BIN="$candidate"
    break
  fi
done
if [ -z "$PYTHON_BIN" ]; then
  echo "Install Python 3.10–3.12, then run Setup.command again."
  exit 1
fi
"$PYTHON_BIN" -m venv .venv
.venv/bin/python -m pip install -r requirements-inference.txt
echo "Setup complete. Open Start.command to review the saved examples."
