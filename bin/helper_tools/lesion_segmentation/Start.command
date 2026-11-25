#!/bin/bash
set -e
cd -- "$(dirname -- "$0")"
if [ ! -x .venv/bin/python ]; then
  echo "Run Setup.command once before starting."
  exit 1
fi
exec .venv/bin/python scripts/start_validation_guidance.py --dataset examples --port 8780
