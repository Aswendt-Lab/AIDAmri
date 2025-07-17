#!/bin/bash

# Usage: ./remove_carets_spaces.sh input.txt > output.txt

while IFS= read -r line; do
    echo "${line//[ ^]/}"
done < "${1:-/dev/stdin}"