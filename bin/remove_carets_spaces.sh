#!/bin/bash

# This script removes carets, spaces, and underscores from the input plain-text file "subject" at the specified path.
# Usage: ./remove_carets_spaces.sh /path/to/subject
input_file="${1}"

if [[ -z "$input_file" ]]; then
    echo "Usage: $0 /path/to/input_file"
    exit 1
fi

subject_dir=$(dirname `realpath "$input_file"`)
echo "Processing file in directory: $subject_dir"

# Check for existence of carets, spaces, underscores in the input file
if [[ ! -s "$input_file" ]]; then
    echo "Input file is empty or does not exist."
    exit 1
fi

# Read from the input file or stdin and remove carets, spaces, and underscores
# Output the modified lines to a subject_clean.txt file
while IFS= read -r line; do
    echo "${line//[ ^]/}"
done < "$input_file" > subject_clean.txt

# Move original file to a backup with name subject_orig.txt
if [[ -n "$input_file" ]]; then
    cp "$input_file" subject_orig.txt
else
    cp "$input_file" subject_orig.txt
    cp  subject_clean.txt subject
fi