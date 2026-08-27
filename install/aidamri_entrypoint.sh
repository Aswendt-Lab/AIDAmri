#!/bin/sh

git_info_source="/aida/build/AIDAmri_git_information.txt"
git_info_target="/aida/DATA/AIDAmri_git_information.txt"

if [ -r "$git_info_source" ] && [ -d "/aida/DATA" ]; then
    cp "$git_info_source" "$git_info_target" 2>/dev/null || true
fi

exec "$@"
