#!/bin/sh

repo_dir="${1:-.}"
output_file="${2:-AIDAmri_git_information.txt}"
unknown="unknown"

sanitize() {
    value="$1"
    value=$(printf "%s" "$value" | tr "\r\n" "  ")
    if [ -n "$value" ]; then
        printf "%s" "$value"
    else
        printf "%s" "$unknown"
    fi
}

git_output() {
    git -C "$repo_dir" "$@" 2>/dev/null || true
}

repo_name_from_remote_url() {
    remote_url=$(printf "%s" "$1" | sed 's:/*$::')
    remote_name=${remote_url##*/}
    remote_name=${remote_name%.git}
    printf "%s" "$remote_name"
}

sanitize_remote_url() {
    remote_url=$(printf "%s" "$1" | tr "\r\n" "  ")
    if printf "%s" "$remote_url" | grep -Eq '^[a-zA-Z][a-zA-Z0-9+.-]*://[^/@]+@'; then
        protocol=${remote_url%%://*}
        remainder=${remote_url#*://}
        printf "%s://%s" "$protocol" "${remainder#*@}"
    else
        printf "%s" "$remote_url"
    fi
}

repo_name_from_worktree() {
    worktree_dir="$1"
    worktree_name=${worktree_dir##*/}
    printf "%s" "$worktree_name"
}

repo_name="$unknown"
remote_url="$unknown"
commit="$unknown"
branch="$unknown"
dirty="$unknown"
commit_author="$unknown"
git_config_user="${AIDAMRI_GIT_CONFIG_USER:-}"

if [ -z "$git_config_user" ] || [ "$git_config_user" = "$unknown" ]; then
    git_config_user="$unknown"
fi

if command -v git >/dev/null 2>&1 && git -C "$repo_dir" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    remote_origin=$(git_output config --get remote.origin.url)
    if [ -n "$remote_origin" ]; then
        remote_url=$(sanitize_remote_url "$remote_origin")
        repo_name=$(repo_name_from_remote_url "$remote_origin")
    fi
    if [ -z "$repo_name" ] || [ "$repo_name" = "$unknown" ]; then
        worktree=$(git_output rev-parse --show-toplevel)
        if [ -n "$worktree" ]; then
            repo_name=$(repo_name_from_worktree "$worktree")
        fi
    fi

    commit=$(git_output rev-parse --verify HEAD)
    branch=$(git_output symbolic-ref --quiet --short HEAD)
    if [ -z "$branch" ] && [ -n "$commit" ]; then
        branch="detached HEAD"
    fi

    status=$(git -C "$repo_dir" status --porcelain --untracked-files=normal 2>/dev/null)
    status_rc=$?
    if [ "$status_rc" -eq 0 ]; then
        if [ -n "$status" ]; then
            dirty="true"
        else
            dirty="false"
        fi
    fi

    commit_author=$(git_output show -s --format=%an HEAD)
    if [ "$git_config_user" = "$unknown" ]; then
        git_config_user=$(git_output config --get user.name)
    fi
fi

mkdir -p "$(dirname "$output_file")"
{
    printf "%s\n" "AIDAmri Git Information"
    printf "%s\n" "-----------------------"
    printf "Git repository          : %s\n" "$(sanitize "$repo_name")"
    printf "Git remote URL          : %s\n" "$(sanitize "$remote_url")"
    printf "Git commit              : %s\n" "$(sanitize "$commit")"
    printf "Git branch              : %s\n" "$(sanitize "$branch")"
    printf "Git uncommitted changes : %s\n" "$(sanitize "$dirty")"
    printf "Commit author           : %s\n" "$(sanitize "$commit_author")"
    printf "Git project user        : %s\n" "$(sanitize "$git_config_user")"
} > "$output_file"
