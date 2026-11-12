#!/bin/sh
# pre-commit hook: bump the patch number of the repo version on every commit, so each commit has
# its own version. The version is the first `version = "x.y.z"` line of pyproject.toml (Python
# packages) or the single line of VERSION (mosaicSource). Only the staged copy and that one line
# of the working copy are changed, so other unstaged edits to the file stay unstaged.
# A version change already staged (a manual major/minor bump) is kept as is.
# Skipped during rebase/merge/cherry-pick, or with SKIP_VERSION_BUMP=1.
# Canonical copy: mosaicSource/scripts/preCommitVersionBump.sh; install with
#   cp preCommitVersionBump.sh <repo>/.git/hooks/pre-commit
GIT_DIR=$(git rev-parse --git-dir) || exit 0
for f in rebase-merge rebase-apply MERGE_HEAD CHERRY_PICK_HEAD; do
    [ -e "$GIT_DIR/$f" ] && exit 0
done
[ "${SKIP_VERSION_BUMP:-0}" = "1" ] && exit 0
cd "$(git rev-parse --show-toplevel)" || exit 0
if [ -f pyproject.toml ]; then
    file=pyproject.toml
    pattern='^version *= *"[0-9][0-9]*[.][0-9][0-9]*[.][0-9][0-9]*"'
elif [ -f VERSION ]; then
    file=VERSION
    pattern='^[0-9][0-9]*[.][0-9][0-9]*[.][0-9][0-9]*$'
else
    exit 0
fi
# manual bump already staged
git diff --cached -U0 -- "$file" | grep -q "^+$(echo "$pattern" | cut -c2-)" && exit 0
bump() {  # stdin -> stdout with the first version line's patch number incremented
    awk -v pat="$pattern" '
        !done && $0 ~ pat {
            match($0, /[0-9]+[.][0-9]+[.][0-9]+/)
            split(substr($0, RSTART, RLENGTH), v, ".")
            $0 = substr($0, 1, RSTART - 1) v[1] "." v[2] "." (v[3] + 1) substr($0, RSTART + RLENGTH)
            done = 1
        }
        { print }'
}
git show ":$file" > /dev/null 2>&1 || exit 0
blob=$(git show ":$file" | bump | git hash-object -w --stdin) || exit 1
git update-index --cacheinfo "100644,$blob,$file" || exit 1
tmp="$file.bump.$$"
bump < "$file" > "$tmp" && cat "$tmp" > "$file"
rm -f "$tmp"
echo "pre-commit: $(git show ":$file" | grep -m 1 -E "$pattern") ($file)"
exit 0
