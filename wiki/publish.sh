#!/usr/bin/env bash
#
# Publishes the pages in this directory to the GitHub wiki.
#
# The GitHub wiki is a separate git repository from the main one, so the
# pages are kept here (where they are reviewed alongside the code they
# describe) and pushed across by this script.
#
# Usage, from the root of the repository:
#
#     ./wiki/publish.sh
#
set -euo pipefail

REMOTE="${WIKI_REMOTE:-https://github.com/johodges/pyfdstools.wiki.git}"
SOURCE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLONE="$(mktemp -d)"
trap 'rm -rf "$CLONE"' EXIT

echo "Cloning $REMOTE"
git clone --quiet "$REMOTE" "$CLONE"

# Copy every page across, then remove the two files that belong to the
# main repository rather than to the wiki.
cp "$SOURCE"/*.md "$CLONE"/
rm -f "$CLONE/README.md" "$CLONE/publish.sh"

cd "$CLONE"
if git diff --quiet HEAD -- . && [ -z "$(git status --porcelain)" ]; then
    echo "The wiki is already up to date."
    exit 0
fi

git add -A
git status --short
git commit --quiet -m "Update wiki from the main repository"
git push --quiet origin HEAD

echo
echo "Published $(ls -1 "$SOURCE"/*.md | grep -cv 'README.md') pages to the wiki."
