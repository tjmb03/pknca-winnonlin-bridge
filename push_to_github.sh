#!/usr/bin/env bash
#
# Push this repository to GitHub.
#
#   ./push_to_github.sh <github-username> [repo-name]
#
# Re-runnable. Verifies the working tree before doing anything, and stops
# rather than force-pushing over an existing remote history.

set -euo pipefail

USER="${1:-}"
REPO="${2:-pknca-winnonlin-bridge}"

if [[ -z "$USER" ]]; then
  echo "usage: ./push_to_github.sh <github-username> [repo-name]" >&2
  exit 1
fi

echo "==> Pre-flight"

# 1. Identity must be real, not the placeholder baked in at packaging time.
NAME=$(git config user.name  || echo "")
MAIL=$(git config user.email || echo "")
if [[ "$NAME" == "REPLACE ME" || "$MAIL" == "replace@example.com" || -z "$NAME" ]]; then
  echo "    git identity is unset or still the placeholder." >&2
  echo "    Run:" >&2
  echo "      git config user.name  'Your Name'" >&2
  echo "      git config user.email 'you@example.com'" >&2
  echo "      git commit --amend --reset-author --no-edit" >&2
  exit 1
fi
echo "    author: $NAME <$MAIL>"

# 2. Nothing uncommitted.
if [[ -n "$(git status --porcelain)" ]]; then
  echo "    working tree is dirty:" >&2
  git status --short >&2
  exit 1
fi
echo "    working tree clean"

# 3. No secrets. The script excludes itself: it contains the search pattern
#    literally, which would otherwise always self-match.
if git grep -nE "api[_-]?key|password|BEGIN [A-Z ]*PRIVATE KEY|AKIA[0-9A-Z]{16}" \
     -- '*.R' '*.py' '*.sh' ':(exclude)push_to_github.sh' 2>/dev/null | grep -q .; then
  echo "    possible credential found in tracked files -- aborting" >&2
  git grep -nE "api[_-]?key|password|BEGIN [A-Z ]*PRIVATE KEY|AKIA[0-9A-Z]{16}" \
     -- '*.R' '*.py' '*.sh' ':(exclude)push_to_github.sh' >&2
  exit 1
fi
echo "    no credentials in tracked files"

# 4. Nothing oversized.
BIG=$(git ls-files -z | xargs -0 du -k 2>/dev/null | awk '$1 > 51200 {print $2}')
if [[ -n "$BIG" ]]; then
  echo "    files over 50MB:" >&2; echo "$BIG" >&2; exit 1
fi
echo "    no oversized files ($(git ls-files | wc -l | tr -d ' ') tracked)"

# 5. Tests, if R is available.
if command -v Rscript >/dev/null 2>&1; then
  echo "==> Tests"
  Rscript -e 'q(status = as.integer(any(as.data.frame(
    testthat::test_dir("tests/testthat", reporter = "silent"))$failed > 0)))' \
    && echo "    passing" || { echo "    FAILING -- aborting" >&2; exit 1; }
fi

# 6. Remote.
echo "==> Remote"
URL="https://github.com/${USER}/${REPO}.git"
if git remote get-url origin >/dev/null 2>&1; then
  echo "    origin already set: $(git remote get-url origin)"
else
  git remote add origin "$URL"
  echo "    added origin: $URL"
fi

# 7. Refuse to clobber an existing remote branch.
if git ls-remote --exit-code --heads origin main >/dev/null 2>&1; then
  echo "    remote 'main' already exists. Pull and merge rather than force-push:" >&2
  echo "      git pull --rebase origin main && git push origin main" >&2
  exit 1
fi

echo "==> Pushing"
echo "    Create the repository first if it does not exist:"
echo "      https://github.com/new  (name: ${REPO}, no README/licence/gitignore)"
echo "    or:  gh repo create ${USER}/${REPO} --public --source=. --remote=origin"
echo
git push -u origin main
echo
echo "Done: https://github.com/${USER}/${REPO}"
