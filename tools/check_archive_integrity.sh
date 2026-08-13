#!/usr/bin/env bash
# Verify every compressed artefact in the repository before committing or submitting.
#
# Added after a 327 kB GitHub HTML error page was committed over
# Fusions/fusions_readsupport_Huh7_discovery.tsv.gz, and two JAFFAL archives were
# committed truncated or empty. All three passed unnoticed because nothing tested them.
#
# Usage:   tools/check_archive_integrity.sh          # scan the whole working tree
#          tools/check_archive_integrity.sh --staged # scan only staged files
#
# To run it automatically, point git at a hooks directory holding a pre-commit that
# calls this script with --staged:
#   git config core.hooksPath .githooks

set -uo pipefail

repo_root=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
cd "$repo_root" || exit 1

if [ "${1:-}" = "--staged" ]; then
  mapfile -t files < <(git diff --cached --name-only --diff-filter=ACM | grep -E '\.(gz|xz|zip|bz2)$' || true)
else
  mapfile -t files < <(find . -path ./.git -prune -o -type f \
    \( -name '*.gz' -o -name '*.xz' -o -name '*.zip' -o -name '*.bz2' \) -print | sort)
fi

failed=0
checked=0

for f in "${files[@]}"; do
  [ -f "$f" ] || continue
  checked=$((checked + 1))
  case "$f" in
    *.gz)  tester=(gzip -t) ;;
    *.xz)  tester=(xz -t) ;;
    *.bz2) tester=(bzip2 -t) ;;
    *.zip) tester=(unzip -qt) ;;
  esac
  if ! "${tester[@]}" "$f" >/dev/null 2>&1; then
    echo "CORRUPT  $f  ($(file -b "$f"))"
    failed=$((failed + 1))
  fi
  # A zero-byte archive tests as neither valid nor obviously corrupt on every tool.
  if [ ! -s "$f" ]; then
    echo "EMPTY    $f"
    failed=$((failed + 1))
  fi
done

if [ "$failed" -gt 0 ]; then
  echo "FAILED: $failed of $checked archives are corrupt or empty."
  exit 1
fi

echo "OK: $checked archives verified."
