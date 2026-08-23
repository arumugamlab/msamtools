#!/bin/sh

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)

if ! MSAMTOOLS="$MSAMTOOLS" "$repo_dir/scripts/update_cli_docs.sh" --check; then
    echo "FAIL: generated CLI documentation is out of date" >&2
    exit 1
fi

echo "PASS: generated CLI documentation matches current --help output"
exit 0
