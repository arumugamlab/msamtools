#!/bin/sh

set -eu

usage()
{
    cat <<'EOF'
Usage:
  scripts/update_cli_docs.sh
  scripts/update_cli_docs.sh --check

Without arguments, update generated command-line help blocks in docs/.
With --check, verify that the committed blocks match the current executable
without modifying any files.

MSAMTOOLS may be set to override the executable path.
EOF
}

mode=update
case "${1-}" in
    "")
        ;;
    --check)
        mode=check
        ;;
    -h|--help)
        usage
        exit 0
        ;;
    *)
        usage >&2
        exit 2
        ;;
esac

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)

MSAMTOOLS=${MSAMTOOLS:-"$repo_dir/msamtools"}

if test ! -x "$MSAMTOOLS"; then
    echo "ERROR: msamtools executable not found or not executable: $MSAMTOOLS" >&2
    echo "Build msamtools first, or set MSAMTOOLS=/path/to/msamtools." >&2
    exit 1
fi

tmpdir=${TMPDIR:-/tmp}/msamtools-cli-docs-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir"

status=0

for command in filter profile coverage summary; do
    doc="$repo_dir/docs/$command.md"
    begin="<!-- BEGIN GENERATED HELP: $command -->"
    end="<!-- END GENERATED HELP: $command -->"
    help="$tmpdir/$command.help"
    expected="$tmpdir/$command.expected"
    current="$tmpdir/$command.current"

    if test ! -f "$doc"; then
        echo "ERROR: documentation file not found: $doc" >&2
        exit 1
    fi

    if ! "$MSAMTOOLS" "$command" --help >"$help"; then
        echo "ERROR: failed to run: $MSAMTOOLS $command --help" >&2
        exit 1
    fi

    {
        printf '%s\n' '```text'
        awk '1' "$help"
        printf '%s\n' '```'
    } >"$expected"

    if ! awk -v begin="$begin" -v end="$end" '
        $0 == begin {
            begin_count++
            inside = 1
            next
        }
        $0 == end {
            end_count++
            inside = 0
            next
        }
        inside {
            print
        }
        END {
            if (begin_count != 1 || end_count != 1 || inside)
                exit 2
        }
    ' "$doc" >"$current"; then
        echo "ERROR: expected exactly one generated-help marker pair in $doc" >&2
        exit 1
    fi

    if cmp -s "$expected" "$current"; then
        printf 'OK: %s\n' "$doc"
        continue
    fi

    if test "$mode" = check; then
        printf 'STALE: %s\n' "$doc" >&2
        status=1
        continue
    fi

    tmp="$doc.tmp.$$"
    trap 'rm -rf "$tmpdir"; rm -f "$tmp"' EXIT HUP INT TERM

    if ! awk -v begin="$begin" -v end="$end" -v replacement="$expected" '
        $0 == begin {
            print
            while ((getline line < replacement) > 0)
                print line
            close(replacement)
            inside = 1
            next
        }
        $0 == end {
            inside = 0
            print
            next
        }
        !inside {
            print
        }
    ' "$doc" >"$tmp"; then
        echo "ERROR: failed to update $doc" >&2
        exit 1
    fi

    mv -f "$tmp" "$doc"
    trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
    printf 'UPDATED: %s\n' "$doc"
done

if test "$status" -ne 0; then
    echo "ERROR: generated CLI documentation is out of date." >&2
    echo "Run scripts/update_cli_docs.sh and commit the resulting changes." >&2
    exit "$status"
fi
