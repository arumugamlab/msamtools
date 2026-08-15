#!/bin/sh
set -u
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

fixture="$script_dir/fixtures/coverage.sam"
tmpdir=${TMPDIR:-/tmp}/msamtools-test-coverage-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir" || exit 1

run_coverage()
{
    name=$1
    shift
    output="$tmpdir/$name.gz"
    stdout="$tmpdir/$name.stdout"
    stderr="$tmpdir/$name.stderr"

    if ! "$MSAMTOOLS" coverage -S -o "$output" "$@" "$fixture" >"$stdout" 2>"$stderr"; then
        cat "$stderr" >&2
        fail "coverage case '$name' should succeed"
    fi
    pass_check "coverage case '$name' should succeed"
    assert_empty "$stdout" "coverage case '$name' should not write to stdout"
    assert_empty "$stderr" "coverage case '$name' should not write to stderr"
}

run_coverage positions -w 4
gzip -cd -- "$tmpdir/positions.gz" >"$tmpdir/positions.txt" ||
    fail "per-position coverage output should be valid gzip"
cat >"$tmpdir/positions.expected" <<'EOF'
>A
1 0 1 1
2 2 1 0
0 1
>B
0 0 0 0
0
>C
0 0 0 0
0 0 0 0
0 1
EOF
if ! cmp -s "$tmpdir/positions.expected" "$tmpdir/positions.txt"; then
    echo "Expected:" >&2; cat "$tmpdir/positions.expected" >&2
    echo "Observed:" >&2; cat "$tmpdir/positions.txt" >&2
    fail "coverage should report exact per-position values and word wrapping"
fi
pass_check "coverage should report exact per-position values and word wrapping"

run_coverage summary --summary
gzip -cd -- "$tmpdir/summary.gz" >"$tmpdir/summary.txt" ||
    fail "coverage summary output should be valid gzip"
cat >"$tmpdir/summary.expected" <<'EOF'
A	0.70000000	0.90
B	0	0
C	0.10000000	0.10
EOF
if ! cmp -s "$tmpdir/summary.expected" "$tmpdir/summary.txt"; then
    echo "Expected:" >&2; cat "$tmpdir/summary.expected" >&2
    echo "Observed:" >&2; cat "$tmpdir/summary.txt" >&2
    fail "coverage summary should include exact covered fractions and mean coverage"
fi
pass_check "coverage summary should include exact covered fractions and mean coverage"

run_coverage skip --summary --skipuncovered
gzip -cd -- "$tmpdir/skip.gz" >"$tmpdir/skip.txt" ||
    fail "skip-uncovered coverage output should be valid gzip"
cat >"$tmpdir/skip.expected" <<'EOF'
A	0.70000000	0.90
C	0.10000000	0.10
EOF
if ! cmp -s "$tmpdir/skip.expected" "$tmpdir/skip.txt"; then
    fail "--skipuncovered should omit only references with zero coverage"
fi
pass_check "--skipuncovered should omit only references with zero coverage"

report_checks
exit 0
