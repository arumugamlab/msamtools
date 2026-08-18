#!/bin/sh

set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

fixture="$script_dir/fixtures/summary.sam"
count_fixture="$script_dir/fixtures/summary_count.sam"
edge_fixture="$script_dir/fixtures/summary_edge.sam"
long_qname_fixture="$script_dir/fixtures/long_qname.sam"
cigar_eqx_fixture="$script_dir/fixtures/cigar_eqx.sam"
tmpdir=${TMPDIR:-/tmp}/msamtools-test-summary-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir" || exit 1

run_summary()
{
    name=$1
    input=$2
    shift 2
    stdout="$tmpdir/$name.stdout"
    stderr="$tmpdir/$name.stderr"

    if ! "$MSAMTOOLS" summary -S "$@" "$input" \
        >"$stdout" 2>"$stderr"; then
        cat "$stderr" >&2
        fail "summary case '$name' should succeed"
    fi
    pass_check "summary case '$name' should succeed"
}

run_summary default "$fixture"
expected_default="$tmpdir/default.expected"
cat >"$expected_default" <<'EOF'
edge	100	A	100	100	100.0
ins	100	A	100	98	98.0
mismatch	100	A	100	99	99.0
perfect	100	A	100	100	100.0
soft	100	A	100	90	90.0
supplementary	100	A	100	100	100.0
EOF
if ! cmp -s "$expected_default" "$tmpdir/default.stdout"; then
    echo "Expected summary output:" >&2
    cat "$expected_default" >&2
    echo "Observed summary output:" >&2
    cat "$tmpdir/default.stdout" >&2
    fail "default summary should report exact read-level statistics"
fi
pass_check "default summary should report exact read-level statistics"
assert_empty "$tmpdir/default.stderr" \
    "default summary should not write to stderr"

run_summary edge10 "$fixture" -e 10
assert_not_contains "$tmpdir/edge10.stdout" "edge	" \
    "--edge should exclude alignments near sequence boundaries"
assert_contains "$tmpdir/edge10.stdout" "perfect	100	A	100	100	100.0" \
    "--edge should retain internal alignments"

run_summary stats_mapped "$fixture" --stats mapped
cat >"$tmpdir/stats_mapped.expected" <<'EOF'
90	1
98	1
99	1
100	3
EOF
if ! cmp -s "$tmpdir/stats_mapped.expected" "$tmpdir/stats_mapped.stdout"; then
    fail "--stats mapped should report the exact mapped-base distribution"
fi
pass_check "--stats mapped should report the exact mapped-base distribution"

run_summary stats_unmapped "$fixture" --stats unmapped
cat >"$tmpdir/stats_unmapped.expected" <<'EOF'
0	3
1	1
2	1
10	1
EOF
if ! cmp -s "$tmpdir/stats_unmapped.expected" "$tmpdir/stats_unmapped.stdout"; then
    fail "--stats unmapped should report the exact unmapped-base distribution"
fi
pass_check "--stats unmapped should report the exact unmapped-base distribution"

run_summary stats_edit "$fixture" --stats edit
if ! cmp -s "$tmpdir/stats_unmapped.expected" "$tmpdir/stats_edit.stdout"; then
    fail "--stats edit should report the exact edit-distance distribution"
fi
pass_check "--stats edit should report the exact edit-distance distribution"

run_summary stats_score "$fixture" --stats score
cat >"$tmpdir/stats_score.expected" <<'EOF'
80	1
96	1
98	1
100	3
EOF
if ! cmp -s "$tmpdir/stats_score.expected" "$tmpdir/stats_score.stdout"; then
    fail "--stats score should report the exact score distribution"
fi
pass_check "--stats score should report the exact score distribution"

run_summary edge_boundary "$edge_fixture" -e 10
cat >"$tmpdir/edge_boundary.expected" <<'EOF'
left_pass	100	A	100	100	100.0
right_pass	100	A	100	100	100.0
EOF
if ! cmp -s "$tmpdir/edge_boundary.expected" "$tmpdir/edge_boundary.stdout"; then
    echo "Expected exact edge-boundary output:" >&2
    cat "$tmpdir/edge_boundary.expected" >&2
    echo "Observed edge-boundary output:" >&2
    cat "$tmpdir/edge_boundary.stdout" >&2
    fail "--edge should retain exactly alignments with at least the requested flanking bases"
fi
pass_check "--edge should retain exactly alignments with at least the requested flanking bases"

count_stdout="$tmpdir/count.stdout"
count_stderr="$tmpdir/count.stderr"
if ! "$MSAMTOOLS" summary -S -c "$count_fixture" \
    >"$count_stdout" 2>"$count_stderr"; then
    cat "$count_stderr" >&2
    fail "summary --count should succeed"
fi
pass_check "summary --count should succeed"
if [ "$(cat "$count_stdout")" != "3" ]; then
    fail "--count should count mapped QNAME groups exactly once"
fi
assert_empty "$count_stderr" \
    "--count should not write to stderr"

long_count_stdout="$tmpdir/long_qname_count.stdout"
long_count_stderr="$tmpdir/long_qname_count.stderr"
if ! "$MSAMTOOLS" summary -S -c "$long_qname_fixture" \
    >"$long_count_stdout" 2>"$long_count_stderr"; then
    cat "$long_count_stderr" >&2
    fail "summary --count should accept full-length legal QNAMEs"
fi
pass_check "summary --count should accept full-length legal QNAMEs"
if test "$(cat "$long_count_stdout")" != "4"; then
    echo "Expected insert count: 4" >&2
    echo -n "Observed insert count: " >&2
    cat "$long_count_stdout" >&2
    fail "summary --count should group complete QNAMEs without truncation"
fi
pass_check "summary --count should group complete QNAMEs without truncation"
assert_empty "$long_count_stderr" \
    "long-QNAME count should not write to stderr"

run_summary cigar_eqx "$cigar_eqx_fixture"
cat >"$tmpdir/cigar_eqx.expected" <<'EOF'
md_eqx0	10	A	10	0	0.0
md_eqx100	100	A	100	100	100.0
md_eqx98	100	A	100	98	98.0
md_m0	10	A	10	0	0.0
md_m100	100	A	100	100	100.0
md_m98	100	A	100	98	98.0
EOF
if ! cmp -s "$tmpdir/cigar_eqx.expected" "$tmpdir/cigar_eqx.stdout"; then
    echo "Expected M/=/X summary output:" >&2
    cat "$tmpdir/cigar_eqx.expected" >&2
    echo "Observed M/=/X summary output:" >&2
    cat "$tmpdir/cigar_eqx.stdout" >&2
    fail "summary should report equivalent statistics for M and =/X CIGARs"
fi
pass_check "summary should report equivalent statistics for M and =/X CIGARs"

report_checks
exit 0
