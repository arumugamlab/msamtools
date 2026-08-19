#!/bin/sh

set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

filter_fixture="$script_dir/fixtures/filter.sam"
profile_fixture="$script_dir/fixtures/profile.sam"
summary_fixture="$script_dir/fixtures/summary.sam"
coverage_fixture="$script_dir/fixtures/coverage.sam"

tmpdir=${TMPDIR:-/tmp}/msamtools-test-errors-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir" || exit 1

# Unknown top-level commands must fail clearly.
if "$MSAMTOOLS" definitely-not-a-command \
    >"$tmpdir/unknown.stdout" 2>"$tmpdir/unknown.stderr"; then
    fail "unknown top-level command should fail"
fi
pass_check "unknown top-level command should fail"
assert_contains "$tmpdir/unknown.stderr" "unrecognized command" \
    "unknown-command error should explain the problem"

# A filter operation is mandatory.  In particular, -h must not cause a SAM
# header to be emitted before this CLI validation fails.
if "$MSAMTOOLS" filter -S -h "$filter_fixture" \
    >"$tmpdir/no_filter.stdout" 2>"$tmpdir/no_filter.stderr"; then
    fail "filter without a filtering operation should fail"
fi
pass_check "filter without a filtering operation should fail"
assert_contains "$tmpdir/no_filter.stdout" "needs -l, -p, --ppt, -z, --besthit or --uniqhit" \
    "missing-filter error should explain the required operation"
assert_not_contains "$tmpdir/no_filter.stdout" "@HD" \
    "filter CLI failure should occur before a SAM header is written"

# Mutually exclusive special-filter modes must be rejected.
if "$MSAMTOOLS" filter -S --besthit --uniqhit "$filter_fixture" \
    >"$tmpdir/exclusive.stdout" 2>"$tmpdir/exclusive.stderr"; then
    fail "--besthit and --uniqhit together should fail"
fi
pass_check "--besthit and --uniqhit together should fail"
assert_contains "$tmpdir/exclusive.stdout" "--besthit cannot be combined with --uniqhit" \
    "mutually exclusive best-hit modes should explain the conflict"

if "$MSAMTOOLS" filter -S -p 95 --ppt 950 "$filter_fixture" \
    >"$tmpdir/exclusive.stdout" 2>"$tmpdir/exclusive.stderr"; then
    fail "-p and --ppt together should fail"
fi
pass_check "-p and --ppt together should fail"
assert_contains "$tmpdir/exclusive.stdout" "-p cannot be combined with --ppt" \
    "mutually exclusive percent identity options should explain the conflict"

# Numeric filter ranges are validated before processing.
if "$MSAMTOOLS" filter -S -p 101 "$filter_fixture" \
    >"$tmpdir/pid.stdout" 2>"$tmpdir/pid.stderr"; then
    fail "out-of-range percent identity should fail"
fi
pass_check "out-of-range percent identity should fail"
assert_contains "$tmpdir/pid.stdout" "-p must be in the range [0,100]" \
    "out-of-range identity error should explain the allowed range"

if "$MSAMTOOLS" filter -S -p -1 "$filter_fixture" \
    >"$tmpdir/pid_negative.stdout" 2>"$tmpdir/pid_negative.stderr"; then
    fail "negative percent identity should fail"
fi
pass_check "negative percent identity should fail"
assert_contains "$tmpdir/pid_negative.stdout" "-p must be in the range [0,100]" \
    "negative identity error should explain the allowed range"

if "$MSAMTOOLS" filter -S -l -1 "$filter_fixture" \
    >"$tmpdir/length_negative.stdout" 2>"$tmpdir/length_negative.stderr"; then
    fail "negative alignment length should fail"
fi
pass_check "negative alignment length should fail"
assert_contains "$tmpdir/length_negative.stdout" "-l must be a non-negative integer" \
    "negative alignment-length error should explain the allowed values"

if "$MSAMTOOLS" filter -S -z 101 "$filter_fixture" \
    >"$tmpdir/qfrac.stdout" 2>"$tmpdir/qfrac.stderr"; then
    fail "out-of-range query fraction should fail"
fi
pass_check "out-of-range query fraction should fail"
assert_contains "$tmpdir/qfrac.stdout" "-z must be in the range [0,100]" \
    "out-of-range query-fraction error should explain the allowed range"

if "$MSAMTOOLS" filter -S -z -1 "$filter_fixture" \
    >"$tmpdir/qfrac_negative.stdout" 2>"$tmpdir/qfrac_negative.stderr"; then
    fail "negative query fraction should fail"
fi
pass_check "negative query fraction should fail"
assert_contains "$tmpdir/qfrac_negative.stdout" "-z must be in the range [0,100]" \
    "negative query-fraction error should explain the allowed range"

# Profile has two mandatory output-identification arguments.
if "$MSAMTOOLS" profile -S "$profile_fixture" \
    >"$tmpdir/profile.stdout" 2>"$tmpdir/profile.stderr"; then
    fail "profile without required output arguments should fail"
fi
pass_check "profile without required output arguments should fail"
assert_contains "$tmpdir/profile.stdout" "Use --help for usage instructions" \
    "profile parse failure should provide usage guidance"

# Profile does not allow --pandas and --no-pandas together.
if "$MSAMTOOLS" profile -S \
    --label conflict \
    --unit ab \
    --nolen \
    --pandas --no-pandas \
    -o "$tmpdir/pandas_conflict.tsv.gz" \
    "$profile_fixture" \
    >"$tmpdir/pandas_conflict.stdout" \
    2>"$tmpdir/pandas_conflict.stderr"; then
    fail "--pandas and --no-pandas together should fail"
fi
pass_check "--pandas and --no-pandas together should fail"
assert_contains "$tmpdir/pandas_conflict.stdout" \
    "--pandas and --no-pandas cannot be used together" \
    "pandas output modes should explain the conflict"

# --mincount must be non-negative
if "$MSAMTOOLS" profile -S \
    --label negative_mincount \
    --mincount -1 \
    -o "$tmpdir/negative_mincount.tsv.gz" \
    "$profile_fixture" \
    >"$tmpdir/mincount_negative.stdout" \
    2>"$tmpdir/mincount_negative.stderr"; then
    fail "negative --mincount should fail"
fi
pass_check "negative --mincount should fail"
assert_contains "$tmpdir/mincount_negative.stdout" \
    "--mincount must be a non-negative integer" \
    "negative --mincount error should explain the allowed values"

if "$MSAMTOOLS" profile -S \
    --label zero_total \
    --total 0 \
    --unit ab \
    --nolen \
    -o "$tmpdir/zero_total.tsv.gz" \
    "$profile_fixture" \
    >"$tmpdir/zero_total.stdout" \
    2>"$tmpdir/zero_total.stderr"; then
    cat "$tmpdir/zero_total.stderr" >&2
    fail "--total 0 should fail"
fi
pass_check "--total 0 should fail"
assert_contains "$tmpdir/zero_total.stdout" \
    "--total must be a positive integer" \
    "negative --total error should explain that a positive value is required"

# Summary rejects unknown statistics modes.
if "$MSAMTOOLS" summary -S --stats nonsense "$summary_fixture" \
    >"$tmpdir/stats.stdout" 2>"$tmpdir/stats.stderr"; then
    fail "unknown summary statistics mode should fail"
fi
pass_check "unknown summary statistics mode should fail"
assert_contains "$tmpdir/stats.stderr" "Do not understand nonsense as mode" \
    "unknown summary mode should identify the invalid value"

# Summary rejects negative edge value.
if "$MSAMTOOLS" summary -e -1 "$summary_fixture" \
    >"$tmpdir/edge.stdout" 2>"$tmpdir/edge.stderr"; then
    fail "negative edge length should fail"
fi
pass_check "negative edge length should fail"
assert_contains "$tmpdir/edge.stdout" "-e must be a positive integer" \
    "negative edge error should explain the allowed values"

# Summary rejects --stats and --count together.
if "$MSAMTOOLS" summary --stats mapped --count "$summary_fixture" \
    >"$tmpdir/stats_count.stdout" 2>"$tmpdir/stats_count.stderr"; then
    fail "--stats and --count together should fail"
fi
pass_check "--stats and --count together should fail"
assert_contains "$tmpdir/stats_count.stdout" \
    "--stats cannot be combined with --count" \
    "summary should explain the stats/count conflict"

# Summary rejects -e and --count together.
if "$MSAMTOOLS" summary -e 10 --count "$summary_fixture" \
    >"$tmpdir/edge_count.stdout" 2>"$tmpdir/edge_count.stderr"; then
    fail "-e and --count together should fail"
fi
pass_check "-e and --count together should fail"
assert_contains "$tmpdir/edge_count.stdout" \
    "-e cannot be combined with --count" \
    "summary should explain the edge/count conflict"

# Coverage rejects negative wordsize value.
if "$MSAMTOOLS" coverage -w -1 -o "$tmpdir/wordsize_negative.gz" "$coverage_fixture" \
    >"$tmpdir/wordsize_negative.stdout" 2>"$tmpdir/wordsize_negative.stderr"; then
    fail "negative word size should fail"
fi
pass_check "negative word size should fail"
assert_contains "$tmpdir/wordsize_negative.stdout" "-w must be a non-zero positive integer" \
    "negative word size should explain the allowed values"

# Coverage rejects zero wordsize value.
if "$MSAMTOOLS" coverage -w 0 -o "$tmpdir/wordsize_zero.gz" "$coverage_fixture" \
    >"$tmpdir/wordsize_zero.stdout" 2>"$tmpdir/wordsize_zero.stderr"; then
    fail "zero word size should fail"
fi
pass_check "zero word size should fail"
assert_contains "$tmpdir/wordsize_zero.stdout" "-w must be a non-zero positive integer" \
    "zero word size should explain the allowed values"

report_checks
exit 0
