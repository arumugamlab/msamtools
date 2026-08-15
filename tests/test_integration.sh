#!/bin/sh

set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

fixture="$script_dir/fixtures/integration.sam"
tmpdir=${TMPDIR:-/tmp}/msamtools-test-integration-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir" || exit 1

run_profile()
{
    name=$1
    input=$2
    output="$tmpdir/$name.tsv.gz"
    stderr="$tmpdir/$name.stderr"

    if ! "$MSAMTOOLS" profile -S \
        --label test \
        --unit ab \
        --nolen \
        --total 4 \
        --multi equal \
        --pandas \
        -o "$output" \
        "$input" \
        >/dev/null 2>"$stderr"; then
        cat "$stderr" >&2
        fail "integration profile case '$name' should succeed"
    fi

    pass_check "integration profile case '$name' should succeed"
}

# Before filtering, filter_to_b and multi each map to both A and B.
run_profile original "$fixture"
assert_profile_contains "$tmpdir/original.tsv.gz" "- Multiple mapped : 2" \
    "unfiltered integration fixture should contain two multimapped inserts"
assert_profile_contains "$tmpdir/original.tsv.gz" "- Uniquely mapped : 2" \
    "unfiltered integration fixture should contain two uniquely mapped inserts"
assert_profile_value "$tmpdir/original.tsv.gz" A 2 1e-9
assert_profile_value "$tmpdir/original.tsv.gz" B 2 1e-9

# Filtering at 95% removes only the poor A alignment of filter_to_b.  That
# insert therefore becomes unique to B, while multi remains ambiguous.
filtered="$tmpdir/filtered.sam"
filter_stderr="$tmpdir/filter.stderr"
if ! "$MSAMTOOLS" filter -S -h -p 95 "$fixture" \
    >"$filtered" 2>"$filter_stderr"; then
    cat "$filter_stderr" >&2
    fail "integration filter should succeed"
fi
pass_check "integration filter should succeed"

assert_sam_records "$filtered" \
    "filter_to_b:256,multi:0,multi:256,uA:0,uB:0" \
    "filter should remove only the low-identity alignment"

# Feed msamtools filter output back into msamtools profile.  The changed
# ambiguity must be reflected exactly in insert counts and equal sharing.
run_profile filtered "$filtered"
assert_profile_contains "$tmpdir/filtered.tsv.gz" "- Multiple mapped : 1" \
    "filtering should reduce the multimapped insert count from two to one"
assert_profile_contains "$tmpdir/filtered.tsv.gz" "- Uniquely mapped : 3" \
    "filtering should increase the uniquely mapped insert count from two to three"
assert_profile_value "$tmpdir/filtered.tsv.gz" A 1.5 1e-9
assert_profile_value "$tmpdir/filtered.tsv.gz" B 2.5 1e-9
assert_profile_contains "$tmpdir/filtered.tsv.gz" \
    "QNAME grouping check: confirmed by input header SO:queryname" \
    "profile should accept QNAME-grouped SAM emitted by filter"

report_checks
exit 0
