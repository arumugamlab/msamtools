#!/bin/sh

set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

fixture="$script_dir/fixtures/profile.sam"
tmpdir=${TMPDIR:-/tmp}/msamtools-test-profile-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir" || exit 1

run_profile()
{
    mode=$1
    output="$tmpdir/$mode.tsv.gz"
    stderr="$tmpdir/$mode.stderr"

    if ! "$MSAMTOOLS" profile -S \
        --label test \
        --unit ab \
        --nolen \
        --total 7 \
        --multi "$mode" \
        --pandas \
        -o "$output" \
        "$fixture" \
        >/dev/null 2>"$stderr"; then
        cat "$stderr" >&2
        fail "profile --multi $mode should succeed"
    fi
    pass_check "profile --multi $mode should succeed"

    assert_profile_contains "$output" \
        "QNAME grouping check: confirmed by input header SO:queryname" \
        "profile should record trusted queryname ordering"
    assert_profile_contains "$output" "Total inserts       : 7" \
        "profile should report seven total inserts"
    assert_profile_contains "$output" "Mapped inserts      : 7" \
        "profile should report seven mapped inserts"
    assert_profile_contains "$output" "- Multiple mapped : 1" \
        "profile should report one multimapped insert"
    assert_profile_contains "$output" "- Uniquely mapped : 6" \
        "profile should report six uniquely mapped inserts"
}

run_profile all
assert_profile_value "$tmpdir/all.tsv.gz" Unknown 0 1e-9
assert_profile_value "$tmpdir/all.tsv.gz" A 6 1e-9
assert_profile_value "$tmpdir/all.tsv.gz" B 2 1e-9

run_profile equal
assert_profile_value "$tmpdir/equal.tsv.gz" Unknown 0 1e-9
assert_profile_value "$tmpdir/equal.tsv.gz" A 5.5 1e-9
assert_profile_value "$tmpdir/equal.tsv.gz" B 1.5 1e-9

run_profile ignore
assert_profile_value "$tmpdir/ignore.tsv.gz" Unknown 1 1e-9
assert_profile_value "$tmpdir/ignore.tsv.gz" A 5 1e-9
assert_profile_value "$tmpdir/ignore.tsv.gz" B 1 1e-9

run_profile proportional
assert_profile_value "$tmpdir/proportional.tsv.gz" Unknown 0 1e-9
assert_profile_value "$tmpdir/proportional.tsv.gz" A 5.833333333333 1e-6
assert_profile_value "$tmpdir/proportional.tsv.gz" B 1.166666666667 1e-6

report_checks
exit 0
