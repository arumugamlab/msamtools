#!/bin/sh

set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

fixture="$script_dir/fixtures/profile.sam"
long_qname_fixture="$script_dir/fixtures/long_qname.sam"
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

long_qname_output="$tmpdir/long_qname.tsv.gz"
long_qname_stderr="$tmpdir/long_qname.stderr"
if ! "$MSAMTOOLS" profile -S \
    --label long_qname \
    --unit ab \
    --nolen \
    --total 7 \
    --multi equal \
    --pandas \
    -o "$long_qname_output" \
    "$long_qname_fixture" \
    >/dev/null 2>"$long_qname_stderr"; then
    cat "$long_qname_stderr" >&2
    fail "profile should accept full-length legal QNAMEs"
fi
pass_check "profile should accept full-length legal QNAMEs"
assert_profile_contains "$long_qname_output" "Mapped inserts      : 4" \
    "profile should count four complete long-QNAME groups"
assert_profile_contains "$long_qname_output" "- Multiple mapped : 4" \
    "profile should classify all long-QNAME groups as multimapped"
assert_profile_contains "$long_qname_output" "- Uniquely mapped : 0" \
    "profile should not split long-QNAME multimappers into unique hits"
assert_profile_value "$long_qname_output" A 2 1e-9
assert_profile_value "$long_qname_output" B 2 1e-9

run_zero_profile()
{
    name=$1
    input=$2
    output="$tmpdir/$name.tsv.gz"
    stderr="$tmpdir/$name.stderr"

    if ! "$MSAMTOOLS" profile -S \
        --label "$name" \
        --unit ab \
        --nolen \
        --multi equal \
        --pandas \
        -o "$output" \
        "$input" \
        >/dev/null 2>"$stderr"; then
        cat "$stderr" >&2
        fail "profile should handle $name input"
    fi
    pass_check "profile should handle $name input"

    assert_profile_contains "$output" "Mapped inserts      :       0" \
        "$name profile should report zero mapped inserts"
    assert_profile_contains "$output" "- Multiple mapped :       0" \
        "$name profile should report zero multimapped inserts"
    assert_profile_contains "$output" "- Uniquely mapped :       0" \
        "$name profile should report zero uniquely mapped inserts"
    assert_profile_contains "$output" "Effective inserts   :       0" \
        "$name profile should report zero effective inserts"

    assert_profile_value "$output" Unknown 0 1e-9
    assert_profile_value "$output" A 0 1e-9
    assert_profile_value "$output" B 0 1e-9
}

run_zero_profile empty "$script_dir/fixtures/profile_empty.sam"
run_zero_profile unmapped "$script_dir/fixtures/profile_unmapped.sam"

report_checks
exit 0
