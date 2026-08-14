#!/bin/sh

set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

tmpdir=${TMPDIR:-/tmp}/msamtools-test-qname-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir" || exit 1

# Explicit SO:queryname is authoritative and should be recorded as confirmed.
queryname_output="$tmpdir/queryname.tsv.gz"
queryname_stderr="$tmpdir/queryname.stderr"
if ! "$MSAMTOOLS" profile -S \
    --label test --unit ab --nolen --total 7 --multi equal --pandas \
    -o "$queryname_output" "$script_dir/fixtures/profile.sam" \
    >/dev/null 2>"$queryname_stderr"; then
    cat "$queryname_stderr" >&2
    fail "SO:queryname input should be accepted"
fi
pass_check "SO:queryname input should be accepted"
assert_profile_contains "$queryname_output" \
    "QNAME grouping check: confirmed by input header SO:queryname" \
    "profile should record SO:queryname confirmation"

# Explicit coordinate sorting is a definite incompatibility and must fail
# before a profile output is created.
coordinate_output="$tmpdir/coordinate.tsv.gz"
coordinate_stdout="$tmpdir/coordinate.stdout"
coordinate_stderr="$tmpdir/coordinate.stderr"
if "$MSAMTOOLS" profile -S \
    --label test --unit ab --nolen --multi equal --pandas \
    -o "$coordinate_output" "$script_dir/fixtures/qname_coordinate.sam" \
    >"$coordinate_stdout" 2>"$coordinate_stderr"; then
    fail "SO:coordinate input should be rejected"
fi
pass_check "SO:coordinate input should be rejected"
assert_contains "$coordinate_stderr" "SO:coordinate" \
    "coordinate-sort error should identify the header declaration"
assert_contains "$coordinate_stderr" "samtools sort -n" \
    "coordinate-sort error should recommend name sorting"
assert_file_absent "$coordinate_output" \
    "coordinate-sort failure should not create a profile"

# A QNAME that reappears after another QNAME is a directly observed grouping
# violation, independent of header metadata.
reopened_output="$tmpdir/reopened.tsv.gz"
reopened_stdout="$tmpdir/reopened.stdout"
reopened_stderr="$tmpdir/reopened.stderr"
if "$MSAMTOOLS" profile -S \
    --label test --unit ab --nolen --multi equal --pandas \
    -o "$reopened_output" "$script_dir/fixtures/qname_reopened.sam" \
    >"$reopened_stdout" 2>"$reopened_stderr"; then
    fail "reopened QNAME should be rejected"
fi
pass_check "reopened QNAME should be rejected"
assert_contains "$reopened_stderr" "not grouped by QNAME" \
    "reopened-QNAME error should explain the violation"
assert_contains "$reopened_stderr" "readA" \
    "reopened-QNAME error should identify an example read"
assert_file_absent "$reopened_output" \
    "reopened-QNAME failure should not create a profile"

report_checks
exit 0
