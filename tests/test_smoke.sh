#!/bin/sh

set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

tmpdir=${TMPDIR:-/tmp}/msamtools-test-smoke-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir" || exit 1

stdout="$tmpdir/stdout"
stderr="$tmpdir/stderr"

if ! "$MSAMTOOLS" help >"$stdout" 2>"$stderr"; then
    fail "msamtools help should succeed"
fi
pass_check "msamtools help should succeed"

assert_contains "$stdout" "Program: msamtools" \
    "help output should identify msamtools"
assert_contains "$stdout" "profile" \
    "help output should list the profile command"
assert_contains "$stdout" "filter" \
    "help output should list the filter command"
assert_empty "$stderr" \
    "msamtools help should not write to stderr"

if "$MSAMTOOLS" >"$stdout" 2>"$stderr"; then
    fail "msamtools without a command should fail"
fi
pass_check "msamtools without a command should fail"

assert_empty "$stdout" \
    "missing-command error should not write to stdout"
assert_contains "$stderr" "Usage:" \
    "missing-command error should print usage information"

report_checks
exit 0
