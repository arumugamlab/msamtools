#!/bin/sh

: "${MSAMTOOLS:?MSAMTOOLS must point to the msamtools executable}"

TEST_COUNT=0

pass_check()
{
    TEST_COUNT=$((TEST_COUNT + 1))
    echo "PASS [$TEST_COUNT]: $*"
}

report_checks()
{
    echo "Completed $TEST_COUNT checks"
}

fail()
{
    echo "FAIL: $*" >&2
    exit 1
}

assert_contains()
{
    file=$1
    expected=$2
    description=$3

    if ! grep -F -- "$expected" "$file" >/dev/null 2>&1; then
        echo "Expected to find:" >&2
        echo "  $expected" >&2
        echo "in $file, whose contents were:" >&2
        cat "$file" >&2
        fail "$description"
    fi

    pass_check "$description"
}

assert_empty()
{
    file=$1
    description=$2

    if test -s "$file"; then
        echo "Expected empty file: $file" >&2
        cat "$file" >&2
        fail "$description"
    fi

    pass_check "$description"
}

assert_failure()
{
    description=$1
    shift

    if "$@"; then
        fail "$description (command unexpectedly succeeded)"
    fi

    pass_check "$description"
}

assert_file_absent()
{
    file=$1
    description=$2

    if test -e "$file"; then
        fail "$description (unexpected file: $file)"
    fi

    pass_check "$description"
}

profile_value()
{
    profile=$1
    feature=$2

    gzip -cd -- "$profile" |
        awk -F '\t' -v feature="$feature" '$1 == feature { print $2; exit }'
}

assert_profile_value()
{
    profile=$1
    feature=$2
    expected=$3
    tolerance=$4

    observed=$(profile_value "$profile" "$feature")
    if test -z "$observed"; then
        fail "profile is missing feature '$feature'"
    fi

    if ! awk -v e="$expected" -v o="$observed" -v t="$tolerance" 'BEGIN {
        d = e - o;
        if (d < 0) d = -d;
        exit !(d <= t);
    }'; then
        fail "feature '$feature' abundance differs from expected value (expected $expected, observed $observed, tolerance $tolerance)"
    fi

    pass_check "feature '$feature' abundance matches expected value"
}

assert_profile_contains()
{
    profile=$1
    expected=$2
    description=$3

    if ! gzip -cd -- "$profile" | grep -F -- "$expected" >/dev/null 2>&1; then
        echo "Expected to find in decompressed profile:" >&2
        echo "  $expected" >&2
        gzip -cd -- "$profile" >&2
        fail "$description"
    fi

    pass_check "$description"
}
