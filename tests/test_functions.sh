#!/bin/sh

: "${MSAMTOOLS:?MSAMTOOLS must point to the msamtools executable}"

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
}
