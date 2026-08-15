#!/bin/sh

set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

tmpdir=${TMPDIR:-/tmp}/msamtools-test-streaming-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir" || exit 1

# Generate 100,001 records without committing a large fixture.  The first
# boundary alignment is record 100,000 and is therefore the last record read
# into the QNAME-preflight replay buffer.  Its second alignment is record
# 100,001 and must be read live from the underlying SAM stream.
generate_boundary_sam()
{
    awk 'BEGIN {
        seq = "";
        qual = "";
        for (j = 0; j < 100; j++) {
            seq = seq "A";
            qual = qual "I";
        }

        print "@HD\tVN:1.6";
        print "@SQ\tSN:A\tLN:1000";
        print "@SQ\tSN:B\tLN:1000";

        # Alternating positions make the prefix definitively non-coordinate
        # ordered, so this test exercises replay rather than the coordinate
        # ordering warning heuristic.
        for (i = 1; i <= 99999; i++) {
            pos = (i % 2 == 0) ? 100 : 200;
            printf "q%06d\t0\tA\t%d\t60\t100M\t*\t0\t0\t%s\t%s\tAS:i:50\tNM:i:0\n",
                   i, pos, seq, qual;
        }

        # Records 100000 and 100001 form one QNAME group across the
        # replay-buffer/live-input transition.  Only the B hit is best.
        printf "boundary\t0\tA\t300\t60\t100M\t*\t0\t0\t%s\t%s\tAS:i:10\tNM:i:0\n",
               seq, qual;
        printf "boundary\t256\tB\t400\t60\t100M\t*\t0\t0\t%s\t%s\tAS:i:20\tNM:i:0\n",
               seq, qual;
    }'
}

output="$tmpdir/boundary.sam"
stderr="$tmpdir/boundary.stderr"

if ! generate_boundary_sam |
    "$MSAMTOOLS" filter -S -h --besthit - \
        >"$output" 2>"$stderr"; then
    cat "$stderr" >&2
    fail "best-hit filtering across the replay boundary should succeed"
fi
pass_check "best-hit filtering across the replay boundary should succeed"

assert_contains "$output" \
    "QNAME grouping check: no QNAME grouping violation detected in first 10000 records" \
    "streaming test should exercise sampled QNAME preflight"

record_count=$(awk -F '\t' '$1 !~ /^@/ { n++ } END { print n + 0 }' "$output")
if test "$record_count" -ne 100000; then
    fail "replay-boundary best-hit output should contain 100000 records (observed $record_count)"
fi
pass_check "replay-boundary best-hit output should contain 100000 records"

boundary_records=$(awk -F '\t' '
    $1 == "boundary" {
        if (n > 0) printf ",";
        printf "%s:%s:%s:%s", $1, $2, $3, $4;
        n++;
    }
    END { print "" }
' "$output")

if test "$boundary_records" != "boundary:256:B:400"; then
    echo "Observed boundary records: $boundary_records" >&2
    fail "QNAME group spanning replay and live input should have one best hit"
fi
pass_check "QNAME group spanning replay and live input should have one best hit"

report_checks
exit 0
