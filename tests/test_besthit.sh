#!/bin/sh

set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

fixture="$script_dir/fixtures/besthit.sam"
rescore_fixture="$script_dir/fixtures/besthit_rescore.sam"
long_qname_fixture="$script_dir/fixtures/long_qname.sam"

tmpdir=${TMPDIR:-/tmp}/msamtools-test-besthit-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir" || exit 1

run_filter()
{
    name=$1
    shift
    output="$tmpdir/$name.sam"
    stderr="$tmpdir/$name.stderr"

    if ! "$MSAMTOOLS" filter -S -h "$@" >"$output" 2>"$stderr"; then
        cat "$stderr" >&2
        fail "$name should succeed"
    fi
    pass_check "$name should succeed"
}

# --besthit keeps every alignment tied for the highest AS score in each
# read/mate pool.
run_filter besthit --besthit "$fixture"
assert_sam_records "$tmpdir/besthit.sam" \
    "filterwin:0,paired:321,paired:129,same_ref:256,single:0,tie2:0,tie2:256,tie3:0,tie3:256,unique2:256,unique3:256" \
    "--besthit should retain all and only highest-scoring hits"
assert_contains "$tmpdir/besthit.sam" \
    "QNAME grouping check: confirmed by input header SO:queryname" \
    "best-hit SAM header should record QNAME-order confirmation"

# --uniqhit retains the best alignment only when the highest score is unique.
run_filter uniqhit --uniqhit "$fixture"
assert_sam_records "$tmpdir/uniqhit.sam" \
    "filterwin:0,paired:321,paired:129,same_ref:256,single:0,unique2:256,unique3:256" \
    "--uniqhit should discard groups tied for the best score"

# Per-alignment filtering happens before best-hit selection.  filterwin's
# higher-AS A alignment is only 90% identical and is removed by -p 95.
run_filter filtered_besthit -p 95 --besthit "$fixture"
assert_sam_records "$tmpdir/filtered_besthit.sam" \
    "filterwin:256,paired:321,paired:129,same_ref:256,single:0,tie2:0,tie2:256,tie3:0,tie3:256,unique2:256,unique3:256" \
    "alignment filtering should precede best-hit selection"

# --rescore recomputes AS before best-hit selection when the normal filtering
# path is active.  Existing AS favors A; edit-distance rescoring favors B.
run_filter rescore_original --besthit "$rescore_fixture"
assert_sam_records "$tmpdir/rescore_original.sam" \
    "rescore:0" \
    "best-hit selection should use existing AS without --rescore"

run_filter rescore_recomputed -l 1 --rescore --besthit "$rescore_fixture"
assert_sam_records "$tmpdir/rescore_recomputed.sam" \
    "rescore:256" \
    "--rescore should change the winning best hit when recalculated AS differs"
assert_sam_record_contains "$tmpdir/rescore_recomputed.sam" \
    "rescore" "256" "AS:i:100" \
    "rescored winning alignment should contain the recomputed AS value"

run_filter long_qname_besthit --besthit "$long_qname_fixture"
awk -F '\t' '
    $1 !~ /^@/ {
        as = "";
        for (i = 12; i <= NF; i++) {
            if ($i ~ /^AS:i:/) {
                as = $i;
                break;
            }
        }
        print length($1) "\t" $2 "\t" $3 "\t" as;
    }
' "$tmpdir/long_qname_besthit.sam" >"$tmpdir/long_qname_besthit.records"
cat >"$tmpdir/long_qname_besthit.expected" <<'EOF'
127	0	A	AS:i:100
128	0	A	AS:i:100
254	0	A	AS:i:100
254	0	B	AS:i:100
EOF
if ! cmp -s "$tmpdir/long_qname_besthit.expected" \
    "$tmpdir/long_qname_besthit.records"; then
    echo "Expected long-QNAME best-hit records:" >&2
    cat "$tmpdir/long_qname_besthit.expected" >&2
    echo "Observed long-QNAME best-hit records:" >&2
    cat "$tmpdir/long_qname_besthit.records" >&2
    fail "--besthit should keep one winner per complete QNAME"
fi
pass_check "--besthit should keep one winner per complete QNAME"

run_filter long_qname_filtered_besthit -p 90 --besthit "$long_qname_fixture"
awk -F '\t' '
    $1 !~ /^@/ {
        as = "";
        for (i = 12; i <= NF; i++) {
            if ($i ~ /^AS:i:/) {
                as = $i;
                break;
            }
        }
        print length($1) "\t" $2 "\t" $3 "\t" as;
    }
' "$tmpdir/long_qname_filtered_besthit.sam" \
    >"$tmpdir/long_qname_filtered_besthit.records"
if ! cmp -s "$tmpdir/long_qname_besthit.expected" \
    "$tmpdir/long_qname_filtered_besthit.records"; then
    echo "Expected filtered long-QNAME best-hit records:" >&2
    cat "$tmpdir/long_qname_besthit.expected" >&2
    echo "Observed filtered long-QNAME best-hit records:" >&2
    cat "$tmpdir/long_qname_filtered_besthit.records" >&2
    fail "filtered --besthit should keep one winner per complete QNAME"
fi
pass_check "filtered --besthit should keep one winner per complete QNAME"

report_checks
exit 0
