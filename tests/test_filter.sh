#!/bin/sh

set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

fixture="$script_dir/fixtures/filter.sam"
tmpdir=${TMPDIR:-/tmp}/msamtools-test-filter-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir" || exit 1

ALL_MAPPED="clip80:0,clip90:0,del99:0,id97:0,id98:0,ins99:0,len80:0,md99:0,md_precedence:0,pairX:65,pairX:129,perfect100:0,secondary98:256,supp98:2048"

run_filter()
{
    name=$1
    shift

    output="$tmpdir/$name.sam"
    stderr="$tmpdir/$name.stderr"

    if ! "$MSAMTOOLS" filter -S "$@" "$fixture" \
        >"$output" 2>"$stderr"; then
        cat "$stderr" >&2
        fail "filter case '$name' should succeed"
    fi

    pass_check "filter case '$name' should succeed"
}

# Alignment-length boundaries are inclusive.
run_filter length_below -l 79
assert_sam_records "$tmpdir/length_below.sam" "$ALL_MAPPED" \
    "-l 79 should retain every mapped alignment"

run_filter length_equal -l 80
assert_sam_records "$tmpdir/length_equal.sam" "$ALL_MAPPED" \
    "-l 80 should retain alignments exactly 80 bases long"

run_filter length_above -l 81
assert_sam_records "$tmpdir/length_above.sam" \
    "clip90:0,del99:0,id97:0,id98:0,ins99:0,md99:0,md_precedence:0,pairX:65,pairX:129,perfect100:0,secondary98:256,supp98:2048" \
    "-l 81 should remove 80-base alignments"

run_filter length_100 -l 100
assert_sam_records "$tmpdir/length_100.sam" \
    "del99:0,id97:0,id98:0,ins99:0,md99:0,md_precedence:0,pairX:65,pairX:129,perfect100:0,secondary98:256,supp98:2048" \
    "-l should count insertions and deletions in alignment length"

# Percent-identity thresholds are inclusive.
run_filter pid_97 -p 97
assert_sam_records "$tmpdir/pid_97.sam" "$ALL_MAPPED" \
    "-p 97 should retain alignments at exactly 97 percent identity"

run_filter pid_98 -p 98
assert_sam_records "$tmpdir/pid_98.sam" \
    "clip80:0,clip90:0,del99:0,id98:0,ins99:0,len80:0,md99:0,md_precedence:0,pairX:65,perfect100:0,secondary98:256,supp98:2048" \
    "-p 98 should remove alignments below 98 percent identity"

run_filter pid_99 -p 99
assert_sam_records "$tmpdir/pid_99.sam" \
    "clip80:0,clip90:0,del99:0,ins99:0,len80:0,md99:0,md_precedence:0,pairX:65,perfect100:0" \
    "-p 99 should retain alignments at exactly 99 percent identity"

run_filter pid_100 -p 100
assert_sam_records "$tmpdir/pid_100.sam" \
    "clip80:0,clip90:0,len80:0,md_precedence:0,pairX:65,perfect100:0" \
    "-p 100 should retain only perfect alignments"

# Parts-per-thousand mode provides the same lower-bound semantics for
# positive thresholds and an inclusive upper bound for negative thresholds.
run_filter ppt_positive --ppt 980
assert_sam_records "$tmpdir/ppt_positive.sam" \
    "clip80:0,clip90:0,del99:0,id98:0,ins99:0,len80:0,md99:0,md_precedence:0,pairX:65,perfect100:0,secondary98:256,supp98:2048" \
    "--ppt 980 should match -p 98"

run_filter ppt_negative --ppt -980
assert_sam_records "$tmpdir/ppt_negative.sam" \
    "id97:0,id98:0,pairX:129,secondary98:256,supp98:2048" \
    "--ppt -980 should retain identity at or below 98 percent"

# Query-aligned fraction (-z) is inclusive at the boundary.
run_filter qfrac_80 -z 80
assert_sam_records "$tmpdir/qfrac_80.sam" "$ALL_MAPPED" \
    "-z 80 should retain an alignment covering exactly 80 percent of the query"

run_filter qfrac_90 -z 90
assert_sam_records "$tmpdir/qfrac_90.sam" \
    "clip90:0,del99:0,id97:0,id98:0,ins99:0,len80:0,md99:0,md_precedence:0,pairX:65,pairX:129,perfect100:0,secondary98:256,supp98:2048" \
    "-z 90 should retain an alignment covering exactly 90 percent of the query"

run_filter qfrac_91 -z 91
assert_sam_records "$tmpdir/qfrac_91.sam" \
    "del99:0,id97:0,id98:0,ins99:0,len80:0,md99:0,md_precedence:0,pairX:65,pairX:129,perfect100:0,secondary98:256,supp98:2048" \
    "-z 91 should remove 90-percent and 80-percent query coverage"

# Exercise every combined filter dispatch path: -lp, -lz, -pz and -lpz.
run_filter length_pid -l 100 -p 99
assert_sam_records "$tmpdir/length_pid.sam" \
    "del99:0,ins99:0,md99:0,md_precedence:0,pairX:65,perfect100:0" \
    "combined -l/-p filtering should require both thresholds"

run_filter length_qfrac -l 90 -z 91
assert_sam_records "$tmpdir/length_qfrac.sam" \
    "del99:0,id97:0,id98:0,ins99:0,md99:0,md_precedence:0,pairX:65,pairX:129,perfect100:0,secondary98:256,supp98:2048" \
    "combined -l/-z filtering should require both thresholds"

run_filter pid_qfrac -p 99 -z 91
assert_sam_records "$tmpdir/pid_qfrac.sam" \
    "del99:0,ins99:0,len80:0,md99:0,md_precedence:0,pairX:65,perfect100:0" \
    "combined -p/-z filtering should require both thresholds"

run_filter all_three -l 100 -p 99 -z 91
assert_sam_records "$tmpdir/all_three.sam" \
    "del99:0,ins99:0,md99:0,md_precedence:0,pairX:65,perfect100:0" \
    "combined -l/-p/-z filtering should require all thresholds"

# MD takes precedence over NM when both tags are present.  md_precedence has
# NM:i:10 but MD:Z:100 and therefore behaves as a perfect alignment.
assert_sam_record_contains "$tmpdir/pid_100.sam" md_precedence 0 "NM:i:10" \
    "MD should take precedence over NM for identity filtering"

# Secondary and supplementary records are filtered using the same identity
# semantics as primary alignments.
assert_sam_record_contains "$tmpdir/pid_98.sam" secondary98 256 "NM:i:2" \
    "secondary alignments should participate in ordinary filtering"
assert_sam_record_contains "$tmpdir/pid_98.sam" supp98 2048 "NM:i:2" \
    "supplementary alignments should participate in ordinary filtering"

# Paired mates are filtered independently in ordinary filtering.
assert_sam_record_contains "$tmpdir/pid_98.sam" pairX 65 "NM:i:0" \
    "passing mate should be retained"
if awk -F '\t' '$1 == "pairX" && $2 == 129 { found = 1 } END { exit !found }' \
    "$tmpdir/pid_98.sam"; then
    fail "failing mate should be removed independently"
fi
pass_check "failing mate should be removed independently"

# Inversion returns the complement among mapped records.  Unmapped records are
# retained only when explicitly requested with --keep_unmapped.
run_filter inverted -p 99 -v
assert_sam_records "$tmpdir/inverted.sam" \
    "id97:0,id98:0,pairX:129,secondary98:256,supp98:2048" \
    "-v should return the mapped complement of the filter"

run_filter inverted_keep_unmapped -p 99 -v -k
assert_sam_records "$tmpdir/inverted_keep_unmapped.sam" \
    "id97:0,id98:0,pairX:129,secondary98:256,supp98:2048,unmapped:4" \
    "--keep_unmapped should retain unmapped reads with an inverted lower-bound filter"

# --rescore replaces/creates AS using the same alignment summary used for
# filtering: matches score +1 and edits score -1.
run_filter rescored -p 98 --rescore
assert_sam_records "$tmpdir/rescored.sam" \
    "clip80:0,clip90:0,del99:0,id98:0,ins99:0,len80:0,md99:0,md_precedence:0,pairX:65,perfect100:0,secondary98:256,supp98:2048" \
    "--rescore should not change which alignments pass -p 98"
assert_sam_record_contains "$tmpdir/rescored.sam" id98 0 "AS:i:96" \
    "--rescore should assign the expected score to a 98-percent identity alignment"
assert_sam_record_contains "$tmpdir/rescored.sam" md_precedence 0 "AS:i:100" \
    "--rescore should use MD-precedence edit counts"

report_checks
exit 0
