#!/bin/sh

set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/test_functions.sh"

fixture="$script_dir/fixtures/profile.sam"
fractional_mincount_fixture="$script_dir/fixtures/profile_fractional_mincount.sam"
long_qname_fixture="$script_dir/fixtures/long_qname.sam"
tmpdir=${TMPDIR:-/tmp}/msamtools-test-profile-$$
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
mkdir -p "$tmpdir" || exit 1

assert_profile_values_nonnegative()
{
    file=$1
    message=$2

    if gzip -cd -- "$file" |
        awk '
            !/^#/ && NR > 1 {
                value = $NF
                if (value ~ /^[+-]?(nan|inf)$/ ||
                    value + 0 < 0) {
                    exit 1
                }
            }
        '
    then
        pass_check "$message"
    else
        fail "$message"
    fi
}

assert_profile_stats_nonnegative()
{
    file=$1
    message=$2

    if gzip -cd -- "$file" |
        awk '
            /^# (Total inserts|Mapped inserts|  - Multiple mapped|  - Uniquely mapped|Purged insert-equivalents|Effective insert-equivalents)/ {
                line = $0
                sub(/^[^:]*:[[:space:]]*/, "", line)
                split(line, x, /[[:space:]]+/)

                if (x[1] != "NA" && x[1] + 0 < 0)
                    exit 1
            }
        '
    then
        pass_check "$message"
    else
        fail "$message"
    fi
}

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
    assert_profile_values_nonnegative "$output" \
        "profile --multi $mode should contain only finite non-negative values"
    assert_profile_stats_nonnegative "$output" \
        "profile --multi $mode should report only non-negative statistics"
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

all_mincount_output="$tmpdir/all_mincount.tsv.gz"
all_mincount_stderr="$tmpdir/all_mincount.stderr"

if ! "$MSAMTOOLS" profile -S \
    --label all_mincount \
    --unit ab \
    --nolen \
    --total 7 \
    --multi all \
    --mincount 10 \
    --pandas \
    -o "$all_mincount_output" \
    "$fixture" \
    >/dev/null 2>"$all_mincount_stderr"; then
    cat "$all_mincount_stderr" >&2
    fail "profile --multi all with --mincount should succeed"
fi
pass_check "profile --multi all with --mincount should succeed"

assert_profile_contains "$all_mincount_output" \
    "Mapped inserts      : 7" \
    "all/mincount profile should report seven mapped inserts"

assert_profile_contains "$all_mincount_output" \
    "Purged insert-equivalents    :       8.00" \
    "all/mincount profile should report eight purged insert-equivalents"

assert_profile_contains "$all_mincount_output" \
    "Effective insert-equivalents :       0.00" \
    "all/mincount profile should report zero effective insert-equivalents"

assert_profile_value "$all_mincount_output" Unknown 8 1e-9
assert_profile_value "$all_mincount_output" A 0 1e-9
assert_profile_value "$all_mincount_output" B 0 1e-9

# Genome-definition tests
genome_def="$tmpdir/genome.tsv"
cat >"$genome_def" <<'EOF'
Genome1	A
Genome2	B
EOF

genome_output="$tmpdir/genome.tsv.gz"
genome_stderr="$tmpdir/genome.stderr"

if ! "$MSAMTOOLS" profile -S \
    --label genome \
    --unit ab \
    --nolen \
    --total 7 \
    --multi equal \
    --genome "$genome_def" \
    --pandas \
    -o "$genome_output" \
    "$fixture" \
    >/dev/null 2>"$genome_stderr"; then
    cat "$genome_stderr" >&2
    fail "profile with valid genome definition should succeed"
fi
pass_check "profile with valid genome definition should succeed"

assert_profile_value "$genome_output" Genome1 5.5 1e-9
assert_profile_value "$genome_output" Genome2 1.5 1e-9

missing_genome_def="$tmpdir/genome_missing.tsv"
cat >"$missing_genome_def" <<'EOF'
Genome1	A
EOF

missing_stderr="$tmpdir/genome_missing.stderr"
if "$MSAMTOOLS" profile -S \
    --label genome_missing \
    --unit ab \
    --nolen \
    --multi equal \
    --genome "$missing_genome_def" \
    -o "$tmpdir/genome_missing.tsv.gz" \
    "$fixture" \
    >/dev/null 2>"$missing_stderr"; then
    fail "profile should reject genome definitions missing a BAM target"
fi
pass_check "profile should reject genome definitions missing a BAM target"

if ! grep -F "B" "$missing_stderr" >/dev/null; then
    cat "$missing_stderr" >&2
    fail "missing-target genome-definition error should mention B"
fi
pass_check "missing-target genome-definition error should mention B"

duplicate_genome_def="$tmpdir/genome_duplicate.tsv"
cat >"$duplicate_genome_def" <<'EOF'
Genome1	A
Genome2	A
Genome2	B
EOF

duplicate_stderr="$tmpdir/genome_duplicate.stderr"
if "$MSAMTOOLS" profile -S \
    --label genome_duplicate \
    --unit ab \
    --nolen \
    --multi equal \
    --genome "$duplicate_genome_def" \
    -o "$tmpdir/genome_duplicate.tsv.gz" \
    "$fixture" \
    >/dev/null 2>"$duplicate_stderr"; then
    fail "profile should reject duplicate sequence names in genome definitions"
fi
pass_check "profile should reject duplicate sequence names in genome definitions"

if ! grep -F "A" "$duplicate_stderr" >/dev/null; then
    cat "$duplicate_stderr" >&2
    fail "duplicate-sequence genome-definition error should mention A"
fi
pass_check "duplicate-sequence genome-definition error should mention A"

unknown_genome_def="$tmpdir/genome_unknown.tsv"
cat >"$unknown_genome_def" <<'EOF'
Genome1	A
Genome2	B
Genome3	C
EOF

unknown_stderr="$tmpdir/genome_unknown.stderr"
if "$MSAMTOOLS" profile -S \
    --label genome_unknown \
    --unit ab \
    --nolen \
    --multi equal \
    --genome "$unknown_genome_def" \
    -o "$tmpdir/genome_unknown.tsv.gz" \
    "$fixture" \
    >/dev/null 2>"$unknown_stderr"; then
    fail "profile should reject genome definitions containing unknown BAM targets"
fi
pass_check "profile should reject genome definitions containing unknown BAM targets"

if ! grep -F "C" "$unknown_stderr" >/dev/null; then
    cat "$unknown_stderr" >&2
    fail "unknown-target genome-definition error should mention C"
fi
pass_check "unknown-target genome-definition error should mention C"

# Genome-definition hash-expansion regression test.
#
# Create 3000 reference features grouped into 300 genomes:
#   f1-f10       -> g1
#   f11-f20      -> g2
#   ...
#   f2991-f3000  -> g300
#
# Only a small number of references receive alignments. Distinct nonzero
# counts are distributed across the genome-ID range so that any mismatch
# between stored genome IDs and hash key order becomes visible.
genome_hash_sam="$tmpdir/genome_hash.sam"
genome_hash_def="$tmpdir/genome_hash.tsv"
genome_hash_output="$tmpdir/genome_hash.profile.tsv.gz"
genome_hash_stderr="$tmpdir/genome_hash.stderr"

{
    printf '@HD\tVN:1.6\tSO:queryname\n'

    i=1
    while test "$i" -le 3000; do
        printf '@SQ\tSN:f%d\tLN:1000\n' "$i"
        i=$((i + 1))
    done

    q=1

    # Helper pattern used below:
    # genome gN is represented by feature f((N-1)*10 + 1).

    # Give selected genomes distinct unique-insert counts.
    # g1 -> 1
    i=1
    while test "$i" -le 1; do
        printf 'q%d\t0\tf1\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
        q=$((q + 1))
        i=$((i + 1))
    done

    # g17 -> 2
    i=1
    while test "$i" -le 2; do
        printf 'q%d\t0\tf161\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
        q=$((q + 1))
        i=$((i + 1))
    done

    # g33 -> 3
    i=1
    while test "$i" -le 3; do
        printf 'q%d\t0\tf321\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
        q=$((q + 1))
        i=$((i + 1))
    done

    # g48 -> 4
    i=1
    while test "$i" -le 4; do
        printf 'q%d\t0\tf471\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
        q=$((q + 1))
        i=$((i + 1))
    done

    # g77 -> 5
    i=1
    while test "$i" -le 5; do
        printf 'q%d\t0\tf761\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
        q=$((q + 1))
        i=$((i + 1))
    done

    # g108 -> 6
    i=1
    while test "$i" -le 6; do
        printf 'q%d\t0\tf1071\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
        q=$((q + 1))
        i=$((i + 1))
    done

    # g134 -> 7
    i=1
    while test "$i" -le 7; do
        printf 'q%d\t0\tf1331\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
        q=$((q + 1))
        i=$((i + 1))
    done

    # g151 -> 8
    i=1
    while test "$i" -le 8; do
        printf 'q%d\t0\tf1501\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
        q=$((q + 1))
        i=$((i + 1))
    done

    # g201 -> 9
    i=1
    while test "$i" -le 9; do
        printf 'q%d\t0\tf2001\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
        q=$((q + 1))
        i=$((i + 1))
    done

    # g250 -> 10
    i=1
    while test "$i" -le 10; do
        printf 'q%d\t0\tf2491\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
        q=$((q + 1))
        i=$((i + 1))
    done

    # g300 -> 11
    i=1
    while test "$i" -le 11; do
        printf 'q%d\t0\tf2991\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
        q=$((q + 1))
        i=$((i + 1))
    done

    # One insert maps to two different features in the same genome.
    # Both f1001 and f1002 belong to g101, so after feature-to-genome
    # collapsing this must contribute exactly one insert to g101.
    printf 'q%d\t0\tf1001\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"
    printf 'q%d\t0\tf1002\t1\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$q"

} >"$genome_hash_sam"

i=1
while test "$i" -le 3000; do
    genome=$(( (i - 1) / 10 + 1 ))
    printf 'g%d\tf%d\n' "$genome" "$i"
    i=$((i + 1))
done >"$genome_hash_def"

if ! "$MSAMTOOLS" profile -S \
    --label genome_hash \
    --multi all \
    --unit ab \
    --nolen \
    --genome "$genome_hash_def" \
    --pandas \
    -o "$genome_hash_output" \
    "$genome_hash_sam" \
    >/dev/null 2>"$genome_hash_stderr"; then
    cat "$genome_hash_stderr" >&2
    fail "profile genome mapping should remain correct after hash expansion"
fi
pass_check "profile genome mapping should remain correct after hash expansion"

# Distinct counts distributed across the genome-ID range. These should all
# fail if genome names become detached from their stored genome IDs.
assert_profile_value "$genome_hash_output" g1   1  1e-9
assert_profile_value "$genome_hash_output" g17  2  1e-9
assert_profile_value "$genome_hash_output" g33  3  1e-9
assert_profile_value "$genome_hash_output" g48  4  1e-9
assert_profile_value "$genome_hash_output" g77  5  1e-9
assert_profile_value "$genome_hash_output" g108 6  1e-9
assert_profile_value "$genome_hash_output" g134 7  1e-9
assert_profile_value "$genome_hash_output" g151 8  1e-9
assert_profile_value "$genome_hash_output" g201 9  1e-9
assert_profile_value "$genome_hash_output" g250 10 1e-9
assert_profile_value "$genome_hash_output" g300 11 1e-9

# Two alignments from the same insert map to different features of the same
# genome and must collapse to one genome-level hit.
assert_profile_value "$genome_hash_output" g101 1 1e-9

# Long QNAME tests

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
    unit=$3
    output="$tmpdir/$name.tsv.gz"
    stderr="$tmpdir/$name.stderr"

    if ! "$MSAMTOOLS" profile -S \
        --label "$name" \
        --unit $unit \
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

    assert_profile_values_nonnegative "$output" \
        "profile --multi $mode should contain only finite non-negative values"
    assert_profile_stats_nonnegative "$output" \
        "profile --multi $mode should report only non-negative statistics"

    assert_profile_contains "$output" "Mapped inserts      :       0" \
        "$name profile should report zero mapped inserts"
    assert_profile_contains "$output" "- Multiple mapped :       0" \
        "$name profile should report zero multimapped inserts"
    assert_profile_contains "$output" "- Uniquely mapped :       0" \
        "$name profile should report zero uniquely mapped inserts"
    assert_profile_contains "$output" "Effective insert-equivalents :       0.00" \
        "$name profile should report zero effective inserts"

    assert_profile_value "$output" Unknown 0 1e-9
    assert_profile_value "$output" A 0 1e-9
    assert_profile_value "$output" B 0 1e-9
}

run_zero_profile empty_ab "$script_dir/fixtures/profile_empty.sam" "ab"
run_zero_profile empty_rel "$script_dir/fixtures/profile_empty.sam" "rel"
run_zero_profile empty_tpm "$script_dir/fixtures/profile_empty.sam" "tpm"
run_zero_profile empty_fpkm "$script_dir/fixtures/profile_empty.sam" "fpkm"
run_zero_profile unmapped_ab "$script_dir/fixtures/profile_unmapped.sam" "ab"
run_zero_profile unmapped_rel "$script_dir/fixtures/profile_unmapped.sam" "rel"
run_zero_profile unmapped_tpm "$script_dir/fixtures/profile_unmapped.sam" "tpm"
run_zero_profile unmapped_fpkm "$script_dir/fixtures/profile_unmapped.sam" "fpkm"

fractional_output="$tmpdir/fractional.tsv.gz"
fractional_stderr="$tmpdir/fractional.stderr"
if ! "$MSAMTOOLS" profile -S \
    --label fractional \
    --unit ab \
    --nolen \
    --total 3 \
    --multi equal \
    --pandas \
    -o "$fractional_output" \
    "$fractional_mincount_fixture" \
    >/dev/null 2>"$fractional_stderr"; then
    cat "$fractional_stderr" >&2
    fail "fractional equal-sharing profile should succeed"
fi
pass_check "fractional equal-sharing profile should succeed"
assert_profile_value "$fractional_output" Unknown 0 1e-9
assert_profile_value "$fractional_output" A 1.333333333333 1e-6
assert_profile_value "$fractional_output" B 1.333333333333 1e-6
assert_profile_value "$fractional_output" C 0.333333333333 1e-6

mincount_output="$tmpdir/fractional_mincount.tsv.gz"
mincount_stderr="$tmpdir/fractional_mincount.stderr"
if ! "$MSAMTOOLS" profile -S \
    --label fractional_mincount \
    --unit ab \
    --nolen \
    --total 3 \
    --multi equal \
    --mincount 1 \
    --pandas \
    -o "$mincount_output" \
    "$fractional_mincount_fixture" \
    >/dev/null 2>"$mincount_stderr"; then
    cat "$mincount_stderr" >&2
    fail "fractional --mincount profile should succeed"
fi
pass_check "fractional --mincount profile should succeed"
assert_profile_value "$mincount_output" Unknown 0.333333333333 1e-6
assert_profile_value "$mincount_output" A 1.333333333333 1e-6
assert_profile_value "$mincount_output" B 1.333333333333 1e-6
assert_profile_value "$mincount_output" C 0 1e-9

# Default output should be pandas-style.
default_output="$tmpdir/default_output.tsv.gz"
if ! "$MSAMTOOLS" profile -S \
    --label default_output \
    --unit ab \
    --nolen \
    --total 7 \
    --multi equal \
    -o "$default_output" \
    "$fixture" \
    >/dev/null 2>"$tmpdir/default_output.stderr"; then
    cat "$tmpdir/default_output.stderr" >&2
    fail "profile should succeed with default output format"
fi
pass_check "profile should succeed with default output format"

default_header=$(gzip -cd -- "$default_output" |
    awk '!/^#/ { print; exit }')
if [ "$default_header" != "ID	default_output" ]; then
    fail "profile should use pandas-style output by default"
fi
pass_check "profile should use pandas-style output by default"


# Explicit --pandas remains accepted and produces the default format.
pandas_output="$tmpdir/pandas_output.tsv.gz"
if ! "$MSAMTOOLS" profile -S \
    --label pandas_output \
    --unit ab \
    --nolen \
    --total 7 \
    --multi equal \
    --pandas \
    -o "$pandas_output" \
    "$fixture" \
    >/dev/null 2>"$tmpdir/pandas_output.stderr"; then
    cat "$tmpdir/pandas_output.stderr" >&2
    fail "profile --pandas should remain accepted"
fi
pass_check "profile --pandas should remain accepted"

pandas_header=$(gzip -cd -- "$pandas_output" |
    awk '!/^#/ { print; exit }')
if [ "$pandas_header" != "ID	pandas_output" ]; then
    fail "--pandas should produce pandas-style output"
fi
pass_check "--pandas should produce pandas-style output"


# --no-pandas restores the legacy header format.
legacy_output="$tmpdir/legacy_output.tsv.gz"
if ! "$MSAMTOOLS" profile -S \
    --label legacy_output \
    --unit ab \
    --nolen \
    --total 7 \
    --multi equal \
    --no-pandas \
    -o "$legacy_output" \
    "$fixture" \
    >/dev/null 2>"$tmpdir/legacy_output.stderr"; then
    cat "$tmpdir/legacy_output.stderr" >&2
    fail "profile --no-pandas should succeed"
fi
pass_check "profile --no-pandas should succeed"

legacy_header=$(gzip -cd -- "$legacy_output" |
    awk '!/^#/ { print; exit }')
if [ "$legacy_header" != "legacy_output" ]; then
    fail "--no-pandas should restore legacy profile output"
fi
pass_check "--no-pandas should restore legacy profile output"

report_checks
exit 0
