#!/usr/bin/env python3
"""Run the Tier-2 msamtools profile validation grid.

This script orchestrates reproducible validation of ``msamtools profile`` using
synthetic paired-end alignments generated from complete bacterial chromosomes.

Default validation design
-------------------------
Main ambiguity experiment:
* every ``*.tsv`` community manifest in the supplied communities directory;
* 1,000, 10,000, and 100,000 simulated inserts;
* RNG seeds 13579, 24680, and 97531;
* multi-mapper modes: all, equal, ignore, and prop;
* generator setting ``--shared-fraction 0.20``.

If the communities directory contains four manifests, this yields
4 x 3 x 3 = 36 simulations.

No-sharing control:
* one designated control community (default: ``baseline``);
* the same three insert counts and three seeds;
* generator settings
  ``--shared-fraction 0 --exclude-cross-genome-multimappers --exclude-within-genome-repeats``.

This adds 9 simulations, for 45 total synthetic datasets and
45 x 4 = 180 profile evaluations.

Primary outputs
---------------
``run_metrics.tsv``
    One row per simulation x multi-mapper mode, including error metrics,
    Unknown abundance, realized ambiguity statistics, and E. coli strain
    estimates/errors.

``feature_estimates.tsv``
    Long-form truth, estimate, and signed/absolute error for every feature,
    including ``Unknown`` with true abundance zero.

``aggregate_metrics.tsv``
    Mean and sample standard deviation across seeds for each
    validation_set x community x insert-count x multi-mode combination.

``plots/``
    Graphical summaries in PNG format:
    * ``bray_curtis_summary.png``
    * ``ecoli_k12_summary.png``
    * ``no_sharing_control_bray_curtis.png``

``docs/validation/profile_validation.md``
    Release-facing Markdown summary generated automatically by the validation
    workflow. The three summary figures are copied to
    ``docs/validation/figures/``. The destination can be changed with
    ``--report-dir``.

The run-level directories retain the compact simulator outputs and the four
profile outputs. ``alignments.sam`` is deleted after successful profiling
unless ``--keep-alignments`` is supplied.

Usage
-----
Run the default grid::

    python3 validate_profiles.py \
        prepared_genomes/validated_genomes.tsv \
        communities \
        --msamtools ../msamtools

Use a specific generator and output directory::

    python3 validate_profiles.py \
        prepared_genomes/validated_genomes.tsv \
        communities \
        --msamtools ../msamtools \
        --generator generate_synthetic_alignments.py \
        --output-dir results/profile_validation

Retain the generated SAM files::

    python3 validate_profiles.py \
        prepared_genomes/validated_genomes.tsv \
        communities \
        --msamtools ../msamtools \
        --keep-alignments

Run a smaller development grid::

    python3 validate_profiles.py \
        prepared_genomes/validated_genomes.tsv \
        communities \
        --msamtools ../msamtools \
        --insert-counts 1000,10000 \
        --seeds 13579

Requirements
------------
* Python 3.9 or newer.
* No third-party Python packages.
* ``generate_synthetic_alignments.py``.
* A built ``msamtools`` executable.
* Chromosomes prepared by ``prepare_genomes.py``.
* Community TSVs accepted by the synthetic alignment generator.

The script exits non-zero immediately if generation, profiling, parsing, or
validation fails. Existing run directories are reused only when ``--resume`` is
supplied and all expected compact outputs/profile files are present.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


MULTI_MODES = ("all", "equal", "ignore", "prop")
DEFAULT_INSERT_COUNTS = (1000, 10000, 100000)
DEFAULT_SEEDS = (13579, 24680, 97531)

ECOLI_K12_ASSEMBLY = "GCF_000005845.2"
ECOLI_SAKAI_ASSEMBLY = "GCF_000008865.2"

PROFILE_MULTI_RE = re.compile(r"Multiple mapped\s*:\s*([0-9]+)")
PROFILE_TOTAL_RE = re.compile(r"Total inserts\s*:\s*([0-9]+)")
PROFILE_MAPPED_RE = re.compile(r"Mapped inserts\s*:\s*([0-9]+)")
PROFILE_VERSION_RE = re.compile(r"^# msamtools version\s+(.+?)\s*$", re.MULTILINE)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Generate and evaluate the Tier-2 msamtools profile validation grid."
        )
    )
    parser.add_argument(
        "validated_genomes",
        type=Path,
        help="validated_genomes.tsv produced by prepare_genomes.py",
    )
    parser.add_argument(
        "communities",
        type=Path,
        help="directory containing community *.tsv abundance manifests",
    )
    parser.add_argument(
        "--msamtools",
        type=Path,
        required=True,
        help="path to the msamtools executable to validate",
    )
    parser.add_argument(
        "--generator",
        type=Path,
        default=Path(__file__).with_name("generate_synthetic_alignments.py"),
        help="path to generate_synthetic_alignments.py (default: alongside this script)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("profile_validation_results"),
        help="validation output directory (default: profile_validation_results)",
    )
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "docs" / "validation",
        help=(
            "directory for the generated Markdown validation report and "
            "release-facing figures (default: <repository>/docs/validation)"
        ),
    )
    parser.add_argument(
        "--insert-counts",
        default=",".join(str(x) for x in DEFAULT_INSERT_COUNTS),
        help="comma-separated insert counts (default: 1000,10000,100000)",
    )
    parser.add_argument(
        "--seeds",
        default=",".join(str(x) for x in DEFAULT_SEEDS),
        help="comma-separated integer RNG seeds (default: 13579,24680,97531)",
    )
    parser.add_argument(
        "--shared-fraction",
        type=float,
        default=0.20,
        help="shared-locus enrichment passed to the generator (default: 0.20)",
    )
    parser.add_argument(
        "--single-mate-fraction",
        type=float,
        default=0.05,
        help="single-mate fraction passed to the generator (default: 0.05)",
    )
    parser.add_argument(
        "--candidate-multiplier",
        type=float,
        default=5.0,
        help="candidate multiplier passed to the generator (default: 5.0)",
    )
    parser.add_argument(
        "--mismatch-weights",
        default="0:0.50,1:0.30,2:0.10,3:0.10",
        help="mismatch weights passed to the generator",
    )
    parser.add_argument(
        "--control-community",
        default="baseline",
        help="stem of the community TSV used for the strict no-sharing control",
    )
    parser.add_argument(
        "--keep-alignments",
        action="store_true",
        help="retain alignments.sam after all four profiles succeed",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="reuse complete existing run outputs instead of regenerating them",
    )
    return parser.parse_args(argv)


def parse_int_list(text: str, label: str) -> tuple[int, ...]:
    """Parse a comma-separated list of positive integers."""

    try:
        values = tuple(int(token.strip()) for token in text.split(",") if token.strip())
    except ValueError as exc:
        raise ValueError(f"Invalid {label}: {text!r}") from exc
    if not values or any(value <= 0 for value in values):
        raise ValueError(f"{label} must contain positive integers")
    if len(set(values)) != len(values):
        raise ValueError(f"{label} contains duplicate values")
    return values


def parse_seed_list(text: str) -> tuple[int, ...]:
    """Parse a comma-separated list of integer RNG seeds."""

    try:
        values = tuple(int(token.strip()) for token in text.split(",") if token.strip())
    except ValueError as exc:
        raise ValueError(f"Invalid seeds: {text!r}") from exc
    if not values:
        raise ValueError("At least one seed is required")
    if len(set(values)) != len(values):
        raise ValueError("Seeds contain duplicate values")
    return values


def discover_communities(directory: Path) -> list[Path]:
    """Return sorted community TSV manifests from a directory."""

    if not directory.is_dir():
        raise ValueError(f"Community directory not found: {directory}")
    files = sorted(path for path in directory.glob("*.tsv") if path.is_file())
    if not files:
        raise ValueError(f"No *.tsv community manifests found in {directory}")
    return files


def run_command(command: Sequence[str], *, log_path: Path | None = None) -> str:
    """Run a subprocess, returning combined stdout/stderr and failing clearly."""

    result = subprocess.run(
        list(command),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    output = result.stdout or ""
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(output, encoding="utf-8")
    if result.returncode != 0:
        quoted = " ".join(command)
        raise RuntimeError(
            f"Command failed with exit status {result.returncode}:\n"
            f"{quoted}\n\n{output}"
        )
    return output


def read_key_value_tsv(path: Path) -> dict[str, str]:
    """Read a two-column ``metric/value`` TSV file."""

    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 2:
            raise ValueError(f"Malformed key/value TSV: {path}")
        for row in reader:
            if len(row) >= 2:
                values[row[0]] = row[1]
    return values


def read_reference_metadata(path: Path) -> dict[str, str]:
    """Map assembly accessions to strain feature IDs."""

    assembly_to_feature: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"assembly_accession", "strain_feature_id"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(f"Missing required columns in {path}")
        for row in reader:
            assembly = row["assembly_accession"].strip()
            feature = row["strain_feature_id"].strip()
            if assembly in assembly_to_feature:
                raise ValueError(f"Duplicate assembly in {path}: {assembly}")
            assembly_to_feature[assembly] = feature
    return assembly_to_feature


def read_truth(
    community_truth_path: Path,
    assembly_to_feature: dict[str, str],
) -> dict[str, float]:
    """Read true cell-relative abundances using strain feature IDs."""

    truth: dict[str, float] = {"Unknown": 0.0}
    with community_truth_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"assembly_accession", "true_cell_relative_abundance"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(f"Missing required columns in {community_truth_path}")
        for row in reader:
            assembly = row["assembly_accession"].strip()
            if assembly not in assembly_to_feature:
                raise ValueError(
                    f"Truth assembly {assembly} absent from reference metadata"
                )
            truth[assembly_to_feature[assembly]] = float(
                row["true_cell_relative_abundance"]
            )

    total = sum(truth.values())
    if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-8):
        raise ValueError(
            f"Truth abundances in {community_truth_path} sum to {total}, not 1"
        )
    return truth


def open_text_auto(path: Path):
    """Open plain text or gzip-compressed text by inspecting file contents.

    ``msamtools profile`` writes gzip-compressed output regardless of the
    filename extension, so suffix-based detection is not reliable here.
    Gzip files are identified by the standard 0x1f 0x8b magic bytes.
    """

    with path.open("rb") as handle:
        magic = handle.read(2)

    if magic == b"\x1f\x8b" or path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_profile(path: Path) -> dict[str, float]:
    """Read a two-column pandas-style msamtools profile."""

    values: dict[str, float] = {}
    with open_text_auto(path) as handle:
        reader = csv.reader(
            (line for line in handle if line.strip() and not line.startswith("#")),
            delimiter="\t",
        )
        header = next(reader, None)
        if header is None or len(header) < 2 or header[0] != "ID":
            raise ValueError(f"Unexpected msamtools profile header in {path}")
        for row in reader:
            if len(row) < 2:
                continue
            feature = row[0]
            value = float(row[1])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"Invalid abundance for {feature!r} in {path}: {value}")
            if feature in values:
                raise ValueError(f"Duplicate feature {feature!r} in {path}")
            values[feature] = value
    return values


def parse_profile_counts(path: Path) -> dict[str, int]:
    """Extract insert counts from the msamtools profile header."""

    with open_text_auto(path) as handle:
        header_text = "".join(line for line in handle if line.startswith("#"))
    result: dict[str, int] = {}
    for key, regex in (
        ("reported_total_inserts", PROFILE_TOTAL_RE),
        ("reported_mapped_inserts", PROFILE_MAPPED_RE),
        ("reported_multimapped_inserts", PROFILE_MULTI_RE),
    ):
        match = regex.search(header_text)
        if match is None:
            raise ValueError(f"Could not parse {key} from {path}")
        result[key] = int(match.group(1))
    return result


def read_expected_insert_counts(
    path: Path,
    assembly_to_feature: dict[str, str],
) -> dict[str, float]:
    """Count simulated source inserts per profiled feature."""

    counts = {
        feature: 0.0
        for feature in assembly_to_feature.values()
    }
    counts["Unknown"] = 0.0

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if (
            reader.fieldnames is None
            or "source_assembly" not in reader.fieldnames
        ):
            raise ValueError(
                f"Missing source_assembly column in {path}"
            )

        for row in reader:
            assembly = row["source_assembly"].strip()

            if assembly not in assembly_to_feature:
                raise ValueError(
                    f"Source assembly {assembly} absent from reference metadata"
                )

            counts[assembly_to_feature[assembly]] += 1.0

    return counts


def validate_exact_recovery(
    run_dir: Path,
    assembly_to_feature: dict[str, str],
    num_inserts: int,
) -> dict[str, dict[str, object]]:
    """Assert exact insert-count recovery in the no-sharing control."""

    expected = read_expected_insert_counts(
        run_dir / "insert_sources.tsv",
        assembly_to_feature,
    )

    if sum(expected.values()) != num_inserts:
        raise ValueError(
            f"Expected source counts sum to {sum(expected.values())}, "
            f"not {num_inserts}"
        )

    observed_by_mode: dict[str, dict[str, float]] = {}
    results: dict[str, dict[str, object]] = {}

    for mode in MULTI_MODES:
        profile_path = (
            run_dir
            / "profiles"
            / f"exact_{mode}.tsv.gz"
        )

        observed = read_profile(profile_path)

        if set(observed) != set(expected):
            missing = sorted(set(expected) - set(observed))
            extra = sorted(set(observed) - set(expected))

            raise ValueError(
                f"Exact-recovery feature mismatch for mode={mode}: "
                f"missing={missing}, unexpected={extra}"
            )

        errors = {
            feature: observed[feature] - expected[feature]
            for feature in expected
        }

        max_error = max(abs(value) for value in errors.values())
        exact_match = max_error == 0

        if not exact_match:
            mismatches = [
                f"{feature}: expected={expected[feature]}, "
                f"observed={observed[feature]}"
                for feature in expected
                if errors[feature] != 0
            ]

            raise ValueError(
                f"Exact insert-count recovery failed for mode={mode}: "
                + "; ".join(mismatches)
            )

        observed_by_mode[mode] = observed
        results[mode] = {
            "exact_count_match": 1,
            "max_count_error": max_error,
        }

    first_mode = MULTI_MODES[0]
    all_modes_identical = all(
        observed_by_mode[mode] == observed_by_mode[first_mode]
        for mode in MULTI_MODES[1:]
    )

    if not all_modes_identical:
        raise ValueError(
            "No-sharing exact-count profiles differ between --multi modes"
        )

    for mode in MULTI_MODES:
        results[mode]["all_modes_identical_counts"] = 1

    return results


def rankdata(values: Sequence[float]) -> list[float]:
    """Return average ranks (1-based) with correct handling of ties."""

    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        average_rank = ((i + 1) + j) / 2.0
        for k in range(i, j):
            ranks[indexed[k][0]] = average_rank
        i = j
    return ranks


def pearson(x: Sequence[float], y: Sequence[float]) -> float:
    """Compute Pearson correlation, returning NaN for zero-variance vectors."""

    if len(x) != len(y) or not x:
        raise ValueError("Pearson vectors must have equal non-zero length")
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    dx = [value - mean_x for value in x]
    dy = [value - mean_y for value in y]
    denom = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    if denom == 0:
        return math.nan
    return sum(a * b for a, b in zip(dx, dy)) / denom


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    """Compute Spearman rank correlation using average ranks for ties."""

    return pearson(rankdata(x), rankdata(y))


def calculate_metrics(truth: dict[str, float], estimate: dict[str, float]) -> dict[str, float]:
    """Calculate abundance-recovery metrics over the complete feature set."""

    if set(estimate) != set(truth):
        missing = sorted(set(truth) - set(estimate))
        extra = sorted(set(estimate) - set(truth))
        parts = []
        if missing:
            parts.append(f"missing={missing}")
        if extra:
            parts.append(f"unexpected={extra}")
        raise ValueError("Profile/truth feature mismatch: " + "; ".join(parts))

    features = sorted(truth)
    t = [truth[f] for f in features]
    e = [estimate[f] for f in features]
    errors = [estimate[f] - truth[f] for f in features]
    abs_errors = [abs(v) for v in errors]

    mae = sum(abs_errors) / len(abs_errors)
    rmse = math.sqrt(sum(v * v for v in errors) / len(errors))
    l1 = sum(abs_errors)
    denominator = sum(t) + sum(e)

    return {
        "mae": mae,
        "rmse": rmse,
        "tvd": 0.5 * l1,
        "bray_curtis": l1 / denominator if denominator > 0 else 0.0,
        "max_abs_error": max(abs_errors),
        "pearson": pearson(t, e),
        "spearman": spearman(t, e),
    }


def run_is_complete(
    run_dir: Path,
    exact_recovery: bool = False,
) -> bool:
    """Return whether all required simulation/profile outputs are present."""

    required = [
        run_dir / "community_truth.tsv",
        run_dir / "simulation_summary.tsv",
        run_dir / "reference_metadata.tsv",
        run_dir / "reference_to_strain.tsv",
        run_dir / "insert_sources.tsv",
    ]

    required.extend(
        run_dir / "profiles" / f"{mode}.tsv.gz"
        for mode in MULTI_MODES
    )
    required.extend(
        run_dir / "profiles" / f"{mode}.log"
        for mode in MULTI_MODES
    )

    if exact_recovery:
        required.extend(
            run_dir / "profiles" / f"exact_{mode}.tsv.gz"
            for mode in MULTI_MODES
        )
        required.extend(
            run_dir / "profiles" / f"exact_{mode}.log"
            for mode in MULTI_MODES
        )

    return all(path.exists() for path in required)


def generate_run(
    args: argparse.Namespace,
    community_path: Path,
    num_inserts: int,
    seed: int,
    run_dir: Path,
    *,
    shared_fraction: float,
    exclude_within_genome_repeats: bool,
    exclude_cross_genome_multimappers: bool,
    exact_recovery: bool,
) -> None:
    """Generate one synthetic SAM dataset unless a complete run is resumed."""

    if args.resume and run_is_complete(run_dir, exact_recovery=exact_recovery):
        print(f"  reusing complete run: {run_dir}", file=sys.stderr)
        return

    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    command = [
        sys.executable,
        str(args.generator.resolve()),
        str(args.validated_genomes.resolve()),
        str(community_path.resolve()),
        "--output-dir",
        str(run_dir.resolve()),
        "--num-inserts",
        str(num_inserts),
        "--seed",
        str(seed),
        "--shared-fraction",
        str(shared_fraction),
        "--single-mate-fraction",
        str(args.single_mate_fraction),
        "--candidate-multiplier",
        str(args.candidate_multiplier),
        "--mismatch-weights",
        args.mismatch_weights,
    ]
    if exclude_within_genome_repeats:
        command.append("--exclude-within-genome-repeats")
    if exclude_cross_genome_multimappers:
        command.append("--exclude-cross-genome-multimappers")
    print(f"  generating {community_path.stem}: n={num_inserts}, seed={seed}", file=sys.stderr)
    run_command(command, log_path=run_dir / "generator.log")


def run_profiles(
    args: argparse.Namespace,
    run_dir: Path,
    num_inserts: int,
    *,
    exact_recovery: bool,
) -> None:
    """Run relative profiles and optional exact-count recovery profiles."""

    if (
        args.resume
        and run_is_complete(
            run_dir,
            exact_recovery=exact_recovery,
        )
        and not (run_dir / "alignments.sam").exists()
    ):
        return

    sam_path = run_dir / "alignments.sam"
    genome_map = run_dir / "reference_to_strain.tsv"

    if not sam_path.exists():
        raise FileNotFoundError(
            f"Missing generated SAM: {sam_path}"
        )

    profile_dir = run_dir / "profiles"
    profile_dir.mkdir(parents=True, exist_ok=True)

    # Standard relative-abundance validation.
    for mode in MULTI_MODES:
        output_path = profile_dir / f"{mode}.tsv.gz"
        log_path = profile_dir / f"{mode}.log"

        command = [
            str(args.msamtools.resolve()),
            "profile",
            "--unit=rel",
            "--pandas",
            "--label",
            "test",
            "--genome",
            str(genome_map.resolve()),
            "--total",
            str(num_inserts),
            "--multi",
            mode,
            "-o",
            str(output_path.resolve()),
            str(sam_path.resolve()),
        ]

        print(
            f"    profile --multi {mode}",
            file=sys.stderr,
        )

        run_command(
            command,
            log_path=log_path,
        )

    # Strict no-sharing control:
    # raw insert counts must be recovered exactly.
    if exact_recovery:
        for mode in MULTI_MODES:
            output_path = (
                profile_dir / f"exact_{mode}.tsv.gz"
            )
            log_path = (
                profile_dir / f"exact_{mode}.log"
            )

            command = [
                str(args.msamtools.resolve()),
                "profile",
                "--unit=ab",
                "--nolen",
                "--pandas",
                "--label",
                "test",
                "--genome",
                str(genome_map.resolve()),
                "--total",
                str(num_inserts),
                "--multi",
                mode,
                "-o",
                str(output_path.resolve()),
                str(sam_path.resolve()),
            ]

            print(
                f"    exact recovery --multi {mode}",
                file=sys.stderr,
            )

            run_command(
                command,
                log_path=log_path,
            )

    if not args.keep_alignments:
        sam_path.unlink()


def evaluate_run(
    validation_set: str,
    community_name: str,
    num_inserts: int,
    seed: int,
    run_dir: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Evaluate all four profile modes for one synthetic dataset."""

    simulation = read_key_value_tsv(run_dir / "simulation_summary.tsv")
    assembly_to_feature = read_reference_metadata(run_dir / "reference_metadata.tsv")
    truth = read_truth(run_dir / "community_truth.tsv", assembly_to_feature)

    ecoli_k12_feature = assembly_to_feature.get(ECOLI_K12_ASSEMBLY)
    ecoli_sakai_feature = assembly_to_feature.get(ECOLI_SAKAI_ASSEMBLY)

    expected_multi = int(simulation["multiple_genomes_inserts"])
    realized_multi_fraction = float(simulation["realized_cross_genome_multimapper_fraction"])
    repeated_within = int(simulation["multiple_loci_single_genome_inserts"])
    single_mate_fraction = float(simulation["realized_single_mate_fraction"])
    sam_records = int(simulation["sam_alignment_records"])

    exact_recovery_results = None

    if validation_set == "no_sharing_control":
        if expected_multi != 0:
            raise ValueError(
                f"No-sharing control contains {expected_multi} "
                "cross-genome multi-mappers"
            )

        if repeated_within != 0:
            raise ValueError(
                f"No-sharing control contains {repeated_within} "
                "within-genome repeated-locus inserts"
            )

        exact_recovery_results = validate_exact_recovery(
            run_dir,
            assembly_to_feature,
            num_inserts,
        )

    run_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []

    for mode in MULTI_MODES:
        profile_path = run_dir / "profiles" / f"{mode}.tsv.gz"
        estimate = read_profile(profile_path)
        metrics = calculate_metrics(truth, estimate)
        counts = parse_profile_counts(profile_path)

        if counts["reported_total_inserts"] != num_inserts:
            raise ValueError(
                f"{community_name} n={num_inserts} seed={seed} mode={mode}: "
                f"reported total inserts {counts['reported_total_inserts']} "
                f"!= requested {num_inserts}"
            )
        if counts["reported_mapped_inserts"] != num_inserts:
            raise ValueError(
                f"{community_name} n={num_inserts} seed={seed} mode={mode}: "
                f"reported mapped inserts {counts['reported_mapped_inserts']} "
                f"!= generated {num_inserts}"
            )
        if counts["reported_multimapped_inserts"] != expected_multi:
            raise ValueError(
                f"{community_name} n={num_inserts} seed={seed} mode={mode}: "
                f"reported multi-mappers {counts['reported_multimapped_inserts']} "
                f"!= simulator truth {expected_multi}"
            )

        estimate_sum = sum(estimate.values())
        if not math.isclose(estimate_sum, 1.0, rel_tol=0.0, abs_tol=5e-6):
            raise ValueError(f"profile sums to {estimate_sum}, not 1")

        row: dict[str, object] = {
            "validation_set": validation_set,
            "community": community_name,
            "num_inserts": num_inserts,
            "seed": seed,
            "multi_mode": mode,
            "realized_multimapper_fraction": realized_multi_fraction,
            "realized_multimapper_inserts": expected_multi,
            "within_genome_repeat_inserts": repeated_within,
            "realized_single_mate_fraction": single_mate_fraction,
            "sam_alignment_records": sam_records,
            **counts,
            **metrics,
            "unknown_estimate": estimate["Unknown"],
        }

        if exact_recovery_results is not None:
            row.update(exact_recovery_results[mode])

        if ecoli_k12_feature is not None:
            k12_truth = truth[ecoli_k12_feature]
            k12_estimate = estimate[ecoli_k12_feature]
            row["ecoli_k12_truth"] = k12_truth
            row["ecoli_k12_estimate"] = k12_estimate
            row["ecoli_k12_error"] = k12_estimate - k12_truth
            row["ecoli_k12_false_positive_abundance"] = k12_estimate if k12_truth == 0 else math.nan

        if ecoli_sakai_feature is not None:
            sakai_truth = truth[ecoli_sakai_feature]
            sakai_estimate = estimate[ecoli_sakai_feature]
            row["ecoli_sakai_truth"] = sakai_truth
            row["ecoli_sakai_estimate"] = sakai_estimate
            row["ecoli_sakai_error"] = sakai_estimate - sakai_truth

        run_rows.append(row)

        for feature in sorted(truth):
            signed_error = estimate[feature] - truth[feature]
            feature_rows.append(
                {
                    "validation_set": validation_set,
                    "community": community_name,
                    "num_inserts": num_inserts,
                    "seed": seed,
                    "multi_mode": mode,
                    "feature": feature,
                    "truth": truth[feature],
                    "estimate": estimate[feature],
                    "error": signed_error,
                    "abs_error": abs(signed_error),
                }
            )

    return run_rows, feature_rows


def write_dict_rows(path: Path, rows: Sequence[dict[str, object]]) -> None:
    """Write dictionaries to TSV using stable first-seen column ordering."""

    if not rows:
        raise ValueError(f"No rows to write: {path}")

    fields = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fields.append(key)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def mean_sd(values: Sequence[float]) -> tuple[float, float]:
    """Return mean and sample SD; SD is zero for one replicate."""

    return (
        statistics.fmean(values),
        statistics.stdev(values) if len(values) > 1 else 0.0,
    )


def aggregate_run_metrics(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    """Aggregate numerical performance metrics across seeds."""

    groups: dict[tuple[str, str, int, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = (
            str(row["validation_set"]),
            str(row["community"]),
            int(row["num_inserts"]),
            str(row["multi_mode"]),
        )
        groups[key].append(row)

    metric_names = [
        "mae",
        "rmse",
        "tvd",
        "bray_curtis",
        "max_abs_error",
        "pearson",
        "spearman",
        "unknown_estimate",
        "ecoli_k12_truth",
        "ecoli_k12_estimate",
        "ecoli_k12_error",
        "ecoli_k12_false_positive_abundance",
        "ecoli_sakai_truth",
        "ecoli_sakai_estimate",
        "ecoli_sakai_error",
        "realized_multimapper_fraction",
    ]

    aggregated: list[dict[str, object]] = []
    for (validation_set, community, num_inserts, mode), group in sorted(groups.items()):
        row: dict[str, object] = {
            "validation_set": validation_set,
            "community": community,
            "num_inserts": num_inserts,
            "multi_mode": mode,
            "n_seeds": len(group),
        }
        for metric in metric_names:
            values = [float(item[metric]) for item in group if metric in item and math.isfinite(float(item[metric]))]
            if values:
                mean, sd = mean_sd(values)
                row[f"mean_{metric}"] = mean
                row[f"sd_{metric}"] = sd
        aggregated.append(row)
    return aggregated


def _select_aggregate_row(aggregate_rows, validation_set, community, num_inserts, multi_mode):
    for row in aggregate_rows:
        if (
            row["validation_set"] == validation_set
            and row["community"] == community
            and int(row["num_inserts"]) == num_inserts
            and row["multi_mode"] == multi_mode
        ):
            return row
    return None


def make_plots(
    output_dir: Path,
    aggregate_rows,
    main_communities,
    control_community,
    insert_counts,
):
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # 1. Bray-Curtis distance from truth
    #    One facet per community, x-axis = insert count.
    # ------------------------------------------------------------
    fig, axes = plt.subplots(
        1,
        len(main_communities),
        figsize=(5 * len(main_communities), 5),
        sharey=True,
    )

    if len(main_communities) == 1:
        axes = [axes]

    for ax, community in zip(axes, main_communities):
        for mode in MULTI_MODES:
            y = []
            yerr = []

            for n in insert_counts:
                row = _select_aggregate_row(
                    aggregate_rows,
                    "main",
                    community,
                    n,
                    mode,
                )

                y.append(float(row["mean_bray_curtis"]))
                yerr.append(float(row["sd_bray_curtis"]))

            ax.errorbar(
                insert_counts,
                y,
                yerr=yerr,
                marker="o",
                label=mode,
            )

        ax.set_xscale("log")
        ax.set_xticks(insert_counts)
        ax.set_xticklabels([str(n) for n in insert_counts])
        ax.set_xlabel("Insert count")
        ax.set_title(community)

    axes[0].set_ylabel("Bray-Curtis dissimilarity")

    fig.suptitle(
        "Bray-Curtis distance from the true community composition"
    )

    axes[-1].legend(title="--multi")

    fig.tight_layout()
    fig.savefig(
        plot_dir / "bray_curtis_summary.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------
    # 2. E. coli K-12 abundance
    #    One facet per community, x-axis = insert count.
    # ------------------------------------------------------------
    fig, axes = plt.subplots(
        1,
        len(main_communities),
        figsize=(5 * len(main_communities), 5),
        sharey=True,
    )

    if len(main_communities) == 1:
        axes = [axes]

    markers = {
        "all": "o",
        "equal": "s",
        "ignore": "^",
        "prop": "D",
    }

    linestyles = {
        "all": "-",
        "equal": "-",
        "ignore": "--",
        "prop": "-.",
    }

    if mode == "ignore":
        markerfacecolor = "none"
    else:
        markerfacecolor = None

    for ax, community in zip(axes, main_communities):
        for mode in MULTI_MODES:
            y = []
            yerr = []

            for n in insert_counts:
                row = _select_aggregate_row(
                    aggregate_rows,
                    "main",
                    community,
                    n,
                    mode,
                )

                y.append(
                    float(row["mean_ecoli_k12_estimate"])
                )
                yerr.append(
                    float(row["sd_ecoli_k12_estimate"])
                )

            ax.errorbar(
                insert_counts,
                y,
                yerr=yerr,
                marker=markers[mode],
                linestyle=linestyles[mode],
                alpha=0.8,
                markersize=6,
                linewidth=1.5,
                markerfacecolor=markerfacecolor,
                label=mode,
            )

        # Truth is identical across seeds, depths, and profiling modes
        # within a given community.
        truth_row = _select_aggregate_row(
            aggregate_rows,
            "main",
            community,
            insert_counts[0],
            MULTI_MODES[0],
        )
        truth = float(
            truth_row["mean_ecoli_k12_truth"]
        )

        ax.axhline(
            truth,
            linestyle="--",
            linewidth=2.0,
            color="black",
            label="truth",
            zorder=10,
        )

        ax.set_xscale("log")
        ax.set_xticks(insert_counts)
        ax.set_xticklabels([str(n) for n in insert_counts])
        ax.set_xlabel("Insert count")
        ax.set_title(community)

    axes[0].set_ylabel(
        "E. coli K-12 relative abundance"
    )

    fig.suptitle(
        "Estimated E. coli K-12 abundance"
    )

    axes[-1].legend(title="Series")

    fig.tight_layout()
    fig.savefig(
        plot_dir / "ecoli_k12_summary.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------
    # 3. Strict no-sharing control
    #    One facet per --multi mode to demonstrate equivalence.
    # ------------------------------------------------------------
    fig, axes = plt.subplots(
        1,
        len(MULTI_MODES),
        figsize=(4.5 * len(MULTI_MODES), 4.5),
        sharey=True,
    )

    if len(MULTI_MODES) == 1:
        axes = [axes]

    for ax, mode in zip(axes, MULTI_MODES):
        y = []
        yerr = []

        for n in insert_counts:
            row = _select_aggregate_row(
                aggregate_rows,
                "no_sharing_control",
                control_community,
                n,
                mode,
            )

            y.append(
                float(row["mean_bray_curtis"])
            )
            yerr.append(
                float(row["sd_bray_curtis"])
            )

        ax.errorbar(
            insert_counts,
            y,
            yerr=yerr,
            marker="o",
        )

        ax.set_xscale("log")
        ax.set_xticks(insert_counts)
        ax.set_xticklabels([str(n) for n in insert_counts])
        ax.set_xlabel("Insert count")
        ax.set_title(mode)

    axes[0].set_ylabel(
        "Bray-Curtis dissimilarity"
    )

    fig.suptitle(
        "No-sharing control: Bray-Curtis distance from the "
        f"true community composition\n{control_community}"
    )

    fig.tight_layout()
    fig.savefig(
        plot_dir / "no_sharing_control_bray_curtis.png",
        dpi=200,
    )
    plt.close(fig)



def read_msamtools_version(path: Path) -> str:
    """Read the msamtools version recorded in a profile header."""

    with open_text_auto(path) as handle:
        header_text = "".join(
            line for line in handle
            if line.startswith("#")
        )

    match = PROFILE_VERSION_RE.search(header_text)
    if match is None:
        raise ValueError(
            f"Could not parse msamtools version from {path}"
        )

    return match.group(1).strip()


def get_repository_commit() -> str:
    """Return the Git commit containing this validation script, if available."""

    script_dir = Path(__file__).resolve().parent
    result = subprocess.run(
        ["git", "-C", str(script_dir), "rev-parse", "HEAD"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )

    if result.returncode != 0:
        return "unknown (not a Git checkout)"

    commit = result.stdout.strip()
    return commit or "unknown"


def format_metric(value: float) -> str:
    """Format report metrics compactly without hiding small values."""

    if value == 0:
        return "0"
    if abs(value) < 1e-4:
        return f"{value:.3e}"
    return f"{value:.6f}".rstrip("0").rstrip(".")


def write_validation_report(
    report_dir: Path,
    output_dir: Path,
    run_specs,
    run_rows,
    community_paths,
    insert_counts,
    seeds,
    shared_fraction: float,
) -> Path:
    """Write the release-facing Markdown report and copy its figures."""

    report_dir = report_dir.resolve()
    figure_dir = report_dir / "figures"
    report_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    figure_names = (
        "bray_curtis_summary.png",
        "ecoli_k12_summary.png",
        "no_sharing_control_bray_curtis.png",
    )

    for name in figure_names:
        source = output_dir / "plots" / name
        if not source.is_file():
            raise FileNotFoundError(
                f"Expected validation figure not found: {source}"
            )
        shutil.copy2(source, figure_dir / name)

    versions = {
        read_msamtools_version(
            Path(spec["run_dir"]) / "profiles" / "prop.tsv.gz"
        )
        for spec in run_specs
    }
    if len(versions) != 1:
        raise ValueError(
            "Validation runs contain multiple msamtools versions: "
            + ", ".join(sorted(versions))
        )
    msamtools_version = next(iter(versions))

    commit = get_repository_commit()
    validation_date = date.today().isoformat()

    main_rows = [
        row for row in run_rows
        if row["validation_set"] == "main"
    ]
    control_rows = [
        row for row in run_rows
        if row["validation_set"] == "no_sharing_control"
    ]

    if not main_rows or not control_rows:
        raise ValueError(
            "Both main and no-sharing validation results are required "
            "to generate the release-facing report"
        )

    exact_count_pass = all(
        int(row.get("exact_count_match", 0)) == 1
        for row in control_rows
    )
    identical_modes_pass = all(
        int(row.get("all_modes_identical_counts", 0)) == 1
        for row in control_rows
    )
    max_count_error = max(
        float(row.get("max_count_error", math.inf))
        for row in control_rows
    )

    if (
        not exact_count_pass
        or not identical_modes_pass
        or max_count_error != 0
    ):
        raise ValueError(
            "Exact-recovery validation did not satisfy all release-report invariants"
        )

    mode_stats = {}
    for mode in MULTI_MODES:
        rows = [
            row for row in main_rows
            if row["multi_mode"] == mode
        ]
        distances = [float(row["bray_curtis"]) for row in rows]
        mode_stats[mode] = {
            "mean": statistics.fmean(distances),
            "max": max(distances),
        }

    prop_rows = [
        row for row in main_rows
        if row["multi_mode"] == "prop"
    ]

    absent_rows = [
        row for row in prop_rows
        if math.isclose(
            float(row["ecoli_k12_truth"]),
            0.0,
            rel_tol=0.0,
            abs_tol=1e-15,
        )
    ]

    absent_community = None
    absent_max_estimate = None
    if absent_rows:
        absent_community = ", ".join(
            sorted({str(row["community"]) for row in absent_rows})
        )
        absent_max_estimate = max(
            float(row["ecoli_k12_estimate"])
            for row in absent_rows
        )

    positive_truths = sorted({
        float(row["ecoli_k12_truth"])
        for row in prop_rows
        if float(row["ecoli_k12_truth"]) > 0
    })

    rare_truth = positive_truths[0] if positive_truths else None
    rare_rows = []
    if rare_truth is not None:
        rare_rows = [
            row for row in prop_rows
            if math.isclose(
                float(row["ecoli_k12_truth"]),
                rare_truth,
                rel_tol=0.0,
                abs_tol=1e-15,
            )
        ]

    rare_community = None
    rare_mean_estimate = None
    rare_max_estimate = None
    if rare_rows:
        rare_community = ", ".join(
            sorted({str(row["community"]) for row in rare_rows})
        )
        rare_estimates = [
            float(row["ecoli_k12_estimate"])
            for row in rare_rows
        ]
        rare_mean_estimate = statistics.fmean(rare_estimates)
        rare_max_estimate = max(rare_estimates)

    control_bray_max = max(
        float(row["bray_curtis"])
        for row in control_rows
    )

    main_run_count = sum(
        spec["validation_set"] == "main"
        for spec in run_specs
    )
    control_run_count = sum(
        spec["validation_set"] == "no_sharing_control"
        for spec in run_specs
    )
    relative_profile_count = len(run_specs) * len(MULTI_MODES)
    exact_profile_count = control_run_count * len(MULTI_MODES)

    first_run_dir = Path(run_specs[0]["run_dir"])
    n_references = len(
        read_reference_metadata(
            first_run_dir / "reference_metadata.tsv"
        )
    )

    communities_text = ", ".join(
        f"`{path.stem}`" for path in community_paths
    )
    depths_text = "; ".join(f"{n:,}" for n in insert_counts)
    seeds_text = ", ".join(str(seed) for seed in seeds)

    lines = [
        "# msamtools profile validation",
        "",
        "> **Generated automatically by `validation/validate_profiles.py`. "
        "Do not edit numerical results manually.**",
        "",
        f"- **msamtools version:** `{msamtools_version}`",
        f"- **Git commit:** `{commit}`",
        f"- **Validation date:** {validation_date}",
        "",
        "## Validation design",
        "",
        f"- {n_references} complete bacterial chromosome references",
        f"- Main communities: {communities_text}",
        f"- Insert counts: {depths_text}",
        f"- RNG seeds: {seeds_text}",
        f"- Multi-mapper modes: {', '.join(f'`{mode}`' for mode in MULTI_MODES)}",
        f"- Main shared-locus target fraction: {shared_fraction:g}",
        f"- Main multimapping simulations: {main_run_count}",
        f"- Strict no-sharing controls: {control_run_count}",
        f"- Relative-abundance profile evaluations: {relative_profile_count}",
        f"- Additional exact-count profile evaluations: {exact_profile_count}",
        "",
        "The strict no-sharing control excludes both cross-genome multimappers "
        "and within-genome repeated loci.",
        "",
        "## Exact recovery without ambiguous mapping",
        "",
        "**PASS**",
        "",
        "- Exact per-reference insert recovery: **PASS**",
        f"- Maximum insert-count error: **{format_metric(max_count_error)}**",
        "- All `--multi` modes produced identical counts: **PASS**",
        "",
        "With ambiguous mappings excluded, `msamtools profile --unit=ab --nolen` "
        "recovered the exact simulated source insert count for every profiled "
        "reference under every multi-mapper mode.",
        "",
        "## Composition recovery with multimapping",
        "",
        "Bray-Curtis dissimilarity was calculated between each estimated "
        "relative-abundance profile and its known generating composition.",
        "",
        "| `--multi` mode | Mean Bray-Curtis | Maximum Bray-Curtis |",
        "|---|---:|---:|",
    ]

    for mode in MULTI_MODES:
        lines.append(
            f"| `{mode}` | {format_metric(mode_stats[mode]['mean'])} | "
            f"{format_metric(mode_stats[mode]['max'])} |"
        )

    lines.extend([
        "",
        "![Bray-Curtis distance from the true community composition]"
        "(figures/bray_curtis_summary.png)",
        "",
        "## Closely related strain stress tests",
        "",
        "The two *Escherichia coli* reference strains provide a natural "
        "multi-mapping challenge in which one strain can be made rare or "
        "completely absent while the related strain remains abundant.",
        "",
        "### Absent strain",
        "",
    ])

    if absent_rows:
        lines.extend([
            f"- Community: `{absent_community}`",
            "- True *E. coli* K-12 relative abundance: **0**",
            "- Maximum abundance estimated with proportional sharing: "
            f"**{format_metric(absent_max_estimate)}**",
        ])
        if absent_max_estimate == 0:
            lines.append(
                "- False-positive abundance with proportional sharing: **none observed**"
            )
    else:
        lines.append(
            "No zero-abundance *E. coli* K-12 community was present in this run."
        )

    lines.extend(["", "### Rare strain", ""])

    if rare_rows:
        lines.extend([
            f"- Community: `{rare_community}`",
            "- True *E. coli* K-12 relative abundance: "
            f"**{format_metric(rare_truth)}**",
            "- Mean abundance estimated with proportional sharing: "
            f"**{format_metric(rare_mean_estimate)}**",
            "- Maximum abundance estimated with proportional sharing: "
            f"**{format_metric(rare_max_estimate)}**",
        ])
    else:
        lines.append(
            "No positive rare-strain community was available in this run."
        )

    lines.extend([
        "",
        "![Estimated E. coli K-12 abundance](figures/ecoli_k12_summary.png)",
        "",
        "## No-sharing control",
        "",
        "In the strict no-sharing control, all profile modes recover the same "
        "exact insert counts. Relative-abundance estimates differ from the "
        "generating composition only through numerical normalization and output precision.",
        "",
        "- Maximum Bray-Curtis distance from truth across all no-sharing runs "
        f"and modes: **{format_metric(control_bray_max)}**",
        "",
        "![No-sharing control](figures/no_sharing_control_bray_curtis.png)",
        "",
        "## Reproducibility",
        "",
        "The validation workflow, genome manifest, community definitions, and "
        "simulation code are maintained under `validation/`. Detailed run-level "
        "TSV files are generated during validation and serve as machine-readable "
        "audit outputs; this Markdown file and the figures above are the "
        "release-facing summary.",
        "",
    ])

    report_path = report_dir / "profile_validation.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path

def validate_paths(args: argparse.Namespace) -> None:
    """Validate required executable/script/input paths."""

    if not args.validated_genomes.is_file():
        raise ValueError(f"Validated genome file not found: {args.validated_genomes}")
    if not args.generator.is_file():
        raise ValueError(f"Generator not found: {args.generator}")
    if not args.msamtools.is_file():
        raise ValueError(f"msamtools executable not found: {args.msamtools}")
    if not (0 <= args.shared_fraction <= 1):
        raise ValueError("--shared-fraction must be between 0 and 1")
    if not (0 <= args.single_mate_fraction <= 1):
        raise ValueError("--single-mate-fraction must be between 0 and 1")
    if args.candidate_multiplier < 1:
        raise ValueError("--candidate-multiplier must be >= 1")

def main(argv: list[str] | None = None) -> int:
    """Run the complete Tier-2 profile validation grid."""

    args = parse_args(argv)
    try:
        validate_paths(args)
        insert_counts = parse_int_list(args.insert_counts, "insert counts")
        seeds = parse_seed_list(args.seeds)
        community_paths = discover_communities(args.communities)
        community_by_stem = {path.stem: path for path in community_paths}
        if args.control_community not in community_by_stem:
            raise ValueError(f"Control community {args.control_community!r} not found in {args.communities}")

        args.output_dir.mkdir(parents=True, exist_ok=True)

        run_specs = []
        for community_path in community_paths:
            for num_inserts in insert_counts:
                for seed in seeds:
                    run_specs.append(
                        {
                            "validation_set": "main",
                            "community_path": community_path,
                            "community_name": community_path.stem,
                            "num_inserts": num_inserts,
                            "seed": seed,
                            "shared_fraction": args.shared_fraction,
                            "exclude_within_genome_repeats": False,
                            "exclude_cross_genome_multimappers": False,
                            "exact_recovery": False,
                            "run_dir": args.output_dir / "main" / community_path.stem / f"n{num_inserts}_seed{seed}",
                        }
                    )

        control_path = community_by_stem[args.control_community]
        for num_inserts in insert_counts:
            for seed in seeds:
                run_specs.append(
                    {
                        "validation_set": "no_sharing_control",
                        "community_path": control_path,
                        "community_name": control_path.stem,
                        "num_inserts": num_inserts,
                        "seed": seed,
                        "shared_fraction": 0.0,
                        "exclude_within_genome_repeats": True,
                        "exclude_cross_genome_multimappers": True,
                        "exact_recovery": True,
                        "run_dir": args.output_dir / "no_sharing_control" / control_path.stem / f"n{num_inserts}_seed{seed}",
                    }
                )

        all_run_rows = []
        all_feature_rows = []

        total_runs = len(run_specs)
        for index, spec in enumerate(run_specs, start=1):
            print(
                f"\n[{index}/{total_runs}] {spec['validation_set']}: {spec['community_name']}, "
                f"n={spec['num_inserts']}, seed={spec['seed']}",
                file=sys.stderr,
            )

            generate_run(
                args,
                Path(spec["community_path"]),
                int(spec["num_inserts"]),
                int(spec["seed"]),
                Path(spec["run_dir"]),
                shared_fraction=float(spec["shared_fraction"]),
                exclude_within_genome_repeats=bool(spec["exclude_within_genome_repeats"]),
                exclude_cross_genome_multimappers=bool(spec["exclude_cross_genome_multimappers"]),
                exact_recovery=bool(spec["exact_recovery"]),
            )
            run_profiles(args, Path(spec["run_dir"]), int(spec["num_inserts"]), exact_recovery=bool(spec["exact_recovery"]))
            run_rows, feature_rows = evaluate_run(
                str(spec["validation_set"]),
                str(spec["community_name"]),
                int(spec["num_inserts"]),
                int(spec["seed"]),
                Path(spec["run_dir"]),
            )
            all_run_rows.extend(run_rows)
            all_feature_rows.extend(feature_rows)

        write_dict_rows(args.output_dir / "run_metrics.tsv", all_run_rows)
        write_dict_rows(args.output_dir / "feature_estimates.tsv", all_feature_rows)
        aggregate_rows = aggregate_run_metrics(all_run_rows)
        write_dict_rows(args.output_dir / "aggregate_metrics.tsv", aggregate_rows)

        make_plots(
            args.output_dir,
            aggregate_rows,
            [path.stem for path in community_paths],
            args.control_community,
            insert_counts,
        )
        report_path = write_validation_report(
            args.report_dir,
            args.output_dir,
            run_specs,
            all_run_rows,
            community_paths,
            insert_counts,
            seeds,
            args.shared_fraction,
        )

        print("\nValidation complete", file=sys.stderr)
        print("-------------------", file=sys.stderr)
        print(f"Synthetic datasets: {total_runs}", file=sys.stderr)
        print(f"Profile evaluations: {total_runs * len(MULTI_MODES)}", file=sys.stderr)
        print(f"Results: {args.output_dir / 'run_metrics.tsv'}", file=sys.stderr)
        print(f"Summary: {args.output_dir / 'aggregate_metrics.tsv'}", file=sys.stderr)
        print(f"Plots:   {args.output_dir / 'plots'}", file=sys.stderr)
        print(f"Report:  {report_path}", file=sys.stderr)
        return 0

    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
