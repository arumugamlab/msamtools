#!/usr/bin/env python3
"""Download and validate chromosome-only reference genomes for msamtools tests.

This script is the first stage of the reproducible msamtools validation-data
workflow. It reads a TSV manifest of NCBI assembly accessions, downloads the
NCBI ``datasets`` command-line program into a local cache, downloads each
assembly, extracts only the chromosome, and checks that each selected assembly
satisfies the assumptions required by the downstream simulation workflow.

Validation criteria
-------------------
For each manifest entry, the script requires:

1. The requested assembly accession is returned by NCBI exactly.
2. The assembly level is ``Complete Genome``.
3. Exactly one sequence is classified by NCBI as a chromosome.
4. The chromosome is marked ``circular`` in its GenBank/RefSeq GBFF record.
5. The extracted chromosome FASTA contains only A/C/G/T bases (case ignored).
6. The FASTA sequence length agrees with the NCBI sequence report.

Assemblies may contain plasmids or other non-chromosomal replicons; these are
reported and ignored. The script validates all manifest entries before exiting,
so a single bad assembly does not prevent problems in other entries from being
reported. The final exit status is non-zero if any assembly fails validation.

Usage
-----
Typical use::

    python3 prepare_genomes.py genomes.tsv

Choose a different output directory::

    python3 prepare_genomes.py genomes.tsv --output-dir prepared_genomes

Force a fresh download of the NCBI ``datasets`` executable::

    python3 prepare_genomes.py genomes.tsv --refresh-datasets

Keep the downloaded NCBI data packages for inspection/debugging::

    python3 prepare_genomes.py genomes.tsv --keep-packages

Input manifest
--------------
The input is a tab-separated file with a header and these required columns::

    assembly_accession    species    strain

Additional columns are allowed and ignored by this script.

Outputs
-------
By default, ``prepared_genomes/`` contains:

``chromosomes/``
    One chromosome-only FASTA file per validated assembly.
``validated_genomes.tsv``
    Machine-readable validation report and metadata for downstream scripts.
``.tools/datasets``
    Locally cached NCBI Datasets executable (``datasets.exe`` on Windows).
``packages/``
    Only retained when ``--keep-packages`` is used.

Requirements
------------
* Python 3.9 or newer.
* Internet access to ``ftp.ncbi.nlm.nih.gov`` and NCBI Datasets services.
* No third-party Python packages are required.

Notes
-----
The NCBI genome sequence report identifies chromosome records but does not
currently expose circular topology. Therefore the script also downloads GBFF
and verifies ``circular`` from the LOCUS line of the selected chromosome record.
The NCBI ``datasets`` executable is downloaded from the official NCBI FTP site
for the current operating system/architecture and cached locally.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator


NCBI_DATASETS_BASE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2"
)
REQUIRED_MANIFEST_COLUMNS = {"assembly_accession", "species", "strain"}
VALID_DNA_BASES = frozenset("ACGT")


@dataclass(frozen=True)
class GenomeRequest:
    """One requested genome from the input manifest."""

    assembly_accession: str
    species: str
    strain: str


@dataclass
class ValidationResult:
    """Validation metadata and status for one requested assembly."""

    requested_accession: str
    expected_species: str
    expected_strain: str
    status: str = "FAIL"
    error: str = ""
    assembly_accession: str = ""
    assembly_level: str = ""
    observed_species: str = ""
    observed_strain: str = ""
    chromosome_accession: str = ""
    chromosome_length: int | None = None
    circular: bool | None = None
    n_count: int | None = None
    non_acgt_count: int | None = None
    ignored_nonchromosomal_sequences: int | None = None
    sha256: str = ""
    fasta_path: str = ""
    notes: list[str] = field(default_factory=list)

    def as_row(self) -> dict[str, str]:
        """Return a flat dictionary suitable for TSV output."""

        return {
            "requested_accession": self.requested_accession,
            "assembly_accession": self.assembly_accession,
            "expected_species": self.expected_species,
            "observed_species": self.observed_species,
            "expected_strain": self.expected_strain,
            "observed_strain": self.observed_strain,
            "assembly_level": self.assembly_level,
            "chromosome_accession": self.chromosome_accession,
            "chromosome_length": _fmt_optional(self.chromosome_length),
            "circular": _fmt_bool(self.circular),
            "n_count": _fmt_optional(self.n_count),
            "non_acgt_count": _fmt_optional(self.non_acgt_count),
            "ignored_nonchromosomal_sequences": _fmt_optional(
                self.ignored_nonchromosomal_sequences
            ),
            "sha256": self.sha256,
            "fasta_path": self.fasta_path,
            "status": self.status,
            "error": self.error,
            "notes": "; ".join(self.notes),
        }


def _fmt_optional(value: object | None) -> str:
    """Format an optional scalar for TSV output."""

    return "" if value is None else str(value)


def _fmt_bool(value: bool | None) -> str:
    """Format an optional boolean as yes/no for human-readable TSV output."""

    if value is None:
        return ""
    return "yes" if value else "no"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Download NCBI assemblies, extract chromosomes, and validate that "
            "they are single, complete, circular, unambiguous genomes."
        )
    )
    parser.add_argument(
        "manifest",
        type=Path,
        help="TSV containing assembly_accession, species, and strain columns",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("prepared_genomes"),
        help="output directory (default: prepared_genomes)",
    )
    parser.add_argument(
        "--refresh-datasets",
        action="store_true",
        help="redownload the cached NCBI datasets executable",
    )
    parser.add_argument(
        "--keep-packages",
        action="store_true",
        help="keep downloaded/unpacked NCBI packages for inspection",
    )
    return parser.parse_args(argv)


def read_manifest(path: Path) -> list[GenomeRequest]:
    """Read and validate the genome manifest.

    Parameters
    ----------
    path
        Path to the tab-separated manifest.

    Returns
    -------
    list of GenomeRequest
        Manifest records in input order.

    Raises
    ------
    ValueError
        If required columns are missing, values are blank, or accessions are
        duplicated.
    """

    requests: list[GenomeRequest] = []
    seen: set[str] = set()

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Manifest has no header: {path}")

        missing = REQUIRED_MANIFEST_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(
                "Manifest is missing required column(s): " + ", ".join(sorted(missing))
            )

        for line_number, row in enumerate(reader, start=2):
            accession = row["assembly_accession"].strip()
            species = row["species"].strip()
            strain = row["strain"].strip()
            if not accession or not species or not strain:
                raise ValueError(
                    f"Blank required value in manifest line {line_number}"
                )
            if accession in seen:
                raise ValueError(f"Duplicate assembly accession: {accession}")
            seen.add(accession)
            requests.append(GenomeRequest(accession, species, strain))

    if not requests:
        raise ValueError("Manifest contains no genome records")
    return requests


def datasets_download_url() -> tuple[str, str]:
    """Return the official NCBI Datasets binary URL and local file name.

    Returns
    -------
    (url, filename)
        Platform-specific NCBI download URL and executable file name.

    Raises
    ------
    RuntimeError
        If the operating system/architecture is not supported by NCBI's
        published standalone binaries.
    """

    system = platform.system().lower()
    machine = platform.machine().lower()

    if system == "linux":
        if machine in {"x86_64", "amd64"}:
            platform_dir = "linux-amd64"
        elif machine in {"aarch64", "arm64"}:
            platform_dir = "linux-arm64"
        elif machine.startswith("arm"):
            platform_dir = "linux-arm"
        else:
            raise RuntimeError(f"Unsupported Linux architecture: {machine}")
        return f"{NCBI_DATASETS_BASE_URL}/{platform_dir}/datasets", "datasets"

    if system == "darwin":
        # NCBI currently publishes a universal macOS binary.
        return f"{NCBI_DATASETS_BASE_URL}/mac/datasets", "datasets"

    if system == "windows":
        if machine not in {"x86_64", "amd64", "x86-64"}:
            raise RuntimeError(f"Unsupported Windows architecture: {machine}")
        return f"{NCBI_DATASETS_BASE_URL}/win64/datasets.exe", "datasets.exe"

    raise RuntimeError(f"Unsupported operating system: {platform.system()}")


def ensure_datasets(tool_dir: Path, refresh: bool = False) -> Path:
    """Download and cache the NCBI ``datasets`` executable if necessary.

    Parameters
    ----------
    tool_dir
        Directory in which to cache the executable.
    refresh
        If True, replace any cached executable with a fresh download.

    Returns
    -------
    pathlib.Path
        Path to the working ``datasets`` executable.
    """

    url, filename = datasets_download_url()
    tool_dir.mkdir(parents=True, exist_ok=True)
    executable = tool_dir / filename

    if refresh and executable.exists():
        executable.unlink()

    if not executable.exists():
        print(f"Downloading NCBI datasets from {url}")
        temporary = executable.with_suffix(executable.suffix + ".download")
        try:
            with urllib.request.urlopen(url, timeout=120) as response, temporary.open(
                "wb"
            ) as out:
                shutil.copyfileobj(response, out)
        except (urllib.error.URLError, OSError) as exc:
            temporary.unlink(missing_ok=True)
            raise RuntimeError(f"Could not download NCBI datasets: {exc}") from exc
        temporary.replace(executable)

    if os.name != "nt":
        executable.chmod(executable.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)

    # A version invocation is a cheap integrity/executability check and gives
    # future maintainers a useful provenance line in logs.
    completed = subprocess.run(
        [str(executable), "--version"],
        check=False,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Downloaded datasets executable could not run: {completed.stderr.strip()}"
        )
    version = (completed.stdout or completed.stderr).strip()
    print(f"Using {version}")
    return executable


def run_datasets_download(
    datasets: Path, accession: str, zip_path: Path
) -> None:
    """Download one NCBI genome data package.

    Genome FASTA, sequence metadata, and GBFF are requested. GBFF is needed
    specifically because the sequence report identifies chromosomes but does
    not currently contain replicon topology.
    """

    command = [
        str(datasets),
        "download",
        "genome",
        "accession",
        accession,
        "--include",
        "genome,seq-report,gbff",
        "--filename",
        str(zip_path),
        "--no-progressbar",
    ]
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"NCBI datasets download failed: {message}")
    if not zip_path.exists() or zip_path.stat().st_size == 0:
        raise RuntimeError("NCBI datasets reported success but produced no ZIP file")


def read_jsonl(path: Path) -> list[dict]:
    """Read a JSON Lines file into a list of dictionaries."""

    records: list[dict] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSON in {path}:{line_number}: {exc}") from exc
    return records


def find_single_file(directory: Path, pattern: str) -> Path:
    """Return exactly one file matching a glob pattern beneath a directory."""

    matches = list(directory.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one {pattern!r} in {directory}, found {len(matches)}"
        )
    return matches[0]


def nested_get(mapping: dict, *keys: str, default: str = "") -> object:
    """Safely retrieve a value from nested dictionaries."""

    current: object = mapping
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def parse_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    """Yield ``(sequence_id, description, sequence)`` records from FASTA.

    The parser is intentionally small and dependency-free because the validation
    workflow should not require Biopython merely to extract one chromosome.
    """

    sequence_id: str | None = None
    description = ""
    chunks: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                if sequence_id is not None:
                    yield sequence_id, description, "".join(chunks)
                description = line[1:].strip()
                sequence_id = description.split()[0]
                chunks = []
            elif line.strip():
                if sequence_id is None:
                    raise ValueError(f"FASTA sequence encountered before header in {path}")
                chunks.append(line.strip())

    if sequence_id is not None:
        yield sequence_id, description, "".join(chunks)


def select_chromosome_sequence(
    fasta_path: Path, chromosome_record: dict
) -> tuple[str, str, str]:
    """Extract the FASTA record corresponding to an NCBI chromosome report."""

    candidate_ids = {
        str(chromosome_record.get("refseqAccession", "")),
        str(chromosome_record.get("genbankAccession", "")),
    }
    candidate_ids.discard("")

    matches = [record for record in parse_fasta(fasta_path) if record[0] in candidate_ids]
    if len(matches) != 1:
        raise ValueError(
            "Could not uniquely match chromosome sequence report to FASTA; "
            f"candidate accessions={sorted(candidate_ids)}, matches={len(matches)}"
        )
    return matches[0]


def iter_gbff_records(path: Path) -> Iterator[tuple[set[str], bool]]:
    """Yield accession aliases and circular topology from GBFF records.

    This deliberately parses only LOCUS/ACCESSION/VERSION, avoiding an otherwise
    unnecessary Biopython dependency. A record is considered circular only when
    the LOCUS line explicitly contains the ``circular`` topology token.
    """

    aliases: set[str] = set()
    circular = False
    in_record = False

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            if raw_line.startswith("LOCUS"):
                if in_record:
                    yield aliases, circular
                in_record = True
                aliases = set()
                circular = " circular " in f" {raw_line.lower()} "
            elif in_record and raw_line.startswith("ACCESSION"):
                aliases.update(raw_line.split()[1:])
            elif in_record and raw_line.startswith("VERSION"):
                parts = raw_line.split()
                if len(parts) >= 2:
                    aliases.add(parts[1])
            elif in_record and raw_line.startswith("//"):
                yield aliases, circular
                in_record = False
                aliases = set()
                circular = False

    if in_record:
        yield aliases, circular


def gbff_is_circular(gbff_path: Path, candidate_accessions: Iterable[str]) -> bool:
    """Return the circular topology for the requested chromosome GBFF record."""

    candidates = {value for value in candidate_accessions if value}
    matches: list[bool] = []
    for aliases, circular in iter_gbff_records(gbff_path):
        if aliases & candidates:
            matches.append(circular)
    if len(matches) != 1:
        raise ValueError(
            "Could not uniquely identify chromosome in GBFF; "
            f"candidate accessions={sorted(candidates)}, matches={len(matches)}"
        )
    return matches[0]


def write_fasta(path: Path, sequence_id: str, description: str, sequence: str) -> None:
    """Write a FASTA record using deterministic 80-base line wrapping."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{sequence_id} {description}\n" if description else f">{sequence_id}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")


def sequence_sha256(sequence: str) -> str:
    """Return a SHA-256 digest of the uppercase chromosome sequence only."""

    return hashlib.sha256(sequence.upper().encode("ascii")).hexdigest()


def validate_assembly(
    request: GenomeRequest,
    datasets: Path,
    output_dir: Path,
    package_dir: Path,
) -> ValidationResult:
    """Download, extract, and validate one requested assembly.

    Any validation exception is converted into a failed ``ValidationResult`` so
    the caller can continue checking the rest of the manifest.
    """

    result = ValidationResult(
        requested_accession=request.assembly_accession,
        expected_species=request.species,
        expected_strain=request.strain,
    )

    try:
        zip_path = package_dir / f"{request.assembly_accession}.zip"
        unpack_dir = package_dir / request.assembly_accession
        if unpack_dir.exists():
            shutil.rmtree(unpack_dir)

        run_datasets_download(datasets, request.assembly_accession, zip_path)
        with zipfile.ZipFile(zip_path) as archive:
            archive.extractall(unpack_dir)

        data_root = unpack_dir / "ncbi_dataset" / "data"
        assembly_report_path = data_root / "assembly_data_report.jsonl"
        assembly_records = read_jsonl(assembly_report_path)
        if len(assembly_records) != 1:
            raise ValueError(
                f"Expected one assembly metadata record, found {len(assembly_records)}"
            )
        assembly = assembly_records[0]

        observed_accession = str(assembly.get("accession", ""))
        result.assembly_accession = observed_accession
        if observed_accession != request.assembly_accession:
            raise ValueError(
                f"Requested {request.assembly_accession} but NCBI returned "
                f"{observed_accession or '<missing accession>'}"
            )

        result.assembly_level = str(
            nested_get(assembly, "assemblyInfo", "assemblyLevel")
        )
        result.observed_species = str(nested_get(assembly, "organism", "organismName"))
        result.observed_strain = str(
            nested_get(assembly, "organism", "infraspecificNames", "strain")
        )

        if result.assembly_level != "Complete Genome":
            raise ValueError(
                f"Assembly level is {result.assembly_level!r}, not 'Complete Genome'"
            )

        # NCBI organismName normally includes the strain designation, e.g.
        # "Bacteroides thetaiotaomicron VPI-5482", whereas the manifest stores
        # species and strain separately. Require the NCBI name to begin with the
        # exact manifest species name followed by either end-of-string or whitespace.
        if (
            result.observed_species
            and result.observed_species != request.species
            and not result.observed_species.startswith(request.species + " ")
        ):
            raise ValueError(
                f"NCBI organism differs from manifest species: "
                f"{result.observed_species}"
            )

        if result.observed_strain and result.observed_strain != request.strain:
            result.notes.append(
                f"NCBI strain differs from manifest: {result.observed_strain}"
            )

        assembly_dir = data_root / request.assembly_accession
        sequence_report_path = assembly_dir / "sequence_report.jsonl"
        sequence_records = read_jsonl(sequence_report_path)
        chromosomes = [
            record
            for record in sequence_records
            if str(record.get("assignedMoleculeLocationType", "")).lower()
            == "chromosome"
        ]
        result.ignored_nonchromosomal_sequences = len(sequence_records) - len(chromosomes)

        if len(chromosomes) != 1:
            raise ValueError(
                f"Expected exactly one chromosome, found {len(chromosomes)}"
            )
        chromosome = chromosomes[0]

        chromosome_accession = str(
            chromosome.get("refseqAccession")
            or chromosome.get("genbankAccession")
            or ""
        )
        if not chromosome_accession:
            raise ValueError("Chromosome report contains no RefSeq/GenBank accession")
        result.chromosome_accession = chromosome_accession
        result.chromosome_length = int(chromosome["length"])

        genomic_fasta = find_single_file(assembly_dir, "*_genomic.fna")
        sequence_id, description, sequence = select_chromosome_sequence(
            genomic_fasta, chromosome
        )
        sequence = sequence.upper()

        if len(sequence) != result.chromosome_length:
            raise ValueError(
                "Chromosome FASTA length does not match NCBI sequence report: "
                f"{len(sequence)} != {result.chromosome_length}"
            )

        result.n_count = sequence.count("N")
        result.non_acgt_count = sum(base not in VALID_DNA_BASES for base in sequence)
        if result.non_acgt_count:
            raise ValueError(
                "Chromosome contains ambiguous/non-ACGT bases: "
                f"N={result.n_count}, total non-ACGT={result.non_acgt_count}"
            )

        gbff_path = find_single_file(assembly_dir, "*.gbff")
        result.circular = gbff_is_circular(
            gbff_path,
            [
                str(chromosome.get("refseqAccession", "")),
                str(chromosome.get("genbankAccession", "")),
            ],
        )
        if not result.circular:
            raise ValueError("Chromosome is not marked circular in the GBFF LOCUS record")

        chromosome_dir = output_dir / "chromosomes"
        fasta_output = chromosome_dir / f"{request.assembly_accession}.fna"
        # Keep the NCBI sequence accession as the FASTA ID; downstream code can
        # use validated_genomes.tsv to connect it to assembly/species metadata.
        write_fasta(fasta_output, sequence_id, description, sequence)

        result.sha256 = sequence_sha256(sequence)
        result.fasta_path = str(fasta_output)
        result.status = "PASS"
        return result

    except Exception as exc:  # Continue through the full panel to report all failures.
        result.error = str(exc)
        return result


def write_validation_report(path: Path, results: list[ValidationResult]) -> None:
    """Write the complete machine-readable validation report as TSV."""

    rows = [result.as_row() for result in results]
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def print_summary(results: list[ValidationResult]) -> None:
    """Print a compact human-readable validation summary."""

    print()
    print(
        f"{'Assembly':<18} {'Species':<34} {'Chr':>3} "
        f"{'Circular':>8} {'Non-ACGT':>9} {'Length':>11} {'Status':>7}"
    )
    print("-" * 100)
    for result in results:
        chromosome_count = "1" if result.chromosome_accession else "-"
        circular = _fmt_bool(result.circular) or "-"
        non_acgt = _fmt_optional(result.non_acgt_count) or "-"
        length = _fmt_optional(result.chromosome_length) or "-"
        print(
            f"{result.requested_accession:<18} "
            f"{result.expected_species:<34.34} "
            f"{chromosome_count:>3} {circular:>8} {non_acgt:>9} "
            f"{length:>11} {result.status:>7}"
        )
        if result.error:
            print(f"  ERROR: {result.error}")
        for note in result.notes:
            print(f"  NOTE:  {note}")


def main(argv: list[str] | None = None) -> int:
    """Run the genome acquisition and validation workflow."""

    args = parse_args(argv)
    try:
        requests = read_manifest(args.manifest)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    tool_dir = output_dir / ".tools"

    try:
        datasets = ensure_datasets(tool_dir, refresh=args.refresh_datasets)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    # Package downloads are normally disposable: the extracted chromosome FASTA
    # plus validation manifest are the durable outputs used by later stages.
    if args.keep_packages:
        package_dir = output_dir / "packages"
        package_dir.mkdir(parents=True, exist_ok=True)
        temp_context = None
    else:
        temp_context = tempfile.TemporaryDirectory(
            prefix="msamtools_genomes_", dir=output_dir
        )
        package_dir = Path(temp_context.name)

    try:
        results: list[ValidationResult] = []
        for index, request in enumerate(requests, start=1):
            print(
                f"[{index}/{len(requests)}] Preparing "
                f"{request.assembly_accession} ({request.species} {request.strain})"
            )
            result = validate_assembly(request, datasets, output_dir, package_dir)
            results.append(result)
            if result.status == "PASS":
                print("  PASS")
            else:
                print(f"  FAIL: {result.error}")

        report_path = output_dir / "validated_genomes.tsv"
        write_validation_report(report_path, results)
        print_summary(results)
        print(f"\nValidation report: {report_path}")

        failures = sum(result.status != "PASS" for result in results)
        if failures:
            print(
                f"\nFAILED: {failures} of {len(results)} assemblies did not satisfy "
                "the validation criteria.",
                file=sys.stderr,
            )
            return 1

        print(f"\nPASS: all {len(results)} assemblies are validated single circular chromosomes.")
        return 0
    finally:
        if temp_context is not None:
            temp_context.cleanup()


if __name__ == "__main__":
    raise SystemExit(main())
