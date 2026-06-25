#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


IUPAC = {
    frozenset(("A", "G")): "R",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S",
    frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K",
    frozenset(("A", "C")): "M",
}


def read_fasta(path):
    records = {}
    name = None
    sequence = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(sequence).upper()
                name = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line.strip())
    if name is not None:
        records[name] = "".join(sequence).upper()
    return records


def interval_is_masked(sequence, start, end):
    return all(base == "N" for base in sequence[start:end])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_directory = Path(args.input)
    output_directory = Path(args.output)
    expected = json.loads((input_directory / "expected.json").read_text())

    consensus_path = (
        output_directory / "03.consensus" / "frog" / "frog.consensus.fa"
    )
    convergence_path = output_directory / "07.qc" / "all_samples.convergence.tsv"
    mask_path = (
        output_directory
        / "03.consensus"
        / "frog"
        / "frog.noncallable_mask.reasons.bed"
    )

    missing = [
        path
        for path in (consensus_path, convergence_path, mask_path)
        if not path.exists()
    ]
    if missing:
        raise SystemExit(
            "Missing expected outputs:\n" + "\n".join(str(path) for path in missing)
        )

    consensus = read_fasta(consensus_path)
    failures = []

    for position_text, variant in expected["homozygous"].items():
        position = int(position_text)
        observed = consensus["chr1"][position - 1]
        if observed != variant["alt"]:
            failures.append(
                f"chr1:{position} homozygous ALT expected {variant['alt']}, observed {observed}"
            )

    for position_text, variant in expected["heterozygous"].items():
        position = int(position_text)
        expected_code = IUPAC[frozenset((variant["ref"], variant["alt"]))]
        observed = consensus["chr1"][position - 1]
        if observed != expected_code:
            failures.append(
                f"chr1:{position} IUPAC expected {expected_code}, observed {observed}"
            )

    gap_start, gap_end = expected["coverage_gap_bed"]
    if not interval_is_masked(consensus["chr1"], gap_start, gap_end):
        failures.append(
            f"Coverage desert chr1:{gap_start}-{gap_end} was not completely masked"
        )

    zero_ploidy_chromosome = expected["ploidy_zero_chromosome"]
    if not interval_is_masked(
        consensus[zero_ploidy_chromosome],
        0,
        len(consensus[zero_ploidy_chromosome]),
    ):
        failures.append(f"{zero_ploidy_chromosome} was not completely masked")

    convergence_lines = convergence_path.read_text().splitlines()
    if not any(
        line.startswith("frog\t2\t") and line.endswith("\ttrue")
        for line in convergence_lines[1:]
    ):
        failures.append("Sample did not report convergence in round 2")

    mask_text = mask_path.read_text()
    for reason in ("LOW_DEPTH", "PLOIDY_ZERO"):
        if reason not in mask_text:
            failures.append(f"Final reason-coded mask is missing {reason}")

    if failures:
        raise SystemExit("Synthetic validation failed:\n- " + "\n- ".join(failures))

    print("Synthetic nf-FROG validation: PASS")


if __name__ == "__main__":
    main()
