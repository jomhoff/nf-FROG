#!/usr/bin/env python3

import argparse
import json
from decimal import Decimal
from pathlib import Path


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--previous", required=True)
    parser.add_argument("--current", required=True)
    parser.add_argument("--new", required=True)
    parser.add_argument("--round", type=int, required=True)
    parser.add_argument("--minimum-rounds", type=int, default=2)
    parser.add_argument("--maximum-changes", type=int, default=0)
    parser.add_argument("--maximum-change-fraction", type=str)
    parser.add_argument("--policy-file", type=Path)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--status", required=True)
    args = parser.parse_args()

    # An explicit work-root policy allows an already-running Nextflow workflow
    # to adopt an approved rule at its next comparison, without restarting it.
    # Only discover it in the standard work/xx/hash task layout.
    policy_path = args.policy_file
    cwd = Path.cwd()
    if policy_path is None and (cwd / ".command.sh").is_file():
        candidate = cwd.parent.parent / ".frog-convergence-policy.json"
        if candidate.is_file():
            policy_path = candidate
    policy = {}
    if policy_path is not None:
        policy = json.loads(policy_path.read_text())
        if set(policy) - {"maximum_change_fraction", "consecutive_rounds", "description"}:
            raise SystemExit("Unrecognized convergence policy fields")
        if policy.get("consecutive_rounds", 2) != 2:
            raise SystemExit("Practical convergence requires exactly two consecutive rounds")
    fraction_text = args.maximum_change_fraction
    if fraction_text is None:
        fraction_text = policy.get("maximum_change_fraction")
    fraction = Decimal(str(fraction_text)) if fraction_text is not None else None
    if fraction is not None and (not fraction.is_finite() or not 0 < fraction <= 1):
        raise SystemExit("Maximum change fraction must be in (0, 1]")
    if args.maximum_changes < 0:
        raise SystemExit("Maximum changes must be nonnegative")

    previous = read_fasta(args.previous)
    current = read_fasta(args.current)
    new = read_fasta(args.new)
    if current.keys() != new.keys() or current.keys() != previous.keys():
        raise SystemExit("Reference contig sets changed between iterations")

    changes = 0
    acgt_changes = 0
    oscillations = 0
    previous_changes = 0
    total_bases = sum(len(sequence) for sequence in current.values())
    if not total_bases:
        raise SystemExit("Reference contains no bases")
    for chrom in current:
        if len(current[chrom]) != len(new[chrom]) or len(current[chrom]) != len(previous[chrom]):
            raise SystemExit(f"Reference length changed for {chrom}; iterative references must remain SNP-only")
        for old_base, new_base, previous_base in zip(current[chrom], new[chrom], previous[chrom]):
            previous_changes += old_base != previous_base
            if old_base != new_base:
                changes += 1
                if old_base in "ACGT" and new_base in "ACGT":
                    acgt_changes += 1
                if new_base == previous_base:
                    oscillations += 1

    exact_converged = args.round >= args.minimum_rounds and changes == 0 and oscillations == 0
    converged = (
        args.round >= args.minimum_rounds
        and changes <= args.maximum_changes
        and oscillations == 0
    )
    mode = "absolute_changes"
    if fraction is not None:
        mode = "practical"
        limit = fraction * total_bases
        converged = (
            args.round >= max(2, args.minimum_rounds)
            and Decimal(changes) < limit
            and Decimal(previous_changes) < limit
        )
    with open(args.metrics, "w") as output:
        output.write("round\tchanged_bases\tchanged_acgt_bases\toscillating_bases\tconverged\n")
        output.write(f"{args.round}\t{changes}\t{acgt_changes}\t{oscillations}\t{str(converged).lower()}\n")
    with open(args.status, "w") as output:
        output.write("converged\n" if converged else "active\n")
    # Keep the legacy TSV columns stable for running workflows and consumers.
    details_path = Path(args.metrics).with_suffix(".details.json")
    details_path.write_text(json.dumps({
        "round": args.round, "mode": mode, "converged": converged,
        "exact_zero_change_convergence": exact_converged,
        "practical_convergence": converged if mode == "practical" else None,
        "reference_bases": total_bases, "changed_bases": changes,
        "previous_round_changed_bases": previous_changes,
        "changed_fraction": changes / total_bases,
        "previous_round_changed_fraction": previous_changes / total_bases,
        "oscillating_bases": oscillations,
        "maximum_change_fraction": str(fraction) if fraction is not None else None,
        "required_consecutive_rounds": 2 if mode == "practical" else 1,
        "policy_file": str(policy_path.resolve()) if policy_path is not None else None,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
