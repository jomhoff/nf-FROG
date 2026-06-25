#!/usr/bin/env python3

import argparse


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
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--status", required=True)
    args = parser.parse_args()

    previous = read_fasta(args.previous)
    current = read_fasta(args.current)
    new = read_fasta(args.new)
    if current.keys() != new.keys() or current.keys() != previous.keys():
        raise SystemExit("Reference contig sets changed between iterations")

    changes = 0
    acgt_changes = 0
    oscillations = 0
    for chrom in current:
        if len(current[chrom]) != len(new[chrom]) or len(current[chrom]) != len(previous[chrom]):
            raise SystemExit(f"Reference length changed for {chrom}; iterative references must remain SNP-only")
        for old_base, new_base, previous_base in zip(current[chrom], new[chrom], previous[chrom]):
            if old_base != new_base:
                changes += 1
                if old_base in "ACGT" and new_base in "ACGT":
                    acgt_changes += 1
                if new_base == previous_base:
                    oscillations += 1

    converged = (
        args.round >= args.minimum_rounds
        and changes <= args.maximum_changes
        and oscillations == 0
    )
    with open(args.metrics, "w") as output:
        output.write("round\tchanged_bases\tchanged_acgt_bases\toscillating_bases\tconverged\n")
        output.write(f"{args.round}\t{changes}\t{acgt_changes}\t{oscillations}\t{str(converged).lower()}\n")
    with open(args.status, "w") as output:
        output.write("converged\n" if converged else "active\n")


if __name__ == "__main__":
    main()
