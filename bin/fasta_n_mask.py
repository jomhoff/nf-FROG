#!/usr/bin/env python3

import argparse


def fasta_records(path):
    name = None
    sequence = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence)
                name = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line.strip())
    if name is not None:
        yield name, "".join(sequence)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--chromosome")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.output, "w") as output:
        for chrom, sequence in fasta_records(args.fasta):
            if args.chromosome and chrom != args.chromosome:
                continue
            start = None
            for index, base in enumerate(sequence.upper()):
                is_masked = base not in "ACGT"
                if is_masked and start is None:
                    start = index
                elif not is_masked and start is not None:
                    output.write(f"{chrom}\t{start}\t{index}\tREFERENCE_AMBIGUOUS\n")
                    start = None
            if start is not None:
                output.write(f"{chrom}\t{start}\t{len(sequence)}\tREFERENCE_AMBIGUOUS\n")


if __name__ == "__main__":
    main()
