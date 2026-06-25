#!/usr/bin/env python3

import argparse
import json
import random
from pathlib import Path


DNA = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def different_base(base):
    return DNA[(DNA.index(base) + 1) % len(DNA)]


def write_fasta(path, records, width=80):
    with open(path, "w") as output:
        for name, sequence in records.items():
            output.write(f">{name}\n")
            for start in range(0, len(sequence), width):
                output.write(sequence[start : start + width] + "\n")


def write_fastq_pair(path1, path2, haplotypes, read_length, insert_size, step, gap):
    quality = "I" * read_length
    count = 0
    with open(path1, "w") as read1_output, open(path2, "w") as read2_output:
        chromosome = "chr1"
        chromosome_length = len(haplotypes[0][chromosome])
        for start in range(0, chromosome_length - insert_size + 1, step):
            fragment_end = start + insert_size
            if start < gap[1] and fragment_end > gap[0]:
                continue
            haplotype_index = count % 2
            sequence = haplotypes[haplotype_index][chromosome]
            read1 = sequence[start : start + read_length]
            read2 = reverse_complement(
                sequence[fragment_end - read_length : fragment_end]
            )
            name = f"synthetic_{count:06d}_{chromosome}_{start + 1}_hap{haplotype_index + 1}"
            read1_output.write(f"@{name}/1\n{read1}\n+\n{quality}\n")
            read2_output.write(f"@{name}/2\n{read2}\n+\n{quality}\n")
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--seed", type=int, default=1731)
    args = parser.parse_args()

    output_directory = Path(args.output).resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    reference = {
        "chr1": "".join(rng.choice(DNA) for _ in range(20_000)),
        "chrZ": "".join(rng.choice(DNA) for _ in range(10_000)),
    }
    haplotype1 = dict(reference)
    haplotype2 = dict(reference)

    homozygous_positions = [4_000, 12_000]
    heterozygous_position = 7_000
    expected = {"homozygous": {}, "heterozygous": {}}

    for position in homozygous_positions:
        ref = reference["chr1"][position - 1]
        alt = different_base(ref)
        for haplotype in (haplotype1, haplotype2):
            sequence = list(haplotype["chr1"])
            sequence[position - 1] = alt
            haplotype["chr1"] = "".join(sequence)
        expected["homozygous"][str(position)] = {"ref": ref, "alt": alt}

    ref = reference["chr1"][heterozygous_position - 1]
    alt = different_base(ref)
    sequence = list(haplotype2["chr1"])
    sequence[heterozygous_position - 1] = alt
    haplotype2["chr1"] = "".join(sequence)
    expected["heterozygous"][str(heterozygous_position)] = {
        "ref": ref,
        "alt": alt,
    }

    coverage_gap = [8_800, 9_700]
    expected["coverage_gap_bed"] = coverage_gap
    expected["ploidy_zero_chromosome"] = "chrZ"

    reference_path = output_directory / "reference.fa"
    write_fasta(reference_path, reference)
    write_fasta(
        output_directory / "truth_haplotype1.fa",
        haplotype1,
    )
    write_fasta(
        output_directory / "truth_haplotype2.fa",
        haplotype2,
    )

    read1_path = output_directory / "frog_R1.fastq"
    read2_path = output_directory / "frog_R2.fastq"
    read_count = write_fastq_pair(
        read1_path,
        read2_path,
        [haplotype1, haplotype2],
        read_length=150,
        insert_size=350,
        step=10,
        gap=coverage_gap,
    )
    expected["read_pairs"] = read_count

    (output_directory / "chromosomes.txt").write_text("chr1\nchrZ\n")
    (output_directory / "samplesheet.tsv").write_text(
        "individual\tread1\tread2\tread_group\tlibrary\tkaryotype\t"
        "default_ploidy\tpcr_free\n"
        f"frog\t{read1_path}\t{read2_path}\tfrog_rg1\tfrog_lib\tXY\t2\ttrue\n"
    )
    (output_directory / "ploidy.tsv").write_text(
        "karyotype\tchromosome\tstart\tend\tploidy\n"
        "XY\tchrZ\t1\t*\t0\n"
    )
    (output_directory / "expected.json").write_text(
        json.dumps(expected, indent=2, sort_keys=True) + "\n"
    )

    print(f"Generated {read_count} read pairs in {output_directory}")


if __name__ == "__main__":
    main()
