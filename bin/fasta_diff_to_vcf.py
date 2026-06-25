#!/usr/bin/env python3

import argparse


IUPAC = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"},
    "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--stats", required=True)
    args = parser.parse_args()

    reference = read_fasta(args.reference)
    consensus = read_fasta(args.consensus)
    if reference.keys() != consensus.keys():
        raise SystemExit("Reference and consensus contig sets differ")

    variant_count = 0
    ambiguous_count = 0
    masked_count = 0
    with open(args.output, "w") as output:
        output.write("##fileformat=VCFv4.2\n")
        for chrom, sequence in reference.items():
            output.write(f"##contig=<ID={chrom},length={len(sequence)}>\n")
        output.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype inferred from consensus sequence">\n')
        output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + args.sample + "\n")

        for chrom, ref_sequence in reference.items():
            cons_sequence = consensus[chrom]
            if len(ref_sequence) != len(cons_sequence):
                raise SystemExit(f"Length changed for {chrom}; original-coordinate VCF requires SNP-only consensus")
            for index, (ref_base, cons_code) in enumerate(zip(ref_sequence, cons_sequence), start=1):
                if cons_code == "N" or cons_code not in IUPAC or ref_base not in "ACGT":
                    masked_count += 1
                    continue
                alleles = sorted(IUPAC[cons_code])
                if alleles == [ref_base]:
                    continue
                alt_alleles = [base for base in alleles if base != ref_base]
                if not alt_alleles:
                    continue
                allele_index = {ref_base: 0}
                for alt_index, alt in enumerate(alt_alleles, start=1):
                    allele_index[alt] = alt_index
                genotype = "/".join(str(allele_index[base]) for base in alleles)
                if len(alleles) > 1:
                    ambiguous_count += 1
                variant_count += 1
                output.write(
                    f"{chrom}\t{index}\t.\t{ref_base}\t{','.join(alt_alleles)}\t.\tPASS\t.\tGT\t{genotype}\n"
                )

    with open(args.stats, "w") as output:
        output.write("sample\tvariants\theterozygous_iupac_sites\tmasked_or_ambiguous_sites\n")
        output.write(f"{args.sample}\t{variant_count}\t{ambiguous_count}\t{masked_count}\n")


if __name__ == "__main__":
    main()
