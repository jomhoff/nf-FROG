#!/usr/bin/env python3

import argparse
import gzip
from collections import Counter


TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    counts = Counter()
    het_ref_fractions = []
    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[0] != args.chromosome:
                continue
            ref = fields[3]
            alts = fields[4].split(",")
            sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
            gt = sample.get("GT", ".").replace("|", "/")
            if "." in gt:
                counts["no_call_variants"] += 1
                continue
            alleles = [int(value) for value in gt.split("/")]
            if len(set(alleles)) > 1:
                counts["heterozygous_variants"] += 1
                ad = sample.get("AD")
                if ad and ad != ".":
                    depths = [int(value) if value not in ("", ".") else 0 for value in ad.split(",")]
                elif sample.get("RO") not in (None, "", "."):
                    depths = [int(sample["RO"])] + [
                        int(value) if value not in ("", ".") else 0
                        for value in sample.get("AO", "").split(",")
                    ]
                else:
                    depths = []
                if depths:
                    total = sum(depths)
                    if total:
                        het_ref_fractions.append(depths[0] / total)
            elif alleles and alleles[0] > 0:
                counts["homozygous_alt_variants"] += 1

            if len(ref) == 1 and len(alts) == 1 and len(alts[0]) == 1:
                counts["snps"] += 1
                if (ref.upper(), alts[0].upper()) in TRANSITIONS:
                    counts["transitions"] += 1
                else:
                    counts["transversions"] += 1
            else:
                counts["indels_or_complex"] += 1

    transitions = counts["transitions"]
    transversions = counts["transversions"]
    titv = transitions / transversions if transversions else float("nan")
    mean_ref_fraction = (
        sum(het_ref_fractions) / len(het_ref_fractions)
        if het_ref_fractions
        else float("nan")
    )
    with open(args.output, "w") as output:
        output.write(
            "chromosome\tsnps\tindels_or_complex\theterozygous_variants\t"
            "homozygous_alt_variants\ttransitions\ttransversions\tti_tv\t"
            "mean_het_reference_allele_fraction\n"
        )
        output.write(
            f"{args.chromosome}\t{counts['snps']}\t{counts['indels_or_complex']}\t"
            f"{counts['heterozygous_variants']}\t{counts['homozygous_alt_variants']}\t"
            f"{transitions}\t{transversions}\t{titv}\t{mean_ref_fraction}\n"
        )


if __name__ == "__main__":
    main()
