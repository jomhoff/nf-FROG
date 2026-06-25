#!/usr/bin/env python3

import argparse
import gzip
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter a single-sample VCF, emit a consensus VCF, and mask every rejected call."
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--het-mode", choices=("retain", "major", "iupac", "mask"), required=True)
    parser.add_argument("--minimum-qual", type=float, default=30)
    parser.add_argument("--minimum-mqm", type=float, default=20.0)
    parser.add_argument("--minimum-alt-mean-quality", type=float, default=20.0)
    parser.add_argument("--balanced-alt-count", type=int, default=4)
    parser.add_argument("--minimum-alt-strand-observations", type=int, default=1)
    parser.add_argument("--minimum-alt-placement-observations", type=int, default=1)
    parser.add_argument("--minimum-het-allele-fraction", type=float, default=0.2)
    parser.add_argument("--maximum-het-major-fraction", type=float, default=0.8)
    parser.add_argument("--minimum-hom-alt-fraction", type=float, default=0.9)
    parser.add_argument("--minimum-allele-depth", type=int, default=3)
    parser.add_argument("--include-indels", action="store_true")
    parser.add_argument("--exclude-bed")
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--reject-bed", required=True)
    parser.add_argument("--stats", required=True)
    return parser.parse_args()


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_thresholds(path, chromosome):
    with open(path) as handle:
        header = handle.readline().rstrip().split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip().split("\t")))
            if row["chromosome"] == chromosome:
                return int(row["min_depth"]), int(row["max_depth"])
    raise SystemExit(f"No depth thresholds found for {chromosome}")


def number(value, kind=float):
    if value in (None, "", "."):
        return None
    try:
        return kind(value)
    except ValueError:
        return None


def parse_info(value):
    info = {}
    for item in value.split(";"):
        if "=" in item:
            key, item_value = item.split("=", 1)
            info[key] = item_value
    return info


def parse_sample(format_field, sample_field):
    keys = format_field.split(":")
    values = sample_field.split(":")
    return keys, values, dict(zip(keys, values))


def genotype_alleles(gt):
    return [int(value) for value in gt.replace("|", "/").split("/")]


def resolve_major_genotype(gt, depths):
    alleles = genotype_alleles(gt)
    candidates = sorted(set(alleles))
    winner = max(
        candidates,
        key=lambda allele: (depths[allele] if allele < len(depths) else 0, -allele),
    )
    return "/".join([str(winner)] * len(alleles))


def numeric_list(value, kind=float):
    if value in (None, "", "."):
        return []
    return [number(item, kind) for item in value.split(",")]


def allele_depths(sample):
    if sample.get("AD") not in (None, "", "."):
        return [number(value, int) or 0 for value in sample["AD"].split(",")]
    reference_depth = number(sample.get("RO"), int)
    alternate_depths = numeric_list(sample.get("AO"), int)
    if reference_depth is None and not alternate_depths:
        return []
    return [reference_depth or 0] + [value or 0 for value in alternate_depths]


def called_alt_indexes(alleles):
    return sorted({allele - 1 for allele in alleles if allele > 0})


def add_freebayes_evidence_failures(reasons, info, alt_indexes, args):
    mqms = numeric_list(info.get("MQM"))
    alt_counts = numeric_list(info.get("AO"), int)
    alt_quality_sums = numeric_list(info.get("QA"))
    strand_forward = numeric_list(info.get("SAF"), int)
    strand_reverse = numeric_list(info.get("SAR"), int)
    placement_left = numeric_list(info.get("RPL"), int)
    placement_right = numeric_list(info.get("RPR"), int)

    for index in alt_indexes:
        mqm = mqms[index] if index < len(mqms) else None
        if mqm is None or mqm < args.minimum_mqm:
            reasons.append("LOW_ALT_MAPPING_QUALITY")

        count = alt_counts[index] if index < len(alt_counts) else None
        quality_sum = alt_quality_sums[index] if index < len(alt_quality_sums) else None
        if count is None or count <= 0 or quality_sum is None:
            reasons.append("MISSING_ALT_QUALITY")
        elif quality_sum / count < args.minimum_alt_mean_quality:
            reasons.append("LOW_ALT_BASE_QUALITY")

        if count is not None and count >= args.balanced_alt_count:
            forward = strand_forward[index] if index < len(strand_forward) else None
            reverse = strand_reverse[index] if index < len(strand_reverse) else None
            if (
                forward is None
                or reverse is None
                or min(forward, reverse) < args.minimum_alt_strand_observations
            ):
                reasons.append("ALT_STRAND_BIAS")

            left = placement_left[index] if index < len(placement_left) else None
            right = placement_right[index] if index < len(placement_right) else None
            if (
                left is None
                or right is None
                or min(left, right) < args.minimum_alt_placement_observations
            ):
                reasons.append("ALT_READ_POSITION_BIAS")


def load_exclusions(path, chromosome):
    intervals = []
    if not path:
        return intervals
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[0] == chromosome:
                intervals.append((int(fields[1]), int(fields[2]), fields[3] if len(fields) > 3 else "EXCLUDED_INTERVAL"))
    return sorted(intervals)


def main():
    args = parse_args()
    min_depth, max_depth = load_thresholds(args.thresholds, args.chromosome)
    exclusions = load_exclusions(args.exclude_bed, args.chromosome)
    counts = Counter()
    rejected = []

    with open_text(args.vcf) as source, open(args.output_vcf, "w") as output:
        for line in source:
            if line.startswith("#"):
                output.write(line)
                continue

            fields = line.rstrip().split("\t")
            if fields[0] != args.chromosome:
                continue

            position = int(fields[1])
            ref = fields[3]
            alts = fields[4].split(",")
            if fields[4] in (".", "<NON_REF>"):
                continue

            reasons = []
            variant_start = position - 1
            variant_end = variant_start + len(ref)
            for start, end, reason in exclusions:
                if start >= variant_end:
                    break
                if end > variant_start and start < variant_end:
                    reasons.append(reason)
            symbolic = any(alt.startswith("<") or alt == "*" for alt in alts)
            is_snp = len(ref) == 1 and all(len(alt) == 1 for alt in alts) and not symbolic
            is_indel = not is_snp
            if symbolic:
                reasons.append("SYMBOLIC_ALLELE")
            if is_indel and not args.include_indels:
                reasons.append("INDEL_EXCLUDED")

            qual = number(fields[5])
            if qual is None or qual < args.minimum_qual:
                reasons.append("LOW_QUAL")
            if fields[6] not in (".", "PASS"):
                reasons.append("SITE_FILTER")

            info = parse_info(fields[7])

            format_keys, sample_values, sample = parse_sample(fields[8], fields[9])
            gt = sample.get("GT")
            dp = number(sample.get("DP"), int)
            if gt is None or "." in gt:
                reasons.append("GENOTYPE_NO_CALL")
            if dp is None:
                reasons.append("GENOTYPE_NO_DEPTH")
            elif dp < min_depth:
                reasons.append("LOW_DEPTH")
            elif dp > max_depth:
                reasons.append("HIGH_DEPTH")
            alleles = [] if gt is None or "." in gt else genotype_alleles(gt)
            depths = allele_depths(sample)
            add_freebayes_evidence_failures(
                reasons, info, called_alt_indexes(alleles), args
            )
            if alleles and not depths:
                reasons.append("MISSING_ALLELE_DEPTH")
            elif alleles and depths:
                total_ad = sum(depths)
                unique_alleles = sorted(set(alleles))
                if total_ad == 0:
                    reasons.append("ZERO_ALLELE_DEPTH")
                elif len(unique_alleles) > 1:
                    called_fractions = [
                        depths[allele] / total_ad if allele < len(depths) else 0
                        for allele in unique_alleles
                    ]
                    if min(called_fractions) < args.minimum_het_allele_fraction:
                        reasons.append("HET_ALLELE_BALANCE")
                    if max(called_fractions) > args.maximum_het_major_fraction:
                        reasons.append("HET_ALLELE_BALANCE")
                    if min(depths[allele] if allele < len(depths) else 0 for allele in unique_alleles) < args.minimum_allele_depth:
                        reasons.append("LOW_ALLELE_DEPTH")
                elif unique_alleles[0] > 0:
                    allele = unique_alleles[0]
                    fraction = depths[allele] / total_ad if allele < len(depths) else 0
                    if fraction < args.minimum_hom_alt_fraction:
                        reasons.append("HOM_ALT_ALLELE_BALANCE")
                    if (depths[allele] if allele < len(depths) else 0) < args.minimum_allele_depth:
                        reasons.append("LOW_ALLELE_DEPTH")

            if reasons:
                for reason in sorted(set(reasons)):
                    rejected.append((fields[0], position - 1, position - 1 + len(ref), reason))
                    counts[f"REJECT_{reason}"] += 1
                continue

            if not alleles or all(allele == 0 for allele in alleles):
                counts["PASS_HOM_REF"] += 1
                continue

            is_het = len(set(alleles)) > 1
            if is_het and args.het_mode == "retain":
                counts["PASS_HET_RETAINED"] += 1
                continue
            if is_het and args.het_mode == "mask":
                rejected.append((fields[0], position - 1, position - 1 + len(ref), "HET_MASKED"))
                counts["REJECT_HET_MASKED"] += 1
                continue
            if is_het and args.het_mode == "major":
                resolved = resolve_major_genotype(gt, depths)
                sample_values[format_keys.index("GT")] = resolved
                if all(allele == 0 for allele in genotype_alleles(resolved)):
                    counts["PASS_HET_MAJOR_REF"] += 1
                    continue
                counts["PASS_HET_MAJOR_ALT"] += 1
            elif is_het:
                counts["PASS_HET_IUPAC"] += 1
            else:
                counts["PASS_HOM_ALT"] += 1

            fields[6] = "PASS"
            fields[9] = ":".join(sample_values)
            output.write("\t".join(fields) + "\n")

    with open(args.reject_bed, "w") as output:
        for row in sorted(set(rejected)):
            output.write("\t".join(map(str, row)) + "\n")

    with open(args.stats, "w") as output:
        output.write("chromosome\tmetric\tcount\n")
        for metric, count in sorted(counts.items()):
            output.write(f"{args.chromosome}\t{metric}\t{count}\n")


if __name__ == "__main__":
    main()
