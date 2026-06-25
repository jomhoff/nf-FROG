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
    parser.add_argument("--minimum-gq", type=int, default=20)
    parser.add_argument("--minimum-qual", type=float, default=30)
    parser.add_argument("--minimum-qd", type=float, default=2.0)
    parser.add_argument("--minimum-mq", type=float, default=40.0)
    parser.add_argument("--maximum-snp-fs", type=float, default=60.0)
    parser.add_argument("--maximum-indel-fs", type=float, default=200.0)
    parser.add_argument("--maximum-snp-sor", type=float, default=3.0)
    parser.add_argument("--maximum-indel-sor", type=float, default=10.0)
    parser.add_argument("--minimum-read-pos-rank-sum", type=float, default=-8.0)
    parser.add_argument("--minimum-mq-rank-sum", type=float, default=-12.5)
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


def add_threshold_failure(reasons, info, key, threshold, comparison, reason):
    value = number(info.get(key))
    if value is not None and comparison(value, threshold):
        reasons.append(reason)


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
            add_threshold_failure(reasons, info, "QD", args.minimum_qd, lambda value, limit: value < limit, "LOW_QD")
            add_threshold_failure(reasons, info, "MQ", args.minimum_mq, lambda value, limit: value < limit, "LOW_MQ")
            add_threshold_failure(
                reasons,
                info,
                "FS",
                args.maximum_snp_fs if is_snp else args.maximum_indel_fs,
                lambda value, limit: value > limit,
                "STRAND_BIAS_FS",
            )
            add_threshold_failure(
                reasons,
                info,
                "SOR",
                args.maximum_snp_sor if is_snp else args.maximum_indel_sor,
                lambda value, limit: value > limit,
                "STRAND_BIAS_SOR",
            )
            add_threshold_failure(
                reasons,
                info,
                "ReadPosRankSum",
                args.minimum_read_pos_rank_sum,
                lambda value, limit: value < limit,
                "READ_POSITION_BIAS",
            )
            add_threshold_failure(
                reasons,
                info,
                "MQRankSum",
                args.minimum_mq_rank_sum,
                lambda value, limit: value < limit,
                "MAPPING_QUALITY_BIAS",
            )

            format_keys, sample_values, sample = parse_sample(fields[8], fields[9])
            gt = sample.get("GT")
            dp = number(sample.get("DP"), int)
            gq = number(sample.get("GQ"), int)
            if gt is None or "." in gt:
                reasons.append("GENOTYPE_NO_CALL")
            if dp is None:
                reasons.append("GENOTYPE_NO_DEPTH")
            elif dp < min_depth:
                reasons.append("LOW_DEPTH")
            elif dp > max_depth:
                reasons.append("HIGH_DEPTH")
            if gq is None or gq < args.minimum_gq:
                reasons.append("LOW_GQ")

            alleles = [] if gt is None or "." in gt else genotype_alleles(gt)
            depths = []
            if sample.get("AD") not in (None, "", "."):
                depths = [number(value, int) or 0 for value in sample["AD"].split(",")]
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
