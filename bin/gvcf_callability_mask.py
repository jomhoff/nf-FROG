#!/usr/bin/env python3

import argparse
import gzip
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a reason-coded BED mask from GVCF reference-confidence records."
    )
    parser.add_argument("--gvcf", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--minimum-gq", type=int, default=20)
    parser.add_argument("--output", required=True)
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


def chromosome_length(path, chromosome):
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if fields[0] == chromosome:
                return int(fields[1])
    raise SystemExit(f"Chromosome {chromosome} is absent from {path}")


def parse_info(value):
    result = {}
    for item in value.split(";"):
        if "=" in item:
            key, item_value = item.split("=", 1)
            result[key] = item_value
        elif item:
            result[item] = True
    return result


def parse_sample(format_field, sample_field):
    return dict(zip(format_field.split(":"), sample_field.split(":")))


def numeric(value):
    if value in (None, "", "."):
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def merge(intervals):
    merged = []
    for chrom, start, end, reason in sorted(intervals, key=lambda item: (item[0], item[3], item[1], item[2])):
        if (
            merged
            and merged[-1][0] == chrom
            and merged[-1][3] == reason
            and start <= merged[-1][2]
        ):
            merged[-1] = (chrom, merged[-1][1], max(end, merged[-1][2]), reason)
        else:
            merged.append((chrom, start, end, reason))
    return merged


def main():
    args = parse_args()
    min_depth, max_depth = load_thresholds(args.thresholds, args.chromosome)
    chrom_length = chromosome_length(args.fai, args.chromosome)

    masks = []
    covered = []
    counts = Counter()
    with open_text(args.gvcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[0] != args.chromosome:
                continue
            position = int(fields[1])
            ref = fields[3]
            info = parse_info(fields[7])
            sample = parse_sample(fields[8], fields[9])
            end_1_based = int(info.get("END", position + len(ref) - 1))
            start = position - 1
            end = min(chrom_length, end_1_based)
            covered.append((start, end))

            depth = numeric(sample.get("MIN_DP"))
            if depth is None:
                depth = numeric(sample.get("DP"))
            gq = numeric(sample.get("GQ"))
            gt = sample.get("GT")

            reasons = []
            if gt is None or "." in gt:
                reasons.append("GVCF_NO_CALL")
            if depth is None:
                reasons.append("GVCF_NO_DEPTH")
            elif depth < min_depth:
                reasons.append("LOW_DEPTH")
            elif depth > max_depth:
                reasons.append("HIGH_DEPTH")
            if args.minimum_gq > 0 and (gq is None or gq < args.minimum_gq):
                reasons.append("LOW_GQ")

            for reason in sorted(set(reasons)):
                masks.append((args.chromosome, start, end, reason))
                counts[reason] += end - start

    cursor = 0
    for start, end in sorted(covered):
        if start > cursor:
            masks.append((args.chromosome, cursor, start, "GVCF_ABSENT"))
            counts["GVCF_ABSENT"] += start - cursor
        cursor = max(cursor, end)
    if cursor < chrom_length:
        masks.append((args.chromosome, cursor, chrom_length, "GVCF_ABSENT"))
        counts["GVCF_ABSENT"] += chrom_length - cursor

    with open(args.output, "w") as output:
        for row in merge(masks):
            output.write("\t".join(map(str, row)) + "\n")

    with open(args.stats, "w") as output:
        output.write("chromosome\treason\tmasked_bases\n")
        for reason, bases in sorted(counts.items()):
            output.write(f"{args.chromosome}\t{reason}\t{bases}\n")


if __name__ == "__main__":
    main()
