#!/usr/bin/env python3

import argparse
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert per-base SAMtools depth into a reason-coded BED mask."
    )
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--stats", required=True)
    return parser.parse_args()


def load_thresholds(path, chromosome):
    with open(path) as handle:
        header = handle.readline().rstrip().split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip().split("\t")))
            if row["chromosome"] == chromosome:
                return int(row["min_depth"]), int(row["max_depth"])
    raise SystemExit(f"No depth thresholds found for {chromosome}")


def main():
    args = parse_args()
    minimum, maximum = load_thresholds(args.thresholds, args.chromosome)
    counts = Counter()
    current = None

    with open(args.output, "w") as output:
        for line in __import__("sys").stdin:
            chromosome, position_text, depth_text = line.rstrip().split("\t")[:3]
            if chromosome != args.chromosome:
                continue
            position = int(position_text)
            depth = int(depth_text)
            reason = "LOW_DEPTH" if depth < minimum else "HIGH_DEPTH" if depth > maximum else None
            counts["TOTAL_BASES"] += 1
            if reason:
                counts[reason] += 1

            start = position - 1
            if current and (reason != current[3] or start != current[2]):
                output.write("\t".join(map(str, current)) + "\n")
                current = None
            if reason:
                if current:
                    current[2] = position
                else:
                    current = [chromosome, start, position, reason]

        if current:
            output.write("\t".join(map(str, current)) + "\n")

    with open(args.stats, "w") as output:
        output.write("chromosome\treason\tbases\n")
        for reason in ("LOW_DEPTH", "HIGH_DEPTH", "TOTAL_BASES"):
            output.write(f"{args.chromosome}\t{reason}\t{counts[reason]}\n")


if __name__ == "__main__":
    main()
