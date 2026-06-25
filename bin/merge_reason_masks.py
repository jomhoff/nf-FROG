#!/usr/bin/env python3

import argparse
import re
from collections import Counter, defaultdict


def sanitize(reason):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", reason)


def merge_intervals(intervals):
    merged = []
    for chrom, start, end in sorted(intervals):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((chrom, start, end))
    return merged


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mask", action="append", default=[], help="LABEL=BED")
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    by_reason = defaultdict(list)
    for specification in args.mask:
        label, path = specification.split("=", 1)
        if not path:
            continue
        with open(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) < 3:
                    continue
                if fields[0] != args.chromosome:
                    continue
                reason = fields[3] if len(fields) >= 4 and fields[3] else label
                by_reason[reason].append((fields[0], int(fields[1]), int(fields[2])))

    events = defaultdict(lambda: defaultdict(lambda: {"add": Counter(), "remove": Counter()}))
    for reason, intervals in by_reason.items():
        merged = merge_intervals(intervals)
        by_reason[reason] = merged
        with open(f"{args.prefix}.reason_{sanitize(reason)}.bed", "w") as output:
            for row in merged:
                output.write("\t".join(map(str, row)) + "\n")
        for chrom, start, end in merged:
            events[chrom][start]["add"][reason] += 1
            events[chrom][end]["remove"][reason] += 1

    if not by_reason:
        open(f"{args.prefix}.reason_NONE.bed", "w").close()

    reason_rows = []
    merged_rows = []
    for chrom, chrom_events in events.items():
        active = Counter()
        positions = sorted(chrom_events)
        for index, position in enumerate(positions):
            for reason, count in chrom_events[position]["remove"].items():
                active[reason] -= count
                if active[reason] <= 0:
                    del active[reason]
            for reason, count in chrom_events[position]["add"].items():
                active[reason] += count
            if index + 1 == len(positions):
                continue
            end = positions[index + 1]
            if active and position < end:
                reason_rows.append((chrom, position, end, ",".join(sorted(active))))
                merged_rows.append((chrom, position, end))

    with open(f"{args.prefix}.reasons.bed", "w") as output:
        for row in reason_rows:
            output.write("\t".join(map(str, row)) + "\n")
    with open(f"{args.prefix}.bed", "w") as output:
        for row in merge_intervals(merged_rows):
            output.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
