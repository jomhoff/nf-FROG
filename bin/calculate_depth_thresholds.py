#!/usr/bin/env python3

import argparse
import math
import sys
from collections import Counter, defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate robust chromosome-specific depth thresholds from samtools depth."
    )
    parser.add_argument("--fai", required=True)
    parser.add_argument("--min-multiplier", type=float, default=0.5)
    parser.add_argument("--max-multiplier", type=float, default=2.0)
    parser.add_argument("--mad-multiplier", type=float, default=6.0)
    parser.add_argument("--upper-quantile", type=float, default=0.99)
    parser.add_argument("--absolute-minimum", type=int, default=5)
    parser.add_argument("--minimum-positive-sites", type=int, default=1000)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def weighted_quantile(histogram, quantile):
    total = sum(histogram.values())
    if total == 0:
        return 0
    target = max(1, math.ceil(total * quantile))
    cumulative = 0
    for value in sorted(histogram):
        cumulative += histogram[value]
        if cumulative >= target:
            return value
    return max(histogram)


def median_absolute_deviation(histogram, median):
    deviations = Counter()
    for depth, count in histogram.items():
        deviations[abs(depth - median)] += count
    return weighted_quantile(deviations, 0.5)


def summarize(histogram, absolute_minimum, min_multiplier, max_multiplier, mad_multiplier, upper_quantile):
    positive_sites = sum(histogram.values())
    if positive_sites == 0:
        return None
    total_depth = sum(depth * count for depth, count in histogram.items())
    mean = total_depth / positive_sites
    median = weighted_quantile(histogram, 0.5)
    q05 = weighted_quantile(histogram, 0.05)
    q95 = weighted_quantile(histogram, 0.95)
    q_upper = weighted_quantile(histogram, upper_quantile)
    mad = median_absolute_deviation(histogram, median)

    minimum = max(absolute_minimum, math.floor(median * min_multiplier))
    multiplier_maximum = math.ceil(median * max_multiplier)
    robust_maximum = math.ceil(median + mad_multiplier * mad) if mad > 0 else multiplier_maximum
    maximum = max(minimum, min(multiplier_maximum, robust_maximum, q_upper))

    return {
        "min_depth": minimum,
        "max_depth": maximum,
        "mean_depth": mean,
        "median_depth": median,
        "mad_depth": mad,
        "q05_depth": q05,
        "q95_depth": q95,
        "upper_quantile_depth": q_upper,
        "positive_sites": positive_sites,
    }


def main():
    args = parse_args()
    lengths = {}
    with open(args.fai) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            lengths[fields[0]] = int(fields[1])

    histograms = defaultdict(Counter)
    genome_histogram = Counter()
    for line in sys.stdin:
        chrom, _position, depth_text = line.rstrip().split("\t")[:3]
        depth = int(depth_text)
        if depth > 0:
            histograms[chrom][depth] += 1
            genome_histogram[depth] += 1

    genome_summary = summarize(
        genome_histogram,
        args.absolute_minimum,
        args.min_multiplier,
        args.max_multiplier,
        args.mad_multiplier,
        args.upper_quantile,
    )
    if genome_summary is None:
        raise SystemExit("No positive-depth observations were available")

    columns = [
        "chromosome",
        "min_depth",
        "max_depth",
        "mean_depth",
        "median_depth",
        "mad_depth",
        "q05_depth",
        "q95_depth",
        "upper_quantile_depth",
        "positive_sites",
        "chromosome_length",
        "threshold_source",
    ]
    with open(args.output, "w") as output:
        output.write("\t".join(columns) + "\n")
        for chrom, length in lengths.items():
            summary = summarize(
                histograms[chrom],
                args.absolute_minimum,
                args.min_multiplier,
                args.max_multiplier,
                args.mad_multiplier,
                args.upper_quantile,
            )
            source = "chromosome"
            if summary is None or summary["positive_sites"] < args.minimum_positive_sites:
                summary = genome_summary
                source = "genome_fallback"
            values = [
                chrom,
                str(summary["min_depth"]),
                str(summary["max_depth"]),
                f"{summary['mean_depth']:.6f}",
                str(summary["median_depth"]),
                str(summary["mad_depth"]),
                str(summary["q05_depth"]),
                str(summary["q95_depth"]),
                str(summary["upper_quantile_depth"]),
                str(summary["positive_sites"]),
                str(length),
                source,
            ]
            output.write("\t".join(values) + "\n")


if __name__ == "__main__":
    main()
