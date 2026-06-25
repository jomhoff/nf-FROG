#!/usr/bin/env python3

import gzip
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class HelperTests(unittest.TestCase):
    def run_tool(self, name, *arguments, input_text=None):
        return subprocess.run(
            ["python3", str(ROOT / "bin" / name), *map(str, arguments)],
            input=input_text,
            text=True,
            capture_output=True,
            check=True,
        )

    def test_filter_masks_and_reference_comparison(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            fai = work / "reference.fa.fai"
            fai.write_text("chr1\t10\t0\t0\t0\n")
            thresholds = work / "thresholds.tsv"
            depth = "".join(f"chr1\t{position}\t10\n" for position in range(1, 11))
            self.run_tool(
                "calculate_depth_thresholds.py",
                "--fai",
                fai,
                "--minimum-positive-sites",
                "1",
                "--output",
                thresholds,
                input_text=depth,
            )

            gvcf = work / "sample.g.vcf.gz"
            with gzip.open(gvcf, "wt") as output:
                output.write(
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                    "chr1\t1\t.\tA\t<*>\t0\t.\tEND=3\tDP:MIN_DP:RO:AO\t10:10:10:0\n"
                    "chr1\t4\t.\tC\t<*>\t0\t.\tEND=5\tDP:MIN_DP:RO:AO\t1:1:1:0\n"
                    "chr1\t7\t.\tG\tA\t50\t.\tAO=5;QA=200;MQM=60;SAF=3;SAR=2;RPL=2;RPR=3\tGT:DP:RO:AO\t0/1:10:5:5\n"
                    "chr1\t8\t.\tT\t<*>\t0\t.\tEND=10\tDP:MIN_DP:RO:AO\t10:10:10:0\n"
                )
            gvcf_mask = work / "gvcf.bed"
            self.run_tool(
                "gvcf_callability_mask.py",
                "--gvcf",
                gvcf,
                "--fai",
                fai,
                "--chromosome",
                "chr1",
                "--output",
                gvcf_mask,
                "--stats",
                work / "gvcf.tsv",
            )
            self.assertIn("GVCF_ABSENT", gvcf_mask.read_text())

            depth_mask = work / "depth.bed"
            depth_rows = "".join(
                f"chr1\t{position}\t{1 if position in (4, 5) else 10}\n"
                for position in range(1, 11)
            )
            self.run_tool(
                "depth_mask.py",
                "--chromosome",
                "chr1",
                "--thresholds",
                thresholds,
                "--output",
                depth_mask,
                "--stats",
                work / "depth.tsv",
                input_text=depth_rows,
            )
            self.assertEqual(depth_mask.read_text(), "chr1\t3\t5\tLOW_DEPTH\n")

            vcf = work / "sample.vcf.gz"
            with gzip.open(vcf, "wt") as output:
                output.write(
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                    "chr1\t1\t.\tA\tG\t50\tPASS\tAO=10;QA=400;MQM=60;SAF=5;SAR=5;RPL=5;RPR=5\tGT:DP:RO:AO\t1/1:10:0:10\n"
                    "chr1\t2\t.\tC\tT\t50\tPASS\tAO=5;QA=200;MQM=60;SAF=3;SAR=2;RPL=2;RPR=3\tGT:DP:RO:AO\t0/1:10:5:5\n"
                    "chr1\t3\t.\tG\tA\t10\tPASS\tAO=10;QA=100;MQM=10;SAF=10;SAR=0;RPL=10;RPR=0\tGT:DP:RO:AO\t1/1:10:0:10\n"
                    "chr1\t4\t.\tT\tTA\t50\tPASS\tAO=10;QA=400;MQM=60;SAF=5;SAR=5;RPL=5;RPR=5\tGT:DP:RO:AO\t1/1:10:0:10\n"
                )
            output_vcf = work / "iterative.vcf"
            reject_bed = work / "reject.bed"
            self.run_tool(
                "prepare_consensus_vcf.py",
                "--vcf",
                vcf,
                "--chromosome",
                "chr1",
                "--thresholds",
                thresholds,
                "--het-mode",
                "retain",
                "--output-vcf",
                output_vcf,
                "--reject-bed",
                reject_bed,
                "--stats",
                work / "filter.tsv",
            )
            positions = [
                line.split("\t")[1]
                for line in output_vcf.read_text().splitlines()
                if not line.startswith("#")
            ]
            self.assertEqual(positions, ["1"])
            self.assertIn("LOW_QUAL", reject_bed.read_text())
            self.assertIn("INDEL_EXCLUDED", reject_bed.read_text())

            previous = work / "previous.fa"
            current = work / "current.fa"
            new = work / "new.fa"
            previous.write_text(">chr1\nACGT\n")
            current.write_text(">chr1\nAGGT\n")
            new.write_text(">chr1\nACGT\n")
            status = work / "status.txt"
            self.run_tool(
                "compare_references.py",
                "--previous",
                previous,
                "--current",
                current,
                "--new",
                new,
                "--round",
                "2",
                "--metrics",
                work / "convergence.tsv",
                "--status",
                status,
            )
            self.assertEqual(status.read_text().strip(), "active")


if __name__ == "__main__":
    unittest.main()
