#!/usr/bin/env python3

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter


def sanitize(reason):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", reason)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge reason-coded BED masks with bounded Python memory use."
    )
    parser.add_argument("--mask", action="append", default=[], help="LABEL=BED")
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--prefix", required=True)
    return parser.parse_args()


def sorted_lines(sort_path, path, *keys):
    command = [sort_path, *keys, path]
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        text=True,
        env={**os.environ, "LC_ALL": "C"},
    )
    assert process.stdout is not None
    for line in process.stdout:
        yield line
    return_code = process.wait()
    if return_code:
        raise SystemExit(f"sort failed with exit status {return_code}: {' '.join(command)}")


def stream_merged_intervals(sort_path, path):
    current = None
    for line in sorted_lines(sort_path, path, "-k1,1", "-k2,2n", "-k3,3n"):
        fields = line.rstrip().split("\t")
        if len(fields) < 3:
            continue
        chrom, start_text, end_text = fields[:3]
        start = int(start_text)
        end = int(end_text)
        if end <= start:
            continue
        if current and current[0] == chrom and start <= current[2]:
            current[2] = max(current[2], end)
        else:
            if current:
                yield tuple(current)
            current = [chrom, start, end]
    if current:
        yield tuple(current)


def write_union_segment(output, chrom, start, end, state):
    if end <= start:
        return state
    union_start, union_end = state
    if union_start is None:
        return start, end
    if start <= union_end:
        return union_start, max(union_end, end)
    output.write(f"{chrom}\t{union_start}\t{union_end}\n")
    return start, end


def close_union(output, chrom, state):
    union_start, union_end = state
    if union_start is not None:
        output.write(f"{chrom}\t{union_start}\t{union_end}\n")


def write_reason_segment(output, chrom, start, end, reasons, state):
    if end <= start:
        return state
    current_start, current_end, current_reasons = state
    if current_start is None:
        return start, end, reasons
    if start == current_end and reasons == current_reasons:
        return current_start, end, current_reasons
    output.write(f"{chrom}\t{current_start}\t{current_end}\t{current_reasons}\n")
    return start, end, reasons


def close_reason_segment(output, chrom, state):
    current_start, current_end, current_reasons = state
    if current_start is not None:
        output.write(f"{chrom}\t{current_start}\t{current_end}\t{current_reasons}\n")


def main():
    args = parse_args()
    sort_path = shutil.which("sort")
    if sort_path is None:
        raise SystemExit("merge_reason_masks.py requires POSIX sort in PATH")

    with tempfile.TemporaryDirectory(prefix="merge_reason_masks.", dir=".") as tempdir:
        reason_handles = {}
        reason_paths = {}
        event_path = os.path.join(tempdir, "events.tsv")
        interval_count = 0

        with open(event_path, "w") as events:
            for specification in args.mask:
                label, path = specification.split("=", 1)
                if not path:
                    continue
                with open(path) as handle:
                    for line in handle:
                        if not line.strip() or line.startswith("#"):
                            continue
                        fields = line.rstrip().split("\t")
                        if len(fields) < 3 or fields[0] != args.chromosome:
                            continue
                        start = int(fields[1])
                        end = int(fields[2])
                        if end <= start:
                            continue
                        reason = fields[3] if len(fields) >= 4 and fields[3] else label
                        if reason not in reason_handles:
                            reason_paths[reason] = os.path.join(
                                tempdir, f"reason_{len(reason_paths)}.bed"
                            )
                            reason_handles[reason] = open(reason_paths[reason], "w")
                        reason_handles[reason].write(f"{fields[0]}\t{start}\t{end}\n")
                        events.write(f"{start}\t1\t{reason}\n")
                        events.write(f"{end}\t-1\t{reason}\n")
                        interval_count += 1

        for handle in reason_handles.values():
            handle.close()

        if interval_count == 0:
            open(f"{args.prefix}.reason_NONE.bed", "w").close()
            open(f"{args.prefix}.reasons.bed", "w").close()
            open(f"{args.prefix}.bed", "w").close()
            return

        used_filenames = Counter()
        for reason, path in reason_paths.items():
            base = sanitize(reason) or "UNLABELED"
            used_filenames[base] += 1
            suffix = "" if used_filenames[base] == 1 else f"_{used_filenames[base]}"
            output_path = f"{args.prefix}.reason_{base}{suffix}.bed"
            with open(output_path, "w") as output:
                for chrom, start, end in stream_merged_intervals(sort_path, path):
                    output.write(f"{chrom}\t{start}\t{end}\n")

        active = Counter()
        previous_position = None
        pending_position = None
        pending_deltas = []
        union_state = (None, None)
        reason_state = (None, None, None)

        with open(f"{args.prefix}.reasons.bed", "w") as reasons_output, open(
            f"{args.prefix}.bed", "w"
        ) as union_output:
            def flush_position(position, deltas):
                nonlocal previous_position, union_state, reason_state
                if previous_position is not None and previous_position < position and active:
                    reason_text = ",".join(sorted(active))
                    reason_state = write_reason_segment(
                        reasons_output,
                        args.chromosome,
                        previous_position,
                        position,
                        reason_text,
                        reason_state,
                    )
                    union_state = write_union_segment(
                        union_output, args.chromosome, previous_position, position, union_state
                    )

                for delta, reason in deltas:
                    active[reason] += delta
                    if active[reason] <= 0:
                        del active[reason]
                previous_position = position

            for line in sorted_lines(sort_path, event_path, "-k1,1n"):
                position_text, delta_text, reason = line.rstrip().split("\t", 2)
                position = int(position_text)
                delta = int(delta_text)
                if pending_position is None:
                    pending_position = position
                if position != pending_position:
                    flush_position(pending_position, pending_deltas)
                    pending_position = position
                    pending_deltas = []
                pending_deltas.append((delta, reason))

            if pending_position is not None:
                flush_position(pending_position, pending_deltas)
            close_reason_segment(reasons_output, args.chromosome, reason_state)
            close_union(union_output, args.chromosome, union_state)


if __name__ == "__main__":
    main()
