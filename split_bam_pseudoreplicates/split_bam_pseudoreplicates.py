#!/usr/bin/env python3
"""Split a BAM into two complementary, reproducible IDR pseudoreplicates."""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import tempfile
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Assign each query name to exactly one of two pseudoreplicates, "
            "then coordinate-sort and index both BAM outputs."
        )
    )
    parser.add_argument("--input", required=True, help="Input BAM")
    parser.add_argument("--output-a", required=True, help="First output BAM")
    parser.add_argument("--output-b", required=True, help="Second output BAM")
    parser.add_argument("--index-a", required=True, help="Index for first BAM")
    parser.add_argument("--index-b", required=True, help="Index for second BAM")
    parser.add_argument("--report", required=True, help="TSV report")
    parser.add_argument("--seed", type=int, default=12345, help="Non-negative seed")
    parser.add_argument("--threads", type=int, default=1, help="Samtools threads")
    return parser.parse_args()


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def assignment(query_name: str, seed: int) -> int:
    """Return 0 or 1 using a stable keyed hash of the complete query name."""
    seed_key = seed.to_bytes(16, byteorder="little", signed=False)
    digest = hashlib.blake2b(
        query_name.encode("utf-8", errors="surrogateescape"),
        digest_size=8,
        key=seed_key,
        person=b"galaxy-idr-split",
    ).digest()
    return digest[0] & 1


def split_to_sam(input_bam: str, sam_a: Path, sam_b: Path, seed: int) -> tuple[int, int, int]:
    counts = [0, 0]
    header_lines = 0
    view = subprocess.Popen(
        ["samtools", "view", "-h", input_bam],
        stdout=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="surrogateescape",
    )
    assert view.stdout is not None

    try:
        with sam_a.open("w", encoding="utf-8", errors="surrogateescape") as out_a, sam_b.open(
            "w", encoding="utf-8", errors="surrogateescape"
        ) as out_b:
            for line in view.stdout:
                if line.startswith("@"):
                    out_a.write(line)
                    out_b.write(line)
                    header_lines += 1
                    continue

                query_name = line.split("\t", 1)[0]
                target = assignment(query_name, seed)
                (out_a if target == 0 else out_b).write(line)
                counts[target] += 1
    except Exception:
        view.kill()
        view.wait()
        raise

    return_code = view.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, view.args)
    if sum(counts) == 0:
        raise ValueError("The input BAM contains no alignment records")
    if 0 in counts:
        raise ValueError(
            "One pseudoreplicate is empty. The input contains too few distinct query names "
            "for this seed; choose another seed or provide a larger BAM."
        )
    return counts[0], counts[1], header_lines


def sam_to_sorted_bam(sam_path: Path, output_bam: str, output_index: str, threads: int) -> None:
    threads_text = str(max(1, threads))
    run(["samtools", "sort", "-@", threads_text, "-O", "BAM", "-o", output_bam, str(sam_path)])
    run(["samtools", "index", "-@", threads_text, "-o", output_index, output_bam])
    run(["samtools", "quickcheck", "-v", output_bam])


def main() -> None:
    args = parse_args()
    if args.seed < 0:
        raise ValueError("--seed must be non-negative")
    if args.seed >= 2**128:
        raise ValueError("--seed must be smaller than 2^128")

    for output in (args.output_a, args.output_b, args.index_a, args.index_b, args.report):
        Path(output).parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="idr_pseudoreps_") as temp_dir:
        temp_path = Path(temp_dir)
        sam_a = temp_path / "pseudoreplicate_a.sam"
        sam_b = temp_path / "pseudoreplicate_b.sam"
        count_a, count_b, header_lines = split_to_sam(args.input, sam_a, sam_b, args.seed)
        sam_to_sorted_bam(sam_a, args.output_a, args.index_a, args.threads)
        sam_to_sorted_bam(sam_b, args.output_b, args.index_b, args.threads)

    total = count_a + count_b
    with open(args.report, "w", encoding="utf-8") as report:
        report.write("metric\tvalue\n")
        report.write(f"seed\t{args.seed}\n")
        report.write("assignment_unit\tquery_name_all_records\n")
        report.write("assignment_method\tblake2b_keyed_hash\n")
        report.write(f"header_lines\t{header_lines}\n")
        report.write(f"input_alignment_records\t{total}\n")
        report.write(f"pseudoreplicate_1_alignment_records\t{count_a}\n")
        report.write(f"pseudoreplicate_2_alignment_records\t{count_b}\n")
        report.write(f"pseudoreplicate_1_fraction\t{count_a / total:.8f}\n")
        report.write(f"pseudoreplicate_2_fraction\t{count_b / total:.8f}\n")
        report.write("complementary_partitions\ttrue\n")
        report.write("coordinate_sorted\ttrue\n")
        report.write("indexed\ttrue\n")


if __name__ == "__main__":
    main()
