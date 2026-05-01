#!/usr/bin/env python3
import re
import sys
import gzip
from pathlib import Path
from typing import TextIO


SPLICE_CIGAR = re.compile(r"\d+N")


def spliced_read_names(sam_path: Path) -> set[str]:
    names: set[str] = set()
    with sam_path.open() as sam:
        for line in sam:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 5 and SPLICE_CIGAR.search(fields[5]):
                names.add(fields[0])
    return names


def fastq_name(header: str) -> str:
    return header[1:].split(None, 1)[0]


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def write_selected_fastq(
    in_fastq: Path, out_fastq: Path, wanted: set[str], max_reads: int
) -> int:
    written = 0
    with open_text(in_fastq) as inp, out_fastq.open("w") as out:
        while True:
            header = inp.readline()
            if not header:
                break
            seq = inp.readline()
            plus = inp.readline()
            qual = inp.readline()
            if not qual:
                raise RuntimeError(f"truncated FASTQ record after {header.rstrip()}")
            if fastq_name(header) in wanted:
                out.write(header)
                out.write(seq)
                out.write(plus)
                out.write(qual)
                written += 1
                if written >= max_reads:
                    break
    return written


def main() -> int:
    if len(sys.argv) != 5:
        print(
            "usage: extract_fastq_by_sam_cigar.py <in.sam> <in.fastq> <out.fastq> <max_reads>",
            file=sys.stderr,
        )
        return 2
    sam_path, in_fastq, out_fastq = map(Path, sys.argv[1:4])
    max_reads = int(sys.argv[4])
    wanted = spliced_read_names(sam_path)
    written = write_selected_fastq(in_fastq, out_fastq, wanted, max_reads)
    print(f"wrote {written} spliced FASTQ records to {out_fastq}")
    return 0 if written > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
