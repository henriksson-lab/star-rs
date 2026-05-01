#!/usr/bin/env python3
import sys
from pathlib import Path


def keep_fasta_chrom(in_fasta: Path, out_fasta: Path, chrom: str) -> None:
    keep = False
    with in_fasta.open() as inp, out_fasta.open("w") as out:
        for line in inp:
            if line.startswith(">"):
                name = line[1:].split()[0]
                keep = name == chrom
            if keep:
                out.write(line)


def keep_gtf_chrom(in_gtf: Path, out_gtf: Path, chrom: str) -> None:
    with in_gtf.open() as inp, out_gtf.open("w") as out:
        for line in inp:
            if line.startswith("#") or line.split("\t", 1)[0] == chrom:
                out.write(line)


def main() -> int:
    if len(sys.argv) != 6:
        print(
            "usage: make_yeast_subset.py <in.fa> <in.gtf> <out.fa> <out.gtf> <chrom>",
            file=sys.stderr,
        )
        return 2
    in_fasta, in_gtf, out_fasta, out_gtf = map(Path, sys.argv[1:5])
    chrom = sys.argv[5]
    keep_fasta_chrom(in_fasta, out_fasta, chrom)
    keep_gtf_chrom(in_gtf, out_gtf, chrom)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
