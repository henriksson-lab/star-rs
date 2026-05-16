#!/usr/bin/env python3
"""Benchmark genome indexing wall time and peak RSS for original STAR vs star-rs.

Runs --runMode genomeGenerate on a small reference and measures:
- wall-clock time (seconds)
- maximum resident set size (KB, from /usr/bin/time -v)
- exit status

Designed for triaging the slow indexing seen on the full GRCh38 input.
Yeast is small enough that both binaries should finish in seconds-to-minutes.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import time


def run_timed(cmd: list[str], log_path: pathlib.Path) -> dict:
    """Run cmd under /usr/bin/time -v, capture wall time and peak RSS."""
    time_log = log_path.with_suffix(".time")
    full_cmd = ["/usr/bin/time", "-v", "-o", str(time_log)] + cmd
    start = time.perf_counter()
    with open(log_path, "wb") as f:
        proc = subprocess.run(full_cmd, stdout=f, stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - start

    rss_kb = None
    wall_str = None
    user_s = None
    sys_s = None
    if time_log.exists():
        txt = time_log.read_text()
        m = re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", txt)
        if m:
            rss_kb = int(m.group(1))
        m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s+(\S+)", txt)
        if m:
            wall_str = m.group(1)
        m = re.search(r"User time \(seconds\):\s+(\S+)", txt)
        if m:
            user_s = float(m.group(1))
        m = re.search(r"System time \(seconds\):\s+(\S+)", txt)
        if m:
            sys_s = float(m.group(1))

    return {
        "elapsed_s": elapsed,
        "rss_kb": rss_kb,
        "wall_str": wall_str,
        "user_s": user_s,
        "sys_s": sys_s,
        "returncode": proc.returncode,
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--star", default="STAR/bin/Linux_x86_64_static/STAR")
    p.add_argument("--star-rs", default="target/release/star-rs")
    p.add_argument("--fasta", required=True, help="FASTA file to index")
    p.add_argument("--gtf", default=None, help="optional GTF file")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--sa-index-nbases", type=int, default=10,
                   help="--genomeSAindexNbases (small genomes need lower than 14)")
    p.add_argument("--chr-bin-nbits", type=int, default=18,
                   help="--genomeChrBinNbits")
    p.add_argument("--sjdb-overhang", type=int, default=100)
    p.add_argument("--limit-genome-ram", type=int, default=32_000_000_000,
                   help="--limitGenomeGenerateRAM in bytes")
    p.add_argument("--out-dir", default=None,
                   help="working directory; created tempdir if omitted")
    p.add_argument("--keep", action="store_true")
    p.add_argument("--label", default=None,
                   help="optional label printed with results")
    p.add_argument("--skip", choices=["none", "cpp", "rust"], default="none")
    return p.parse_args()


def index_cmd(binary: pathlib.Path, genome_dir: pathlib.Path, args: argparse.Namespace,
              fasta: pathlib.Path, gtf: pathlib.Path | None) -> list[str]:
    cmd = [
        str(binary),
        "--runMode", "genomeGenerate",
        "--runThreadN", str(args.threads),
        "--genomeDir", str(genome_dir),
        "--genomeFastaFiles", str(fasta),
        "--genomeSAindexNbases", str(args.sa_index_nbases),
        "--genomeChrBinNbits", str(args.chr_bin_nbits),
        "--limitGenomeGenerateRAM", str(args.limit_genome_ram),
    ]
    if gtf is not None:
        cmd += ["--sjdbGTFfile", str(gtf), "--sjdbOverhang", str(args.sjdb_overhang)]
    return cmd


def fmt_rss(rss_kb: int | None) -> str:
    if rss_kb is None:
        return "?"
    mb = rss_kb / 1024.0
    if mb < 1024:
        return f"{mb:.1f} MB"
    return f"{mb / 1024.0:.2f} GB"


def main() -> int:
    args = parse_args()
    root = pathlib.Path.cwd()
    star = (root / args.star).resolve()
    star_rs = (root / args.star_rs).resolve()
    fasta = pathlib.Path(args.fasta).resolve()
    gtf = pathlib.Path(args.gtf).resolve() if args.gtf else None

    for path in [star, star_rs, fasta]:
        if not path.exists():
            print(f"missing: {path}", file=sys.stderr)
            return 2
    if gtf and not gtf.exists():
        print(f"missing gtf: {gtf}", file=sys.stderr)
        return 2

    if args.out_dir:
        work = pathlib.Path(args.out_dir).resolve()
        work.mkdir(parents=True, exist_ok=True)
    else:
        work = pathlib.Path(tempfile.mkdtemp(prefix="star_index_bench_"))

    cpp_dir = work / "cpp_index"
    rust_dir = work / "rust_index"
    for d in (cpp_dir, rust_dir):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

    label = args.label or fasta.name
    fasta_bytes = fasta.stat().st_size
    print(f"# bench label: {label}")
    print(f"# fasta: {fasta} ({fasta_bytes / 1e6:.2f} MB)")
    if gtf:
        print(f"# gtf:   {gtf} ({gtf.stat().st_size / 1e6:.2f} MB)")
    print(f"# threads={args.threads} saIndexNbases={args.sa_index_nbases} "
          f"chrBinNbits={args.chr_bin_nbits} sjdbOverhang={args.sjdb_overhang}")
    print(f"# work_dir: {work}")
    print()

    results = {}

    if args.skip != "cpp":
        print("# running original STAR ...", flush=True)
        cmd = index_cmd(star, cpp_dir, args, fasta, gtf)
        # original STAR cares about CWD for Log.out etc
        cwd_cpp = work / "cpp_run"
        cwd_cpp.mkdir(exist_ok=True)
        # Run inside cwd
        old = pathlib.Path.cwd()
        try:
            import os
            os.chdir(cwd_cpp)
            results["cpp"] = run_timed(cmd, work / "cpp.log")
        finally:
            os.chdir(old)

    if args.skip != "rust":
        print("# running star-rs ...", flush=True)
        cmd = index_cmd(star_rs, rust_dir, args, fasta, gtf)
        cwd_rust = work / "rust_run"
        cwd_rust.mkdir(exist_ok=True)
        old = pathlib.Path.cwd()
        try:
            import os
            os.chdir(cwd_rust)
            results["rust"] = run_timed(cmd, work / "rust.log")
        finally:
            os.chdir(old)

    print()
    print(f"{'binary':<12} {'wall_s':>10} {'user_s':>10} {'sys_s':>10} {'rss':>10} {'rc':>4}")
    for key, r in results.items():
        rss = fmt_rss(r["rss_kb"])
        wall = f"{r['elapsed_s']:.2f}"
        user = f"{r['user_s']:.2f}" if r['user_s'] is not None else "?"
        sysv = f"{r['sys_s']:.2f}" if r['sys_s'] is not None else "?"
        print(f"{key:<12} {wall:>10} {user:>10} {sysv:>10} {rss:>10} {r['returncode']:>4}")

    if "cpp" in results and "rust" in results:
        c, r = results["cpp"], results["rust"]
        if c["returncode"] == 0 and r["returncode"] == 0 and c["elapsed_s"] > 0:
            print()
            print(f"ratio_rust_over_cpp_wall\t{r['elapsed_s'] / c['elapsed_s']:.3f}")
            if c["rss_kb"] and r["rss_kb"]:
                print(f"ratio_rust_over_cpp_rss\t{r['rss_kb'] / c['rss_kb']:.3f}")

    if not args.keep and not args.out_dir:
        shutil.rmtree(work, ignore_errors=True)
    else:
        print(f"\n# kept {work}")

    rcs = [r["returncode"] for r in results.values()]
    return 0 if all(rc == 0 for rc in rcs) else 1


if __name__ == "__main__":
    raise SystemExit(main())
