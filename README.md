# star-rs

Rust translation of the STAR RNA-seq aligner 

**more auditing needed; uint mistranslated**

* 2026-08-01: CI added. More audit pending
* 2026-05-16: Restructured code, fixed STAR index bug + multithreading
* 2026-05-13: Automatically detects gzip files if provided by name (not piped, CLI). Binary is now an optional feature
* 2026-05-02: Passes basic tests but a larger battery is needed before we can assert that this translation is fully functional


## This is an LLM-mediated faithful (hopefully) translation, not the original code! 

Most users should probably first see if the existing original code works for them, unless they have reason otherwise. The original source
may have newer features and it has had more love in terms of fixing bugs. In fact, we aim to replicate bugs if they are present, for the
sake of reproducibility! (but then we might have added a few more in the process)

There are however cases when you might prefer this Rust version. We generally agree with [this manifesto](https://rewrites.bio/) but more specifically:
* We have had many issues with ensuring that our software works using existing containers (Docker, PodMan, Singularity). One size does not fit all and it eats our resources trying to keep up with every way of delivering software
* Common package managers do not work well. It was great when we had a few Linux distributions with stable procedures, but now there are just too many ecosystems (Homebrew, Conda). Conda has an NP-complete resolver which does not scale. Homebrew is only so-stable. And our dependencies in Python still break. These can no longer be considered professional serious options. Meanwhile, Cargo enables multiple versions of packages to be available, even within the same program(!)
* The future is the web. We deploy software in the web browser, and until now that has meant Javascript. This is a language where even the == operator is broken. Typescript is one step up, but a game changer is the ability to compile Rust code into webassembly, enabling performance and sharing of code with the backend. Translating code to Rust enables new ways of deployment and running code in the browser has especial benefits for science - researchers do not have deep pockets to run servers, so pushing compute to the user enables deployment that otherwise would be impossible
* Old CLI-based utilities are bad for the environment(!). A large amount of compute resources are spent creating and communicating via small files, which we can bypass by using code as libraries. Even better, we can avoid frequent reloading of databases by hoisting this stage, with up to 100x speedups in some cases. Less compute means faster compute and less electricity wasted
* LLM-mediated translations may actually be safer to use than the original code. This article shows that [running the same code on different operating systems can give somewhat different answers](https://doi.org/10.1038/nbt.3820). This is a gap that Rust+Cargo can reduce. Typesafe interfaces also reduce coding mistakes and error handling, as opposed to typical command-line scripting

But:

* **This approach should still be considered experimental**. The LLM technology is immature and has sharp corners. But there are opportunities to reap, and the genie is not going back into the bottle. This translation is as much aimed to learn how to improve the technology and get feedback on the results.
* Translations are not endorsed by the original authors unless otherwise noted. **Do not send bug reports to the original developers**. Use our Github issues page instead.
* **Do not trust the benchmarks on this page**. They are used to help evaluate the translation. If you want improved performance, you generally have to use this code as a library, and use the additional tricks it offers. We generally accept performance losses in order to reduce our dependency issues
* **Check the original Github pages for information about the package**. This README is kept sparse on purpose. It is not meant to be the primary source of information
* **If you are the author of the original code and wish to move to Rust, you can obtain ownership of this repository and crate**. Until then, our commitment is to offer an as-faithful-as-possible translation of a snapshot of your code. If we find serious bugs, we will report them to you. Otherwise we will just replicate them, to ensure comparability across studies that claim to use package XYZ v.666. Think of this like a fancy Ubuntu .deb-package of your software - that is how we treat it

This blurb might be out of date. Go to [this page](https://github.com/henriksson-lab/rustification) for the latest information and further information about how we approach translation

## Building

Build the optional CLI binary with:

```sh
cargo build --release --features binary
```

The CLI binary is named `star-rs`. The STAR-style CLI implementation is also
available as the library module `star_rs::cli`.

## Running

The CLI accepts STAR-style arguments. For example:

```sh
cargo run --features binary -- \
  --genomeDir path/to/genome \
  --readFilesIn reads.fq \
  --outSAMtype SAM \
  --outFileNamePrefix out/
```

Genome generation uses the same CLI entry point:

```sh
cargo run --features binary -- \
  --runMode genomeGenerate \
  --genomeDir path/to/genome \
  --genomeFastaFiles genome.fa
```

## Behavior Differences From Original STAR

This translation aims to preserve original STAR behavior, but it also includes
some Rust-specific CLI/filesystem glue where that makes use more practical.

- Named input files are automatically detected as gzip-compressed by their magic
  bytes and decompressed when loaded. This applies to read files, GTF files, and
  genome index files loaded from `--genomeDir`.
- Explicit streaming or command-based input remains under user control. In
  particular, `--readFilesCommand` still runs the command requested by the user,
  and data provided through such commands is not re-detected or re-decompressed
  by star-rs.

## Testing

Run the full Rust test suite:

```sh
cargo test
```

Run focused test groups while working on the aligner core:

```sh
cargo test --test core_leaf
cargo test --test cli
cargo test --test cpp_parity
```

Several parity tests compare Rust behavior with the vendored C++ STAR source or
binary. Prefer real-world or externally derived data when adding broader tests;
synthetic fixtures are useful for isolating failures but should not be the only
evidence for parity.

The env-gated real-world parity test can be prepared with a splice-enriched
yeast RNA-seq fixture:

```sh
tools/prepare_yeast_conformance.sh
source .tmp/yeast_conformance/env.sh
cargo test real_world_conformance_from_env_matches_original_star_core_sam_fields --test cpp_parity -- --nocapture
```

The preparation script uses original STAR to select real reads whose alignments
contain splice-junction CIGARs, and the test requires at least one splice event
unless `STAR_RS_REAL_MIN_SPLICES=0` is set.

## Benchmarking

Original benchmark baseline: vendored upstream STAR commit `b1edc1208d91` (`2.7.11b`).

The included benchmark helpers are smoke checks for translation work, not
published performance claims. Build the release binary first:

```sh
cargo build --release --features binary
```

Run the small read-mapping baseline with:

```sh
tools/perf_baseline.py
```

The output includes both wall-time and peak RSS comparisons:

```text
ratio_star_rs_over_original_wall    ...
ratio_star_rs_over_original_rss     ...
```

Current smoke result from this fixture, measured on 2026-07-07 with original
STAR 2.7.11b and `cargo build --release --features binary` on Linux 6.8
(`x86_64`, Intel Xeon Gold 6138):

| benchmark | original STAR | star-rs | star-rs / original |
| --- | ---: | ---: | ---: |
| read mapping wall time, median of 5 | 0.125 s | 0.068 s | 0.537 |
| read mapping peak RSS | 156,160 KB | 79,360 KB | 0.508 |

The wall-time ratio varied from 0.387 to 0.580 across the five runs, so treat
the numbers as a regression smoke signal rather than a stable throughput
benchmark.

Run a genome-generation comparison on a chosen FASTA with:

```sh
tools/bench_genome_generate.py --fasta path/to/genome.fa
```

When both binaries finish successfully, it reports:

```text
ratio_rust_over_cpp_wall    ...
ratio_rust_over_cpp_rss     ...
```

## Citing

Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras, STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, January 2013, Pages 15–21, https://doi.org/10.1093/bioinformatics/bts635


## License

MIT license
