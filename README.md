# star-rs

This is a Rust translation of the STAR RNA-seq aligner 

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

## Layout

- [`STAR/`](STAR/) contains the original C++ STAR source and documentation used as
  the translation source.
- [`src/generated/functions.rs`](src/generated/functions.rs) contains translated
  functions. Function names include source file and line information to help
  audit against the original implementation.
- [`src/generated/structs.rs`](src/generated/structs.rs) contains translated data
  structures.
- [`src/cli.rs`](src/cli.rs) contains the Rust CLI wrapper and filesystem-facing
  glue.
- [`src/direct.rs`](src/direct.rs) contains intentional Rust helpers for zero-copy
  direct access to the aligner core. These are user-approved deviations from pure
  CCC one-to-one translation where they make integration practical.
- [`ccc_mapping.toml`](ccc_mapping.toml) maps Rust functions back to their C++
  counterparts.
- [`ccc/`](ccc/) contains code-complexity-comparator output and porting order
  artifacts.
- [`tests/`](tests/) contains focused Rust tests, CLI tests, and C++ parity tests.

## Building

```sh
cargo build
```

Build an optimized binary with:

```sh
cargo build --release
```

The CLI binary is named `star-rs`.

## Running

The CLI accepts STAR-style arguments. For example:

```sh
cargo run -- \
  --genomeDir path/to/genome \
  --readFilesIn reads.fq \
  --outSAMtype SAM \
  --outFileNamePrefix out/
```

Genome generation uses the same CLI entry point:

```sh
cargo run -- \
  --runMode genomeGenerate \
  --genomeDir path/to/genome \
  --genomeFastaFiles genome.fa
```

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

## Translation Workflow

1. Keep the C++ source in [`STAR/source`](STAR/source/) as the ground truth.
2. Use [`ccc_mapping.toml`](ccc_mapping.toml) to preserve traceability from Rust
   functions back to original C++ functions.
3. Use code-complexity-comparator output in [`ccc/`](ccc/) to choose bottom-up
   implementation order.
4. Translate complete logic in the first pass for each function where feasible.
5. Prefer preserving original control flow and data layout over idiomatic Rust
   refactors when auditability and output parity would otherwise suffer.
6. Add parity tests around real command-line behavior before relying on larger
   refactors.

Helper functions are allowed in this repository when they serve the approved
direct-access aligner interface or isolate Rust-only filesystem/CLI glue. For
generated translation code, new helpers should still be treated skeptically unless
they are necessary for faithful behavior or safe Rust ownership.

## Optional Tracehash Feature

The crate has an optional `tracehash` feature for trace-oriented translation
verification:

```sh
cargo test --features tracehash
```



## Citing

Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras, STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, January 2013, Pages 15–21, https://doi.org/10.1093/bioinformatics/bts635


## License

MIT license
