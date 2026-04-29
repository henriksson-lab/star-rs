use star_rs::cli::run_cli;
use std::path::{Path, PathBuf};
use std::process::Command;

fn unique_temp_dir(label: &str) -> PathBuf {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("star_rs_{label}_{}_{}", std::process::id(), unique))
}

fn run_original_star(star_bin: &Path, args: &[String]) {
    let run_dir = unique_temp_dir("original_star_cwd");
    std::fs::create_dir_all(&run_dir).unwrap();
    let output = Command::new(star_bin)
        .args(args)
        .current_dir(&run_dir)
        .output()
        .expect("failed to run original STAR binary");
    let _ = std::fs::remove_dir_all(&run_dir);
    assert!(
        output.status.success(),
        "original STAR failed\nstatus: {}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn sam_records(sam: &str) -> Vec<String> {
    sam.lines()
        .filter(|line| !line.starts_with('@'))
        .filter(|line| !line.as_bytes().iter().all(|b| *b == 0))
        .map(ToString::to_string)
        .collect()
}

#[test]
fn rust_matches_original_star_on_cpp_generated_unique_read() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let star_bin = manifest_dir.join("STAR/bin/Linux_x86_64_static/STAR");
    if !star_bin.exists() {
        eprintln!(
            "skipping original STAR parity test: {} is absent",
            star_bin.display()
        );
        return;
    }

    let dir = unique_temp_dir("cpp_parity_unique");
    let cpp_genome_dir = dir.join("cpp_genome");
    let rust_genome_dir = dir.join("rust_genome");
    let cpp_prefix = dir.join("cpp/");
    let rust_prefix = dir.join("rust/");
    std::fs::create_dir_all(&cpp_genome_dir).unwrap();
    std::fs::create_dir_all(&rust_genome_dir).unwrap();

    let seq = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA";
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, format!(">chr1\n{seq}\n")).unwrap();
    std::fs::write(
        &reads,
        format!("@r1\n{seq}\n+\n{}\n", "F".repeat(seq.len())),
    )
    .unwrap();

    run_original_star(
        &star_bin,
        &[
            "--runMode".to_string(),
            "genomeGenerate".to_string(),
            "--genomeDir".to_string(),
            cpp_genome_dir.to_string_lossy().to_string(),
            "--genomeFastaFiles".to_string(),
            fasta.to_string_lossy().to_string(),
            "--genomeSAindexNbases".to_string(),
            "2".to_string(),
            "--genomeChrBinNbits".to_string(),
            "4".to_string(),
            "--limitGenomeGenerateRAM".to_string(),
            "1000000".to_string(),
        ],
    );
    run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        rust_genome_dir.to_string_lossy().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_string_lossy().to_string(),
        "--genomeSAindexNbases".to_string(),
        "2".to_string(),
        "--genomeChrBinNbits".to_string(),
        "4".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();

    for index_file in [
        "Genome",
        "SA",
        "SAindex",
        "chrName.txt",
        "chrLength.txt",
        "chrStart.txt",
    ] {
        assert_eq!(
            std::fs::read(cpp_genome_dir.join(index_file)).unwrap(),
            std::fs::read(rust_genome_dir.join(index_file)).unwrap(),
            "{index_file} differs between original STAR and Rust genomeGenerate"
        );
    }

    let cpp_align_args = [
        "--genomeDir",
        cpp_genome_dir.to_str().unwrap(),
        "--readFilesIn",
        reads.to_str().unwrap(),
        "--outSAMtype",
        "SAM",
        "--seedSplitMin",
        "12",
        "--seedMapMin",
        "5",
        "--outFilterMatchNmin",
        "20",
        "--outFilterScoreMin",
        "20",
    ];
    let rust_align_args = [
        "--genomeDir",
        rust_genome_dir.to_str().unwrap(),
        "--readFilesIn",
        reads.to_str().unwrap(),
        "--outSAMtype",
        "SAM",
        "--seedSplitMin",
        "12",
        "--seedMapMin",
        "5",
        "--outFilterMatchNmin",
        "20",
        "--outFilterScoreMin",
        "20",
    ];

    run_original_star(
        &star_bin,
        &[
            cpp_align_args
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>(),
            vec![
                "--outFileNamePrefix".to_string(),
                cpp_prefix.to_string_lossy().to_string(),
            ],
        ]
        .concat(),
    );

    let rust_result = run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--outFileNamePrefix".to_string(),
                rust_prefix.to_string_lossy().to_string(),
            ],
            rust_align_args
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>(),
        ]
        .concat(),
    )
    .unwrap();
    assert_eq!(rust_result.exit_code, 0);
    assert_eq!(rust_result.stats_all.read_n, 1);
    assert_eq!(rust_result.stats_all.mapped_reads_u, 1);

    let cpp_sam = std::fs::read_to_string(dir.join("cpp/Aligned.out.sam")).unwrap();
    let rust_sam = std::fs::read_to_string(dir.join("rust/Aligned.out.sam")).unwrap();
    assert_eq!(sam_records(&cpp_sam), sam_records(&rust_sam));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn rust_matches_original_star_too_many_paired_repeated_read() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let star_bin = manifest_dir.join("STAR/bin/Linux_x86_64_static/STAR");
    if !star_bin.exists() {
        eprintln!(
            "skipping original STAR parity test: {} is absent",
            star_bin.display()
        );
        return;
    }

    let dir = unique_temp_dir("cpp_parity_paired_too_many");
    let genome_dir = dir.join("genome");
    let cpp_prefix = dir.join("cpp/");
    let rust_prefix = dir.join("rust/");
    std::fs::create_dir_all(&genome_dir).unwrap();

    let fasta = dir.join("tiny.fa");
    let reads_1 = dir.join("reads_1.fq");
    let reads_2 = dir.join("reads_2.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads_1, "@pair1/1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();
    std::fs::write(&reads_2, "@pair1/2\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    run_original_star(
        &star_bin,
        &[
            "--runMode".to_string(),
            "genomeGenerate".to_string(),
            "--genomeDir".to_string(),
            genome_dir.to_string_lossy().to_string(),
            "--genomeFastaFiles".to_string(),
            fasta.to_string_lossy().to_string(),
            "--genomeSAindexNbases".to_string(),
            "1".to_string(),
            "--genomeChrBinNbits".to_string(),
            "2".to_string(),
            "--limitGenomeGenerateRAM".to_string(),
            "1000000".to_string(),
        ],
    );

    let align_args = [
        "--genomeDir",
        genome_dir.to_str().unwrap(),
        "--readFilesIn",
        reads_1.to_str().unwrap(),
        reads_2.to_str().unwrap(),
        "--outSAMtype",
        "SAM",
        "--seedSplitMin",
        "1",
        "--seedMapMin",
        "0",
        "--outFilterMatchNmin",
        "0",
        "--outFilterScoreMin",
        "0",
    ];

    run_original_star(
        &star_bin,
        &[
            align_args.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
            vec![
                "--outFileNamePrefix".to_string(),
                cpp_prefix.to_string_lossy().to_string(),
            ],
        ]
        .concat(),
    );

    let rust_result = run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--outFileNamePrefix".to_string(),
                rust_prefix.to_string_lossy().to_string(),
            ],
            align_args.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
        ]
        .concat(),
    )
    .unwrap();
    assert_eq!(rust_result.exit_code, 0);
    assert_eq!(rust_result.stats_all.read_n, 1);
    assert_eq!(rust_result.stats_all.mapped_reads_m, 0);
    assert_eq!(rust_result.stats_all.unmapped_multi, 1);

    let cpp_sam = std::fs::read_to_string(dir.join("cpp/Aligned.out.sam")).unwrap();
    let rust_sam = std::fs::read_to_string(dir.join("rust/Aligned.out.sam")).unwrap();
    assert_eq!(sam_records(&cpp_sam), sam_records(&rust_sam));
    assert!(sam_records(&rust_sam).is_empty());

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn rust_matches_original_star_high_nmax_paired_repeated_read() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let star_bin = manifest_dir.join("STAR/bin/Linux_x86_64_static/STAR");
    if !star_bin.exists() {
        eprintln!(
            "skipping original STAR parity test: {} is absent",
            star_bin.display()
        );
        return;
    }

    let dir = unique_temp_dir("cpp_parity_paired_high_nmax");
    let genome_dir = dir.join("genome");
    let cpp_prefix = dir.join("cpp/");
    let rust_prefix = dir.join("rust/");
    std::fs::create_dir_all(&genome_dir).unwrap();

    let fasta = dir.join("tiny.fa");
    let reads_1 = dir.join("reads_1.fq");
    let reads_2 = dir.join("reads_2.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads_1, "@pair1/1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();
    std::fs::write(&reads_2, "@pair1/2\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    run_original_star(
        &star_bin,
        &[
            "--runMode".to_string(),
            "genomeGenerate".to_string(),
            "--genomeDir".to_string(),
            genome_dir.to_string_lossy().to_string(),
            "--genomeFastaFiles".to_string(),
            fasta.to_string_lossy().to_string(),
            "--genomeSAindexNbases".to_string(),
            "1".to_string(),
            "--genomeChrBinNbits".to_string(),
            "2".to_string(),
            "--limitGenomeGenerateRAM".to_string(),
            "1000000".to_string(),
        ],
    );

    let align_args = [
        "--genomeDir",
        genome_dir.to_str().unwrap(),
        "--readFilesIn",
        reads_1.to_str().unwrap(),
        reads_2.to_str().unwrap(),
        "--outSAMtype",
        "SAM",
        "--seedSplitMin",
        "1",
        "--seedMapMin",
        "0",
        "--outFilterMatchNmin",
        "0",
        "--outFilterScoreMin",
        "0",
        "--outFilterMultimapNmax",
        "100",
        "--outSAMmultNmax",
        "100",
    ];

    run_original_star(
        &star_bin,
        &[
            align_args.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
            vec![
                "--outFileNamePrefix".to_string(),
                cpp_prefix.to_string_lossy().to_string(),
            ],
        ]
        .concat(),
    );

    let rust_result = run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--outFileNamePrefix".to_string(),
                rust_prefix.to_string_lossy().to_string(),
            ],
            align_args.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
        ]
        .concat(),
    )
    .unwrap();
    assert_eq!(rust_result.exit_code, 0);
    assert_eq!(rust_result.stats_all.read_n, 1);
    assert_eq!(rust_result.stats_all.mapped_reads_m, 1);
    assert_eq!(rust_result.stats_all.unmapped_multi, 0);

    let cpp_sam = std::fs::read_to_string(dir.join("cpp/Aligned.out.sam")).unwrap();
    let rust_sam = std::fs::read_to_string(dir.join("rust/Aligned.out.sam")).unwrap();
    let cpp_records = sam_records(&cpp_sam);
    let rust_records = sam_records(&rust_sam);
    assert_eq!(cpp_records.len(), 112);
    assert_eq!(rust_records.len(), 112);
    assert_eq!(cpp_records, rust_records);

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn rust_matches_original_star_annotated_spliced_read() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let star_bin = manifest_dir.join("STAR/bin/Linux_x86_64_static/STAR");
    if !star_bin.exists() {
        eprintln!(
            "skipping original STAR parity test: {} is absent",
            star_bin.display()
        );
        return;
    }

    let dir = unique_temp_dir("cpp_parity_annotated_splice");
    let cpp_genome_dir = dir.join("cpp_genome");
    let rust_genome_dir = dir.join("rust_genome");
    let cpp_prefix = dir.join("cpp/");
    let rust_prefix = dir.join("rust/");
    std::fs::create_dir_all(&cpp_genome_dir).unwrap();
    std::fs::create_dir_all(&rust_genome_dir).unwrap();

    let exon1 = "ACGTACGTAA";
    let intron = "C".repeat(90);
    let exon2 = "TTGGAACCTT";
    let fasta = dir.join("spliced.fa");
    let gtf = dir.join("spliced.gtf");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, format!(">chr1\n{exon1}{intron}{exon2}\n")).unwrap();
    std::fs::write(
        &gtf,
        concat!(
            "chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n",
            "chr1\ttest\texon\t101\t110\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &reads,
        format!("@splice1\n{exon1}{exon2}\n+\n{}\n", "F".repeat(20)),
    )
    .unwrap();

    let generate_args = [
        "--runMode",
        "genomeGenerate",
        "--genomeFastaFiles",
        fasta.to_str().unwrap(),
        "--sjdbGTFfile",
        gtf.to_str().unwrap(),
        "--sjdbOverhang",
        "9",
        "--genomeSAindexNbases",
        "2",
        "--genomeChrBinNbits",
        "4",
        "--limitGenomeGenerateRAM",
        "1000000",
    ];

    run_original_star(
        &star_bin,
        &[
            generate_args
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>(),
            vec![
                "--genomeDir".to_string(),
                cpp_genome_dir.to_string_lossy().to_string(),
            ],
        ]
        .concat(),
    );
    run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                rust_genome_dir.to_string_lossy().to_string(),
            ],
            generate_args
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>(),
        ]
        .concat(),
    )
    .unwrap();

    let align_args = [
        "--readFilesIn",
        reads.to_str().unwrap(),
        "--outSAMtype",
        "SAM",
        "--seedSplitMin",
        "1",
        "--seedMapMin",
        "0",
        "--outFilterMatchNmin",
        "0",
        "--outFilterScoreMin",
        "0",
    ];

    run_original_star(
        &star_bin,
        &[
            vec![
                "--genomeDir".to_string(),
                cpp_genome_dir.to_string_lossy().to_string(),
                "--outFileNamePrefix".to_string(),
                cpp_prefix.to_string_lossy().to_string(),
            ],
            align_args.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
        ]
        .concat(),
    );

    let rust_result = run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                rust_genome_dir.to_string_lossy().to_string(),
                "--outFileNamePrefix".to_string(),
                rust_prefix.to_string_lossy().to_string(),
            ],
            align_args.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
        ]
        .concat(),
    )
    .unwrap();
    assert_eq!(rust_result.exit_code, 0);

    let cpp_sam = std::fs::read_to_string(dir.join("cpp/Aligned.out.sam")).unwrap();
    let rust_sam = std::fs::read_to_string(dir.join("rust/Aligned.out.sam")).unwrap();
    assert_eq!(sam_records(&cpp_sam), sam_records(&rust_sam));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn rust_matches_original_star_novel_spliced_read() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let star_bin = manifest_dir.join("STAR/bin/Linux_x86_64_static/STAR");
    if !star_bin.exists() {
        eprintln!(
            "skipping original STAR parity test: {} is absent",
            star_bin.display()
        );
        return;
    }

    let dir = unique_temp_dir("cpp_parity_novel_splice");
    let cpp_genome_dir = dir.join("cpp_genome");
    let rust_genome_dir = dir.join("rust_genome");
    let cpp_prefix = dir.join("cpp/");
    let rust_prefix = dir.join("rust/");
    std::fs::create_dir_all(&cpp_genome_dir).unwrap();
    std::fs::create_dir_all(&rust_genome_dir).unwrap();

    let exon1 = "ACGTACGTAA";
    let intron = format!("GT{}AG", "C".repeat(86));
    let exon2 = "TTGGAACCTT";
    let fasta = dir.join("novel_spliced.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, format!(">chr1\n{exon1}{intron}{exon2}\n")).unwrap();
    std::fs::write(
        &reads,
        format!("@splice1\n{exon1}{exon2}\n+\n{}\n", "F".repeat(20)),
    )
    .unwrap();

    let generate_args = [
        "--runMode",
        "genomeGenerate",
        "--genomeFastaFiles",
        fasta.to_str().unwrap(),
        "--genomeSAindexNbases",
        "2",
        "--genomeChrBinNbits",
        "4",
        "--limitGenomeGenerateRAM",
        "1000000",
    ];

    run_original_star(
        &star_bin,
        &[
            generate_args
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>(),
            vec![
                "--genomeDir".to_string(),
                cpp_genome_dir.to_string_lossy().to_string(),
            ],
        ]
        .concat(),
    );
    run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                rust_genome_dir.to_string_lossy().to_string(),
            ],
            generate_args
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>(),
        ]
        .concat(),
    )
    .unwrap();

    let align_args = [
        "--readFilesIn",
        reads.to_str().unwrap(),
        "--outSAMtype",
        "SAM",
        "--seedSplitMin",
        "1",
        "--seedMapMin",
        "0",
        "--outFilterMatchNmin",
        "0",
        "--outFilterScoreMin",
        "0",
    ];

    run_original_star(
        &star_bin,
        &[
            vec![
                "--genomeDir".to_string(),
                cpp_genome_dir.to_string_lossy().to_string(),
                "--outFileNamePrefix".to_string(),
                cpp_prefix.to_string_lossy().to_string(),
            ],
            align_args.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
        ]
        .concat(),
    );

    let rust_result = run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                rust_genome_dir.to_string_lossy().to_string(),
                "--outFileNamePrefix".to_string(),
                rust_prefix.to_string_lossy().to_string(),
            ],
            align_args.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
        ]
        .concat(),
    )
    .unwrap();
    assert_eq!(rust_result.exit_code, 0);
    assert_eq!(rust_result.stats_all.read_n, 1);
    assert_eq!(rust_result.stats_all.mapped_reads_u, 1);
    assert_eq!(rust_result.stats_all.splices_n.iter().sum::<u32>(), 1);
    assert_eq!(
        rust_result.stats_all.splices_n[1] + rust_result.stats_all.splices_n[2],
        1
    );

    let cpp_sam = std::fs::read_to_string(dir.join("cpp/Aligned.out.sam")).unwrap();
    let rust_sam = std::fs::read_to_string(dir.join("rust/Aligned.out.sam")).unwrap();
    assert_eq!(sam_records(&cpp_sam), sam_records(&rust_sam));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn rust_matches_original_star_chimeric_split_read() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let star_bin = manifest_dir.join("STAR/bin/Linux_x86_64_static/STAR");
    if !star_bin.exists() {
        eprintln!(
            "skipping original STAR parity test: {} is absent",
            star_bin.display()
        );
        return;
    }

    let dir = unique_temp_dir("cpp_parity_chimeric");
    let cpp_genome_dir = dir.join("cpp_genome");
    let rust_genome_dir = dir.join("rust_genome");
    let cpp_prefix = dir.join("cpp/");
    let rust_prefix = dir.join("rust/");
    std::fs::create_dir_all(&cpp_genome_dir).unwrap();
    std::fs::create_dir_all(&rust_genome_dir).unwrap();

    let chr_a = "AAAGCGGCACTTGTGAAGTGTTCCCCACGCCGCTTGGGTCTTCTGTGTTGTTCGCGTGGTGCTGAGACAAAGCACGCCATAAGGCCAAAAAAAGGCCCATACCAAGAGGTAGTAGTCTCAGAATCTTGCGGGTACAGACCCATCACCTAGACGGTGACATTCAACAAACCACATTGTCCTTAATCATGAAGGGGATAAGCATATTTCAAGAGGACTCAGTTCGTAGAAAGTCAATATGGTCGGTTTTGTCCTGTAAAGCCTAAACGTCGTCGACTAGCGCCTCTGCTTATCTATGTGTTG";
    let chr_b = "GACCTTAGTTCAATCTCATCGCTCATTGCTCAGATATGTGTAAGCTGCACTTTGCAGTAGATTCGTCTGAGGGGGTACTCAGACTCGAAATGCGGAGTGCTTGTCTCGGCACTCGCGCCCGTTGGGTGAGGTTCGGTTACGTCAAGCGATAGCTGTCGGCTACCGGCTGGAGCCCAGGACCATTGCGAGTCATTTGATTTCTTTAATCACATGTAGAGCCACTAGTATCATCACAACAGCCGTACACATCACTGTCACCCTCGGTCTCTGGAATGGTGCTCAACCCTACAGTACCGACAC";
    let read = format!("{}{}", &chr_a[70..105], &chr_b[140..175]);
    let fasta = dir.join("chim.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, format!(">chrA\n{chr_a}\n>chrB\n{chr_b}\n")).unwrap();
    std::fs::write(
        &reads,
        format!("@chim1\n{read}\n+\n{}\n", "F".repeat(read.len())),
    )
    .unwrap();

    let generate_args = [
        "--runMode",
        "genomeGenerate",
        "--genomeFastaFiles",
        fasta.to_str().unwrap(),
        "--genomeSAindexNbases",
        "3",
        "--genomeChrBinNbits",
        "5",
        "--limitGenomeGenerateRAM",
        "1000000",
    ];

    run_original_star(
        &star_bin,
        &[
            generate_args
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>(),
            vec![
                "--genomeDir".to_string(),
                cpp_genome_dir.to_string_lossy().to_string(),
            ],
        ]
        .concat(),
    );
    run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                rust_genome_dir.to_string_lossy().to_string(),
            ],
            generate_args
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>(),
        ]
        .concat(),
    )
    .unwrap();

    let align_args = [
        "--readFilesIn",
        reads.to_str().unwrap(),
        "--outSAMtype",
        "SAM",
        "--seedSplitMin",
        "8",
        "--seedMapMin",
        "5",
        "--seedPerWindowNmax",
        "1000",
        "--outFilterMatchNmin",
        "0",
        "--outFilterScoreMin",
        "0",
        "--chimSegmentMin",
        "12",
        "--chimJunctionOverhangMin",
        "12",
        "--chimOutType",
        "Junctions",
        "SeparateSAMold",
        "--chimOutJunctionFormat",
        "1",
        "--chimScoreMin",
        "1",
        "--chimScoreDropMax",
        "1000",
        "--chimNonchimScoreDropMin",
        "0",
    ];

    run_original_star(
        &star_bin,
        &[
            vec![
                "--genomeDir".to_string(),
                cpp_genome_dir.to_string_lossy().to_string(),
                "--outFileNamePrefix".to_string(),
                cpp_prefix.to_string_lossy().to_string(),
            ],
            align_args.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
        ]
        .concat(),
    );

    let rust_result = run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                rust_genome_dir.to_string_lossy().to_string(),
                "--outFileNamePrefix".to_string(),
                rust_prefix.to_string_lossy().to_string(),
            ],
            align_args.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
        ]
        .concat(),
    )
    .unwrap();
    assert_eq!(rust_result.exit_code, 0);
    assert_eq!(rust_result.stats_all.read_n, 1);
    assert_eq!(rust_result.stats_all.chimeric_all, 1);

    let cpp_junction = std::fs::read_to_string(dir.join("cpp/Chimeric.out.junction")).unwrap();
    let rust_junction = std::fs::read_to_string(dir.join("rust/Chimeric.out.junction")).unwrap();
    assert_eq!(cpp_junction.lines().next(), rust_junction.lines().next());

    let cpp_sam = std::fs::read_to_string(dir.join("cpp/Chimeric.out.sam")).unwrap();
    let rust_sam = std::fs::read_to_string(dir.join("rust/Chimeric.out.sam")).unwrap();
    assert_eq!(sam_records(&cpp_sam), sam_records(&rust_sam));

    std::fs::remove_dir_all(dir).unwrap();
}
