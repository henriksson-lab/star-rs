use star_rs::cli::run_cli;
use star_rs::direct::{DirectReadPair, DirectStarRun};

#[test]
fn direct_star_run_maps_in_memory_read_pair() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_direct_smoke_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("direct/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();

    run_cli(&[
        "STAR".to_string(),
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
    ])
    .unwrap();

    let mut run = DirectStarRun::new(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_string_lossy().to_string(),
        "--readFilesIn".to_string(),
        "in-memory-r1.fq".to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_string_lossy().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
    ])
    .unwrap();
    run.clear_chunk_input();
    run.append_read_pair(&DirectReadPair {
        name: "r1",
        r1: b"ACGTACGT",
        q1: b"FFFFFFFF",
        r2: b"",
        q2: b"",
    });
    run.finalize_and_map_chunk().unwrap();

    let result = run.finish();
    assert_eq!(result.exit_code, 0);
    assert_eq!(result.process_chunks.len(), 1);
    assert_eq!(result.process_chunks[0].map_chunks.len(), 1);
    let sam = String::from_utf8_lossy(&result.process_chunks[0].map_chunks[0].direct_sam_output);
    assert!(sam.contains("r1\t0\tchr1\t1\t"), "{sam}");
    assert_eq!(result.stats_all.read_n, 1);
    assert_eq!(result.stats_all.mapped_reads_u, 1);

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn direct_star_run_maps_in_memory_paired_read() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_direct_paired_smoke_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("direct_paired/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();

    run_cli(&[
        "STAR".to_string(),
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
    ])
    .unwrap();

    let mut run = DirectStarRun::new(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_string_lossy().to_string(),
        "--readFilesIn".to_string(),
        "in-memory-r1.fq".to_string(),
        "in-memory-r2.fq".to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_string_lossy().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
    ])
    .unwrap();
    run.clear_chunk_input();
    run.append_read_pair(&DirectReadPair {
        name: "pair1",
        r1: b"ACGTACGT",
        q1: b"FFFFFFFF",
        r2: b"ACGTACGT",
        q2: b"FFFFFFFF",
    });
    run.finalize_and_map_chunk().unwrap();

    let result = run.finish();
    assert_eq!(result.exit_code, 0);
    assert_eq!(result.parameters.read_nends, 2);
    assert_eq!(result.parameters.read_files_in.len(), 2);
    assert_eq!(result.process_chunks.len(), 1);
    assert_eq!(result.process_chunks[0].map_chunks.len(), 1);
    let sam = String::from_utf8_lossy(&result.process_chunks[0].map_chunks[0].direct_sam_output);
    let pair_lines = sam
        .lines()
        .filter(|line| line.starts_with("pair1\t"))
        .collect::<Vec<_>>();
    assert_eq!(pair_lines.len(), 0, "{sam}");
    assert_eq!(result.stats_all.read_n, 1);
    assert_eq!(result.stats_all.mapped_reads_u, 0);
    assert_eq!(result.stats_all.mapped_reads_m, 0);
    assert_eq!(result.stats_all.unmapped_multi, 1);

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn direct_star_run_maps_nonrepetitive_paired_read_without_suffix_index_panic() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_direct_paired_unique_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("direct_paired_unique/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    std::fs::write(&fasta, ">chr1\nACGTACGTTGCAAGTC\n").unwrap();

    run_cli(&[
        "STAR".to_string(),
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
    ])
    .unwrap();

    let mut run = DirectStarRun::new(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_string_lossy().to_string(),
        "--readFilesIn".to_string(),
        "in-memory-r1.fq".to_string(),
        "in-memory-r2.fq".to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_string_lossy().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
    ])
    .unwrap();
    run.clear_chunk_input();
    run.append_read_pair(&DirectReadPair {
        name: "pair_unique",
        r1: b"ACGTACGT",
        q1: b"FFFFFFFF",
        r2: b"TGCAAGTC",
        q2: b"HHHHHHHH",
    });
    run.finalize_and_map_chunk().unwrap();

    let result = run.finish();
    assert_eq!(result.exit_code, 0);
    assert_eq!(result.stats_all.read_n, 1);
    assert_eq!(result.process_chunks.len(), 1);
    assert_eq!(result.process_chunks[0].map_chunks.len(), 1);

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn direct_star_run_matches_cli_for_paired_multimapping_records() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_direct_cli_parity_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let cli_prefix = dir.join("cli/");
    let direct_prefix = dir.join("direct/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads_1 = dir.join("reads_1.fq");
    let reads_2 = dir.join("reads_2.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads_1, "@pair1/1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();
    std::fs::write(&reads_2, "@pair1/2\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    run_cli(&[
        "STAR".to_string(),
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
    ])
    .unwrap();

    let cli_result = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_string_lossy().to_string(),
        "--readFilesIn".to_string(),
        reads_1.to_string_lossy().to_string(),
        reads_2.to_string_lossy().to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        cli_prefix.to_string_lossy().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
    ])
    .unwrap();

    let mut direct_run = DirectStarRun::new(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_string_lossy().to_string(),
        "--readFilesIn".to_string(),
        "in-memory-r1.fq".to_string(),
        "in-memory-r2.fq".to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        direct_prefix.to_string_lossy().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
    ])
    .unwrap();
    direct_run.clear_chunk_input();
    direct_run.append_read_pair(&DirectReadPair {
        name: "pair1",
        r1: b"ACGTACGT",
        q1: b"FFFFFFFF",
        r2: b"ACGTACGT",
        q2: b"FFFFFFFF",
    });
    direct_run.finalize_and_map_chunk().unwrap();
    let direct_result = direct_run.finish();

    assert_eq!(cli_result.exit_code, direct_result.exit_code);
    assert_eq!(cli_result.stats_all.read_n, direct_result.stats_all.read_n);
    assert_eq!(
        cli_result.stats_all.mapped_reads_u,
        direct_result.stats_all.mapped_reads_u
    );
    assert_eq!(
        cli_result.stats_all.mapped_reads_m,
        direct_result.stats_all.mapped_reads_m
    );

    let cli_sam = std::fs::read_to_string(dir.join("cli/Aligned.out.sam")).unwrap();
    let direct_sam =
        String::from_utf8_lossy(&direct_result.process_chunks[0].map_chunks[0].direct_sam_output);
    assert_eq!(sam_records(&cli_sam), sam_records(&direct_sam));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn direct_star_run_matches_cli_for_mixed_paired_fixture() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_direct_cli_mixed_parity_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let cli_prefix = dir.join("cli/");
    let direct_prefix = dir.join("direct/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("mixed.fa");
    let reads_1 = dir.join("reads_1.fq");
    let reads_2 = dir.join("reads_2.fq");
    std::fs::write(
        &fasta,
        ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTNNNNTTGGAACCTGACTGAC\n",
    )
    .unwrap();
    std::fs::write(
        &reads_1,
        "@multi/1\nACGTACGT\n+\nFFFFFFFF\n@unique/1\nTTGGAACC\n+\nHHHHHHHH\n@nohit/1\nCCCCCCCC\n+\nIIIIIIII\n",
    )
    .unwrap();
    std::fs::write(
        &reads_2,
        "@multi/2\nACGTACGT\n+\nFFFFFFFF\n@unique/2\nTGACTGAC\n+\nHHHHHHHH\n@nohit/2\nCCCCCCCC\n+\nIIIIIIII\n",
    )
    .unwrap();

    run_cli(&[
        "STAR".to_string(),
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
    ])
    .unwrap();

    let common_args = [
        "--outSAMtype",
        "SAM",
        "--seedSplitMin",
        "1",
        "--seedMapMin",
        "0",
        "--outFilterMatchNmin",
        "8",
        "--outFilterScoreMin",
        "0",
        "--outFilterMismatchNmax",
        "0",
    ];

    let cli_result = run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                genome_dir.to_string_lossy().to_string(),
                "--readFilesIn".to_string(),
                reads_1.to_string_lossy().to_string(),
                reads_2.to_string_lossy().to_string(),
                "--outFileNamePrefix".to_string(),
                cli_prefix.to_string_lossy().to_string(),
            ],
            common_args.iter().map(|s| s.to_string()).collect(),
        ]
        .concat(),
    )
    .unwrap();

    let mut direct_run = DirectStarRun::new(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                genome_dir.to_string_lossy().to_string(),
                "--readFilesIn".to_string(),
                "in-memory-r1.fq".to_string(),
                "in-memory-r2.fq".to_string(),
                "--outFileNamePrefix".to_string(),
                direct_prefix.to_string_lossy().to_string(),
            ],
            common_args.iter().map(|s| s.to_string()).collect(),
        ]
        .concat(),
    )
    .unwrap();
    direct_run.clear_chunk_input();
    for read in [
        DirectReadPair {
            name: "multi",
            r1: b"ACGTACGT",
            q1: b"FFFFFFFF",
            r2: b"ACGTACGT",
            q2: b"FFFFFFFF",
        },
        DirectReadPair {
            name: "unique",
            r1: b"TTGGAACC",
            q1: b"HHHHHHHH",
            r2: b"TGACTGAC",
            q2: b"HHHHHHHH",
        },
        DirectReadPair {
            name: "nohit",
            r1: b"CCCCCCCC",
            q1: b"IIIIIIII",
            r2: b"CCCCCCCC",
            q2: b"IIIIIIII",
        },
    ] {
        direct_run.append_read_pair(&read);
    }
    direct_run.finalize_and_map_chunk().unwrap();
    let direct_result = direct_run.finish();

    assert_eq!(cli_result.exit_code, direct_result.exit_code);
    assert_eq!(cli_result.stats_all.read_n, direct_result.stats_all.read_n);
    assert_eq!(
        cli_result.stats_all.mapped_reads_u,
        direct_result.stats_all.mapped_reads_u
    );
    assert_eq!(
        cli_result.stats_all.mapped_reads_m,
        direct_result.stats_all.mapped_reads_m
    );

    let cli_sam = std::fs::read_to_string(dir.join("cli/Aligned.out.sam")).unwrap();
    let direct_sam =
        String::from_utf8_lossy(&direct_result.process_chunks[0].map_chunks[0].direct_sam_output);
    let cli_records = sam_records(&cli_sam);
    let direct_records = sam_records(&direct_sam);
    assert_eq!(cli_records, direct_records);
    assert!(!cli_records.iter().any(|line| line.starts_with("multi\t")));
    assert!(!cli_records.iter().any(|line| line.starts_with("nohit\t")));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn direct_star_run_matches_cli_for_unmapped_within_sam() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_direct_cli_unmapped_within_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let cli_prefix = dir.join("cli/");
    let direct_prefix = dir.join("direct/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("mixed.fa");
    let reads_1 = dir.join("reads_1.fq");
    let reads_2 = dir.join("reads_2.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    std::fs::write(
        &reads_1,
        "@multi/1\nACGTACGT\n+\nFFFFFFFF\n@nohit/1\nTTTTTTTT\n+\nIIIIIIII\n",
    )
    .unwrap();
    std::fs::write(
        &reads_2,
        "@multi/2\nACGTACGT\n+\nFFFFFFFF\n@nohit/2\nTTTTTTTT\n+\nIIIIIIII\n",
    )
    .unwrap();

    run_cli(&[
        "STAR".to_string(),
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
    ])
    .unwrap();

    let common_args = [
        "--outSAMtype",
        "SAM",
        "--outSAMunmapped",
        "Within",
        "--seedSplitMin",
        "1",
        "--seedMapMin",
        "0",
        "--outFilterMatchNmin",
        "8",
        "--outFilterScoreMin",
        "0",
        "--outFilterMismatchNmax",
        "0",
    ];

    let cli_result = run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                genome_dir.to_string_lossy().to_string(),
                "--readFilesIn".to_string(),
                reads_1.to_string_lossy().to_string(),
                reads_2.to_string_lossy().to_string(),
                "--outFileNamePrefix".to_string(),
                cli_prefix.to_string_lossy().to_string(),
            ],
            common_args.iter().map(|s| s.to_string()).collect(),
        ]
        .concat(),
    )
    .unwrap();

    let mut direct_run = DirectStarRun::new(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                genome_dir.to_string_lossy().to_string(),
                "--readFilesIn".to_string(),
                "in-memory-r1.fq".to_string(),
                "in-memory-r2.fq".to_string(),
                "--outFileNamePrefix".to_string(),
                direct_prefix.to_string_lossy().to_string(),
            ],
            common_args.iter().map(|s| s.to_string()).collect(),
        ]
        .concat(),
    )
    .unwrap();
    direct_run.clear_chunk_input();
    for read in [
        DirectReadPair {
            name: "multi",
            r1: b"ACGTACGT",
            q1: b"FFFFFFFF",
            r2: b"ACGTACGT",
            q2: b"FFFFFFFF",
        },
        DirectReadPair {
            name: "nohit",
            r1: b"TTTTTTTT",
            q1: b"IIIIIIII",
            r2: b"TTTTTTTT",
            q2: b"IIIIIIII",
        },
    ] {
        direct_run.append_read_pair(&read);
    }
    direct_run.finalize_and_map_chunk().unwrap();
    let direct_result = direct_run.finish();

    assert_eq!(cli_result.exit_code, direct_result.exit_code);
    assert_eq!(cli_result.stats_all.read_n, direct_result.stats_all.read_n);
    assert_eq!(
        cli_result.stats_all.mapped_reads_u,
        direct_result.stats_all.mapped_reads_u
    );
    assert_eq!(
        cli_result.stats_all.mapped_reads_m,
        direct_result.stats_all.mapped_reads_m
    );
    assert_eq!(
        cli_result.stats_all.unmapped_other,
        direct_result.stats_all.unmapped_other
    );
    assert_eq!(
        cli_result.stats_all.unmapped_short,
        direct_result.stats_all.unmapped_short
    );

    let cli_sam = std::fs::read_to_string(dir.join("cli/Aligned.out.sam")).unwrap();
    let direct_sam =
        String::from_utf8_lossy(&direct_result.process_chunks[0].map_chunks[0].direct_sam_output);
    let cli_records = sam_records(&cli_sam);
    let direct_records = sam_records(&direct_sam);
    assert_eq!(cli_records, direct_records);
    assert!(cli_records.iter().any(|line| line.starts_with("multi\t")));
    assert!(cli_records.iter().any(|line| line.starts_with("nohit\t")));
    assert!(
        cli_records
            .iter()
            .any(|line| line.starts_with("nohit\t") && line.split('\t').nth(2) == Some("*"))
    );

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn direct_star_run_split_chunks_match_cli_record_order() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_direct_cli_split_chunks_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let cli_prefix = dir.join("cli/");
    let direct_prefix = dir.join("direct/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads_1 = dir.join("reads_1.fq");
    let reads_2 = dir.join("reads_2.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    std::fs::write(
        &reads_1,
        "@pair1/1\nACGTACGT\n+\nFFFFFFFF\n@pair2/1\nACGTACGT\n+\nHHHHHHHH\n",
    )
    .unwrap();
    std::fs::write(
        &reads_2,
        "@pair1/2\nACGTACGT\n+\nFFFFFFFF\n@pair2/2\nACGTACGT\n+\nHHHHHHHH\n",
    )
    .unwrap();

    run_cli(&[
        "STAR".to_string(),
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
    ])
    .unwrap();

    let common_args = [
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

    let cli_result = run_cli(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                genome_dir.to_string_lossy().to_string(),
                "--readFilesIn".to_string(),
                reads_1.to_string_lossy().to_string(),
                reads_2.to_string_lossy().to_string(),
                "--outFileNamePrefix".to_string(),
                cli_prefix.to_string_lossy().to_string(),
            ],
            common_args.iter().map(|s| s.to_string()).collect(),
        ]
        .concat(),
    )
    .unwrap();

    let mut direct_run = DirectStarRun::new(
        &[
            vec![
                "STAR".to_string(),
                "--genomeDir".to_string(),
                genome_dir.to_string_lossy().to_string(),
                "--readFilesIn".to_string(),
                "in-memory-r1.fq".to_string(),
                "in-memory-r2.fq".to_string(),
                "--outFileNamePrefix".to_string(),
                direct_prefix.to_string_lossy().to_string(),
            ],
            common_args.iter().map(|s| s.to_string()).collect(),
        ]
        .concat(),
    )
    .unwrap();

    direct_run.clear_chunk_input();
    direct_run.append_read_pair(&DirectReadPair {
        name: "pair1",
        r1: b"ACGTACGT",
        q1: b"FFFFFFFF",
        r2: b"ACGTACGT",
        q2: b"FFFFFFFF",
    });
    direct_run.finalize_and_map_chunk().unwrap();

    direct_run.clear_chunk_input();
    direct_run.append_read_pair(&DirectReadPair {
        name: "pair2",
        r1: b"ACGTACGT",
        q1: b"HHHHHHHH",
        r2: b"ACGTACGT",
        q2: b"HHHHHHHH",
    });
    direct_run.finalize_and_map_chunk().unwrap();

    let direct_result = direct_run.finish();
    assert_eq!(direct_result.process_chunks.len(), 1);
    assert_eq!(direct_result.process_chunks[0].chunks_read, 2);
    assert_eq!(direct_result.process_chunks[0].map_chunks.len(), 2);
    assert_eq!(cli_result.exit_code, direct_result.exit_code);
    assert_eq!(cli_result.stats_all.read_n, direct_result.stats_all.read_n);
    assert_eq!(
        cli_result.stats_all.mapped_reads_u,
        direct_result.stats_all.mapped_reads_u
    );
    assert_eq!(
        cli_result.stats_all.mapped_reads_m,
        direct_result.stats_all.mapped_reads_m
    );

    let cli_sam = std::fs::read_to_string(dir.join("cli/Aligned.out.sam")).unwrap();
    let direct_sam = direct_sam_output(&direct_result);
    let cli_records = sam_records(&cli_sam);
    let direct_records = sam_records(&direct_sam);
    assert_eq!(cli_records, direct_records);
    assert!(
        !direct_records
            .iter()
            .any(|line| line.starts_with("pair1\t"))
    );
    assert!(
        !direct_records
            .iter()
            .any(|line| line.starts_with("pair2\t"))
    );

    std::fs::remove_dir_all(dir).unwrap();
}

fn sam_records(sam: &str) -> Vec<&str> {
    sam.lines()
        .filter(|line| !line.trim_matches('\0').is_empty() && !line.starts_with('@'))
        .collect()
}

fn direct_sam_output(result: &star_rs::generated::structs::StarMainResult) -> String {
    let mut sam = String::new();
    for process in &result.process_chunks {
        for map_chunk in &process.map_chunks {
            sam.push_str(&String::from_utf8_lossy(&map_chunk.direct_sam_output));
        }
    }
    sam
}
