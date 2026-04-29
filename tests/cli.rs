use star_rs::cli::{existing_read_files_from_args, run_cli};
use star_rs::generated::structs::Parameters;

#[test]
fn cli_help_uses_translated_star_usage() {
    let result = run_cli(&["STAR".to_string(), "--help".to_string()]).unwrap();

    assert_eq!(result.exit_code, 0);
    assert!(result.usage.starts_with("Usage: STAR"));
    assert!(result.usage.contains("runMode"));
}

#[test]
fn cli_rejects_invalid_out_std_like_star() {
    let err = run_cli(&[
        "STAR".to_string(),
        "--outStd".to_string(),
        "BadMode".to_string(),
    ])
    .unwrap_err();

    assert!(err.contains("outStd=BadMode is not a valid value"));
}

#[test]
fn cli_rejects_out_sam_filter_without_added_references_like_star() {
    let err = run_cli(&[
        "STAR".to_string(),
        "--outSAMfilter".to_string(),
        "KeepOnlyAddedReferences".to_string(),
    ])
    .unwrap_err();

    assert!(err.contains("--outSAMfilter KeepOnlyAddedReferences OR KeepAllAddedReferences"));
    assert!(err.contains("--genomeFastaFiles"));
}

#[test]
fn cli_rejects_unknown_out_sam_filter_like_star() {
    let err = run_cli(&[
        "STAR".to_string(),
        "--outSAMfilter".to_string(),
        "BadFilter".to_string(),
    ])
    .unwrap_err();

    assert!(err.contains("unknown/unimplemented value for --outSAMfilter: BadFilter"));
}

#[test]
fn cli_rejects_unknown_out_multimapper_order_like_star() {
    let err = run_cli(&[
        "STAR".to_string(),
        "--outMultimapperOrder".to_string(),
        "BadOrder".to_string(),
    ])
    .unwrap_err();

    assert!(err.contains("unknown/unimplemented value for --outMultimapperOrder: BadOrder"));
}

#[test]
fn cli_rejects_paired_keep_input_order_without_sam_like_star() {
    let err = run_cli(&[
        "STAR".to_string(),
        "--outSAMtype".to_string(),
        "BAM".to_string(),
        "Unsorted".to_string(),
        "--outSAMorder".to_string(),
        "PairedKeepInputOrder".to_string(),
    ])
    .unwrap_err();

    assert!(err.contains("--outSAMorder=PairedKeepInputOrder"));
    assert!(err.contains("only compatible with SAM output"));
}

#[test]
fn cli_rejects_bysjout_with_paired_keep_input_order_like_star() {
    let err = run_cli(&[
        "STAR".to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outSAMorder".to_string(),
        "PairedKeepInputOrder".to_string(),
        "--outFilterType".to_string(),
        "BySJout".to_string(),
        "--outSJtype".to_string(),
        "Standard".to_string(),
    ])
    .unwrap_err();

    assert!(err.contains("--outFilterType=BySJout"));
    assert!(err.contains("--outSAMorder=PairedKeepInputOrder"));
}

#[test]
fn cli_rejects_twopass_reads_without_mode_like_star() {
    let err = run_cli(&[
        "STAR".to_string(),
        "--twopass1readsN".to_string(),
        "10".to_string(),
    ])
    .unwrap_err();

    assert!(err.contains("--twopass1readsN is defined"));
    assert!(err.contains("--twoPassMode is not defined"));
}

#[test]
fn cli_rejects_twopass_mode_outside_align_reads_like_star() {
    let err = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--twopassMode".to_string(),
        "Basic".to_string(),
    ])
    .unwrap_err();

    assert!(err.contains("2-pass mapping option"));
    assert!(err.contains("--runMode alignReads"));
}

#[test]
fn cli_rejects_unknown_twopass_mode_like_star() {
    let err = run_cli(&[
        "STAR".to_string(),
        "--twopassMode".to_string(),
        "BadMode".to_string(),
    ])
    .unwrap_err();

    assert!(err.contains("unrecognized value of --twopassMode=BadMode"));
}

#[test]
fn cli_rejects_zero_twopass_reads_like_star() {
    let err = run_cli(&[
        "STAR".to_string(),
        "--twopassMode".to_string(),
        "Basic".to_string(),
        "--twopass1readsN".to_string(),
        "0".to_string(),
    ])
    .unwrap_err();

    assert!(err.contains("--twopass1readsN = 0"));
}

#[test]
fn cli_rejects_twopass_with_shared_genome_like_star() {
    let err = run_cli(&[
        "STAR".to_string(),
        "--twopassMode".to_string(),
        "Basic".to_string(),
        "--genomeLoad".to_string(),
        "LoadAndKeep".to_string(),
    ])
    .unwrap_err();

    assert!(err.contains("2-pass method is not compatible with --genomeLoad LoadAndKeep"));
}

#[test]
fn cli_genome_generate_applies_star_style_arguments() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_tmp = dir.join("tmp");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();

    let result = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "0".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
        "--outTmpDir".to_string(),
        out_tmp.to_str().unwrap().to_string(),
    ])
    .unwrap();

    assert_eq!(result.exit_code, 0);
    assert!(result.log_main.contains("DONE: Genome generation, EXITING"));
    assert!(!result.parameters.out_sam_filter_yes);
    assert_eq!(
        std::fs::read_to_string(genome_dir.join("chrName.txt")).unwrap(),
        "chr1\n"
    );

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_accepts_out_sam_filter_with_added_references() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_out_sam_filter_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();

    let result = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "0".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
        "--outSAMfilter".to_string(),
        "KeepAllAddedReferences".to_string(),
    ])
    .unwrap();

    assert!(result.parameters.out_sam_filter_yes);
    assert!(result.parameters.out_sam_filter_keep_all_added_references);
    assert!(!result.parameters.out_sam_filter_keep_only_added_references);

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_reads_solo_whitelist_before_parameter_initialization() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_solo_wl_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let whitelist = dir.join("whitelist.txt");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&whitelist, "ACGTACGTACGTACGT\n").unwrap();

    let result = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
        "--soloType".to_string(),
        "CB_UMI_Simple".to_string(),
        "--soloCBwhitelist".to_string(),
        whitelist.to_str().unwrap().to_string(),
    ])
    .unwrap();

    assert_eq!(result.exit_code, 0);
    assert!(result.log_main.contains("DONE: Genome generation, EXITING"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_genome_generate_reads_sjdb_gtf_file() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_gtf_genome_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_tmp = dir.join("tmp");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let gtf = dir.join("genes.gtf");
    std::fs::write(&fasta, ">chr1\nACCCGTAGTTTTCCCC\n").unwrap();
    std::fs::write(
        &gtf,
        concat!(
            "chr1\tsrc\texon\t1\t4\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
            "chr1\tsrc\texon\t9\t12\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
        ),
    )
    .unwrap();

    let result = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
        "--limitSjdbInsertNsj".to_string(),
        "1000".to_string(),
        "--sjdbOverhang".to_string(),
        "1".to_string(),
        "--sjdbGTFfile".to_string(),
        gtf.to_str().unwrap().to_string(),
        "--outTmpDir".to_string(),
        out_tmp.to_str().unwrap().to_string(),
    ])
    .unwrap();

    assert_eq!(result.exit_code, 0);
    assert!(result.log_main.contains("processing annotations GTF"));
    assert!(result.log_main.contains("1 transcripts"));
    assert!(result.log_main.contains("2 exons"));
    assert!(result.log_main.contains("1 collapsed junctions"));
    assert_eq!(
        result
            .genome_generate
            .last()
            .and_then(|generation| generation.gtf.as_ref())
            .map(|gtf| gtf.sjdb_list_from_gtf_out_tab.as_str()),
        Some("chr1\t5\t8\t+\t1\n")
    );
    assert_eq!(
        std::fs::read_to_string(genome_dir.join("geneInfo.tab")).unwrap(),
        "1\nG1\tG1\tMissingGeneType\n"
    );
    assert_eq!(
        std::fs::read_to_string(genome_dir.join("exonGeTrInfo.tab")).unwrap(),
        "2\n0\t3\t1\t0\t0\n8\t11\t1\t0\t0\n"
    );

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_loads_genome_and_processes_fastq_read() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_align_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("align_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
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

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned.log_final_out.contains("Number of input reads |\t1"));
    assert!(
        aligned
            .log_final_out
            .contains("Uniquely mapped reads number |\t1")
    );
    assert!(aligned.log_main.contains("Completed: thread #0"));
    let sam = std::fs::read_to_string(dir.join("align_out/Aligned.out.sam")).unwrap();
    assert!(sam.starts_with("@HD\tVN:1.4\n"));
    assert!(sam.contains("@SQ\tSN:chr1\tLN:16\n"));
    assert!(sam.contains("r1\t0\tchr1\t1\t"));
    assert!(
        std::fs::read_to_string(dir.join("align_out/Log.final.out"))
            .unwrap()
            .contains("Number of input reads |\t1")
    );

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_applies_multimapper_output_controls() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_multimapper_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("multi_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("duplicate.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n>chr2\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
        "--outMultimapperOrder".to_string(),
        "Random".to_string(),
        "--outSAMprimaryFlag".to_string(),
        "AllBestScore".to_string(),
    ])
    .unwrap();

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned.parameters.out_multimapper_order_random);
    let sam = std::fs::read_to_string(dir.join("multi_out/Aligned.out.sam")).unwrap();
    let records: Vec<&str> = sam
        .lines()
        .filter(|line| line.starts_with("r1\t"))
        .collect();
    assert!(records.len() > 1);
    assert!(records.iter().all(|record| {
        record
            .split('\t')
            .nth(1)
            .and_then(|flag| flag.parse::<u16>().ok())
            .map(|flag| flag & 0x100 == 0)
            .unwrap_or(false)
    }));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_twopass_basic_initializes_pass1_outputs() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_twopass_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("twopass_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
        "--twopassMode".to_string(),
        "Basic".to_string(),
        "--twopass1readsN".to_string(),
        "1".to_string(),
    ])
    .unwrap();

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned.parameters.two_pass_yes);
    assert!(aligned.parameters.two_pass_pass2);
    assert_eq!(
        aligned.parameters.two_pass_dir,
        format!("{}_STARpass1/", out_prefix.to_str().unwrap())
    );
    assert!(std::path::Path::new(&aligned.parameters.two_pass_pass1sj_file).exists());

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_marks_limited_multimapper_output() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_multimapper_limit_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("multi_limit_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("duplicate.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n>chr2\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
        "--outSAMmultNmax".to_string(),
        "1".to_string(),
    ])
    .unwrap();

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned.parameters.out_sam_mult_nmax_is_limited);
    let sam = std::fs::read_to_string(dir.join("multi_limit_out/Aligned.out.sam")).unwrap();
    let records: Vec<&str> = sam
        .lines()
        .filter(|line| line.starts_with("r1\t"))
        .collect();
    assert_eq!(records.len(), 1);

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_writes_chimeric_output_files_when_requested() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_chimeric_out_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("chim_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
        "--chimSegmentMin".to_string(),
        "1".to_string(),
        "--chimOutType".to_string(),
        "Junctions".to_string(),
        "SeparateSAMold".to_string(),
        "--chimOutJunctionFormat".to_string(),
        "1".to_string(),
    ])
    .unwrap();

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned.parameters.p_ch.out_chim_junction_opened);
    assert!(aligned.parameters.p_ch.out_chim_sam_opened);
    assert!(dir.join("chim_out/Chimeric.out.junction").exists());
    let chim_junction =
        std::fs::read_to_string(dir.join("chim_out/Chimeric.out.junction")).unwrap();
    assert!(chim_junction.contains("# Nreads "));
    let chim_sam = std::fs::read_to_string(dir.join("chim_out/Chimeric.out.sam")).unwrap();
    assert!(chim_sam.starts_with("@HD\tVN:1.4\n"));
    assert!(chim_sam.contains("@SQ\tSN:chr1\tLN:16\n"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_accepts_chimeric_within_bam_output() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_chimeric_bam_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("chim_bam_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "BAM".to_string(),
        "Unsorted".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
        "--chimSegmentMin".to_string(),
        "1".to_string(),
        "--chimOutType".to_string(),
        "WithinBAM".to_string(),
    ])
    .unwrap();

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned.parameters.p_ch.out_bam);
    let bam = std::fs::read(dir.join("chim_bam_out/Aligned.out.bam")).unwrap();
    assert!(bam.starts_with(b"BAM\x01"));
    assert!(bam.len() > aligned.parameters.out_bam_unsorted_header.len());

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_applies_custom_sam_headers() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_sam_header_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("header_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    let comments = dir.join("comments.sam");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();
    std::fs::write(&comments, "\n@CO\tfrom-file\n   \n@CO\tsecond\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outSAMheaderHD".to_string(),
        "@HD".to_string(),
        "VN:1.6".to_string(),
        "--outSAMheaderPG".to_string(),
        "@PG".to_string(),
        "ID:custom".to_string(),
        "PN:tool".to_string(),
        "--outSAMheaderCommentFile".to_string(),
        comments.to_str().unwrap().to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
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

    assert_eq!(aligned.exit_code, 0);
    let sam = std::fs::read_to_string(dir.join("header_out/Aligned.out.sam")).unwrap();
    assert!(sam.starts_with("@HD\tVN:1.6\n"));
    assert!(sam.contains("@PG\tID:custom\tPN:tool\n"));
    assert!(sam.contains("@CO\tfrom-file\n"));
    assert!(sam.contains("@CO\tsecond\n"));
    assert!(!sam.contains("\n   \n"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_out_std_sam_keeps_alignment_on_stdout_path() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_out_std_sam_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("stdout_sam_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outStd".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
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

    assert_eq!(aligned.exit_code, 0);
    assert_eq!(aligned.parameters.out_std, "SAM");
    assert!(
        aligned
            .parameters
            .out_sam_contents
            .starts_with("@HD\tVN:1.4\n")
    );
    assert!(
        aligned
            .parameters
            .out_sam_contents
            .contains("r1\t0\tchr1\t1\t")
    );
    assert!(!dir.join("stdout_sam_out/Aligned.out.sam").exists());
    assert!(dir.join("stdout_sam_out/Log.std.out").exists());

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_writes_unmapped_fastx_output() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_unmapped_fastx_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("unmapped_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r_un\nTTTTTTTT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outReadsUnmapped".to_string(),
        "Fastx".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "4".to_string(),
    ])
    .unwrap();

    assert_eq!(aligned.exit_code, 0);
    assert_eq!(aligned.parameters.out_reads_unmapped, "Fastx");
    let unmapped = std::fs::read_to_string(dir.join("unmapped_out/Unmapped.out.mate1")).unwrap();
    assert!(unmapped.starts_with("@r_un 0:"));
    assert!(unmapped.contains("\nTTTTTTTT\n+\nFFFFFFFF\n"));
    assert!(
        aligned
            .log_final_out
            .contains("Number of reads unmapped: too short |\t1")
    );

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_writes_unsorted_bam_output() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_bam_align_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("bam_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "BAM".to_string(),
        "Unsorted".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
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

    assert_eq!(aligned.exit_code, 0);
    let bam = std::fs::read(dir.join("bam_out/Aligned.out.bam")).unwrap();
    assert!(bam.starts_with(b"BAM\x01"));
    assert!(bam.len() > aligned.parameters.out_bam_unsorted_header.len());
    assert!(aligned.log_final_out.contains("Number of input reads |\t1"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_writes_sorted_bam_output() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_sorted_bam_align_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("sorted_bam_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "BAM".to_string(),
        "SortedByCoordinate".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
        "--limitBAMsortRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();

    assert_eq!(aligned.exit_code, 0);
    let bam = std::fs::read(dir.join("sorted_bam_out/Aligned.sortedByCoord.out.bam")).unwrap();
    assert!(bam.starts_with(b"BAM\x01"));
    assert!(bam.len() > aligned.parameters.out_bam_unsorted_header.len());
    assert!(
        aligned
            .bam_sort
            .as_ref()
            .map(|sort| !sort.output_bam.is_empty())
            .unwrap_or(false)
    );
    assert!(aligned.log_final_out.contains("Number of input reads |\t1"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_writes_wiggle_signal_from_sorted_bam_records() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_signal_align_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("signal_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "BAM".to_string(),
        "SortedByCoordinate".to_string(),
        "--outWigType".to_string(),
        "bedGraph".to_string(),
        "--outWigReferencesPrefix".to_string(),
        "chr".to_string(),
        "--outWigStrand".to_string(),
        "Unstranded".to_string(),
        "--outWigNorm".to_string(),
        "None".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
        "--limitBAMsortRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned.signal.as_ref().unwrap().n_unique == 0.0);
    let unique_signal =
        std::fs::read_to_string(dir.join("signal_out/Signal.Unique.str1.out.bg")).unwrap();
    let unique_multiple_signal =
        std::fs::read_to_string(dir.join("signal_out/Signal.UniqueMultiple.str1.out.bg")).unwrap();
    assert_eq!(unique_signal, "chr1\t0\t8\t1\n");
    assert_eq!(unique_multiple_signal, "chr1\t0\t8\t1\n");

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_input_alignments_from_bam_writes_signal_from_existing_bam() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_input_bam_signal_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let align_prefix = dir.join("align/");
    let signal_prefix = dir.join("signal/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "BAM".to_string(),
        "Unsorted".to_string(),
        "--outFileNamePrefix".to_string(),
        align_prefix.to_str().unwrap().to_string(),
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
    assert_eq!(aligned.exit_code, 0);
    let bam_path = dir.join("align/Aligned.out.bam");
    assert!(bam_path.exists());
    let bgzf_bam_path = dir.join("align/Aligned.bgzf.bam");
    let raw_bam = std::fs::read(&bam_path).unwrap();
    let mut encoder =
        flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
    std::io::Write::write_all(&mut encoder, &raw_bam).unwrap();
    let compressed_bam = encoder.finish().unwrap();
    let block_size = 18 + compressed_bam.len() + 8;
    let mut bgzf_bam = Vec::new();
    bgzf_bam.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 255]);
    bgzf_bam.extend_from_slice(&6u16.to_le_bytes());
    bgzf_bam.extend_from_slice(b"BC");
    bgzf_bam.extend_from_slice(&2u16.to_le_bytes());
    bgzf_bam.extend_from_slice(&((block_size - 1) as u16).to_le_bytes());
    bgzf_bam.extend_from_slice(&compressed_bam);
    bgzf_bam.extend_from_slice(&0u32.to_le_bytes());
    bgzf_bam.extend_from_slice(&(raw_bam.len() as u32).to_le_bytes());
    std::fs::write(&bgzf_bam_path, &bgzf_bam).unwrap();

    let signal = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "inputAlignmentsFromBAM".to_string(),
        "--inputBAMfile".to_string(),
        bam_path.to_str().unwrap().to_string(),
        "--outWigType".to_string(),
        "bedGraph".to_string(),
        "--outFileNamePrefix".to_string(),
        signal_prefix.to_str().unwrap().to_string(),
    ])
    .unwrap();
    assert_eq!(signal.exit_code, 0);
    assert!(signal.signal.is_some());
    let signal_file = dir.join("signal/Signal.UniqueMultiple.str1.out.bg");
    assert!(signal_file.exists());
    let signal_contents = std::fs::read_to_string(signal_file).unwrap();
    assert!(signal_contents.contains("chr1\t0\t8\t1"));

    let signal_bgzf_prefix = dir.join("signal_bgzf/");
    let signal_bgzf = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "inputAlignmentsFromBAM".to_string(),
        "--inputBAMfile".to_string(),
        bgzf_bam_path.to_str().unwrap().to_string(),
        "--outWigType".to_string(),
        "bedGraph".to_string(),
        "--outFileNamePrefix".to_string(),
        signal_bgzf_prefix.to_str().unwrap().to_string(),
    ])
    .unwrap();
    assert_eq!(signal_bgzf.exit_code, 0);
    let signal_bgzf_file = dir.join("signal_bgzf/Signal.UniqueMultiple.str1.out.bg");
    let signal_bgzf_contents = std::fs::read_to_string(signal_bgzf_file).unwrap();
    assert!(signal_bgzf_contents.contains("chr1\t0\t8\t1"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_input_alignments_from_bam_marks_duplicates_and_writes_processed_bam() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_input_bam_dedup_test_{}_{}",
        std::process::id(),
        unique
    ));
    let out_prefix = dir.join("dedup/");
    std::fs::create_dir_all(&dir).unwrap();
    let bam_path = dir.join("duplicates.bam");

    let mut bam = Vec::new();
    bam.extend_from_slice(b"BAM\x01");
    let sam_header = "@SQ\tSN:chr1\tLN:1000\n";
    bam.extend_from_slice(&(sam_header.len() as i32).to_ne_bytes());
    bam.extend_from_slice(sam_header.as_bytes());
    bam.extend_from_slice(&1i32.to_ne_bytes());
    bam.extend_from_slice(&5i32.to_ne_bytes());
    bam.extend_from_slice(b"chr1\0");
    bam.extend_from_slice(&1000i32.to_ne_bytes());

    for (name, pos, flag, as_score) in [
        ("readA", 100u32, 0u32, 20i32),
        ("readA", 300u32, 0x80u32, 20i32),
        ("readB", 100u32, 0u32, 25i32),
        ("readB", 300u32, 0x80u32, 25i32),
    ] {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NHi");
        aux.extend_from_slice(&1i32.to_ne_bytes());
        aux.extend_from_slice(b"ASi");
        aux.extend_from_slice(&as_score.to_ne_bytes());
        let qname = format!("{}\0", name);
        let cigar = 6u32 << 4;
        let seq = [0x12u8, 0x34, 0x56];
        let qual = [30u8; 6];
        let block_len =
            32 + qname.len() + std::mem::size_of::<u32>() + seq.len() + qual.len() + aux.len();
        bam.extend_from_slice(&(block_len as i32).to_ne_bytes());
        bam.extend_from_slice(&0i32.to_ne_bytes());
        bam.extend_from_slice(&(pos as i32).to_ne_bytes());
        bam.extend_from_slice(&(((255u32) << 8) | qname.len() as u32).to_ne_bytes());
        bam.extend_from_slice(&((flag << 16) | 1u32).to_ne_bytes());
        bam.extend_from_slice(&6i32.to_ne_bytes());
        bam.extend_from_slice(&0i32.to_ne_bytes());
        bam.extend_from_slice(&0i32.to_ne_bytes());
        bam.extend_from_slice(&0i32.to_ne_bytes());
        bam.extend_from_slice(qname.as_bytes());
        bam.extend_from_slice(&cigar.to_ne_bytes());
        bam.extend_from_slice(&seq);
        bam.extend_from_slice(&qual);
        bam.extend_from_slice(&aux);
    }
    std::fs::write(&bam_path, &bam).unwrap();

    let result = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "inputAlignmentsFromBAM".to_string(),
        "--inputBAMfile".to_string(),
        bam_path.to_str().unwrap().to_string(),
        "--bamRemoveDuplicatesType".to_string(),
        "UniqueIdentical".to_string(),
        "--bamRemoveDuplicatesMate2basesN".to_string(),
        "6".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
    ])
    .unwrap();
    assert_eq!(result.exit_code, 0);
    assert!(!result.processed_bam_output.is_empty());
    let processed_path = dir.join("dedup/Processed.out.bam");
    assert!(processed_path.exists());
    let processed = std::fs::read(processed_path).unwrap();

    let mut pos = 4usize;
    let header_len = i32::from_ne_bytes(processed[pos..pos + 4].try_into().unwrap()) as usize;
    pos += 4 + header_len;
    let ref_n = i32::from_ne_bytes(processed[pos..pos + 4].try_into().unwrap());
    pos += 4;
    for _ in 0..ref_n {
        let name_len = i32::from_ne_bytes(processed[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4 + name_len + 4;
    }
    let mut flags = Vec::new();
    while pos < processed.len() {
        let block_len = i32::from_ne_bytes(processed[pos..pos + 4].try_into().unwrap()) as usize;
        let flag_nc = u32::from_ne_bytes(processed[pos + 16..pos + 20].try_into().unwrap());
        flags.push(flag_nc >> 16);
        pos += 4 + block_len;
    }
    assert_eq!(flags.len(), 4);
    assert_ne!(flags[0] & 0x400, 0);
    assert_ne!(flags[1] & 0x400, 0);
    assert_eq!(flags[2] & 0x400, 0);
    assert_eq!(flags[3] & 0x400, 0);

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_quant_mode_gene_counts_writes_reads_per_gene() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_gene_counts_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("counts_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let gtf = dir.join("genes.gtf");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(
        &gtf,
        "chr1\tsrc\texon\t1\t8\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
    )
    .unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
        "--limitSjdbInsertNsj".to_string(),
        "1000".to_string(),
        "--sjdbOverhang".to_string(),
        "1".to_string(),
        "--sjdbGTFfile".to_string(),
        gtf.to_str().unwrap().to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--quantMode".to_string(),
        "GeneCounts".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
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

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned.parameters.quant_ge_count_yes);
    assert_eq!(
        std::fs::read_to_string(dir.join("counts_out/ReadsPerGene.out.tab")).unwrap(),
        concat!(
            "N_unmapped\t0\t0\t0\n",
            "N_multimapping\t0\t0\t0\n",
            "N_noFeature\t0\t0\t1\n",
            "N_ambiguous\t0\t0\t0\n",
            "G1\t1\t1\t0\n",
        )
    );

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_quant_mode_transcriptome_sam_writes_bam() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_transcriptome_bam_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("quant_bam_out/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let gtf = dir.join("genes.gtf");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(
        &gtf,
        "chr1\tsrc\texon\t1\t8\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
    )
    .unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
        "--limitSjdbInsertNsj".to_string(),
        "1000".to_string(),
        "--sjdbOverhang".to_string(),
        "1".to_string(),
        "--sjdbGTFfile".to_string(),
        gtf.to_str().unwrap().to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--quantMode".to_string(),
        "TranscriptomeSAM".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
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

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned.parameters.quant_tr_sam_bam_yes);
    let bam = std::fs::read(dir.join("quant_bam_out/Aligned.toTranscriptome.out.bam")).unwrap();
    assert!(bam.starts_with(b"BAM\x01"));
    assert!(bam.len() > aligned.parameters.out_quant_bam_header.len());

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_uses_read_files_command_for_gzip_fastq() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_gzip_align_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads_plain = dir.join("reads.fq");
    let reads_gz = dir.join("reads.fq.gz");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads_plain, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();
    let gzip_status = std::process::Command::new("gzip")
        .arg("-c")
        .arg(&reads_plain)
        .output()
        .unwrap();
    assert!(gzip_status.status.success());
    std::fs::write(&reads_gz, gzip_status.stdout).unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads_gz.to_str().unwrap().to_string(),
        "--readFilesCommand".to_string(),
        "gzip".to_string(),
        "-cd".to_string(),
        "--outSAMtype".to_string(),
        "None".to_string(),
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

    assert_eq!(aligned.exit_code, 0);
    assert_eq!(aligned.parameters.read_files_command, vec!["gzip", "-cd"]);
    assert!(aligned.log_final_out.contains("Number of input reads |\t1"));
    assert!(
        aligned
            .log_final_out
            .contains("Uniquely mapped reads number |\t1")
    );

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_accepts_paired_read_files_as_star_vector_argument() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_paired_align_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads_1 = dir.join("reads_1.fq");
    let reads_2 = dir.join("reads_2.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads_1, "@r1/1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();
    std::fs::write(&reads_2, "@r1/2\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads_1.to_str().unwrap().to_string(),
        reads_2.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "None".to_string(),
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

    assert_eq!(aligned.exit_code, 0);
    assert_eq!(
        aligned.parameters.read_files_in,
        vec![
            reads_1.to_str().unwrap().to_string(),
            reads_2.to_str().unwrap().to_string()
        ]
    );
    assert!(aligned.log_final_out.contains("Number of input reads |\t1"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_applies_read_limit_numeric_ids_flags_and_quality_shift() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_sam_controls_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("sam_controls/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(
        &reads,
        "@r1\nACGTACGT\n+\nFFFFFFFF\n@r2\nTACGTACG\n+\nFFFFFFFF\n",
    )
    .unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
        "--seedSplitMin".to_string(),
        "1".to_string(),
        "--seedMapMin".to_string(),
        "0".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
        "--readMapNumber".to_string(),
        "1".to_string(),
        "--outSAMreadID".to_string(),
        "Number".to_string(),
        "--outSAMflagOR".to_string(),
        "256".to_string(),
        "--outSAMflagAND".to_string(),
        "65535".to_string(),
        "--outQSconversionAdd".to_string(),
        "2".to_string(),
    ])
    .unwrap();

    assert_eq!(aligned.exit_code, 0);
    assert_eq!(aligned.parameters.i_read_all, 1);
    let sam = std::fs::read_to_string(dir.join("sam_controls/Aligned.out.sam")).unwrap();
    let records: Vec<&str> = sam.lines().filter(|line| !line.starts_with('@')).collect();
    assert_eq!(records.len(), 1);
    let fields: Vec<&str> = records[0].split('\t').collect();
    assert_eq!(fields[0], "1");
    assert_ne!(fields[1].parse::<u16>().unwrap() & 256, 0);
    assert_eq!(fields[10], "HHHHHHHH");
    assert!(!sam.contains("r2"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_applies_out_sam_mode_noqs_and_none() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_out_sam_mode_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let noqs_prefix = dir.join("noqs/");
    let none_prefix = dir.join("none/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let noqs = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outSAMmode".to_string(),
        "NoQS".to_string(),
        "--outFileNamePrefix".to_string(),
        noqs_prefix.to_str().unwrap().to_string(),
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
    assert_eq!(noqs.exit_code, 0);
    let noqs_sam = std::fs::read_to_string(dir.join("noqs/Aligned.out.sam")).unwrap();
    let noqs_record = noqs_sam
        .lines()
        .find(|line| !line.starts_with('@'))
        .unwrap();
    let noqs_fields: Vec<&str> = noqs_record.split('\t').collect();
    assert_eq!(noqs_fields[10], "*");

    let none = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMmode".to_string(),
        "None".to_string(),
        "--outFileNamePrefix".to_string(),
        none_prefix.to_str().unwrap().to_string(),
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
    assert_eq!(none.exit_code, 0);
    assert!(!dir.join("none/Aligned.out.sam").exists());

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_accepts_comma_separated_read_file_chunks() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_multifile_align_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads_1 = dir.join("reads_1.fq");
    let reads_2 = dir.join("reads_2.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads_1, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();
    std::fs::write(&reads_2, "@r2\nACGTACGT\n+\nFFFFFFFF\n").unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        format!(
            "{},{}",
            reads_1.to_str().unwrap(),
            reads_2.to_str().unwrap()
        ),
        "--outSAMtype".to_string(),
        "None".to_string(),
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

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned.log_final_out.contains("Number of input reads |\t2"));
    assert!(aligned.log_main.contains("Starting to map file # 0"));
    assert!(aligned.log_main.contains("Starting to map file # 1"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_accepts_read_files_manifest() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_manifest_align_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    let manifest = dir.join("manifest.tsv");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();
    std::fs::write(
        &manifest,
        format!("{}\t-\tID:rg1\n", reads.to_str().unwrap()),
    )
    .unwrap();

    let generated = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "genomeGenerate".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--genomeFastaFiles".to_string(),
        fasta.to_str().unwrap().to_string(),
        "--genomeSAindexNbases".to_string(),
        "1".to_string(),
        "--genomeChrBinNbits".to_string(),
        "2".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000".to_string(),
    ])
    .unwrap();
    assert_eq!(generated.exit_code, 0);

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesManifest".to_string(),
        manifest.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "None".to_string(),
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

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned.log_final_out.contains("Number of input reads |\t1"));
    assert_eq!(
        aligned.parameters.read_files_names[0][0],
        reads.to_str().unwrap()
    );
    assert_eq!(aligned.parameters.out_sam_attr_rg[0], "rg1");

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_existing_read_files_use_finalized_read_files_prefix() {
    let args = vec![
        "STAR".to_string(),
        "--readFilesIn".to_string(),
        "r1a.fq,r1b.fq,".to_string(),
        "r2a.fq,r2b.fq".to_string(),
    ];
    let parameters = Parameters {
        read_files_prefix_final: "data/".to_string(),
        read_files_in: vec!["r1a.fq,r1b.fq,".to_string(), "r2a.fq,r2b.fq".to_string()],
        ..Default::default()
    };

    let files = existing_read_files_from_args(&args, Some(&parameters));

    assert!(files.contains("data/r1a.fq"));
    assert!(files.contains("data/r1b.fq"));
    assert!(files.contains("data/r2a.fq"));
    assert!(files.contains("data/r2b.fq"));
    assert!(!files.contains("r1a.fq"));
}
