use star_rs::cli::{existing_read_files_from_args, run_cli};
use star_rs::generated::functions::{
    quantifications_l3_quantifications_quantifications, star_l58_main, SOLO_TYPE_CB_SAM_TAG_OUT,
};
use star_rs::generated::structs::{
    Genome, Parameters, ParametersSolo, ReadAlign, ReadAlignChunk, SoloRead, SoloReadBarcode,
    SoloReadBarcodeStats, Transcriptome,
};
use std::collections::BTreeSet;
use std::io::Read;
use std::path::PathBuf;

const UNIQUE_TEST_GENOME: &str = ">chr1\nACGTTGCAAGTCCTGA\n";
const UNIQUE_TEST_READ: &str = "@r1\nACGTTGCA\n+\nFFFFFFFF\n";

fn unique_temp_dir(label: &str) -> PathBuf {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("star_rs_{label}_{}_{}", std::process::id(), unique))
}

fn bam_payload(bytes: &[u8]) -> Vec<u8> {
    if bytes.starts_with(b"BAM\x01") {
        bytes.to_vec()
    } else {
        let mut out = Vec::new();
        let mut decoder = flate2::read::MultiGzDecoder::new(bytes);
        decoder.read_to_end(&mut out).unwrap();
        out
    }
}

fn gzip_bytes(bytes: &[u8]) -> Vec<u8> {
    use std::io::Write;

    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

#[test]
fn cli_help_uses_translated_star_usage() {
    let result = run_cli(&["STAR".to_string(), "--help".to_string()]).unwrap();

    assert_eq!(result.exit_code, 0);
    assert!(result.usage.starts_with("Usage: STAR"));
    assert!(result.usage.contains("runMode"));
}

#[test]
fn cli_malformed_fastq_quality_length_returns_error_without_panic() {
    let dir = unique_temp_dir("malformed_fastq_no_panic");
    let genome_dir = dir.join("genome");
    let out_dir = dir.join("align/");
    std::fs::create_dir_all(&genome_dir).unwrap();

    let fasta = dir.join("genome.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, "@bad\nACGTACGT\n+\nFFFF\n").unwrap();

    run_cli(&[
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

    let result = std::panic::catch_unwind(|| {
        run_cli(&[
            "STAR".to_string(),
            "--genomeDir".to_string(),
            genome_dir.to_str().unwrap().to_string(),
            "--readFilesIn".to_string(),
            reads.to_str().unwrap().to_string(),
            "--outFileNamePrefix".to_string(),
            out_dir.to_str().unwrap().to_string(),
        ])
    });
    let err = result.expect("malformed FASTQ must return an error, not panic");
    assert!(err
        .unwrap_err()
        .contains("quality string length is not equal to sequence length"));

    std::fs::remove_dir_all(dir).unwrap();
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
fn cli_liftover_writes_lifted_gtf_for_each_chain_file() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_liftover_test_{}_{}",
        std::process::id(),
        unique
    ));
    let out_prefix = dir.join("lift/");
    std::fs::create_dir_all(&dir).unwrap();
    let chain1 = dir.join("one.chain");
    let chain2 = dir.join("two.chain");
    let gtf = dir.join("input.gtf");
    let chain_contents = "chain 100 chr1 1000 + 10 80 chrA 2000 + 20 90 1\n\
                          25 3 4\n\
                          30 5 6\n\
                          12\n";
    std::fs::write(&chain1, chain_contents).unwrap();
    std::fs::write(&chain2, chain_contents.replace("chrA", "chrB")).unwrap();
    std::fs::write(
        &gtf,
        "chr1\tsrc\texon\t12\t20\t.\t+\t.\tgene_id \"a\";\n\
         chr1\tsrc\texon\t90\t95\t.\t+\t.\tgene_id \"bad\";\n",
    )
    .unwrap();

    let result = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "liftOver".to_string(),
        "--genomeChainFiles".to_string(),
        chain1.to_str().unwrap().to_string(),
        chain2.to_str().unwrap().to_string(),
        "--sjdbGTFfile".to_string(),
        gtf.to_str().unwrap().to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
    ])
    .unwrap();

    assert_eq!(result.exit_code, 0);
    assert!(result
        .log_main
        .contains("DONE: lift-over of GTF file, EXITING"));
    assert_eq!(
        std::fs::read_to_string(dir.join("lift/GTFliftOver_1.gtf")).unwrap(),
        "chrA\tsrc\texon\t22\t30\t.\t+\t.\tgene_id \"a\";\n"
    );
    assert_eq!(
        std::fs::read_to_string(dir.join("lift/GTFliftOver_2.gtf")).unwrap(),
        "chrB\tsrc\texon\t22\t30\t.\t+\t.\tgene_id \"a\";\n"
    );
    assert_eq!(
        std::fs::read_to_string(dir.join("lift/GTFliftOver_1.gtf.unlifted")).unwrap(),
        "chr1\tsrc\texon\t90\t95\t.\t+\t.\tgene_id \"bad\";\n"
    );

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_liftover_malformed_chain_returns_error_without_panic() {
    let dir = unique_temp_dir("liftover_bad_chain_no_panic");
    let out_prefix = dir.join("lift/");
    std::fs::create_dir_all(&dir).unwrap();
    let chain = dir.join("bad.chain");
    let gtf = dir.join("input.gtf");
    std::fs::write(
        &chain,
        "chain 100 chr1 1000 + bad 80 chrA 2000 + 20 90 1\n25\n",
    )
    .unwrap();
    std::fs::write(&gtf, "chr1\tsrc\texon\t12\t20\t.\t+\t.\tgene_id \"a\";\n").unwrap();

    let result = std::panic::catch_unwind(|| {
        run_cli(&[
            "STAR".to_string(),
            "--runMode".to_string(),
            "liftOver".to_string(),
            "--genomeChainFiles".to_string(),
            chain.to_str().unwrap().to_string(),
            "--sjdbGTFfile".to_string(),
            gtf.to_str().unwrap().to_string(),
            "--outFileNamePrefix".to_string(),
            out_prefix.to_str().unwrap().to_string(),
        ])
    });
    let err = result.expect("malformed chain must return an error, not panic");
    assert!(err.unwrap_err().contains("invalid source start bad"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_liftover_malformed_gtf_coordinate_returns_error_without_panic() {
    let dir = unique_temp_dir("liftover_bad_gtf_no_panic");
    let out_prefix = dir.join("lift/");
    std::fs::create_dir_all(&dir).unwrap();
    let chain = dir.join("one.chain");
    let gtf = dir.join("bad.gtf");
    std::fs::write(
        &chain,
        "chain 100 chr1 1000 + 10 80 chrA 2000 + 20 90 1\n25\n",
    )
    .unwrap();
    std::fs::write(&gtf, "chr1\tsrc\texon\tbad\t20\t.\t+\t.\tgene_id \"a\";\n").unwrap();

    let result = std::panic::catch_unwind(|| {
        run_cli(&[
            "STAR".to_string(),
            "--runMode".to_string(),
            "liftOver".to_string(),
            "--genomeChainFiles".to_string(),
            chain.to_str().unwrap().to_string(),
            "--sjdbGTFfile".to_string(),
            gtf.to_str().unwrap().to_string(),
            "--outFileNamePrefix".to_string(),
            out_prefix.to_str().unwrap().to_string(),
        ])
    });
    let err = result.expect("malformed GTF coordinate must return an error, not panic");
    assert!(err.unwrap_err().contains("has invalid coordinate bad"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_solo_cell_filtering_loads_raw_matrix_and_writes_filtered_output() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_solo_filter_test_{}_{}",
        std::process::id(),
        unique
    ));
    let raw_dir = dir.join("raw");
    let filtered_prefix = dir.join("filtered/");
    let out_prefix = dir.join("star/");
    std::fs::create_dir_all(&raw_dir).unwrap();
    std::fs::write(&raw_dir.join("features.tsv"), "g1\tG1\ng2\tG2\n").unwrap();
    std::fs::write(&raw_dir.join("barcodes.tsv"), "CB1\nCB2\n").unwrap();
    std::fs::write(
        raw_dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n%\n2 2 3\n1 1 5\n2 1 2\n2 2 1\n",
    )
    .unwrap();

    let result = run_cli(&[
        "STAR".to_string(),
        "--runMode".to_string(),
        "soloCellFiltering".to_string(),
        raw_dir.to_str().unwrap().to_string(),
        filtered_prefix.to_str().unwrap().to_string(),
        "--soloCellFilter".to_string(),
        "TopCells".to_string(),
        "1".to_string(),
        "--outFileNamePrefix".to_string(),
        out_prefix.to_str().unwrap().to_string(),
    ])
    .unwrap();

    assert_eq!(result.exit_code, 0);
    assert!(result.solo_cell_filtering.as_ref().unwrap().exited);
    assert!(result.log_stdout.contains("starting SoloCellFiltering"));
    assert_eq!(
        std::fs::read_to_string(filtered_prefix.join("matrix.mtx")).unwrap(),
        "%%MatrixMarket matrix coordinate integer general\n%\n2 1 2\n1 1 5\n2 1 2\n"
    );
    assert_eq!(
        std::fs::read_to_string(filtered_prefix.join("barcodes.tsv")).unwrap(),
        "CB1\n"
    );
    assert!(std::fs::read_to_string(out_prefix.join("Log.out"))
        .unwrap()
        .ends_with("ALL DONE!\n"));

    std::fs::remove_dir_all(dir).unwrap();
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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();

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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();

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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
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
fn cli_genome_generate_auto_detects_gzip_sjdb_gtf_file() {
    let dir = unique_temp_dir("cli_auto_gzip_gtf_genome");
    let genome_dir = dir.join("genome");
    let out_tmp = dir.join("tmp");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let gtf = dir.join("genes.gtf.gz");
    let gtf_contents = concat!(
        "chr1\tsrc\texon\t1\t4\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
        "chr1\tsrc\texon\t9\t12\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
    );
    std::fs::write(&fasta, ">chr1\nACCCGTAGTTTTCCCC\n").unwrap();
    std::fs::write(&gtf, gzip_bytes(gtf_contents.as_bytes())).unwrap();

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
    assert_eq!(
        result
            .genome_generate
            .last()
            .and_then(|generation| generation.gtf.as_ref())
            .map(|gtf| gtf.sjdb_list_from_gtf_out_tab.as_str()),
        Some("chr1\t5\t8\t+\t1\n")
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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    assert!(aligned
        .log_final_out
        .contains("Uniquely mapped reads number |\t1"));
    assert!(aligned.log_main.contains("Completed: thread #0"));
    let sam = std::fs::read_to_string(dir.join("align_out/Aligned.out.sam")).unwrap();
    assert!(sam.starts_with("@HD\tVN:1.4\n"));
    assert!(sam.contains("@SQ\tSN:chr1\tLN:16\n"));
    assert!(sam.contains("r1\t0\tchr1\t1\t"));
    assert!(std::fs::read_to_string(dir.join("align_out/Log.final.out"))
        .unwrap()
        .contains("Number of input reads |\t1"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_paired_keep_input_order_writes_chunk_file_output() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_paired_order_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("paired_order/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

    run_cli(&[
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

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--runThreadN".to_string(),
        "2".to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outSAMorder".to_string(),
        "PairedKeepInputOrder".to_string(),
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
    assert_eq!(aligned.process_chunks.len(), 1);
    let map_chunk = &aligned.process_chunks[0].map_chunks[0];
    let expected_chunk = format!("{}/Aligned.out.sam.chunk0", aligned.parameters.out_file_tmp);
    assert_eq!(
        map_chunk.paired_keep_input_order_final_name.as_deref(),
        Some(expected_chunk.as_str())
    );
    assert!(!dir
        .join("paired_order/_STARtmp/Aligned.out.sam.chunk0")
        .exists());
    let sam = std::fs::read_to_string(dir.join("paired_order/Aligned.out.sam")).unwrap();
    assert!(sam.contains("r1\t0\tchr1\t1\t"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_bysjout_runs_second_mapping_stage() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_bysjout_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("bysjout/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

    run_cli(&[
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

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
        "--outSAMtype".to_string(),
        "SAM".to_string(),
        "--outFilterType".to_string(),
        "BySJout".to_string(),
        "--outSJtype".to_string(),
        "Standard".to_string(),
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
    assert!(aligned
        .log_main
        .contains("Completed stage 1 mapping of outFilterBySJout mapping"));
    assert_eq!(aligned.process_chunks.len(), 2);
    assert!(dir.join("bysjout/SJ.out.tab").exists());

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn star_main_align_reads_calls_solo_process_and_output() {
    let p = Parameters {
        run_mode_in: vec!["alignReads".to_string()],
        run_thread_n: 1,
        out_file_name_prefix: "out/".to_string(),
        p_solo: ParametersSolo {
            solo_type: SOLO_TYPE_CB_SAM_TAG_OUT,
            out_file_names: vec!["Solo.out/".to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    let chunks = vec![ReadAlignChunk {
        ra: ReadAlign {
            solo_read: SoloRead {
                read_bar: Some(SoloReadBarcode {
                    stats: SoloReadBarcodeStats {
                        names: vec!["ok".to_string()],
                        v: vec![7],
                    },
                    ..Default::default()
                }),
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    }];

    let result = star_l58_main(
        &[
            "STAR".to_string(),
            "--runMode".to_string(),
            "alignReads".to_string(),
        ],
        p,
        b"",
        Some(Genome::default()),
        Some(Transcriptome::default()),
        None,
        Some(chunks),
        &BTreeSet::new(),
        &[],
        &[],
    )
    .unwrap();

    let solo = result.solo_process_and_output.unwrap();
    assert!(solo.returned_after_barcode_output);
    assert!(solo.files["out/Solo.out/Barcodes.stats"].contains("yesWLmatchExact"));
}

#[test]
fn star_main_align_reads_sums_quant_gene_counts_from_chunks() {
    let p = Parameters {
        run_mode_in: vec!["alignReads".to_string()],
        run_thread_n: 2,
        quant_ge_count_yes: true,
        ..Default::default()
    };
    let mut q0 = quantifications_l3_quantifications_quantifications(2);
    q0.gene_counts.c_multi = 1;
    q0.gene_counts.c_none = vec![2, 3, 4];
    q0.gene_counts.c_ambig = vec![5, 6, 7];
    q0.gene_counts.g_count = vec![vec![10, 11], vec![12, 13], vec![14, 15]];
    let mut q1 = quantifications_l3_quantifications_quantifications(2);
    q1.gene_counts.c_multi = 20;
    q1.gene_counts.c_none = vec![30, 40, 50];
    q1.gene_counts.c_ambig = vec![60, 70, 80];
    q1.gene_counts.g_count = vec![vec![1, 2], vec![3, 4], vec![5, 6]];

    let chunks = vec![
        ReadAlignChunk {
            chunk_tr: Some(Transcriptome {
                n_ge: 2,
                ge_id: vec!["g0".to_string(), "g1".to_string()],
                quants: q0,
                ..Default::default()
            }),
            ..Default::default()
        },
        ReadAlignChunk {
            chunk_tr: Some(Transcriptome {
                n_ge: 2,
                ge_id: vec!["g0".to_string(), "g1".to_string()],
                quants: q1,
                ..Default::default()
            }),
            ..Default::default()
        },
    ];

    let result = star_l58_main(
        &[
            "STAR".to_string(),
            "--runMode".to_string(),
            "alignReads".to_string(),
        ],
        p,
        b"",
        Some(Genome::default()),
        Some(Transcriptome::default()),
        None,
        Some(chunks),
        &BTreeSet::new(),
        &[],
        &[],
    )
    .unwrap();

    let gc = &result.transcriptome.unwrap().quants.gene_counts;
    assert_eq!(gc.c_multi, 21);
    assert_eq!(gc.c_none, vec![32, 43, 54]);
    assert_eq!(gc.c_ambig, vec![65, 76, 87]);
    assert_eq!(gc.g_count, vec![vec![11, 13], vec![15, 17], vec![19, 21]]);
}

#[test]
fn star_main_align_reads_concatenates_paired_keep_input_order_chunks() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_star_main_chunk_cat_test_{}_{}",
        std::process::id(),
        unique
    ));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("Aligned.out.sam.chunk0"), "first\n").unwrap();
    std::fs::write(dir.join("Aligned.out.sam.chunk1"), "second\n").unwrap();
    let p = Parameters {
        run_mode_in: vec!["alignReads".to_string()],
        run_thread_n: 2,
        out_sam_bool: true,
        out_sam_order: "PairedKeepInputOrder".to_string(),
        out_file_tmp: dir.to_string_lossy().to_string(),
        out_sam_contents: "@HD\n".to_string(),
        out_tmp_keep: "All".to_string(),
        ..Default::default()
    };

    let result = star_l58_main(
        &[
            "STAR".to_string(),
            "--runMode".to_string(),
            "alignReads".to_string(),
        ],
        p,
        b"",
        Some(Genome::default()),
        Some(Transcriptome::default()),
        None,
        Some(vec![ReadAlignChunk::default(), ReadAlignChunk::default()]),
        &BTreeSet::new(),
        &[],
        &[],
    )
    .unwrap();

    assert_eq!(result.parameters.out_sam_contents, "@HD\nfirst\nsecond\n");
    assert!(!dir.join("Aligned.out.sam.chunk0").exists());
    assert!(!dir.join("Aligned.out.sam.chunk1").exists());

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
    std::fs::write(&fasta, ">chr1\nACGTTGCAAGTCCTGA\n>chr2\nACGTTGCAAGTCCTGA\n").unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    std::fs::write(&fasta, ">chr1\nACGTTGCAAGTCCTGA\n>chr2\nACGTTGCAAGTCCTGA\n").unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    assert_eq!(chim_junction.matches("# Nreads ").count(), 1);
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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    assert!(bam.starts_with(&[0x1f, 0x8b]));
    let bam = bam_payload(&bam);
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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();
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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    assert!(aligned
        .parameters
        .out_sam_contents
        .starts_with("@HD\tVN:1.4\n"));
    assert!(aligned
        .parameters
        .out_sam_contents
        .contains("r1\t0\tchr1\t1\t"));
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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
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
    assert!(aligned
        .log_final_out
        .contains("Number of reads unmapped: too short |\t1"));

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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    assert!(bam.starts_with(&[0x1f, 0x8b]));
    let bam = bam_payload(&bam);
    assert!(bam.starts_with(b"BAM\x01"));
    assert!(bam.len() > aligned.parameters.out_bam_unsorted_header.len());
    let samtools_available = std::process::Command::new("samtools")
        .arg("--version")
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false);
    if samtools_available {
        let bam_path = dir.join("bam_out/Aligned.out.bam");
        let quickcheck = std::process::Command::new("samtools")
            .arg("quickcheck")
            .arg("-v")
            .arg(&bam_path)
            .output()
            .unwrap();
        assert!(
            quickcheck.status.success(),
            "{}",
            String::from_utf8_lossy(&quickcheck.stderr)
        );
        let viewed = std::process::Command::new("samtools")
            .arg("view")
            .arg(&bam_path)
            .output()
            .unwrap();
        assert!(
            viewed.status.success(),
            "{}",
            String::from_utf8_lossy(&viewed.stderr)
        );
        let viewed = String::from_utf8_lossy(&viewed.stdout);
        assert!(viewed.contains("r1\t0\tchr1\t1\t255\t8M"));
    }
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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    assert!(bam.starts_with(&[0x1f, 0x8b]));
    let bam = bam_payload(&bam);
    assert!(bam.starts_with(b"BAM\x01"));
    assert!(bam.len() > aligned.parameters.out_bam_unsorted_header.len());
    assert!(aligned
        .bam_sort
        .as_ref()
        .map(|sort| !sort.output_bam.is_empty())
        .unwrap_or(false));
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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    assert!(std::fs::read(&bam_path).unwrap().starts_with(&[0x1f, 0x8b]));

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
        bam_path.to_str().unwrap().to_string(),
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

    let samtools_available = std::process::Command::new("samtools")
        .arg("--version")
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false);
    if samtools_available {
        let external_sam = dir.join("external.sam");
        let external_bam = dir.join("external.bam");
        std::fs::write(
            &external_sam,
            "@HD\tVN:1.4\n@SQ\tSN:chr1\tLN:16\nrExt\t0\tchr1\t1\t255\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF\tNH:i:1\n",
        )
        .unwrap();
        let converted = std::process::Command::new("samtools")
            .arg("view")
            .arg("-b")
            .arg("-o")
            .arg(&external_bam)
            .arg(&external_sam)
            .output()
            .unwrap();
        assert!(
            converted.status.success(),
            "{}",
            String::from_utf8_lossy(&converted.stderr)
        );
        let external_signal_prefix = dir.join("signal_external/");
        let external_signal = run_cli(&[
            "STAR".to_string(),
            "--runMode".to_string(),
            "inputAlignmentsFromBAM".to_string(),
            "--inputBAMfile".to_string(),
            external_bam.to_str().unwrap().to_string(),
            "--outWigType".to_string(),
            "bedGraph".to_string(),
            "--outFileNamePrefix".to_string(),
            external_signal_prefix.to_str().unwrap().to_string(),
        ])
        .unwrap();
        assert_eq!(external_signal.exit_code, 0);
        let external_signal_file = dir.join("signal_external/Signal.UniqueMultiple.str1.out.bg");
        let external_signal_contents = std::fs::read_to_string(external_signal_file).unwrap();
        assert!(external_signal_contents.contains("chr1\t0\t8\t1"));
    }

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
    let processed_file = std::fs::read(processed_path).unwrap();
    assert!(processed_file.starts_with(&[0x1f, 0x8b]));
    let processed = bam_payload(&processed_file);

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
fn cli_input_alignments_from_bam_malformed_header_returns_error_without_panic() {
    let dir = unique_temp_dir("input_bam_bad_header_no_panic");
    let out_prefix = dir.join("bam/");
    std::fs::create_dir_all(&dir).unwrap();
    let bam_path = dir.join("bad.bam");

    let mut bam = Vec::new();
    bam.extend_from_slice(b"BAM\x01");
    bam.extend_from_slice(&(-1i32).to_ne_bytes());
    bam.extend_from_slice(&0i32.to_ne_bytes());
    std::fs::write(&bam_path, &bam).unwrap();

    let result = std::panic::catch_unwind(|| {
        run_cli(&[
            "STAR".to_string(),
            "--runMode".to_string(),
            "inputAlignmentsFromBAM".to_string(),
            "--inputBAMfile".to_string(),
            bam_path.to_str().unwrap().to_string(),
            "--outFileNamePrefix".to_string(),
            out_prefix.to_str().unwrap().to_string(),
        ])
    });
    let err = result.expect("malformed BAM header must return an error, not panic");
    assert!(err.unwrap_err().contains("malformed BAM header"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_input_alignments_from_bam_malformed_record_layout_returns_error_without_panic() {
    let dir = unique_temp_dir("input_bam_bad_record_layout_no_panic");
    let out_prefix = dir.join("bam/");
    std::fs::create_dir_all(&dir).unwrap();
    let bam_path = dir.join("bad_record.bam");

    let mut bam = Vec::new();
    bam.extend_from_slice(b"BAM\x01");
    bam.extend_from_slice(&0i32.to_ne_bytes());
    bam.extend_from_slice(&1i32.to_ne_bytes());
    bam.extend_from_slice(&5i32.to_ne_bytes());
    bam.extend_from_slice(b"chr1\0");
    bam.extend_from_slice(&1000i32.to_ne_bytes());

    bam.extend_from_slice(&32i32.to_ne_bytes());
    bam.extend_from_slice(&0i32.to_ne_bytes());
    bam.extend_from_slice(&0i32.to_ne_bytes());
    bam.extend_from_slice(&250u32.to_ne_bytes());
    bam.extend_from_slice(&0u32.to_ne_bytes());
    bam.extend_from_slice(&0i32.to_ne_bytes());
    bam.extend_from_slice(&(-1i32).to_ne_bytes());
    bam.extend_from_slice(&(-1i32).to_ne_bytes());
    bam.extend_from_slice(&0i32.to_ne_bytes());
    std::fs::write(&bam_path, &bam).unwrap();

    let result = std::panic::catch_unwind(|| {
        run_cli(&[
            "STAR".to_string(),
            "--runMode".to_string(),
            "inputAlignmentsFromBAM".to_string(),
            "--inputBAMfile".to_string(),
            bam_path.to_str().unwrap().to_string(),
            "--outFileNamePrefix".to_string(),
            out_prefix.to_str().unwrap().to_string(),
        ])
    });
    let err = result.expect("malformed BAM record layout must return an error, not panic");
    assert!(err.unwrap_err().contains("malformed BAM record"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_input_alignments_from_bam_malformed_bgzf_returns_error_without_panic() {
    let dir = unique_temp_dir("input_bam_bad_bgzf_no_panic");
    let out_prefix = dir.join("bam/");
    std::fs::create_dir_all(&dir).unwrap();
    let bam_path = dir.join("bad_bgzf.bam");

    let mut bgzf = Vec::new();
    bgzf.extend_from_slice(&[0x1f, 0x8b, 8, 4]);
    bgzf.extend_from_slice(&[0, 0, 0, 0, 0, 255]);
    bgzf.extend_from_slice(&6u16.to_le_bytes());
    bgzf.extend_from_slice(b"BC");
    bgzf.extend_from_slice(&2u16.to_le_bytes());
    bgzf.extend_from_slice(&25u16.to_le_bytes());
    bgzf.extend_from_slice(&[0; 4]);
    std::fs::write(&bam_path, &bgzf).unwrap();

    let result = std::panic::catch_unwind(|| {
        run_cli(&[
            "STAR".to_string(),
            "--runMode".to_string(),
            "inputAlignmentsFromBAM".to_string(),
            "--inputBAMfile".to_string(),
            bam_path.to_str().unwrap().to_string(),
            "--outFileNamePrefix".to_string(),
            out_prefix.to_str().unwrap().to_string(),
        ])
    });
    let err = result.expect("malformed BGZF BAM must return an error, not panic");
    let err = err.unwrap_err();
    assert!(
        err.contains("truncated BAM record")
            || err.contains("malformed BAM record")
            || err.contains("could not read BAM header")
    );

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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(
        &gtf,
        "chr1\tsrc\texon\t1\t8\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
    )
    .unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(
        &gtf,
        "chr1\tsrc\texon\t1\t8\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
    )
    .unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    assert!(bam.starts_with(&[0x1f, 0x8b]));
    let bam = bam_payload(&bam);
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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads_plain, UNIQUE_TEST_READ).unwrap();
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
    assert!(aligned
        .log_final_out
        .contains("Uniquely mapped reads number |\t1"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_auto_detects_gzip_fastq_without_read_files_command() {
    let dir = unique_temp_dir("cli_auto_gzip_align");
    let genome_dir = dir.join("genome");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads_gz = dir.join("reads.fq.gz");
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads_gz, gzip_bytes(UNIQUE_TEST_READ.as_bytes())).unwrap();

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
    assert!(aligned.parameters.read_files_command_string.is_empty());
    assert!(aligned.log_final_out.contains("Number of input reads |\t1"));
    assert!(aligned
        .log_final_out
        .contains("Uniquely mapped reads number |\t1"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_auto_detects_gzip_genome_index_files() {
    let dir = unique_temp_dir("cli_auto_gzip_index");
    let genome_dir = dir.join("genome");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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

    for file_name in ["Genome", "SA", "SAindex"] {
        let path = genome_dir.join(file_name);
        let plain = std::fs::read(&path).unwrap();
        std::fs::write(&path, gzip_bytes(&plain)).unwrap();
    }

    let aligned = run_cli(&[
        "STAR".to_string(),
        "--genomeDir".to_string(),
        genome_dir.to_str().unwrap().to_string(),
        "--readFilesIn".to_string(),
        reads.to_str().unwrap().to_string(),
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
    assert!(aligned
        .log_final_out
        .contains("Uniquely mapped reads number |\t1"));

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
fn cli_align_reads_nonrepetitive_paired_read_does_not_panic_in_suffix_index_search() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_paired_unique_align_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("paired_unique/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("tiny.fa");
    let reads_1 = dir.join("reads_1.fq");
    let reads_2 = dir.join("reads_2.fq");
    std::fs::write(&fasta, ">chr1\nACGTACGTTGCAAGTC\n").unwrap();
    std::fs::write(&reads_1, "@pair_unique/1\nACGTACGT\n+\nFFFFFFFF\n").unwrap();
    std::fs::write(&reads_2, "@pair_unique/2\nTGCAAGTC\n+\nHHHHHHHH\n").unwrap();

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
        std::fs::read_to_string(dir.join("paired_unique/Aligned.out.sam"))
            .unwrap()
            .starts_with("@HD\tVN:1.4\n")
    );

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cli_align_reads_randomized_single_end_does_not_panic_in_stitching() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_cli_stitching_regression_test_{}_{}",
        std::process::id(),
        unique
    ));
    let genome_dir = dir.join("genome");
    let out_prefix = dir.join("stitching/");
    std::fs::create_dir_all(&genome_dir).unwrap();
    let fasta = dir.join("genome.fa");
    let reads = dir.join("reads.fq");
    std::fs::write(
        &fasta,
        concat!(
            ">chr1\n",
            "GGTTCTCCACGTTAAAAGCCGATTAGCACCTGAGACGGCTCCGTAGCCATGATTGAAGACTATTGTAGCCTGGAAAGTTGACGTTCCAGCAGAATGACGCCACACTGAGATTGGTTTCAAAGGCGGCTTCGACGAAAACTTCTAGGCGCCTCATAAAGCATGCATCTATACCCAGCGGAATCTCATCAATCGAGGCGCGTATGATCTGGGCTAAAGATCGTAACGCCCTACATATCTCTTCCGTATTAGGTTTAAAGACTTTACCTCGCCGGACGCGGGTAACCGGTTTGGAGAGCGGGGACAAGGCTCCGGTCTGTTTCGGCCCAAAAAGGGCCTACAGATATGGCCGCACACACCTCCCCCGTGGTGGTCTATAATTGTAATGATAAAAAAAACTCACCTCGTTCAAAGGGCTGGGTACACATACCCTTCGCGGATTCCGGCATTCAGTCTGCGTGAAGTACGTTCAGCCGGGGGTTTGATAACTCAGGCTGCTGGTAGTAGTCTCAGTCGTGGACGTGTACTTGTCTTCCATGCCTTTCGCCCTCAGACGATTTACCTAAAAGTTGCAGTCCCGGCTCGATACATTTATCGTATTCAACAGGAGAAATTGTGCAGGTGGTCCTGAATCTCGAGGAGAGGTCAAAAAATCCCGTCAAACAATAGACATACCTCATGAACGCGCGCCTTAATACCCAAG\n",
            ">chr2\n",
            "GGGACACCATCTAGTCTGCTCGATACTCCTGACATGTCACGGTACGCCACGAGCATTTTTTGCCCGCTGGTAACGCCTACGACGCCCACTTCCGGTATTAGTGTCACCTTGCCGAGGCCATTTGAGTATGTGAGAGCACTGTTCCAGCGCAAAGCCTTTAACTGGAATGTTGACGTTTAGAAGTTCCACGCTCAGGACCGCTCTCCCCAGGCCAGCCTCTTGGGCCGCCCCCCAGGTTTAACTTCCTAAAGGATTGGTCCCAAGGACCCAAAGCCGTTTATAATGGACGTCTGCACTTAGGACTGGGCGGTTAAGCTGGTTTAGTCGAATGACCTAACGGGTCTACGTCGCCCGAGCCCCTTTGCTCTCCGGCAAGCAAGTGCAATCGACGCGCTCCCCAGTACCCAGGATGCATTTACCTGGCCTTTTAAAAGTCTTAGGCCTCGATTGAAAACGGTTGTGACCATCCATGGGCGAACTGTGCCAGTCGGGCTAAGTGCCTACAGGAATCTGCGGTATAACGTGCTCACTCAAGGGGACTTAGTTTCTTACACCTTGCTCAGCGACTGGACCTTTTTCCGAACTCGATCTGGACGGGTGCTTGGCGGTGTAGCAATGACAGACTACAGTACTCTCTGGGGCAGCGCTCT\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &reads,
        concat!(
            "@r0_45\nAAGGACCCAAAGCCGGTTAAAATGGACGTCTGCACTTAGGACTGG\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
            "@r1_36\nCGAACTCGATCTGGACGGGTGCTTGGGGGTGTAGCC\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
            "@r2_36\nTCGCGGATTCCGCCTTTCAGTCTCCGAGAAGTACGT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
            "@r3_36\nAGTCGTGGACGTGTACTTGTCTTCCAAGCCTTTCGC\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
            "@r4_36\nCCAGGTTTAGCTCCCTAAAGGATTGGTCCCAAGGAC\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
            "@r5_75\nCAAAGGGCTGTGTACACATACCCATCGCGGATTCCGGCATTCAGTCTGCGTGAAGTACGTTCAGCCGGGGGTTTG\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
            "@r6_30\nTTATTGTAATGATAAAAAAAACTCACCTCG\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
            "@r7_75\nATTGATAACACAGGCTGCTGGTAGTAGTCTCAGTCGTTGACGTGTAATTGTCTTCCATGCCTTTCGCCCTCAGAC\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
            "@r8_36\nCAGGAGCAATTGTGCAGGTGGTCCTGAATATCGAGG\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
            "@r9_30\nGACCCAAAGCCGTTTATAATGGACGTCTGC\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
            "@r10_30\nGTTTATAATGGACGTCTGCACTTAGGACTG\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
            "@r11_30\nGTGACATGTCACGGTACGCCACGAGCATTT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n",
        ),
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
        "3".to_string(),
        "--genomeChrBinNbits".to_string(),
        "4".to_string(),
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
        "12".to_string(),
        "--seedMapMin".to_string(),
        "5".to_string(),
        "--outFilterMatchNmin".to_string(),
        "0".to_string(),
        "--outFilterScoreMin".to_string(),
        "0".to_string(),
        "--outFilterMultimapNmax".to_string(),
        "20".to_string(),
        "--outSAMmultNmax".to_string(),
        "20".to_string(),
    ])
    .unwrap();

    assert_eq!(aligned.exit_code, 0);
    assert!(aligned
        .log_final_out
        .contains("Number of input reads |\t12"));
    let sam = std::fs::read(dir.join("stitching/Aligned.out.sam")).unwrap();
    assert!(!sam.contains(&0));

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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(
        &reads,
        "@r1\nACGTTGCA\n+\nFFFFFFFF\n@r2\nCGTTGCAA\n+\nFFFFFFFF\n",
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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();

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
    std::fs::write(&fasta, UNIQUE_TEST_GENOME).unwrap();
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
    std::fs::write(&reads, UNIQUE_TEST_READ).unwrap();
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
