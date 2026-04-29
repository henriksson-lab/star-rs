use std::path::Path;

use star_rs::generated::functions::*;
use star_rs::generated::structs::{
    CBMatchWL, ChimericAlign, ChimericDetection, ChimericSegment, ClipMate, GTF, Genome, GenomeOut,
    InOutStreams, JunctionRecord, MultiMappers, OutSJ, OutWigFlags, PackedArray,
    ParameterScanEntry, Parameters, ParametersChimeric, ParametersClip, ParametersGenome,
    ParametersGenomeTransform, ParametersScanState, ParametersSolo, ReadAlign, ReadAlignChunk,
    ReadAlignGenomeTransformResult, ReadAlignPeOverlap, ReadAnnotFeature, ReadAnnotations,
    ReadClipInput, ReadSoloFeatures, SNP, SharedMemoryException, SignalFromBamRecord, SjdbClass,
    Solo, SoloBarcode, SoloCellFilter, SoloFeature, SoloFeatureReadInfo, SoloFilteredCells,
    SoloRead, SoloReadBarcode, SoloReadFeature, SoloReadFlagClass, SoloUmiFiltering, Stats,
    SuperTranscript, SuperTranscriptome, ThreadControl, TrTypeStruct, Transcript, Transcriptome,
    TranscriptomeGeneFull, UMIdedup, Variation, WaspMapOutcome, sjInfo,
};

#[test]
fn simple_uint_and_sj_comparators_match_original_ordering() {
    assert_eq!(outputsj_l7_compareuint(4, 3), 1);
    assert_eq!(outputsj_l7_compareuint(3, 4), -1);
    assert_eq!(outputsj_l7_compareuint(4, 4), 0);

    let mut a = [0u8; 8];
    let mut b = [0u8; 8];
    a[0..4].copy_from_slice(&10u32.to_ne_bytes());
    a[4..8].copy_from_slice(&2u32.to_ne_bytes());
    b[0..4].copy_from_slice(&10u32.to_ne_bytes());
    b[4..8].copy_from_slice(&3u32.to_ne_bytes());
    assert_eq!(outsj_l15_comparesj(&a, &b), -1);
    b[0..4].copy_from_slice(&9u32.to_ne_bytes());
    assert_eq!(outsj_l15_comparesj(&a, &b), 1);
}

#[test]
fn parameters_scan_one_line_matches_comments_levels_and_errors() {
    let mut state = ParametersScanState {
        par_array: vec![
            ParameterScanEntry {
                name_string: "runThreadN".to_string(),
                input_level_allowed: -1,
                input_level: -1,
                value_line: String::new(),
            },
            ParameterScanEntry {
                name_string: "genomeDir".to_string(),
                input_level_allowed: 1,
                input_level: 1,
                value_line: "cmdGenome".to_string(),
            },
        ],
        parameter_input_name: vec![
            "Default".to_string(),
            "Command-Line-Initial".to_string(),
            "Command-Line".to_string(),
        ],
        ..Default::default()
    };

    assert_eq!(
        parameters_l1205_parameters_scanoneline(&mut "".to_string(), 0, -1, &mut state).unwrap(),
        0
    );
    assert_eq!(
        parameters_l1205_parameters_scanoneline(
            &mut "  runThreadN 8".to_string(),
            0,
            -1,
            &mut state
        )
        .unwrap(),
        0
    );
    assert_eq!(
        parameters_l1205_parameters_scanoneline(
            &mut "# runThreadN 8".to_string(),
            1,
            -1,
            &mut state
        )
        .unwrap(),
        0
    );
    assert_eq!(
        parameters_l1205_parameters_scanoneline(
            &mut "runThreadN 8 9".to_string(),
            1,
            2,
            &mut state
        )
        .unwrap(),
        1
    );
    assert_eq!(state.par_array[0].input_level, -1);

    assert_eq!(
        parameters_l1205_parameters_scanoneline(
            &mut "runThreadN 8 9".to_string(),
            1,
            -1,
            &mut state
        )
        .unwrap(),
        1
    );
    assert_eq!(state.par_array[0].input_level, 1);
    assert_eq!(state.par_array[0].value_line, "8 9");
    assert!(state.log_main.contains("runThreadN"));
    assert!(state.log_main.contains("~RE-DEFINED"));

    let duplicate = parameters_l1205_parameters_scanoneline(
        &mut "runThreadN 10".to_string(),
        1,
        -1,
        &mut state,
    )
    .unwrap_err();
    assert!(duplicate.contains("duplicate parameter \"runThreadN\""));

    let too_late = parameters_l1205_parameters_scanoneline(
        &mut "genomeDir later".to_string(),
        2,
        -1,
        &mut state,
    )
    .unwrap_err();
    assert!(too_late.contains("cannot be defined at the input level \"Command-Line\""));

    let empty_unknown_first =
        parameters_l1205_parameters_scanoneline(&mut "unknown".to_string(), 1, -1, &mut state)
            .unwrap_err();
    assert!(empty_unknown_first.contains("empty value for parameter \"unknown\""));

    let unknown = parameters_l1205_parameters_scanoneline(
        &mut "unknown value".to_string(),
        1,
        -1,
        &mut state,
    )
    .unwrap_err();
    assert!(unknown.contains("unrecognized parameter name \"unknown\""));
}

#[test]
fn parameters_scan_all_lines_applies_each_input_line() {
    let mut state = ParametersScanState {
        par_array: vec![
            ParameterScanEntry {
                name_string: "alpha".to_string(),
                input_level_allowed: -1,
                input_level: -1,
                value_line: String::new(),
            },
            ParameterScanEntry {
                name_string: "beta".to_string(),
                input_level_allowed: -1,
                input_level: -1,
                value_line: String::new(),
            },
        ],
        parameter_input_name: vec!["Default".to_string(), "Command-Line".to_string()],
        ..Default::default()
    };

    parameters_l1197_parameters_scanalllines(
        "# comment\nalpha one\n\nbeta two words\n",
        1,
        -1,
        &mut state,
    )
    .unwrap();

    assert_eq!(state.par_array[0].value_line, "one");
    assert_eq!(state.par_array[1].value_line, "two words");
    assert_eq!(state.par_array[0].input_level, 1);
    assert_eq!(state.par_array[1].input_level, 1);
}

#[test]
fn parameters_constructor_registers_star_parameter_scan_table() {
    let mut state = parameters_l19_parameters_parameters();
    assert_eq!(state.par_array.len(), 203);
    assert_eq!(state.parameter_input_name[0], "Default");
    assert_eq!(state.parameter_input_name[1], "Command-Line-Initial");
    assert_eq!(state.parameter_input_name[2], "Command-Line");
    assert_eq!(state.parameter_input_name[3], "genomeParameters.txt");
    assert_eq!(state.par_array[0].name_string, "versionGenome");
    assert_eq!(state.par_array[1].name_string, "parametersFiles");
    assert_eq!(state.par_array[1].input_level_allowed, 2);
    assert_eq!(state.par_array[42].name_string, "outFileNamePrefix");
    assert_eq!(state.par_array[42].input_level_allowed, 2);
    assert_eq!(state.par_array[127].name_string, "seedMapMin");
    assert_eq!(state.par_array[202].name_string, "soloCBtype");

    parameters_l1205_parameters_scanoneline(&mut "runThreadN 16".to_string(), 1, -1, &mut state)
        .unwrap();
    parameters_l1205_parameters_scanoneline(
        &mut "parametersFiles custom.params".to_string(),
        2,
        -1,
        &mut state,
    )
    .unwrap();
    assert_eq!(state.par_array[4].value_line, "16");
    assert_eq!(state.par_array[1].value_line, "custom.params");

    let too_late = parameters_l1205_parameters_scanoneline(
        &mut "parametersFiles other.params".to_string(),
        3,
        -1,
        &mut state,
    )
    .unwrap_err();
    assert!(too_late.contains("cannot be defined at the input level \"genomeParameters.txt\""));
}

#[test]
fn parameters_input_parameters_scans_defaults_files_and_command_line() {
    let default_parameters = include_str!("../STAR/source/parametersDefault");
    let args = vec![
        "STAR".to_string(),
        "--runThreadN".to_string(),
        "8".to_string(),
        "--sysShell".to_string(),
        "/bin/sh".to_string(),
        "--outFileNamePrefix=with space/".to_string(),
        "--outBAMcompression".to_string(),
        "7".to_string(),
        "--scoreGenomicLengthLog2scale".to_string(),
        "0.5".to_string(),
        "--genomeChainFiles".to_string(),
        "chain1".to_string(),
        "chain2".to_string(),
        "--readMatesLengthsIn".to_string(),
        "Equal".to_string(),
        "--inputBAMfile".to_string(),
        "input.bam".to_string(),
        "--bamRemoveDuplicatesType".to_string(),
        "UniqueIdentical".to_string(),
        "--bamRemoveDuplicatesMate2basesN".to_string(),
        "11".to_string(),
        "--limitNreadsSoft".to_string(),
        "123".to_string(),
        "--seedNoneLociPerWindow".to_string(),
        "4".to_string(),
        "--winReadCoverageRelativeMin".to_string(),
        "0.75".to_string(),
        "--winReadCoverageBasesMin".to_string(),
        "9".to_string(),
        "--genomeDir".to_string(),
        "cmdGenome".to_string(),
    ];
    let mut p = Parameters::default();
    let state = parameters_l310_parameters_inputparameters(
        &mut p,
        &args,
        default_parameters,
        &[(
            "user.params".to_string(),
            "genomeDir fileGenome\nreadFilesIn file_1.fq file_2.fq\n".to_string(),
        )],
        None,
        &[],
    )
    .unwrap();

    let run_thread_n = state
        .par_array
        .iter()
        .find(|p| p.name_string == "runThreadN")
        .unwrap();
    assert_eq!(run_thread_n.input_level, 2);
    assert_eq!(run_thread_n.value_line, "8");

    let genome_dir = state
        .par_array
        .iter()
        .find(|p| p.name_string == "genomeDir")
        .unwrap();
    assert_eq!(genome_dir.input_level, 2);
    assert_eq!(genome_dir.value_line, "cmdGenome");

    let read_files = state
        .par_array
        .iter()
        .find(|p| p.name_string == "readFilesIn")
        .unwrap();
    assert_eq!(read_files.input_level, 4);
    assert_eq!(read_files.value_line, "file_1.fq file_2.fq");

    assert!(p.command_line.contains("--outFileNamePrefix=with space/"));
    assert!(p.command_line_full.contains("--runThreadN 8"));
    assert!(p.command_line_full.contains("--genomeDir cmdGenome"));
    assert!(
        p.command_line_full
            .contains("--outFileNamePrefix \"with space/\"")
    );
    assert!(
        p.command_line_full
            .contains("--readFilesIn file_1.fq file_2.fq")
    );
    assert!(
        state
            .log_main
            .contains("##### Final effective command line:")
    );
    assert_eq!(p.run_thread_n, 8);
    assert_eq!(p.sys_shell, "/bin/sh");
    assert_eq!(p.out_bam_compression, 7);
    assert_eq!(p.score_genomic_length_log2scale, 0.5);
    assert_eq!(p.p_ge.g_chain_files, vec!["chain1", "chain2"]);
    assert_eq!(p.read_mates_lengths_in, "Equal");
    assert_eq!(p.input_bam_file, "input.bam");
    assert_eq!(p.bam_remove_duplicates_type, "UniqueIdentical");
    assert!(p.bam_remove_duplicates_yes);
    assert!(p.bam_remove_duplicates_mark_multi);
    assert_eq!(p.bam_remove_duplicates_mate2bases_n, 11);
    assert_eq!(p.limit_nreads_soft, 123);
    assert_eq!(p.seed_none_loci_per_window, 4);
    assert_eq!(p.win_read_coverage_relative_min, 0.75);
    assert_eq!(p.win_read_coverage_bases_min, 9);
    assert_eq!(p.p_ge.g_dir, "cmdGenome/");
    assert_eq!(p.read_files_in, vec!["file_1.fq", "file_2.fq"]);
    assert_eq!(p.read_nends, 2);
    assert_eq!(p.read_files_type_n, 1);
    assert_eq!(p.read_files_names[0], vec!["file_1.fq"]);
    assert_eq!(p.read_files_names[1], vec!["file_2.fq"]);
    assert_eq!(p.out_file_tmp, "with space/_STARtmp/");
    assert!(p.out_sam_bool);
    assert_eq!(p.p_solo.solo_type, SOLO_TYPE_NONE);
    assert!(!p.p_solo.yes);
}

#[test]
fn parameters_input_parameters_rejects_zero_limit_genome_generate_ram_like_star() {
    let default_parameters = include_str!("../STAR/source/parametersDefault");
    let args = vec![
        "STAR".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "0".to_string(),
    ];
    let mut p = Parameters::default();

    let err = parameters_l310_parameters_inputparameters(
        &mut p,
        &args,
        default_parameters,
        &[],
        None,
        &[],
    )
    .unwrap_err();

    assert!(err.contains("limitGenomeGenerateRAM=0"));
    assert!(err.contains("please specify a >0 value"));
}

#[test]
fn parameters_input_parameters_warns_on_huge_limit_genome_generate_ram_like_star() {
    let default_parameters = include_str!("../STAR/source/parametersDefault");
    let args = vec![
        "STAR".to_string(),
        "--limitGenomeGenerateRAM".to_string(),
        "1000000000001".to_string(),
    ];
    let mut p = Parameters::default();

    let state = parameters_l310_parameters_inputparameters(
        &mut p,
        &args,
        default_parameters,
        &[],
        None,
        &[],
    )
    .unwrap();

    assert!(state.log_main.contains(
        "WARNING: specified limitGenomeGenerateRAM=1000000000001 bytes appears to be too large"
    ));
    assert_eq!(p.score_genomic_length_log2scale, -0.25);
}

#[test]
fn parameters_input_parameters_rejects_zero_bam_sort_ram_with_shared_genome_like_star() {
    let default_parameters = include_str!("../STAR/source/parametersDefault");
    let args = vec![
        "STAR".to_string(),
        "--outSAMtype".to_string(),
        "BAM".to_string(),
        "SortedByCoordinate".to_string(),
        "--genomeLoad".to_string(),
        "LoadAndKeep".to_string(),
        "--limitBAMsortRAM".to_string(),
        "0".to_string(),
    ];
    let mut p = Parameters::default();

    let err = parameters_l310_parameters_inputparameters(
        &mut p,
        &args,
        default_parameters,
        &[],
        None,
        &[],
    )
    .unwrap_err();

    assert!(
        err.contains("limitBAMsortRAM=0 (default) cannot be used with --genomeLoad=LoadAndKeep")
    );
}

#[test]
fn parameters_input_parameters_rejects_unknown_bam_duplicate_mode_like_star() {
    let default_parameters = include_str!("../STAR/source/parametersDefault");
    let args = vec![
        "STAR".to_string(),
        "--bamRemoveDuplicatesType".to_string(),
        "Identical".to_string(),
    ];
    let mut p = Parameters::default();

    let err = parameters_l310_parameters_inputparameters(
        &mut p,
        &args,
        default_parameters,
        &[],
        None,
        &[],
    )
    .unwrap_err();

    assert!(err.contains("unrecognized option in of --bamRemoveDuplicatesType=Identical"));
}

#[test]
fn signal_from_bam_emits_bedgraph_full_signal_and_respects_prefix() {
    let p = Parameters {
        out_wig_references_prefix: "chr".to_string(),
        out_wig_flags: OutWigFlags {
            format: 0,
            type_: 0,
            norm: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    let target_names = vec!["chr1".to_string(), "skip".to_string()];
    let target_lens = vec![5, 5];
    let records = vec![
        SignalFromBamRecord {
            tid: 0,
            pos: 1,
            cigar: vec![3 << 4],
            nh: Some(1),
            ..Default::default()
        },
        SignalFromBamRecord {
            tid: 0,
            pos: 2,
            cigar: vec![2 << 4],
            nh: Some(2),
            ..Default::default()
        },
        SignalFromBamRecord {
            tid: 1,
            pos: 0,
            cigar: vec![2 << 4],
            nh: Some(1),
            ..Default::default()
        },
    ];

    let result =
        signalfrombam_l5_signalfrombam("Signal", &p, &target_names, &target_lens, &records)
            .unwrap();

    assert_eq!(result.files["Signal.Unique.str1.out.bg"], "chr1\t1\t4\t1\n");
    assert_eq!(
        result.files["Signal.UniqueMultiple.str1.out.bg"],
        "chr1\t1\t2\t1\nchr1\t2\t4\t1.5\n"
    );
}

#[test]
fn signal_from_bam_emits_stranded_wiggle_read_end_and_rpm() {
    let p = Parameters {
        out_wig_references_prefix: "-".to_string(),
        out_wig_flags: OutWigFlags {
            strand: true,
            type_: 1,
            format: 1,
            norm: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    let target_names = vec!["chr1".to_string()];
    let target_lens = vec![6];
    let records = vec![
        SignalFromBamRecord {
            tid: 0,
            pos: 1,
            flag: 0,
            cigar: vec![3 << 4],
            nh: Some(1),
        },
        SignalFromBamRecord {
            tid: 0,
            pos: 2,
            flag: 0x10,
            cigar: vec![3 << 4],
            nh: Some(2),
        },
        SignalFromBamRecord {
            tid: 0,
            pos: 4,
            flag: 0x80,
            cigar: vec![1 << 4],
            nh: Some(1),
        },
    ];

    let result =
        signalfrombam_l5_signalfrombam("Signal", &p, &target_names, &target_lens, &records)
            .unwrap();

    assert_eq!(result.n_unique, 2.0);
    assert_eq!(result.n_multiple, 0.5);
    assert_eq!(
        result.files["Signal.Unique.str1.out.wig"],
        "variableStep chrom=chr1\n2\t500000.00000\n"
    );
    assert_eq!(
        result.files["Signal.UniqueMultiple.str2.out.wig"],
        "variableStep chrom=chr1\n5\t200000.00000\n"
    );
}

#[test]
fn output_sj_collapses_filters_and_emits_stage_outputs() {
    let gen_out = Genome {
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 10,
            ..Default::default()
        },
        chr_bin: vec![0, 0],
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![1000],
        ..Default::default()
    };
    let chunk0 = ReadAlignChunk {
        chunk_out_sj: OutSJ {
            n: 2,
            junctions: vec![
                JunctionRecord {
                    start: 1000,
                    gap: 50,
                    strand: 1,
                    motif: 1,
                    annot: 0,
                    count_unique: 1,
                    overhang_left: 10,
                    overhang_right: 8,
                    ..Default::default()
                },
                JunctionRecord {
                    start: 1200,
                    gap: 30,
                    strand: 2,
                    motif: 2,
                    annot: 1,
                    count_unique: 0,
                    overhang_left: 0,
                    overhang_right: 0,
                    ..Default::default()
                },
            ],
            ..Default::default()
        },
        chunk_out_sj1: OutSJ {
            n: 1,
            junctions: vec![JunctionRecord {
                start: 1300,
                gap: 40,
                strand: 1,
                motif: 1,
                annot: 0,
                count_unique: 2,
                overhang_left: 12,
                overhang_right: 12,
                ..Default::default()
            }],
            ..Default::default()
        },
        ..Default::default()
    };
    let chunk1 = ReadAlignChunk {
        chunk_out_sj: OutSJ {
            n: 1,
            junctions: vec![JunctionRecord {
                start: 1000,
                gap: 50,
                strand: 1,
                motif: 1,
                annot: 0,
                count_unique: 3,
                count_multiple: 1,
                overhang_left: 7,
                overhang_right: 20,
                ..Default::default()
            }],
            ..Default::default()
        },
        chunk_out_sj1: OutSJ::default(),
        ..Default::default()
    };
    let mut p = Parameters {
        run_thread_n: 2,
        limit_out_sj_collapsed: 4,
        out_filter_by_sjout_stage: 0,
        out_sjfilter_count_unique_min: vec![1; 5],
        out_sjfilter_count_total_min: vec![1; 5],
        out_sjfilter_overhang_min: vec![1; 5],
        out_sjfilter_intron_max_vs_read_n: vec![1_000; 10],
        out_sjfilter_dist_to_other_sj_min: vec![0; 5],
        ..Default::default()
    };

    let result = outputsj_l20_outputsj(&[chunk0.clone(), chunk1], &mut p, &gen_out).unwrap();

    assert_eq!(
        result.sj_out_tab,
        "chr1\t1\t50\t1\t1\t0\t4\t1\t10\nchr1\t201\t230\t2\t2\t1\t0\t0\t0\n"
    );
    assert_eq!(result.sj_start_gap_tsv, "1000\t50\n1200\t30\n");
    assert_eq!(p.sj_all, [vec![1000, 1200], vec![50, 30]]);

    p.out_filter_by_sjout_stage = 1;
    p.run_thread_n = 1;
    p.sj_all = [Vec::new(), Vec::new()];
    let result = outputsj_l20_outputsj(&[chunk0], &mut p, &gen_out).unwrap();
    assert_eq!(result.sj_out_tab, "");
    assert_eq!(p.sj_novel_n, 1);
    assert_eq!(p.sj_novel_start, vec![1300]);
    assert_eq!(p.sj_novel_end, vec![1339]);
    assert!(
        result
            .log_main
            .contains("Detected 1 novel junctions that passed filtering")
    );
}

#[test]
fn two_pass_run_pass1_noops_when_disabled() {
    let mut p = Parameters::default();
    let mut genome = Genome::default();
    let mut sjdb = SjdbClass::default();

    let result = twopassrunpass1_l9_twopassrunpass1(
        &mut p,
        &mut genome,
        None,
        &mut sjdb,
        None,
        &std::collections::BTreeSet::new(),
    )
    .unwrap();

    assert_eq!(result, Default::default());
    assert!(!p.two_pass_pass2);
}

#[test]
fn two_pass_run_pass1_sets_pass2_writes_sj_and_reopens_reads() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir =
        std::env::temp_dir().join(format!("star_rs_twopass_{}_{}", std::process::id(), unique));
    std::fs::create_dir_all(&dir).unwrap();

    let mut p = Parameters {
        two_pass_yes: true,
        two_pass_dir: dir.to_str().unwrap().to_string(),
        two_pass_pass1reads_n: 7,
        read_map_number: 10,
        run_thread_n: 1,
        read_nends: 1,
        out_sam_type: vec!["SAM".to_string()],
        out_sam_bool: true,
        out_bam_unsorted: true,
        out_bam_coord: true,
        p_ch: ParametersChimeric {
            segment_min: 12,
            ..Default::default()
        },
        quant_yes: true,
        quant_tr_sam_yes: true,
        quant_tr_sam_bam_yes: true,
        quant_gene_full_yes: true,
        quant_ge_count_yes: true,
        quant_gene_yes: true,
        out_sam_unmapped_within: true,
        out_filter_by_sjout_stage: 2,
        out_reads_unmapped: "Fastx".to_string(),
        wasp_yes: true,
        wasp_output_mode: "SAMtag".to_string(),
        wasp_sam_tag: true,
        p_solo: ParametersSolo {
            type_str: "CB_UMI_Simple".to_string(),
            solo_type: 1,
            ..Default::default()
        },
        p_ge: ParametersGenome {
            transform: ParametersGenomeTransform {
                out_yes: true,
                out_quant: true,
                out_sam: true,
                out_sj: true,
                ..Default::default()
            },
            ..Default::default()
        },
        limit_sjdb_insert_nsj: 0,
        sjdb_insert_out_dir: dir.to_str().unwrap().to_string(),
        read_files_in: vec!["r1.fq".to_string()],
        read_in_open: vec![true],
        read_files_command_pid: vec![123],
        ..Default::default()
    };
    let mut genome = Genome {
        g: vec![
            0,
            1,
            2,
            3,
            GENOME_SPACING_CHAR,
            GENOME_SPACING_CHAR,
            GENOME_SPACING_CHAR,
            GENOME_SPACING_CHAR,
        ],
        n_genome: 8,
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![0, 8],
        chr_length: vec![4],
        genome_out: GenomeOut {
            conv_yes: true,
            ..Default::default()
        },
        p_ge: ParametersGenome {
            g_dir: dir.to_str().unwrap().to_string(),
            sjdb_file_chr_start_end: vec!["-".to_string()],
            sjdb_gtf_file: "-".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    let mut sjdb = SjdbClass::default();
    let existing = std::collections::BTreeSet::from(["r1.fq".to_string()]);

    let result = twopassrunpass1_l9_twopassrunpass1(
        &mut p,
        &mut genome,
        None,
        &mut sjdb,
        Some(vec![ReadAlignChunk::default()]),
        &existing,
    )
    .unwrap();

    assert!(p.two_pass_pass2);
    assert_eq!(
        p.two_pass_pass1sj_file,
        format!("{}/SJ.out.tab", dir.to_str().unwrap())
    );
    assert_eq!(
        std::fs::read_to_string(&p.two_pass_pass1sj_file).unwrap(),
        ""
    );
    assert_eq!(result.pass1_parameters.out_sam_type[0], "None");
    assert!(!result.pass1_parameters.out_sam_bool);
    assert!(!result.pass1_parameters.out_bam_unsorted);
    assert!(!result.pass1_parameters.out_bam_coord);
    assert_eq!(result.pass1_parameters.p_ch.segment_min, 0);
    assert!(!result.pass1_parameters.quant_yes);
    assert!(!result.pass1_parameters.out_sam_unmapped_within);
    assert_eq!(result.pass1_parameters.out_reads_unmapped, "None");
    assert_eq!(
        result.pass1_parameters.out_file_name_prefix,
        dir.to_str().unwrap()
    );
    assert_eq!(result.pass1_parameters.read_map_number, 7);
    assert!(!result.pass1_parameters.wasp_yes);
    assert_eq!(result.pass1_parameters.wasp_output_mode, "None");
    assert!(!result.pass1_parameters.wasp_sam_tag);
    assert_eq!(result.pass1_parameters.p_solo.solo_type, 0);
    assert!(!result.pass1_parameters.p_ge.transform.out_yes);
    assert!(result.log_progress.contains("Started 1st pass mapping"));
    assert!(result.log_progress.contains("Finished 1st pass mapping"));
    assert_eq!(result.killed_read_command_pids, vec![123]);
    assert_eq!(
        result.reopened_reads.opened_inputs,
        vec!["r1.fq".to_string()]
    );
    assert!(genome.genome_out.conv_yes);
    assert_eq!(genome.sjdb_n, 0);

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn parameters_read_sam_header_consumes_only_leading_extra_file_headers() {
    let mut p = Parameters {
        read_in_0: "@HD\tVN:1.6\n@SQ\tSN:chr1\n@RG\tID:rg1\n@CO\tcomment\nread1\t4\t*\n@RG\tlate\n"
            .to_string(),
        ..Default::default()
    };

    let commands =
        parameters_readsamheader_l6_parameters_readsamheader(&mut p, "", &[], &[]).unwrap();

    assert!(commands.is_empty());
    assert_eq!(p.sam_header_extra, "@RG\tID:rg1\n@CO\tcomment\n");
    assert_eq!(&p.read_in_0[p.read_in_0_pos..], "read1\t4\t*\n@RG\tlate\n");
}

#[test]
fn parameters_read_sam_header_command_mode_records_commands_and_pass2_skip() {
    let names = vec!["a.fq".to_string(), "b.fq".to_string()];
    let outputs = vec![
        "@HD\tVN:1.6\n@RG\tID:a\nSEQ\n".to_string(),
        "@SQ\tSN:chr1\n@CO\tb\n".to_string(),
    ];
    let mut p = Parameters {
        out_file_tmp: "/tmp/star/".to_string(),
        ..Default::default()
    };

    let commands =
        parameters_readsamheader_l6_parameters_readsamheader(&mut p, "zcat", &names, &outputs)
            .unwrap();

    assert_eq!(
        commands,
        vec![
            "zcat   a.fq > /tmp/star/tmp.fifo.header&",
            "zcat   b.fq > /tmp/star/tmp.fifo.header&"
        ]
    );
    assert_eq!(p.sam_header_extra, "@RG\tID:a\n@CO\tb\n");

    let mut p_pass2 = Parameters {
        out_file_tmp: "/tmp/star/".to_string(),
        two_pass_pass2: true,
        ..Default::default()
    };
    parameters_readsamheader_l6_parameters_readsamheader(&mut p_pass2, "zcat", &names, &outputs)
        .unwrap();
    assert!(p_pass2.sam_header_extra.is_empty());
}

#[test]
fn parameters_close_reads_files_closes_open_streams_and_reports_positive_pids() {
    let mut p = Parameters {
        read_files_in: vec![
            "r1.fq".to_string(),
            "r2.fq".to_string(),
            "r3.fq".to_string(),
        ],
        read_in_open: vec![true, false, true],
        read_files_command_pid: vec![123, 0, -7],
        ..Default::default()
    };

    let killed = parameters_closereadsfiles_l5_parameters_closereadsfiles(&mut p);

    assert_eq!(p.read_in_open, vec![false, false, false]);
    assert_eq!(killed, vec![123]);
}

#[test]
fn parameters_close_reads_files_uses_read_files_in_as_loop_bound() {
    let mut p = Parameters {
        read_files_in: vec!["r1.fq".to_string()],
        read_in_open: vec![true, true],
        read_files_command_pid: vec![5, 6],
        ..Default::default()
    };

    let killed = parameters_closereadsfiles_l5_parameters_closereadsfiles(&mut p);

    assert_eq!(p.read_in_open, vec![false, true]);
    assert_eq!(killed, vec![5]);
}

#[test]
fn parameters_read_files_init_expands_fastx_files_rg_and_command() {
    let mut p = Parameters {
        read_files_type: vec!["Fastx".to_string()],
        read_files_prefix: "data/".to_string(),
        read_files_manifest: vec!["-".to_string()],
        read_files_in: vec!["r1a.fq,r1b.fq,".to_string(), "r2a.fq,r2b.fq".to_string()],
        out_sam_attr_rgline: vec![
            "ID:rg1".to_string(),
            "SM:s1".to_string(),
            ",".to_string(),
            "ID:rg2".to_string(),
            "SM:s2".to_string(),
        ],
        read_files_command: vec!["-".to_string()],
        ..Default::default()
    };

    let log = parameters_readfilesinit_l8_parameters_readfilesinit(&mut p, None).unwrap();

    assert_eq!(p.read_files_type_n, 1);
    assert_eq!(p.read_files_prefix_final, "data/");
    assert_eq!(
        p.read_files_names,
        vec![
            vec!["data/r1a.fq".to_string(), "data/r1b.fq".to_string()],
            vec!["data/r2a.fq".to_string(), "data/r2b.fq".to_string()]
        ]
    );
    assert_eq!(p.read_files_n, 2);
    assert_eq!(
        p.out_sam_attr_rgline_split,
        vec!["ID:rg1\tSM:s1".to_string(), "ID:rg2\tSM:s2".to_string()]
    );
    assert_eq!(
        p.out_sam_attr_rg,
        vec!["rg1".to_string(), "rg2".to_string()]
    );
    assert_eq!(p.read_files_command_string, "cat   ");
    assert_eq!(p.read_nends, 2);
    assert_eq!(p.read_nmates, 2);
    assert_eq!(log, "Number of fastq files for each mate = 2\n");
}

#[test]
fn parameters_read_files_init_parses_sam_manifest_and_tag_keep() {
    let mut p = Parameters {
        read_files_type: vec!["SAM".to_string(), "SE".to_string()],
        read_files_prefix: "-".to_string(),
        read_files_manifest: vec!["manifest.tsv".to_string()],
        read_files_sam_attr_keep_in: vec!["NM".to_string(), "MD".to_string()],
        read_files_command: vec!["gzip".to_string(), "-cd".to_string()],
        ..Default::default()
    };
    let manifest = "a.sam\t-\tRG1\tSM:s1\n\nb.sam\t-\tID:RG2\tSM:s2\n";

    let log = parameters_readfilesinit_l8_parameters_readfilesinit(&mut p, Some(manifest)).unwrap();

    assert_eq!(p.read_files_type_n, 10);
    assert!(
        p.read_files_sam_attr_keep
            .contains(&u16::from_ne_bytes(*b"NM"))
    );
    assert!(
        p.read_files_sam_attr_keep
            .contains(&u16::from_ne_bytes(*b"MD"))
    );
    assert_eq!(
        p.read_files_names,
        vec![vec!["a.sam".to_string(), "b.sam".to_string()]]
    );
    assert_eq!(
        p.out_sam_attr_rgline_split,
        vec!["ID:RG1\tSM:s1".to_string(), "ID:RG2\tSM:s2".to_string()]
    );
    assert_eq!(
        p.out_sam_attr_rg,
        vec!["RG1".to_string(), "RG2".to_string()]
    );
    assert_eq!(p.read_files_command_string, "gzip   -cd   ");
    assert_eq!(p.read_nends, 1);
    assert_eq!(p.read_nmates, 1);
    assert!(
        log.contains(
            "Reading input file names and read groups from readFileManifest manifest.tsv\n"
        )
    );
    assert!(log.ends_with("Number of fastq files for each mate = 2\n"));
}

#[test]
fn parameters_open_reads_files_opens_direct_inputs_and_reports_missing_file() {
    let mut p = Parameters {
        read_files_command_string: String::new(),
        read_files_prefix_final: "data/".to_string(),
        read_files_in: vec!["r1.fq".to_string(), "r2.fq".to_string()],
        read_in_open: vec![true, false],
        read_files_command_pid: vec![99, 99],
        ..Default::default()
    };
    let existing =
        std::collections::BTreeSet::from(["data/r1.fq".to_string(), "data/r2.fq".to_string()]);

    let out = parameters_openreadsfiles_l5_parameters_openreadsfiles(&mut p, &existing).unwrap();

    assert_eq!(out.opened_inputs, vec!["data/r1.fq", "data/r2.fq"]);
    assert_eq!(p.read_in_open, vec![true, true]);
    assert_eq!(p.read_files_command_pid, vec![0, 0]);
    assert_eq!(p.read_files_index, 0);

    let mut missing = p.clone();
    missing.read_files_in = vec!["missing.fq".to_string()];
    let err = parameters_openreadsfiles_l5_parameters_openreadsfiles(&mut missing, &existing)
        .unwrap_err();
    assert!(err.contains("could not open readFilesIn=data/missing.fq"));
}

#[test]
fn parameters_open_reads_files_builds_fifo_commands_for_preprocessing() {
    let mut p = Parameters {
        read_files_command_string: "zcat   ".to_string(),
        read_files_names: vec![
            vec!["r1a.fq.gz".to_string(), "r1b.fq.gz".to_string()],
            vec!["r2a.fq.gz".to_string(), "r2b.fq.gz".to_string()],
        ],
        read_files_n: 2,
        out_file_tmp: "/tmp/star/".to_string(),
        sys_shell: "/bin/sh".to_string(),
        ..Default::default()
    };
    let existing = std::collections::BTreeSet::from([
        "r1a.fq.gz".to_string(),
        "r1b.fq.gz".to_string(),
        "r2a.fq.gz".to_string(),
        "r2b.fq.gz".to_string(),
    ]);

    let out = parameters_openreadsfiles_l5_parameters_openreadsfiles(&mut p, &existing).unwrap();

    assert_eq!(
        p.read_files_in_tmp,
        vec![
            "/tmp/star/tmp.fifo.read1".to_string(),
            "/tmp/star/tmp.fifo.read2".to_string()
        ]
    );
    assert_eq!(p.read_files_command_pid, vec![1, 2]);
    assert_eq!(p.read_in_open, vec![true, true]);
    assert_eq!(
        out.reads_command_file_names,
        vec![
            "/tmp/star//readsCommand_read1".to_string(),
            "/tmp/star//readsCommand_read2".to_string()
        ]
    );
    assert_eq!(
        out.reads_command_file_contents[0],
        concat!(
            "#!/bin/sh\n",
            "exec > \"/tmp/star/tmp.fifo.read1\"\n",
            "echo FILE 0\n",
            "zcat      \"r1a.fq.gz\"\n",
            "echo FILE 1\n",
            "zcat      \"r1b.fq.gz\"\n"
        )
    );
    assert!(out.log_main.contains("Input read files for mate 1"));
    assert_eq!(
        out.opened_inputs,
        vec!["/tmp/star/tmp.fifo.read1", "/tmp/star/tmp.fifo.read2"]
    );
}

#[test]
fn chimeric_detection_constructor_preserves_initializer_fields() {
    let p = Parameters {
        run_thread_n: 3,
        ..Default::default()
    };
    let ra = ReadAlign {
        l_read: 101,
        ..Default::default()
    };
    let out_gen = Genome {
        n_genome: 99,
        ..Default::default()
    };
    let tr_all = vec![
        vec![Transcript {
            max_score: 7,
            ..Default::default()
        }],
        vec![Transcript {
            max_score: 3,
            ..Default::default()
        }],
    ];

    let cd = chimericdetection_l3_chimericdetection_chimericdetection(
        p.clone(),
        tr_all.clone(),
        vec![1, 1],
        [vec![0, 1], vec![3, 2]],
        out_gen.clone(),
        true,
        ra.clone(),
    );

    assert_eq!(cd.p, p);
    assert_eq!(cd.ra, Some(ra));
    assert_eq!(cd.tr_all, tr_all);
    assert_eq!(cd.n_w, 2);
    assert_eq!(cd.n_win_tr, vec![1, 1]);
    assert_eq!(cd.read1, [vec![0, 1], vec![3, 2]]);
    assert_eq!(cd.out_gen, out_gen);
    assert!(cd.ostream_chim_junction_attached);
    assert!(cd.chim_aligns.is_empty());
}

#[test]
fn find_char_returns_first_offset_without_boundary_checks() {
    assert_eq!(clipmate_clipchunk_l55_findchar(b"abc\ndef", b'\n'), 3);
    assert_eq!(clipmate_clipchunk_l55_findchar(b"abc\ndef", b'd'), 4);
}

#[test]
fn solo_barcode_position_parser_matches_anchor_assignment() {
    let mut barcode = SoloBarcode::default();
    solobarcode_l37_solobarcode_extractpositionsfromstring(&mut barcode, "1_25_2_30");
    assert_eq!(barcode.anchor_type, [1, 2]);
    assert_eq!(barcode.anchor_dist, [25, 30]);
}

#[test]
fn solo_barcode_extract_barcode_uses_anchor_positions_like_original() {
    let seq = "AACCGGTTAACC";
    let qual = "abcdefghijkl";

    let barcode = SoloBarcode {
        anchor_type: [0, 2],
        anchor_dist: [1, -1],
        adapter_length: 4,
        ..Default::default()
    };
    assert_eq!(
        solobarcode_extractbarcode_l4_solobarcode_extractbarcode(&barcode, seq, qual, 6),
        Some(("ACCGG".to_string(), "bcdef".to_string()))
    );

    let adapter_end = SoloBarcode {
        anchor_type: [3, 1],
        anchor_dist: [0, 0],
        adapter_length: 3,
        ..Default::default()
    };
    assert_eq!(
        solobarcode_extractbarcode_l4_solobarcode_extractbarcode(&adapter_end, seq, qual, 2),
        Some(("GGTTAACC".to_string(), "efghijkl".to_string()))
    );

    let negative_start = SoloBarcode {
        anchor_type: [0, 1],
        anchor_dist: [-1, 0],
        ..Default::default()
    };
    assert_eq!(
        solobarcode_extractbarcode_l4_solobarcode_extractbarcode(&negative_start, seq, qual, 0),
        None
    );

    let reversed = SoloBarcode {
        anchor_type: [1, 0],
        anchor_dist: [0, 0],
        ..Default::default()
    };
    assert_eq!(
        solobarcode_extractbarcode_l4_solobarcode_extractbarcode(&reversed, seq, qual, 0),
        None
    );
}

#[test]
fn parameters_solo_umi_swap_halves_matches_original_bit_operations() {
    let p_solo = ParametersSolo {
        umi_l: 4,
        umi_mask_low: 0x0f,
        ..Default::default()
    };

    let mut umi = 0xab;
    parameterssolo_l496_parameterssolo_umiswaphalves(&p_solo, &mut umi);
    assert_eq!(umi, 0xba);

    let mut wider = 0x1234;
    parameterssolo_l496_parameterssolo_umiswaphalves(
        &ParametersSolo {
            umi_l: 8,
            umi_mask_low: 0x00ff,
            ..Default::default()
        },
        &mut wider,
    );
    assert_eq!(wider, 0x3412);
}

#[test]
fn umi_dedup_and_multimappers_initialize_match_type_tables() {
    let mut umi = UMIdedup {
        types_in: vec!["Exact".to_string(), "1MM_CR".to_string()],
        count_ind_i: [7; 6],
        yes_b: [true; 6],
        ..Default::default()
    };
    parameterssolo_l585_umidedup_initialize(&mut umi, &ParametersSolo::default()).unwrap();
    assert_eq!(umi.types, vec![1, 4]);
    assert_eq!(umi.yes_n, 2);
    assert!(umi.yes_b[1]);
    assert!(umi.yes_b[4]);
    assert!(!umi.yes_b[0]);
    assert_eq!(umi.count_ind_i[1], 1);
    assert_eq!(umi.count_ind_i[4], 2);
    assert_eq!(umi.count_ind_i[0], u32::MAX);
    assert_eq!(umi.type_main, 1);
    assert_eq!(umi.count_ind_main, 1);

    let mut smartseq_bad = UMIdedup {
        types_in: vec!["1MM_All".to_string()],
        ..Default::default()
    };
    assert!(
        parameterssolo_l585_umidedup_initialize(
            &mut smartseq_bad,
            &ParametersSolo {
                solo_type: 4,
                ..Default::default()
            },
        )
        .unwrap_err()
        .contains("SmartSeq")
    );

    let mut bad_umi = UMIdedup {
        types_in: vec!["Bad".to_string()],
        ..Default::default()
    };
    assert!(
        parameterssolo_l585_umidedup_initialize(&mut bad_umi, &ParametersSolo::default())
            .unwrap_err()
            .contains("--soloUMIdedup")
    );

    let p_solo = ParametersSolo {
        umi_dedup: UMIdedup { yes_n: 2, ..umi },
        ..Default::default()
    };
    let mut multi = MultiMappers {
        types_in: vec![
            "Unique".to_string(),
            "Uniform".to_string(),
            "EM".to_string(),
        ],
        count_ind_i: [9; 5],
        yes_b: [true; 5],
        yes_multi: true,
        ..Default::default()
    };
    parameterssolo_l624_multimappers_initialize(&mut multi, &p_solo).unwrap();
    assert_eq!(multi.types, vec![1, 4]);
    assert_eq!(multi.yes_n, 2);
    assert!(multi.yes_b[1]);
    assert!(multi.yes_b[4]);
    assert!(!multi.yes_b[0]);
    assert_eq!(multi.count_ind_i[1], 1);
    assert_eq!(multi.count_ind_i[4], 3);
    assert_eq!(multi.count_ind_i[2], u32::MAX);
    assert_eq!(multi.type_main, 1);
    assert_eq!(multi.count_ind_main, 1);
    assert!(multi.yes_multi);

    let mut unique_only = MultiMappers {
        types_in: vec!["Unique".to_string()],
        yes_multi: true,
        ..Default::default()
    };
    parameterssolo_l624_multimappers_initialize(&mut unique_only, &p_solo).unwrap();
    assert_eq!(unique_only.yes_n, 0);
    assert!(!unique_only.yes_multi);

    let mut bad_multi = MultiMappers {
        types_in: vec!["Bad".to_string()],
        ..Default::default()
    };
    assert!(
        parameterssolo_l624_multimappers_initialize(&mut bad_multi, &p_solo)
            .unwrap_err()
            .contains("--soloMultiMappers")
    );
}

#[test]
fn parameters_solo_init_cbmatchwl_sets_flags_and_rejects_incompatible_modes() {
    let mut exact = ParametersSolo {
        type_str: "CB_UMI_Simple".to_string(),
        cb_match_wl: CBMatchWL {
            type_: "Exact".to_string(),
            mm1: true,
            mm1_multi: true,
            mm1_multi_pc: true,
            mm1_multi_nbase: true,
            edit_dist_2: true,
            ..Default::default()
        },
        ..Default::default()
    };
    parameterssolo_l675_parameterssolo_init_cbmatchwl(&mut exact).unwrap();
    assert!(exact.cb_match_wl.one_exact);
    assert!(!exact.cb_match_wl.mm1);
    assert!(!exact.cb_match_wl.mm1_multi);
    assert!(!exact.cb_match_wl.mm1_multi_pc);
    assert!(!exact.cb_match_wl.mm1_multi_nbase);
    assert!(!exact.cb_match_wl.edit_dist_2);

    let mut one_mm_multi = ParametersSolo {
        type_str: "CB_UMI_Simple".to_string(),
        cb_match_wl: CBMatchWL {
            type_: "1MM_multi".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    parameterssolo_l675_parameterssolo_init_cbmatchwl(&mut one_mm_multi).unwrap();
    assert!(one_mm_multi.cb_match_wl.mm1);
    assert!(one_mm_multi.cb_match_wl.mm1_multi);
    assert!(one_mm_multi.cb_match_wl.one_exact);
    assert!(!one_mm_multi.cb_match_wl.mm1_multi_pc);

    let mut nbase_pc = ParametersSolo {
        type_str: "CB_UMI_Simple".to_string(),
        cb_match_wl: CBMatchWL {
            type_: "1MM_multi_Nbase_pseudocounts".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    parameterssolo_l675_parameterssolo_init_cbmatchwl(&mut nbase_pc).unwrap();
    assert!(nbase_pc.cb_match_wl.mm1);
    assert!(nbase_pc.cb_match_wl.mm1_multi);
    assert!(nbase_pc.cb_match_wl.mm1_multi_pc);
    assert!(nbase_pc.cb_match_wl.mm1_multi_nbase);
    assert!(!nbase_pc.cb_match_wl.one_exact);

    let mut edit_dist = ParametersSolo {
        type_str: "CB_UMI_Complex".to_string(),
        cb_match_wl: CBMatchWL {
            type_: "EditDist_2".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    parameterssolo_l675_parameterssolo_init_cbmatchwl(&mut edit_dist).unwrap();
    assert!(edit_dist.cb_match_wl.edit_dist_2);

    let mut incompatible = ParametersSolo {
        type_str: "CB_samTagOut".to_string(),
        cb_match_wl: CBMatchWL {
            type_: "1MM_multi".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(
        parameterssolo_l675_parameterssolo_init_cbmatchwl(&mut incompatible)
            .unwrap_err()
            .contains("does not work")
    );

    let mut bad = ParametersSolo {
        type_str: "CB_UMI_Simple".to_string(),
        cb_match_wl: CBMatchWL {
            type_: "Bad".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(
        parameterssolo_l675_parameterssolo_init_cbmatchwl(&mut bad)
            .unwrap_err()
            .contains("unrecognized option")
    );
}

#[test]
fn parameters_solo_initialize_simple_reads_whitelist_and_sets_feature_state() {
    let mut p = Parameters {
        run_thread_n: 2,
        read_files_type_n: 1,
        read_nends: 2,
        read_nmates: 2,
        p_solo: ParametersSolo {
            type_str: "CB_UMI_Simple".to_string(),
            barcode_read: 0,
            cb_s: 1,
            cb_l: 4,
            umi_s: 5,
            umi_l: 4,
            b_l: 1,
            solo_cb_whitelist: vec!["wl.txt".to_string()],
            features: vec![SOLO_FEATURE_GENE as u32, SOLO_FEATURE_SJ as u32],
            cell_filter: SoloCellFilter {
                type_: vec!["None".to_string()],
                ..Default::default()
            },
            cb_match_wl: CBMatchWL {
                type_: "Exact".to_string(),
                ..Default::default()
            },
            umi_dedup: UMIdedup {
                types_in: vec!["Exact".to_string()],
                ..Default::default()
            },
            multi_map: MultiMappers {
                types_in: vec!["Unique".to_string()],
                ..Default::default()
            },
            umi_filtering: SoloUmiFiltering {
                type_: vec!["-".to_string()],
                ..Default::default()
            },
            read_stats_type: "Standard".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };

    let log = parameterssolo_l10_parameterssolo_initialize(
        &mut p,
        &[("wl.txt".to_string(), "ACGT\nTGCA\nACGT\nNNNN\n".to_string())],
        0,
    )
    .unwrap();

    assert_eq!(p.p_solo.solo_type, SOLO_TYPE_CB_UMI_SIMPLE);
    assert!(p.p_solo.yes);
    assert_eq!(p.p_solo.redistr_reads_nfiles, 6);
    assert_eq!(p.read_nmates, 1);
    assert_eq!(p.p_solo.barcode_read, 1);
    assert_eq!(p.p_solo.cbumi_l, 8);
    assert_eq!(p.p_solo.b_l, 8);
    assert_eq!(p.p_solo.barcode_start, 0);
    assert_eq!(p.p_solo.barcode_end, 7);
    assert_eq!(p.p_solo.cb_wl_size, 2);
    assert_eq!(
        p.p_solo.cb_wl_str,
        vec!["ACGT".to_string(), "TGCA".to_string()]
    );
    assert!(log.contains("Number of CBs in the whitelist = 2"));
    assert!(log.contains("WARNING: CB whitelist sequence contains non-ACGT base"));
    assert!(p.p_solo.feature_yes[SOLO_FEATURE_GENE as usize]);
    assert!(p.p_solo.feature_yes[SOLO_FEATURE_SJ as usize]);
    assert_eq!(p.p_solo.feature_ind[SOLO_FEATURE_SJ as usize], 0);
    assert_eq!(p.p_solo.feature_ind[SOLO_FEATURE_GENE as usize], 1);
    assert_eq!(p.p_solo.feature_first, SOLO_FEATURE_GENE);
    assert_eq!(p.p_solo.umi_dedup.yes_n, 1);
    assert_eq!(p.p_solo.umi_dedup.count_ind_i[1], 1);
    assert!(!p.p_solo.multi_map.yes_multi);
    assert!(p.p_solo.read_stats_yes[SOLO_FEATURE_GENE as usize]);
    assert!(!p.p_solo.read_stats_yes[SOLO_FEATURE_SJ as usize]);
    assert!(p.quant_gene_yes);
    assert!(p.quant_yes);
}

#[test]
fn parameters_solo_initialize_sam_tags_and_multimap_errors_match_original_guards() {
    let mut sam_p = Parameters {
        read_files_type_n: 10,
        read_nends: 1,
        read_nmates: 1,
        p_solo: ParametersSolo {
            type_str: "CB_samTagOut".to_string(),
            cb_s: 1,
            cb_l: 4,
            umi_s: 5,
            umi_l: 4,
            b_l: 1,
            solo_cb_whitelist: vec!["None".to_string()],
            sam_attr_barcode_seq: vec!["CR".to_string()],
            sam_attr_barcode_qual: vec!["-".to_string()],
            features: vec![SOLO_FEATURE_GENE as u32],
            cell_filter: SoloCellFilter {
                type_: vec!["None".to_string()],
                ..Default::default()
            },
            cb_match_wl: CBMatchWL {
                type_: "Exact".to_string(),
                ..Default::default()
            },
            umi_dedup: UMIdedup {
                types_in: vec!["Exact".to_string()],
                ..Default::default()
            },
            multi_map: MultiMappers {
                types_in: vec!["Unique".to_string()],
                ..Default::default()
            },
            umi_filtering: SoloUmiFiltering {
                type_: vec!["-".to_string()],
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };

    let log = parameterssolo_l10_parameterssolo_initialize(&mut sam_p, &[], 0).unwrap();
    assert_eq!(sam_p.p_solo.solo_type, SOLO_TYPE_CB_SAM_TAG_OUT);
    assert_eq!(sam_p.p_solo.sam_attr_barcode_seq, vec!["\tCR".to_string()]);
    assert!(sam_p.p_solo.sam_attr_barcode_qual.is_empty());
    assert!(log.contains("qualities for barcode read will be replaced with 'H'"));

    let mut bad_multi = Parameters {
        read_files_type_n: 1,
        read_nends: 1,
        read_nmates: 1,
        p_solo: ParametersSolo {
            type_str: "SmartSeq".to_string(),
            features: vec![SOLO_FEATURE_GENE as u32],
            cell_filter: SoloCellFilter {
                type_: vec!["None".to_string()],
                ..Default::default()
            },
            cb_match_wl: CBMatchWL {
                type_: "Exact".to_string(),
                ..Default::default()
            },
            umi_dedup: UMIdedup {
                types_in: vec!["Exact".to_string()],
                ..Default::default()
            },
            multi_map: MultiMappers {
                types_in: vec!["Uniform".to_string()],
                ..Default::default()
            },
            umi_filtering: SoloUmiFiltering {
                type_: vec!["-".to_string()],
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    let err = parameterssolo_l10_parameterssolo_initialize(&mut bad_multi, &[], 0).unwrap_err();
    assert!(err.contains("multimapping options do not work for --soloType SmartSeq"));
}

#[test]
fn parameters_solo_cell_filtering_parses_modes_and_defaults() {
    let mut cr = ParametersSolo {
        cell_filter: star_rs::generated::structs::SoloCellFilter {
            type_: vec!["CellRanger2.2".to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    let log = parameterssolo_l533_parameterssolo_cellfiltering(&mut cr).unwrap();
    assert!(log.contains("using hardcoded filtering parameters"));
    assert_eq!(cr.cell_filter.knee.n_expected_cells, 3000.0);
    assert_eq!(cr.cell_filter.knee.max_percentile, 0.99);
    assert_eq!(cr.cell_filter.knee.max_min_ratio, 10.0);

    let mut ed = ParametersSolo {
        cell_filter: star_rs::generated::structs::SoloCellFilter {
            type_: vec![
                "EmptyDrops_CR".to_string(),
                "4000".to_string(),
                "0.95".to_string(),
                "8".to_string(),
                "100".to_string(),
                "200".to_string(),
                "20".to_string(),
                "0.05".to_string(),
                "300".to_string(),
                "0.02".to_string(),
                "500".to_string(),
            ],
            ..Default::default()
        },
        ..Default::default()
    };
    let log = parameterssolo_l533_parameterssolo_cellfiltering(&mut ed).unwrap();
    assert!(log.contains("EmptyDrops_CR filtering parameters"));
    assert_eq!(ed.cell_filter.knee.n_expected_cells, 4000.0);
    assert_eq!(ed.cell_filter.ed_cr.ind_min, 100);
    assert_eq!(ed.cell_filter.ed_cr.ind_max, 200);
    assert_eq!(ed.cell_filter.ed_cr.umi_min, 20);
    assert_eq!(ed.cell_filter.ed_cr.umi_min_frac_median, 0.05);
    assert_eq!(ed.cell_filter.ed_cr.cand_max_n, 300);
    assert_eq!(ed.cell_filter.ed_cr.fdr, 0.02);
    assert_eq!(ed.cell_filter.ed_cr.sim_n, 500);

    let mut top = ParametersSolo {
        cell_filter: star_rs::generated::structs::SoloCellFilter {
            type_: vec!["TopCells".to_string(), "123".to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    assert_eq!(
        parameterssolo_l533_parameterssolo_cellfiltering(&mut top).unwrap(),
        ""
    );
    assert_eq!(top.cell_filter.top_cells, 123);

    let mut bad = ParametersSolo {
        cell_filter: star_rs::generated::structs::SoloCellFilter {
            type_: vec!["Bad".to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(
        parameterssolo_l533_parameterssolo_cellfiltering(&mut bad)
            .unwrap_err()
            .contains("unrecognized option")
    );
}

#[test]
fn solo_read_barcode_match_cb_to_wl_matches_exact_and_no_whitelist_modes() {
    let no_wl = ParametersSolo {
        cb_wl_yes: false,
        ..Default::default()
    };
    assert_eq!(
        soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(&no_wl, "ACGT", "IIII", &[]),
        (0, vec![27], "27".to_string())
    );
    assert_eq!(
        soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(&no_wl, "ACNT", "IIII", &[]),
        (-2, Vec::new(), String::new())
    );

    let exact = ParametersSolo {
        cb_wl_yes: true,
        cb_match_wl: CBMatchWL {
            mm1: false,
            ..Default::default()
        },
        ..Default::default()
    };
    assert_eq!(
        soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(
            &exact,
            "ACGT",
            "IIII",
            &[2, 27, 33],
        ),
        (0, vec![1], "1".to_string())
    );
    assert_eq!(
        soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(
            &exact,
            "ACGA",
            "IIII",
            &[2, 27, 33],
        ),
        (-1, Vec::new(), String::new())
    );
}

#[test]
fn solo_read_barcode_match_cb_to_wl_handles_nbase_and_mismatch_corrections() {
    let one_mm_multi = ParametersSolo {
        cb_wl_yes: true,
        cb_match_wl: CBMatchWL {
            mm1: true,
            mm1_multi: true,
            mm1_multi_nbase: true,
            ..Default::default()
        },
        ..Default::default()
    };

    assert_eq!(
        soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(
            &one_mm_multi,
            "ACNT",
            "ABCD",
            &[19, 23, 27, 31],
        ),
        (4, vec![0, 1, 2, 3], " 0 C 1 C 2 C 3 C".to_string())
    );

    assert_eq!(
        soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(
            &one_mm_multi,
            "ACGT",
            "WXYZ",
            &[11, 19, 25, 26, 30],
        ),
        (4, vec![3, 2, 1, 0], " 3 Z 2 Z 1 Y 0 X".to_string())
    );
}

#[test]
fn solo_read_barcode_match_cb_to_wl_rejects_disallowed_multiple_matches() {
    let one_mm_single = ParametersSolo {
        cb_wl_yes: true,
        cb_match_wl: CBMatchWL {
            mm1: true,
            mm1_multi: false,
            mm1_multi_nbase: false,
            ..Default::default()
        },
        ..Default::default()
    };

    assert_eq!(
        soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(
            &one_mm_single,
            "ACNT",
            "ABCD",
            &[19, 23],
        ),
        (-3, Vec::new(), String::new())
    );

    assert_eq!(
        soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(
            &one_mm_single,
            "ACGT",
            "WXYZ",
            &[11, 19],
        ),
        (-3, Vec::new(), String::new())
    );
}

#[test]
fn solo_barcode_sort_whitelist_matches_offsets_and_dedup() {
    let mut barcode = SoloBarcode {
        wl: vec![vec![], vec![3, 1, 3], vec![], vec![9, 8, 9]],
        ..Default::default()
    };
    solobarcode_l9_solobarcode_sortwhitelist(&mut barcode, false);
    assert_eq!(barcode.wl[1], vec![1, 3]);
    assert_eq!(barcode.wl[3], vec![8, 9]);
    assert_eq!(barcode.wl_add, vec![0, 0, 2, 2]);
    assert_eq!(barcode.min_len, 1);
    assert_eq!(barcode.total_size, 4);

    let mut barcode_ed = SoloBarcode {
        wl: vec![vec![], vec![0], vec![1]],
        ..Default::default()
    };
    solobarcode_l9_solobarcode_sortwhitelist(&mut barcode_ed, true);
    assert_eq!(barcode_ed.wl_add, vec![0, 0, 1]);
    assert_eq!(barcode_ed.total_size, 2);
    assert_eq!(barcode_ed.wl_ed.len(), 3);
    assert!(!barcode_ed.wl_ed[1].is_empty());
    assert_eq!(barcode_ed.wl_ed[1].len(), barcode_ed.wl_ed_ind[1].len());
}

#[test]
fn whitelist_mismatch_expansion_keeps_unique_matches_only() {
    let mut edited = Vec::new();
    let mut indices = Vec::new();
    solobarcode_l47_wladdmismatches(1, 1, &[0], &mut edited, &mut indices);
    assert_eq!(edited, vec![0, 1, 2, 3]);
    assert_eq!(indices, vec![0, 0, 0, 0]);

    solobarcode_l47_wladdmismatches(1, 1, &[0, 1], &mut edited, &mut indices);
    assert_eq!(edited, vec![0, 1]);
    assert_eq!(indices, vec![0, 1]);
}

#[test]
fn sys_remove_dir_removes_nested_directory_tree() {
    let dir = std::env::temp_dir().join(format!("star-rs-rmdir-{}", std::process::id()));
    let nested = dir.join("a/b");
    std::fs::create_dir_all(&nested).unwrap();
    std::fs::write(nested.join("x.txt"), b"x").unwrap();

    sysremovedir_l25_sysremovedir(Path::new(&dir)).unwrap();
    assert!(!dir.exists());
}

#[test]
fn stream_open_wrappers_match_original_modes_and_error_messages() {
    use std::io::{Read, Seek, SeekFrom, Write};

    let dir = std::env::temp_dir().join(format!("star-rs-stream-{}", std::process::id()));
    let file = dir.join("nested/out.txt");
    let log = streamfuns_l10_createdirectory(file.to_str().unwrap(), 0o700, "--outFileNamePrefix")
        .unwrap();
    assert!(log.contains("directory"));
    assert!(dir.join("nested").is_dir());

    let mut out = streamfuns_l91_ofstropen(file.to_str().unwrap(), "ERR").unwrap();
    out.write_all(b"abcdef").unwrap();
    drop(out);
    let mut out2 = streamfuns_l91_ofstropen(file.to_str().unwrap(), "ERR").unwrap();
    out2.write_all(b"xy").unwrap();
    drop(out2);
    assert_eq!(std::fs::read(&file).unwrap(), b"xy");

    let mut f = streamfuns_l102_fstropen(file.to_str().unwrap(), "ERR", false).unwrap();
    f.seek(SeekFrom::End(0)).unwrap();
    f.write_all(b"z").unwrap();
    drop(f);
    assert_eq!(std::fs::read(&file).unwrap(), b"xyz");

    let mut f_trunc = streamfuns_l102_fstropen(file.to_str().unwrap(), "ERR", true).unwrap();
    f_trunc.write_all(b"q").unwrap();
    drop(f_trunc);
    assert_eq!(std::fs::read(&file).unwrap(), b"q");

    let mut input = streamfuns_l124_ifstropen(file.to_str().unwrap(), "ERR", "").unwrap();
    let mut text = String::new();
    input.read_to_string(&mut text).unwrap();
    assert_eq!(text, "q");

    let genome_dir = dir.join("genome");
    std::fs::create_dir_all(&genome_dir).unwrap();
    std::fs::write(genome_dir.join("Genome"), b"g").unwrap();
    assert!(
        streamfuns_l139_ifstropengenomefile("Genome", "ERR", genome_dir.to_str().unwrap()).is_ok()
    );
    let err = streamfuns_l124_ifstropen("/definitely/not/here", "ERR", "extra").unwrap_err();
    assert!(err.contains("*INPUT FILE*"));
    assert!(err.contains("extra"));

    sysremovedir_l25_sysremovedir(Path::new(&dir)).unwrap();
}

#[test]
fn genome_open_stream_matches_size_paths_and_fatal_errors() {
    use std::io::Read;

    let dir = std::env::temp_dir().join(format!("star-rs-genome-open-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("Genome"), b"abcdef").unwrap();
    std::fs::write(dir.join("SA"), b"").unwrap();
    let genome = Genome {
        p_ge: ParametersGenome {
            g_dir: dir.to_str().unwrap().to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    let p = Parameters {
        run_thread_n: 2,
        ..Default::default()
    };

    let (mut file, size, log) = genome_l48_genome_openstream(&genome, "Genome", 99, &p, 0).unwrap();
    let mut text = String::new();
    file.read_to_string(&mut text).unwrap();
    assert_eq!(text, "abcdef");
    assert_eq!(size, 99);
    assert_eq!(log, "Genome: size given as a parameter = 99\n");

    let (_, detected_size, detected_log) =
        genome_l48_genome_openstream(&genome, "Genome", 0, &p, 0).unwrap();
    assert_eq!(detected_size, 6);
    assert_eq!(
        detected_log,
        "Checking Genome sizefile size: 6 bytes; state: good=1 eof=0 fail=0 bad=0\n"
    );

    let missing = genome_l48_genome_openstream(&genome, "Missing", 0, &p, 0).unwrap_err();
    assert_eq!(missing.error_int, EXIT_CODE_GENOME_FILES);
    assert!(missing.stream_out1.contains("could not open genome file"));
    assert!(missing.thread_mutex_locked);

    let empty = genome_l48_genome_openstream(&genome, "SA", 0, &p, 0).unwrap_err();
    assert_eq!(empty.error_int, 1);
    assert!(
        empty
            .stream_out1
            .contains("failed reading from genome file")
    );

    sysremovedir_l25_sysremovedir(Path::new(&dir)).unwrap();
}

#[test]
fn genome_handle_shared_memory_exception_matches_error_messages_and_exit_codes() {
    let p = Parameters::default();

    let exists = genome_l83_genome_handlesharedmemoryexception(
        &SharedMemoryException {
            error_code: SHM_EEXISTS,
            error_detail: 12,
        },
        345,
        &p,
        0,
    );
    assert_eq!(exists.error_int, EXIT_CODE_SHM);
    assert!(
        exists
            .stream_out1
            .contains("Shared memory error: 5, errno: 12(12)")
    );
    assert!(exists.stream_out1.contains("2000000345 bytes"));
    assert!(exists.stream_out1.contains("ulimit -v 2000000345"));

    let trunc = genome_l83_genome_handlesharedmemoryexception(
        &SharedMemoryException {
            error_code: SHM_EFTRUNCATE,
            error_detail: 0,
        },
        100,
        &p,
        0,
    );
    assert_eq!(trunc.error_int, EXIT_CODE_MEMORY_ALLOCATION);
    assert!(
        trunc
            .stream_out1
            .contains("ftruncate() error shared memory")
    );

    let unknown = genome_l83_genome_handlesharedmemoryexception(
        &SharedMemoryException {
            error_code: 777,
            error_detail: 0,
        },
        100,
        &p,
        0,
    );
    assert_eq!(unknown.error_int, EXIT_CODE_SHM);
    assert!(
        unknown
            .stream_out1
            .contains("There was an issue with the shared memory allocation")
    );
}

#[test]
fn genome_chr_info_load_reads_index_files_logs_and_marks_mito_names() {
    let dir = std::env::temp_dir().join(format!("star-rs-chr-info-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("chrName.txt"), "chr1\nchrM\nchr2\n\nignored\n").unwrap();
    std::fs::write(dir.join("chrLength.txt"), "100\n17\n250\n").unwrap();
    std::fs::write(dir.join("chrStart.txt"), "0\n112\n144\n400\n").unwrap();

    let mut genome = Genome {
        p_ge: ParametersGenome {
            g_dir: dir.to_str().unwrap().to_string(),
            chr_set_mito_strings: vec!["chrM".to_string(), "absent".to_string()],
            ..Default::default()
        },
        chr_name: vec!["old".to_string()],
        chr_start: vec![9],
        chr_length: vec![9],
        ..Default::default()
    };

    let log = genome_l139_genome_chrinfoload(&mut genome, &Parameters::default(), 0).unwrap();

    assert_eq!(genome.n_chr_real, 3);
    assert_eq!(genome.chr_name, vec!["chr1", "chrM", "chr2"]);
    assert_eq!(genome.chr_length, vec![100, 17, 250]);
    assert_eq!(genome.chr_start, vec![0, 112, 144, 400]);
    assert_eq!(genome.chr_name_index.get("chrM"), Some(&1));
    assert!(genome.p_ge.chr_set_mito.contains(&1));
    assert!(genome.p_ge.chr_set_mito.contains(&3));
    assert_eq!(
        log,
        "Number of real (reference) chromosomes= 3\n1\tchr1\t100\t0\n2\tchrM\t17\t112\n3\tchr2\t250\t144\n"
    );

    std::fs::remove_file(dir.join("chrLength.txt")).unwrap();
    let err = genome_l139_genome_chrinfoload(&mut genome, &Parameters::default(), 0).unwrap_err();
    assert_eq!(err.error_int, EXIT_CODE_INPUT_FILES);
    assert!(err.stream_out1.contains("could not open file"));

    sysremovedir_l25_sysremovedir(Path::new(&dir)).unwrap();
}

#[test]
fn genome_sequence_allocate_matches_spacing_initialization_and_ram_guard() {
    let p = Parameters {
        limit_genome_generate_ram: 1_000,
        ..Default::default()
    };
    let (n_alloc, g1, g_offset) = genome_l219_genome_genomesequenceallocate(50, &p, 0).unwrap();
    assert_eq!(n_alloc, 300);
    assert_eq!(g1.len(), 300);
    assert!(g1.iter().all(|&base| base == GENOME_SPACING_CHAR));
    assert_eq!(g_offset, 100);

    let err = genome_l219_genome_genomesequenceallocate(
        50,
        &Parameters {
            limit_genome_generate_ram: 399,
            ..Default::default()
        },
        0,
    )
    .unwrap_err();
    assert_eq!(err.error_int, EXIT_CODE_INPUT_FILES);
    assert!(
        err.stream_out1
            .contains("limitGenomeGenerateRAM=399is too small")
    );
    assert!(err.stream_out1.contains("not less than 400"));
}

#[test]
fn remove_file_or_dir_reports_bad_typeflag() {
    assert_eq!(
        sysremovedir_l7_removefileordir(Path::new("/tmp/nonexistent"), 99),
        -1
    );
}

#[test]
fn thread_control_constructor_initializes_chunk_counters() {
    let tc = threadcontrol_l3_threadcontrol_threadcontrol();
    assert_eq!(tc.chunk_in_n, 0);
    assert_eq!(tc.chunk_out_n, 0);
}

#[test]
fn map_threads_spawn_logs_create_process_join_and_reports_status_errors() {
    let mut processed = 0;
    let log = mapthreadsspawn_l6_mapthreadsspawn(3, &[0, 0, 0], &[0, 0, 0], || {
        processed += 1;
        Ok("process main\n".to_string())
    })
    .unwrap();

    assert_eq!(processed, 1);
    assert_eq!(
        log,
        "Created thread # 1\nCreated thread # 2\nprocess main\nJoined thread # 1\nJoined thread # 2\n"
    );

    let err = mapthreadsspawn_l6_mapthreadsspawn(3, &[0, 5, 0], &[0, 0, 0], || Ok(String::new()))
        .unwrap_err();
    assert_eq!(
        err,
        "EXITING because of FATAL ERROR: phtread error while creating thread # 1, error code: 5"
    );

    let err = mapthreadsspawn_l6_mapthreadsspawn(3, &[0, 0, 0], &[0, 0, 7], || Ok(String::new()))
        .unwrap_err();
    assert_eq!(
        err,
        "EXITING because of FATAL ERROR: phtread error while joining thread # 2, error code: 7"
    );
}

#[test]
fn splicegraph_constructor_and_destructor_allocate_original_buffers() {
    let super_trome = SuperTranscriptome {
        n: 3,
        sj_donor_nmax: 4,
        ..Default::default()
    };
    let ra = ReadAlign {
        l_read: 99,
        ..Default::default()
    };

    let mut graph = splicegraph_l8_splicegraph_splicegraph(super_trome.clone(), Some(ra.clone()));
    assert_eq!(graph.super_trome, super_trome);
    assert_eq!(graph.ra, Some(ra));
    assert_eq!(graph.super_tr_seed_count.len(), 6);
    assert_eq!(graph.scoring_matrix.len(), 6);
    assert!(
        graph
            .scoring_matrix
            .iter()
            .all(|row| row.len() == SPLICEGRAPH_MAX_SEQ_LENGTH)
    );
    assert_eq!(graph.score_two_columns[0].len(), SPLICEGRAPH_MAX_SEQ_LENGTH);
    assert_eq!(graph.score_two_columns[1].len(), SPLICEGRAPH_MAX_SEQ_LENGTH);
    assert_eq!(graph.sj_dindex.len(), 4);
    assert_eq!(graph.gap_penalty, -1);
    assert_eq!(graph.match_score, 1);
    assert_eq!(graph.mismatch_penalty, -1);

    splicegraph_l28_splicegraph_splicegraph(&mut graph);
    assert!(graph.super_tr_seed_count.is_empty());
    assert!(graph.scoring_matrix.is_empty());
    assert!(graph.score_two_columns[0].is_empty());
    assert!(graph.score_two_columns[1].is_empty());
    assert!(graph.direction_matrix.is_empty());
    assert!(graph.sj_dindex.is_empty());
}

#[test]
fn splicegraph_sw_score_spliced_aligns_unspliced_and_spliced_reads() {
    let mut graph = splicegraph_l8_splicegraph_splicegraph(
        SuperTranscriptome {
            n: 1,
            sj_donor_nmax: 1,
            ..Default::default()
        },
        None,
    );

    let mut cigar = Vec::new();
    let score = splicegraph_swscorespliced_l8_splicegraph_swscorespliced(
        &mut graph,
        b"ACG",
        3,
        &SuperTranscript {
            seq: b"ACG".to_vec(),
            length: 3,
            ..Default::default()
        },
        &mut cigar,
    );
    assert_eq!(score, 3);
    assert_eq!(cigar, vec![[BAM_CIGAR_M, 3]]);
    assert_eq!(graph.align_info.n_map, 3);
    assert_eq!(graph.align_info.n_mm, 0);
    assert_eq!(graph.align_info.a_start, [0, 0]);
    assert_eq!(graph.align_info.a_end, [2, 2]);

    let score = splicegraph_swscorespliced_l8_splicegraph_swscorespliced(
        &mut graph,
        b"AT",
        2,
        &SuperTranscript {
            seq: b"ACGT".to_vec(),
            length: 4,
            sj_c: vec![[0, 3, 0]],
            sj_donor: vec![0],
        },
        &mut cigar,
    );
    assert_eq!(score, 2);
    assert_eq!(
        cigar,
        vec![[BAM_CIGAR_M, 1], [BAM_CIGAR_N, 2], [BAM_CIGAR_M, 1]]
    );
    assert_eq!(graph.align_info.n_map, 2);
    assert_eq!(graph.align_info.n_d, 0);
    assert_eq!(graph.align_info.n_sj, 1);
    assert_eq!(graph.align_info.a_start, [0, 0]);
    assert_eq!(graph.align_info.a_end, [1, 3]);
}

#[test]
fn splicegraph_find_super_tr_ranks_seeded_supertranscripts() {
    let mut graph = splicegraph_l8_splicegraph_splicegraph(
        SuperTranscriptome {
            n: 1,
            sj_donor_nmax: 1,
            super_trs: vec![SuperTranscript {
                seq: vec![0, 1, 2],
                length: 3,
                ..Default::default()
            }],
            ..Default::default()
        },
        Some(ReadAlign::default()),
    );
    let map_gen = Genome {
        sa: vec![0, 1, 2],
        sai: vec![0, 1, 2, 3 | (1 << 30)],
        genome_sa_index_start: vec![0, 4],
        n_sa: 3,
        n_genome: 3,
        sj_gstart: u64::MAX,
        gstrand_bit: 28,
        gstrand_mask: (1 << 28) - 1,
        sai_mark_absent_mask_c: 1 << 30,
        sai_mark_nmask: !(1 << 29),
        chr_bin: vec![0, 0, 0],
        chr_start: vec![100],
        p_ge: ParametersGenome {
            g_saindex_nbases: 1,
            g_chr_bin_nbits: 0,
            ..Default::default()
        },
        ..Default::default()
    };

    let out = splicegraph_findsupertr_l5_splicegraph_findsupertr(
        &mut graph,
        &[0, 1, 2],
        &[1, 2, 3],
        3,
        "read1",
        &map_gen,
        10,
    );

    assert_eq!(graph.super_tr_seed_count, vec![3, 0]);
    assert!(
        out.starts_with("read1\t0\t0\t3\t3\t1\t1\t3\t1\t0\t2\t0\t2\n"),
        "{out:?}"
    );
    let ra = graph.ra.as_ref().unwrap();
    assert_eq!(ra.n_w, 1);
    assert_eq!(ra.n_win_tr, vec![1]);
    assert_eq!(ra.tr_all.len(), 1);
    assert_eq!(ra.tr_all[0][0].chr, 0);
    assert_eq!(ra.tr_all[0][0].str_, 0);
    assert_eq!(ra.tr_all[0][0].max_score, 3);
    assert_eq!(ra.tr_all[0][0].g_start, 100);
    assert_eq!(ra.tr_all[0][0].cigar, vec![[BAM_CIGAR_M, 3]]);
    assert_eq!(ra.tr_best.max_score, 3);
}

#[test]
fn read_align_map_one_read_splice_graph_initializes_and_dispatches_to_splice_graph() {
    let mut ra = ReadAlign {
        l_read: 3,
        read_nmates: 1,
        read_length: vec![3],
        read_length_original: vec![3],
        read_length_pair_original: 3,
        read_name: "read1".to_string(),
        i_read_all: 77,
        read1: [vec![0, 1, 2], Vec::new(), vec![1, 2, 3]],
        max_score_mate: vec![9],
        n_w: 9,
        n_tr: 8,
        tr_best: Transcript {
            max_score: -100,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut graph = splicegraph_l8_splicegraph_splicegraph(
        SuperTranscriptome {
            n: 1,
            sj_donor_nmax: 1,
            super_trs: vec![SuperTranscript {
                seq: vec![0, 1, 2],
                length: 3,
                ..Default::default()
            }],
            ..Default::default()
        },
        None,
    );
    let map_gen = Genome {
        sa: vec![0, 1, 2],
        sai: vec![0, 1, 2, 3 | (1 << 30)],
        genome_sa_index_start: vec![0, 4],
        n_sa: 3,
        n_genome: 3,
        sj_gstart: u64::MAX,
        gstrand_bit: 28,
        gstrand_mask: (1 << 28) - 1,
        sai_mark_absent_mask_c: 1 << 30,
        sai_mark_nmask: !(1 << 29),
        chr_bin: vec![0, 0, 0],
        chr_start: vec![100],
        p_ge: ParametersGenome {
            g_saindex_nbases: 1,
            g_chr_bin_nbits: 0,
            ..Default::default()
        },
        ..Default::default()
    };

    let out = readalign_maponereadsplicegraph_l6_readalign_maponereadsplicegraph(
        &mut ra, &mut graph, &map_gen, 2, 10,
    );

    assert!(out.starts_with("read1\t0\t0\t3\t3\t1\t1\t3\t1\t0\t2\t0\t2\n"));
    assert_eq!(ra.map_marker, 0);
    assert_eq!(ra.n_tr, 0);
    assert_eq!(ra.max_score_mate[0], 0);
    assert_eq!(ra.n_w, 1);
    assert_eq!(ra.n_win_tr, vec![1]);
    assert_eq!(ra.tr_all[0][0].max_score, 3);
    assert_eq!(ra.tr_best.max_score, 3);
    let graph_ra = graph.ra.as_ref().unwrap();
    assert_eq!(graph_ra.tr_best.i_read, 77);
    assert_eq!(graph_ra.tr_best.l_read, 3);
    assert_eq!(graph_ra.tr_best.read_name, "read1");
}

#[test]
fn read_align_map_one_read_splice_graph_marks_too_short_without_reset() {
    let mut ra = ReadAlign {
        l_read: 3,
        map_marker: 11,
        n_w: 7,
        n_tr: 5,
        tr_best: Transcript {
            r_length: 99,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut graph = splicegraph_l8_splicegraph_splicegraph(SuperTranscriptome::default(), None);
    let map_gen = Genome::default();

    let out = readalign_maponereadsplicegraph_l6_readalign_maponereadsplicegraph(
        &mut ra, &mut graph, &map_gen, 4, 10,
    );

    assert_eq!(out, "");
    assert_eq!(ra.map_marker, MARKER_READ_TOO_SHORT);
    assert_eq!(ra.tr_best.r_length, 0);
    assert_eq!(ra.n_w, 0);
    assert_eq!(ra.n_tr, 5);
    assert!(graph.ra.is_none());
}

#[test]
fn supertranscriptome_sj_collapse_sorts_deduplicates_and_outputs_tsv() {
    let mut st = SuperTranscriptome {
        seq: vec![vec![0], vec![1], vec![2]],
        sj: vec![
            sjInfo {
                super_: 1,
                start: 20,
                end: 30,
                tr: 7,
            },
            sjInfo {
                super_: 0,
                start: 10,
                end: 20,
                tr: 1,
            },
            sjInfo {
                super_: 1,
                start: 15,
                end: 25,
                tr: 2,
            },
            sjInfo {
                super_: 1,
                start: 15,
                end: 25,
                tr: 3,
            },
        ],
        ..Default::default()
    };

    let (tsv, log) = supertranscriptome_l4_supertranscriptome_sjcollapse(&mut st);

    assert_eq!(tsv, "0\t10\t20\n1\t15\t25\n1\t20\t30\n");
    assert_eq!(
        log,
        "Number of splice junctions in superTranscripts = 4\nNumber of collapsed splice junctions in superTranscripts = 3\n"
    );
    assert_eq!(
        st.sj
            .iter()
            .map(|sj| (sj.super_, sj.start, sj.end, sj.tr))
            .collect::<Vec<_>>(),
        vec![
            (0, 10, 20, 1),
            (1, 15, 25, 2),
            (1, 15, 25, 3),
            (1, 20, 30, 7)
        ]
    );
}

#[test]
fn solo_feature_stats_output_writes_summary_umi_and_cell_reads_tables() {
    let mut read_bar_sum = SoloReadBarcode {
        qual_hist: vec![0; 256],
        ..Default::default()
    };
    read_bar_sum.stats.v = vec![0; SOLO_READ_BARCODE_N_STATS];
    read_bar_sum.stats.v[0] = 7;
    read_bar_sum.qual_hist[20] = 10;
    read_bar_sum.qual_hist[40] = 10;

    let mut read_feat_sum = SoloReadFeature::default();
    read_feat_sum.stats.v = vec![0; SOLO_READ_FEATURE_N_STATS];
    read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_NO_TOO_MANY_WL_MATCHES] = 1;
    read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_NO_MM_TO_WL_WITHOUT_EXACT] = 1;
    read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_WL_MATCH] = 70;
    read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE] = 50;
    read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_UMIS] = 25;

    let mut cb0 = [0_u64; SOLO_READ_FLAG_N_BITS];
    cb0[SOLO_READ_FLAG_CB_MATCH as usize] = 3;
    cb0[SOLO_READ_FLAG_COUNTED_U as usize] = 2;
    let mut cb2 = [0_u64; SOLO_READ_FLAG_N_BITS];
    cb2[SOLO_READ_FLAG_CB_MM_UNIQUE as usize] = 4;
    cb2[SOLO_READ_FLAG_COUNTED_M as usize] = 1;
    let mut flag_counts = std::collections::BTreeMap::new();
    flag_counts.insert(0, cb0);
    flag_counts.insert(2, cb2);

    let solo_feature = SoloFeature {
        feature_type: SOLO_FEATURE_GENE,
        read_bar_sum: Some(read_bar_sum),
        read_feat_sum: Some(read_feat_sum),
        output_prefix: "solo/Gene/".to_string(),
        ind_cb_wl: vec![0, u32::MAX, 1],
        n_umi_per_cb: vec![5, 7],
        n_gene_per_cb: vec![2, 3],
        n_umi_per_cb_multi: vec![1, 2],
        n_gene_per_cb_multi: vec![1, 2],
        n_umi_per_cb_sorted: vec![9, 6, 0, 4],
        read_flag_counts: SoloReadFlagClass {
            flag_counts,
            flag_counts_no_cb: {
                let mut a = [0_u64; SOLO_READ_FLAG_N_BITS];
                a[SOLO_READ_FLAG_CB_MATCH as usize] = 1;
                a
            },
            ..Default::default()
        },
        filtered_cells: SoloFilteredCells {
            n_cells: 2,
            n_read_in_cells_unique: 40,
            mean_read_per_cell_unique: 20,
            median_read_per_cell_unique: 21,
            n_umi_in_cells: 30,
            mean_umi_per_cell: 15,
            median_umi_per_cell: 16,
            mean_gene_per_cell: 3,
            median_gene_per_cell: 4,
            n_gene_detected: 6,
            ..Default::default()
        },
        ..Default::default()
    };

    let mut p_solo = ParametersSolo {
        solo_type: SOLO_TYPE_CB_UMI_SIMPLE,
        cb_wl_str: vec!["CBA".to_string(), "CBB".to_string(), "CBC".to_string()],
        multi_map: MultiMappers {
            yes_multi: true,
            ..Default::default()
        },
        cell_filter: SoloCellFilter {
            type_: vec!["CellRanger2.2".to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    p_solo
        .read_stats_yes
        .resize(SOLO_FEATURE_GENE as usize + 1, false);
    p_solo.read_stats_yes[SOLO_FEATURE_GENE as usize] = true;

    let p = Parameters {
        run_thread_n: 1,
        read_nmates: 1,
        read_quality_score_base: 10,
        ..Default::default()
    };
    let ra_chunks = vec![ReadAlign {
        qual_hist: {
            let mut hist = vec![vec![0_u64; 256]];
            hist[0][20] = 20;
            hist[0][40] = 20;
            hist
        },
        ..Default::default()
    }];

    let out = solofeature_statsoutput_l6_solofeature_statsoutput(
        &solo_feature,
        &p,
        &p_solo,
        &Stats {
            read_n: 100,
            mapped_reads_u: 60,
            mapped_reads_m: 20,
            ..Default::default()
        },
        &ra_chunks,
    );

    assert_eq!(
        out.files["solo/Gene/Summary.csv"],
        concat!(
            "Number of Reads,100\n",
            "Reads With Valid Barcodes,0.91\n",
            "Sequencing Saturation,0.5\n",
            "Q30 Bases in CB+UMI,0.5\n",
            "Q30 Bases in RNA read,0.5\n",
            "Reads Mapped to Genome: Unique+Multiple,0.8\n",
            "Reads Mapped to Genome: Unique,0.6\n",
            "Reads Mapped to Gene: Unique+Multiple Gene,0.7\n",
            "Reads Mapped to Gene: Unique Gene,0.5\n",
            "Estimated Number of Cells,2\n",
            "Unique Reads in Cells Mapped to Gene,40\n",
            "Fraction of Unique Reads in Cells,0.8\n",
            "Mean Reads per Cell,20\n",
            "Median Reads per Cell,21\n",
            "UMIs in Cells,30\n",
            "Mean UMI per Cell,15\n",
            "Median UMI per Cell,16\n",
            "Mean Gene per Cell,3\n",
            "Median Gene per Cell,4\n",
            "Total Gene Detected,6\n"
        )
    );
    assert_eq!(out.files["solo/Gene/UMIperCellSorted.txt"], "9\n6\n");
    assert!(out.files["solo/Gene/CellReads.stats"].starts_with(
        "CB\tcbMatch\tcbPerfect\tcbMMunique\tcbMMmultiple\tgenomeU\tgenomeM\tfeatureU\tfeatureM\t"
    ));
    assert!(
        out.files["solo/Gene/CellReads.stats"]
            .contains("CBnotInPasslist\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n")
    );
    assert!(
        out.files["solo/Gene/CellReads.stats"]
            .contains("CBA\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t2\t0\t5\t2\t1\t1\n")
    );
    assert!(
        out.files["solo/Gene/CellReads.stats"]
            .contains("CBC\t0\t0\t4\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1\t7\t3\t2\t2\n")
    );
}

#[test]
fn supertranscriptome_load_slices_sequences_and_indexes_junctions() {
    let mut st = SuperTranscriptome::default();
    let g = b"AAACCCGGTT".to_vec();
    let log = supertranscriptome_l32_supertranscriptome_load(
        &mut st,
        &g,
        &[0, 3, 6],
        &[3, 3, 4],
        "0\t10\t30\n0\t10\t20\n0\t15\t25\n2\t8\t12\n",
    );

    assert_eq!(st.n, 3);
    assert_eq!(st.super_trs[0].length, 3);
    assert_eq!(st.super_trs[0].seq, b"AAA");
    assert_eq!(st.super_trs[1].seq, b"CCC");
    assert_eq!(st.super_trs[2].seq, b"GGTT");
    assert_eq!(st.super_trs[0].sj_donor, vec![10, 15]);
    assert_eq!(
        st.super_trs[0].sj_c,
        vec![[10, 20, 0], [15, 25, 1], [10, 30, 0]]
    );
    assert_eq!(st.super_trs[1].sj_c, Vec::<[u32; 3]>::new());
    assert_eq!(st.super_trs[2].sj_c, vec![[8, 12, 0]]);
    assert_eq!(st.sj_nmax, 3);
    assert_eq!(st.sj_donor_nmax, 2);
    assert_eq!(
        log,
        "Max number of splice junctions in a superTranscript = 3\nMax number of donor sites in a superTranscript = 2\n"
    );
}

#[test]
fn sam_attr_requires_bam_matches_error_and_warning_behavior() {
    let no_attr = Parameters::default();
    assert_eq!(
        parameters_samattributes_l251_parameters_samattrrequiresbam(&no_attr, false, "CB").unwrap(),
        None
    );

    let err = parameters_samattributes_l251_parameters_samattrrequiresbam(&no_attr, true, "CB")
        .unwrap_err();
    assert!(err.contains("CB tag, which requires BAM output"));

    let bam_only = Parameters {
        out_bam_unsorted: true,
        out_sam_bool: false,
        ..Default::default()
    };
    assert_eq!(
        parameters_samattributes_l251_parameters_samattrrequiresbam(&bam_only, true, "CB").unwrap(),
        None
    );

    let bam_and_sam = Parameters {
        out_bam_coord: true,
        out_sam_bool: true,
        ..Default::default()
    };
    assert_eq!(
        parameters_samattributes_l251_parameters_samattrrequiresbam(&bam_and_sam, true, "vW")
            .unwrap(),
        Some(
            "WARNING: --outSAMattributes contains vW tag. It will be output into BAM file(s), but not SAM file.\n"
                .to_string()
        )
    );
}

#[test]
fn output_transcript_sj_records_junctions_and_updates_duplicate_overhangs() {
    let mut out_sj = OutSJ::default();
    let tr = Transcript {
        exons: vec![[0, 100, 8, 0, 0], [10, 120, 5, 0, 0], [20, 140, 7, 0, 0]],
        canon_sj: vec![1, -1],
        sj_annot: vec![0, 0],
        n_exons: 3,
        ..Default::default()
    };
    readalign_outputtranscriptsj_l4_readalign_outputtranscriptsj(&tr, 1, &mut out_sj, 0);
    assert_eq!(out_sj.n, 1);
    assert_eq!(out_sj.junctions[0].start, 108);
    assert_eq!(out_sj.junctions[0].gap, 12);
    assert_eq!(out_sj.junctions[0].overhang_left, 5);
    assert_eq!(out_sj.junctions[0].overhang_right, 5);
    assert_eq!(out_sj.junctions[0].motif, 1);
    assert_eq!(out_sj.junctions[0].strand, 1);
    assert_eq!(out_sj.junctions[0].annot, 0);
    assert_eq!(out_sj.junctions[0].count_unique, 1);
    assert_eq!(out_sj.junctions[0].count_multiple, 0);

    let tr_dup = Transcript {
        exons: vec![[0, 88, 20, 0, 0], [25, 120, 30, 0, 0]],
        canon_sj: vec![2],
        sj_annot: vec![1],
        n_exons: 2,
        ..Default::default()
    };
    readalign_outputtranscriptsj_l4_readalign_outputtranscriptsj(&tr_dup, 3, &mut out_sj, 0);
    assert_eq!(out_sj.n, 1);
    assert_eq!(out_sj.junctions[0].overhang_left, 20);
    assert_eq!(out_sj.junctions[0].overhang_right, 20);

    readalign_outputtranscriptsj_l4_readalign_outputtranscriptsj(&tr_dup, 3, &mut out_sj, 1);
    assert_eq!(out_sj.n, 2);
    assert_eq!(out_sj.junctions[1].count_unique, 0);
    assert_eq!(out_sj.junctions[1].count_multiple, 1);
    assert_eq!(out_sj.junctions[1].strand, 2);
}

#[test]
fn outsj_collapse_sorts_and_merges_matching_junctions() {
    let mut out_sj = OutSJ {
        n: 4,
        n_store: 4,
        junctions: vec![
            star_rs::generated::structs::JunctionRecord {
                start: 30,
                gap: 5,
                motif: 2,
                annot: 1,
                count_unique: 1,
                count_multiple: 2,
                overhang_left: 3,
                overhang_right: 7,
                ..Default::default()
            },
            star_rs::generated::structs::JunctionRecord {
                start: 10,
                gap: 8,
                motif: 1,
                annot: 0,
                count_unique: 4,
                count_multiple: 0,
                overhang_left: 9,
                overhang_right: 1,
                ..Default::default()
            },
            star_rs::generated::structs::JunctionRecord {
                start: 30,
                gap: 5,
                motif: 2,
                annot: 1,
                count_unique: 5,
                count_multiple: 6,
                overhang_left: 8,
                overhang_right: 4,
                ..Default::default()
            },
            star_rs::generated::structs::JunctionRecord {
                start: 10,
                gap: 2,
                motif: 0,
                annot: 0,
                count_unique: 7,
                ..Default::default()
            },
        ],
    };

    outsj_l36_outsj_collapsesj(&mut out_sj).unwrap();
    assert_eq!(out_sj.n, 3);
    assert_eq!(
        out_sj
            .junctions
            .iter()
            .map(|sj| (sj.start, sj.gap))
            .collect::<Vec<_>>(),
        vec![(10, 2), (10, 8), (30, 5)]
    );
    assert_eq!(out_sj.junctions[2].count_unique, 6);
    assert_eq!(out_sj.junctions[2].count_multiple, 8);
    assert_eq!(out_sj.junctions[2].overhang_left, 8);
    assert_eq!(out_sj.junctions[2].overhang_right, 7);
}

#[test]
fn outsj_constructor_and_storage_growth_match_counter_state() {
    let mut out_sj = outsj_l4_outsj_outsj(3);
    assert_eq!(out_sj.n, 0);
    assert_eq!(out_sj.n_store, 3);
    assert!(out_sj.junctions.is_empty());
    assert!(out_sj.junctions.capacity() >= 3);

    out_sj
        .junctions
        .push(star_rs::generated::structs::JunctionRecord {
            start: 1,
            ..Default::default()
        });
    out_sj.n = 1;
    outsj_l62_outsj_datasizeincrease(&mut out_sj);
    assert_eq!(out_sj.n, 1);
    assert_eq!(out_sj.n_store, 6);
    assert_eq!(out_sj.junctions[0].start, 1);
    assert!(out_sj.junctions.capacity() >= 6);

    assert_eq!(
        outsj_l68_junction_junction(),
        star_rs::generated::structs::Junction::default()
    );
}

#[test]
fn junction_pointer_and_output_stream_match_original_layout_and_line() {
    let genome = Genome {
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 4,
            ..Default::default()
        },
        chr_bin: vec![0, 1, 1, 1],
        chr_name: vec!["chrA".to_string(), "chrB".to_string()],
        chr_start: vec![0, 16],
        ..Default::default()
    };
    let records = vec![
        star_rs::generated::structs::JunctionRecord {
            start: 3,
            gap: 5,
            strand: 1,
            motif: 2,
            annot: 0,
            count_unique: 7,
            count_multiple: 1,
            overhang_left: 11,
            overhang_right: 13,
        },
        star_rs::generated::structs::JunctionRecord {
            start: 40,
            gap: 9,
            strand: 2,
            motif: 5,
            annot: 1,
            count_unique: 3,
            count_multiple: 4,
            overhang_left: 8,
            overhang_right: 10,
        },
    ];
    let mut junction = star_rs::generated::structs::Junction {
        gen_out: genome,
        ..Default::default()
    };

    outsj_l72_junction_junctionpointer(&mut junction, &records, 1).unwrap();
    assert_eq!(junction.record, Some(records[1].clone()));
    assert_eq!(
        outsj_l85_junction_outputstream(&junction).unwrap(),
        "chrB\t25\t33\t2\t5\t1\t3\t4\t8\n"
    );

    assert!(outsj_l72_junction_junctionpointer(&mut junction, &records, 2).is_err());
    let unbound = star_rs::generated::structs::Junction::default();
    assert!(outsj_l85_junction_outputstream(&unbound).is_err());
}

#[test]
fn chain_constructor_loads_chain_blocks_like_original_parser() {
    let path = std::env::temp_dir().join(format!(
        "star_rs_chain_{}_{}.chain",
        std::process::id(),
        "blocks"
    ));
    std::fs::write(
        &path,
        "chain 100 chr1 1000 + 10 80 chrA 2000 + 20 90 1\n\
         25 3 4\n\
         30 5 6\n\
         12\n\
         \n\
         chain 50 chr2 500 + 7 40 chrB 700 + 9 42 2\n\
         11\n",
    )
    .unwrap();

    let chain = chain_l5_chain_chain(path.to_str().unwrap()).unwrap();
    let chr1 = chain.chr_chains.get("chr1").unwrap();
    assert_eq!(chr1.chr1, "chr1");
    assert_eq!(chr1.chr2, "chrA");
    assert_eq!(chr1.b_start1, vec![10, 38, 73]);
    assert_eq!(chr1.b_start2, vec![20, 49, 85]);
    assert_eq!(chr1.b_len, vec![25, 30, 12]);
    assert_eq!(chr1.b_n, 3);

    let chr2 = chain.chr_chains.get("chr2").unwrap();
    assert_eq!(chr2.b_start1, vec![7]);
    assert_eq!(chr2.b_start2, vec![9]);
    assert_eq!(chr2.b_len, vec![11]);
    assert_eq!(chr2.b_n, 1);

    std::fs::remove_file(path).unwrap();
}

#[test]
fn chain_liftover_gtf_writes_lifted_and_unlifted_records() {
    let base = std::env::temp_dir().join(format!("star_rs_chain_liftover_{}", std::process::id()));
    let chain_path = base.with_extension("chain");
    let gtf_path = base.with_extension("gtf");
    let out_path = base.with_extension("out.gtf");
    let unlifted_path = format!("{}.unlifted", out_path.to_str().unwrap());

    std::fs::write(
        &chain_path,
        "chain 100 chr1 1000 + 10 80 chrA 2000 + 20 90 1\n\
         25 3 4\n\
         30 5 6\n\
         12\n",
    )
    .unwrap();
    std::fs::write(
        &gtf_path,
        "#skip me\n\
         \n\
         chr1\tsrc\texon\t12\t20\t.\t+\t.\tgene_id \"a\";\n\
         chr1\tsrc\texon\t5\t11\t.\t+\t.\tgene_id \"left\";\n\
         chr1\tsrc\texon\t90\t95\t.\t+\t.\tgene_id \"bad\";\n",
    )
    .unwrap();

    let chain = chain_l5_chain_chain(chain_path.to_str().unwrap()).unwrap();
    chain_l58_chain_liftovergtf(
        &chain,
        gtf_path.to_str().unwrap(),
        out_path.to_str().unwrap(),
    )
    .unwrap();

    assert_eq!(
        std::fs::read_to_string(&out_path).unwrap(),
        "chrA\tsrc\texon\t22\t30\t.\t+\t.\tgene_id \"a\";\n\
         chrA\tsrc\texon\t20\t21\t.\t+\t.\tgene_id \"left\";\n"
    );
    assert_eq!(
        std::fs::read_to_string(&unlifted_path).unwrap(),
        "chr1\tsrc\texon\t90\t95\t.\t+\t.\tgene_id \"bad\";\n"
    );

    std::fs::remove_file(chain_path).unwrap();
    std::fs::remove_file(gtf_path).unwrap();
    std::fs::remove_file(out_path).unwrap();
    std::fs::remove_file(unlifted_path).unwrap();
}

#[test]
fn junction_collapse_one_sj_rejects_incompatible_duplicates() {
    let mut base = star_rs::generated::structs::JunctionRecord {
        start: 5,
        gap: 9,
        motif: 1,
        annot: 0,
        ..Default::default()
    };
    let motif_mismatch = star_rs::generated::structs::JunctionRecord {
        start: 5,
        gap: 9,
        motif: 2,
        annot: 0,
        ..Default::default()
    };
    assert!(outsj_l92_junction_collapseonesj(&mut base, &motif_mismatch).is_err());

    let annot_mismatch = star_rs::generated::structs::JunctionRecord {
        start: 5,
        gap: 9,
        motif: 1,
        annot: 1,
        ..Default::default()
    };
    assert!(outsj_l92_junction_collapseonesj(&mut base, &annot_mismatch).is_err());
}

#[test]
fn record_sj_applies_output_switch_and_multimapper_filter() {
    let tr1 = Transcript {
        exons: vec![[0, 100, 8, 0, 0], [10, 120, 5, 0, 0]],
        canon_sj: vec![1],
        sj_annot: vec![0],
        n_exons: 2,
        ..Default::default()
    };
    let tr2 = Transcript {
        exons: vec![[0, 200, 9, 0, 0], [10, 240, 6, 0, 0]],
        canon_sj: vec![2],
        sj_annot: vec![1],
        n_exons: 2,
        ..Default::default()
    };

    let mut out = OutSJ::default();
    readalign_outputalignments_l76_readalign_recordsj(
        2,
        &[tr1.clone(), tr2.clone()],
        &mut out,
        false,
        "All",
    );
    assert_eq!(out.n, 0);

    readalign_outputalignments_l76_readalign_recordsj(
        2,
        &[tr1.clone(), tr2.clone()],
        &mut out,
        true,
        "Unique",
    );
    assert_eq!(out.n, 0);

    readalign_outputalignments_l76_readalign_recordsj(
        2,
        &[tr1.clone(), tr2.clone()],
        &mut out,
        true,
        "All",
    );
    assert_eq!(out.n, 2);
    assert_eq!(out.junctions[0].count_unique, 0);
    assert_eq!(out.junctions[0].count_multiple, 1);
    assert_eq!(out.junctions[1].start, 209);

    let mut unique_out = OutSJ::default();
    readalign_outputalignments_l76_readalign_recordsj(1, &[tr1], &mut unique_out, true, "Unique");
    assert_eq!(unique_out.n, 1);
    assert_eq!(unique_out.junctions[0].count_unique, 1);
    assert_eq!(unique_out.junctions[0].count_multiple, 0);
}

#[test]
fn inoutstreams_constructor_and_destructor_state_matches_flush_close_order() {
    let streams = inoutstreams_l4_inoutstreams_inoutstreams();
    assert!(!streams.log_stdout_attached);
    assert!(!streams.out_sam_attached);

    let mut streams = InOutStreams {
        log_stdout_attached: true,
        out_sam_attached: true,
        out_unmapped_reads_open: [true, false],
        ..Default::default()
    };
    inoutstreams_l11_inoutstreams_inoutstreams(&mut streams);
    assert!(streams.log_stdout_flushed);
    assert!(streams.out_sam_flushed);
    assert!(streams.log_stdout_file_flushed);
    assert!(streams.out_sam_file_flushed);
    assert!(streams.out_chim_sam_flushed);
    assert!(streams.out_chim_junction_flushed);
    assert!(streams.log_progress_flushed);
    assert!(streams.log_main_flushed);
    assert!(streams.log_final_flushed);
    assert!(streams.out_local_chains_flushed);
    assert!(streams.out_sam_file_closed);
    assert!(streams.out_chim_sam_closed);
    assert!(streams.out_chim_junction_closed);
    assert!(streams.log_progress_closed);
    assert!(streams.log_final_closed);
    assert!(streams.out_local_chains_closed);
    assert!(!streams.out_unmapped_reads_open[0]);
    assert!(streams.out_unmapped_reads_flushed[0]);
    assert!(streams.out_unmapped_reads_closed[0]);
    assert!(!streams.out_unmapped_reads_flushed[1]);
    assert!(!streams.out_unmapped_reads_closed[1]);
}

#[test]
fn parameters_genome_initialize_validates_transform_and_genome_modes() {
    let mut pg = ParametersGenome {
        g_dir: "/tmp/genome".to_string(),
        transform: ParametersGenomeTransform {
            type_string: "Diploid".to_string(),
            output: vec!["SAM".to_string(), "SJ".to_string(), "Quant".to_string()],
            ..Default::default()
        },
        g_type_string: "Full".to_string(),
        g_load: "NoSharedMemory".to_string(),
        ..Default::default()
    };
    parametersgenome_l5_parametersgenome_initialize(&mut pg).unwrap();
    assert_eq!(pg.g_dir, "/tmp/genome/");
    assert_eq!(pg.transform.type_, 2);
    assert!(pg.transform.out_yes);
    assert!(pg.transform.out_sam);
    assert!(pg.transform.out_sj);
    assert!(pg.transform.out_quant);

    let mut none = ParametersGenome {
        g_dir: "/tmp/genome/".to_string(),
        transform: ParametersGenomeTransform {
            type_string: "None".to_string(),
            output: vec!["None".to_string()],
            ..Default::default()
        },
        g_type_string: "Transcriptome".to_string(),
        g_load: "LoadAndKeep".to_string(),
        ..Default::default()
    };
    parametersgenome_l5_parametersgenome_initialize(&mut none).unwrap();
    assert_eq!(none.transform.type_, 0);
    assert!(!none.transform.out_yes);

    let mut bad = none.clone();
    bad.transform.type_string = "Triploid".to_string();
    assert!(
        parametersgenome_l5_parametersgenome_initialize(&mut bad)
            .unwrap_err()
            .contains("--outTransformType")
    );
}

#[test]
fn parameters_chimeric_initialize_matches_output_flags_and_validations() {
    let mut p = Parameters {
        out_bam_unsorted: true,
        pe_overlap_nbases_min: 0,
        ..Default::default()
    };
    let mut pc = ParametersChimeric {
        segment_min: 12,
        multimap_nmax: 3,
        out_type: vec![
            "WithinBAM".to_string(),
            "Junctions".to_string(),
            "SoftClip".to_string(),
        ],
        filter_string_in: vec!["banGenomicN".to_string()],
        ..Default::default()
    };
    parameterschimeric_initialize_l6_parameterschimeric_initialize(&mut pc, &mut p, "@HD\n")
        .unwrap();
    assert!(pc.out_bam);
    assert!(pc.out_junctions);
    assert!(!pc.out_sam_old);
    assert!(!pc.out_bam_hard_clip);
    assert!(pc.filter_genomic_n);
    assert!(pc.out_chim_junction_opened);
    assert!(pc.out_chim_junction_contents.starts_with("chr_donorA\t"));
    assert_eq!(p.out_sam_attr_order, vec![ATTR_NM]);
    assert!(pc.log_main.contains("WARNING --chimOutType=WithinBAM"));

    let mut disabled = ParametersChimeric {
        segment_min: 0,
        out_type: vec!["bad".to_string()],
        ..Default::default()
    };
    parameterschimeric_initialize_l6_parameterschimeric_initialize(
        &mut disabled,
        &mut Parameters::default(),
        "",
    )
    .unwrap();
    assert!(!disabled.out_bam);
    assert!(disabled.out_bam_hard_clip);

    let mut no_bam_output = ParametersChimeric {
        segment_min: 1,
        out_type: vec!["WithinBAM".to_string()],
        ..Default::default()
    };
    assert!(
        parameterschimeric_initialize_l6_parameterschimeric_initialize(
            &mut no_bam_output,
            &mut Parameters::default(),
            "",
        )
        .unwrap_err()
        .contains("WithinBAM requires BAM output")
    );

    let mut old_sam = ParametersChimeric {
        segment_min: 1,
        multimap_nmax: 2,
        out_type: vec!["SeparateSAMold".to_string()],
        ..Default::default()
    };
    let mut p_old = Parameters::default();
    assert!(
        parameterschimeric_initialize_l6_parameterschimeric_initialize(
            &mut old_sam,
            &mut p_old,
            "@SQ\n",
        )
        .unwrap_err()
        .contains("chimMultimapNmax > 0")
    );
    assert!(old_sam.out_chim_sam_opened);
    assert_eq!(old_sam.out_chim_sam_contents, "@SQ\n");

    let mut pe_overlap = ParametersChimeric {
        segment_min: 1,
        out_type: vec!["Junctions".to_string()],
        ..Default::default()
    };
    let mut p_pe = Parameters {
        pe_overlap_nbases_min: 1,
        ..Default::default()
    };
    assert!(
        parameterschimeric_initialize_l6_parameterschimeric_initialize(
            &mut pe_overlap,
            &mut p_pe,
            "",
        )
        .unwrap_err()
        .contains("peOverlapNbasesMin")
    );

    let mut bad_filter = ParametersChimeric {
        segment_min: 1,
        filter_string_in: vec!["bad".to_string()],
        ..Default::default()
    };
    assert!(
        parameterschimeric_initialize_l6_parameterschimeric_initialize(
            &mut bad_filter,
            &mut Parameters::default(),
            "",
        )
        .unwrap_err()
        .contains("unrecognized value of --chimFilter")
    );
}

#[test]
fn chimeric_align_constructor_orders_segments_and_check_matches_original_guards() {
    let seg_late = ChimericSegment {
        align: Transcript {
            ro_start: 50,
            str_: 0,
            n_exons: 2,
            exons: vec![[0, 0, 5, 1, 0], [10, 10, 12, 1, 0]],
            ..Default::default()
        },
        ..Default::default()
    };
    let seg_early = ChimericSegment {
        align: Transcript {
            ro_start: 10,
            str_: 1,
            n_exons: 1,
            exons: vec![[0, 0, 10, 0, 0]],
            ..Default::default()
        },
        ..Default::default()
    };
    let chim = chimericalign_l3_chimericalign_chimericalign(seg_late, seg_early, 33, 8);
    assert_eq!(chim.chim_score, 33);
    assert!(!chim.stitching_done);
    assert_eq!(chim.al1.ro_start, 10);
    assert_eq!(chim.al2.ro_start, 50);
    assert_eq!(chim.ex1, 0);
    assert_eq!(chim.ex2, 0);
    assert!(chimericalign_l17_chimericalign_chimericcheck(&chim));

    let same_mate_short = chimericalign_l3_chimericalign_chimericalign(
        ChimericSegment {
            align: Transcript {
                ro_start: 1,
                str_: 1,
                n_exons: 1,
                exons: vec![[0, 0, 7, 0, 0]],
                ..Default::default()
            },
            ..Default::default()
        },
        ChimericSegment {
            align: Transcript {
                ro_start: 2,
                str_: 0,
                n_exons: 1,
                exons: vec![[0, 0, 10, 0, 0]],
                ..Default::default()
            },
            ..Default::default()
        },
        1,
        8,
    );
    assert!(!chimericalign_l17_chimericalign_chimericcheck(
        &same_mate_short
    ));
}

#[test]
fn chimeric_align_stitching_handles_mate_bracketed_junction() {
    let mut chim = ChimericAlign {
        seg1: ChimericSegment {
            str_: 1,
            ..Default::default()
        },
        seg2: ChimericSegment {
            str_: 0,
            ..Default::default()
        },
        chim_score: 0,
        al1: Transcript {
            exons: vec![[0, 2, 3, 0, 0]],
            n_exons: 1,
            l_read: 6,
            read_length: vec![3, 2],
            str_: 0,
            ..Default::default()
        },
        al2: Transcript {
            exons: vec![[4, 8, 2, 1, 0]],
            n_exons: 1,
            l_read: 6,
            read_length: vec![3, 2],
            str_: 0,
            ..Default::default()
        },
        ex1: 0,
        ex2: 0,
        ..Default::default()
    };
    let gen_seq = vec![4, 4, 0, 1, 2, 4, 4, 4, 0, 1, 4];
    let read = vec![0, 1, 2, 4, 0, 1];
    let read_rev = read.clone();

    chimericalign_chimericstitching_l3_chimericalign_chimericstitching(
        &mut chim,
        &gen_seq,
        [&read, &read_rev],
        &ParametersChimeric {
            junction_overhang_min: 1,
            score_junction_non_gtag: -5,
            ..Default::default()
        },
        0,
        -1,
        -2,
        -1,
        -2,
        -4,
        -8,
        -2,
        -4,
        0.0,
    );

    assert!(chim.stitching_done);
    assert_eq!(chim.chim_str, 1);
    assert_eq!(chim.chim_motif, -1);
    assert_eq!(chim.chim_j1, 5);
    assert_eq!(chim.chim_j2, 7);
    assert_eq!(chim.chim_repeat1, 0);
    assert_eq!(chim.chim_repeat2, 0);
    assert_eq!(chim.chim_score, 5);
}

#[test]
fn chimeric_detection_mult_stitches_and_outputs_unique_candidate() {
    let read_length = vec![10, 10];
    let left = Transcript {
        exons: vec![[0, 100, 6, 0, 0]],
        n_exons: 1,
        l_read: 21,
        read_length: read_length.clone(),
        read_length_original: read_length.clone(),
        read_length_pair_original: 21,
        read_nmates: 2,
        r_length: 6,
        max_score: 8,
        chr: 0,
        str_: 0,
        ro_start: 0,
        ..Default::default()
    };
    let right = Transcript {
        exons: vec![[11, 200, 6, 1, 0]],
        n_exons: 1,
        l_read: 21,
        read_length: read_length.clone(),
        read_length_original: read_length.clone(),
        read_length_pair_original: 21,
        read_nmates: 2,
        r_length: 6,
        max_score: 9,
        chr: 1,
        str_: 0,
        ro_start: 11,
        ..Default::default()
    };
    let mut chim_det = ChimericDetection {
        p: Parameters {
            p_ch: ParametersChimeric {
                segment_min: 3,
                junction_overhang_min: 2,
                score_min: 10,
                score_drop_max: 100,
                multimap_score_range: 0,
                multimap_nmax: 10,
                out_junctions: true,
                ..Default::default()
            },
            ..Default::default()
        },
        ra: Some(ReadAlign {
            read_name: "@chim1".to_string(),
            read_files_index: 0,
            ..Default::default()
        }),
        tr_all: vec![vec![left], vec![right]],
        n_win_tr: vec![1, 1],
        read1: [vec![0; 21], vec![0; 21]],
        out_gen: Genome {
            chr_name: vec!["chr1".to_string(), "chr2".to_string()],
            chr_start: vec![0, 0],
            g: vec![0; 256],
            ..Default::default()
        },
        ..Default::default()
    };

    let result =
        chimericdetection_chimericdetectionmult_l23_chimericdetection_chimericdetectionmult(
            &mut chim_det,
            2,
            &read_length,
            0,
            None,
            None,
            0,
            -1,
            -2,
            -1,
            -2,
            -4,
            -8,
            -2,
            -4,
            0.0,
        );

    assert!(result.chim_record);
    assert_eq!(result.chim_n, 1);
    assert_eq!(result.chim_score_best, 12);
    assert_eq!(chim_det.chim_aligns.len(), 1);
    assert_eq!(chim_det.chim_aligns[0].chim_motif, -1);
    assert!(result.chim_junction.contains("chr1\t107\t+"));
    assert!(result.chim_junction.contains("chr2\t200\t+"));
    assert!(result.bam_outputs.is_empty());
}

#[test]
fn chimeric_detection_mult_rejects_candidates_below_score_floor() {
    let read_length = vec![10, 10];
    let tr = |r_start, g_start, ifrag, score| Transcript {
        exons: vec![[r_start, g_start, 6, ifrag, 0]],
        n_exons: 1,
        l_read: 21,
        read_length: read_length.clone(),
        r_length: 6,
        max_score: score,
        str_: 0,
        ..Default::default()
    };
    let mut chim_det = ChimericDetection {
        p: Parameters {
            p_ch: ParametersChimeric {
                segment_min: 3,
                junction_overhang_min: 2,
                score_min: 30,
                score_drop_max: 100,
                multimap_nmax: 10,
                ..Default::default()
            },
            ..Default::default()
        },
        tr_all: vec![vec![tr(0, 100, 0, 8)], vec![tr(11, 200, 1, 9)]],
        n_win_tr: vec![1, 1],
        read1: [vec![0; 21], vec![0; 21]],
        out_gen: Genome {
            g: vec![0; 256],
            ..Default::default()
        },
        ..Default::default()
    };

    let result =
        chimericdetection_chimericdetectionmult_l23_chimericdetection_chimericdetectionmult(
            &mut chim_det,
            2,
            &read_length,
            0,
            None,
            None,
            0,
            -1,
            -2,
            -1,
            -2,
            -4,
            -8,
            -2,
            -4,
            0.0,
        );

    assert!(!result.chim_record);
    assert_eq!(result.chim_n, 0);
    assert!(chim_det.chim_aligns.is_empty());
}

#[test]
fn chimeric_junction_output_matches_original_tabular_record() {
    let map_gen = Genome {
        chr_name: vec!["chr1".to_string(), "chr2".to_string()],
        chr_start: vec![100, 1000],
        ..Default::default()
    };
    let p = Parameters {
        out_sam_attr_present: star_rs::generated::structs::SamAttrPresent {
            rg: true,
            ..Default::default()
        },
        out_sam_attr_rg: vec!["RG1".to_string()],
        p_solo: ParametersSolo {
            solo_type: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    let chim = ChimericAlign {
        chim_j1: 150,
        chim_j2: 1030,
        chim_repeat1: 2,
        chim_repeat2: 3,
        chim_motif: 1,
        chim_score: 37,
        al1: Transcript {
            chr: 0,
            str_: 0,
            exons: vec![[0, 120, 10, 0, 0]],
            n_exons: 1,
            read_nmates: 1,
            read_length_original: vec![20],
            read_length_pair_original: 20,
            ..Default::default()
        },
        al2: Transcript {
            chr: 1,
            str_: 1,
            exons: vec![[2, 1010, 8, 0, 0]],
            n_exons: 1,
            read_nmates: 1,
            read_length_original: vec![20],
            read_length_pair_original: 20,
            ..Default::default()
        },
        ..Default::default()
    };
    let solo_bar = SoloReadBarcode {
        cb_seq: "CELL".to_string(),
        umi_seq: "UMI".to_string(),
        ..Default::default()
    };

    let out = chimericalign_chimericjunctionoutput_l4_chimericalign_chimericjunctionoutput(
        &chim,
        &map_gen,
        &p,
        "@readA",
        0,
        4,
        25,
        true,
        40,
        45,
        Some(&solo_bar),
    );

    assert_eq!(
        out,
        "chr1\t51\t+\tchr2\t31\t-\t1\t2\t3\treadA\t21\t10M10S\t11\t2S8M10S\t4\t45\t25\t37\t40\t1\tRG1\tCELL\tUMI\n"
    );
}

#[test]
fn clipcr4_constructor_fill_and_polytail_match_original_defaults() {
    let mut clip = clipcr4_l3_clipcr4_clipcr4();
    assert_eq!(clip.db_n, 64);
    assert_eq!(clip.read_len, 91);
    assert_eq!(clip.alphabet_length, 5);
    assert_eq!(clip.gap_open, 2);
    assert_eq!(clip.gap_ext, 2);
    assert_eq!(clip.score_matrix.len(), 25);
    assert_eq!(&clip.score_matrix[..5], &[1, -2, -2, -2, -2]);
    assert_eq!(clip.db_seq_arr.len(), 64 * 91);
    assert_eq!(clip.db_seqs_len, vec![91; 64]);
    assert_eq!(clip.opal_res.len(), 64);

    clipcr4_l43_clipcr4_opalfilloneseq(&mut clip, 2, b"ACGTX", 5);
    let start = 2 * 91;
    assert_eq!(
        &clip.db_seq_arr[start..start + 8],
        &[0, 1, 2, 3, 4, 4, 4, 4]
    );
    clip.opal_res[0].score = 77;
    clip.opal_res[0].score_set = 1;
    clip.opal_res[0].end_location_query = 12;
    clip.opal_res[0].end_location_target = 13;
    clipcr4_l43_clipcr4_opalfilloneseq(&mut clip, 0, b"ACG", 3);
    clip.db_seqs_len[0] = 3;
    assert_eq!(
        clipcr4_l72_clipcr4_opalalign(&mut clip, &[0, 1, 2], 3, 1),
        0
    );
    assert_eq!(clip.opal_res[0].score, 3);
    assert_eq!(clip.opal_res[0].end_location_query, 2);
    assert_eq!(clip.opal_res[0].end_location_target, 2);
    assert_eq!(clip.opal_res[0].alignment, None);
    assert_eq!(clip.opal_res[0].alignment_length, -1);

    let mut tail = vec![2u8; 8];
    tail.extend(std::iter::repeat_n(0u8, 22));
    assert_eq!(clipcr4_l82_clipcr4_polytail3p(&tail, tail.len() as u32), 22);
    assert_eq!(clipcr4_l82_clipcr4_polytail3p(&tail[..19], 19), 0);
}

#[test]
fn clipmate_clipchunk_records_cr4_clip_lengths_in_fastq_plus_lines() {
    let mut clip_mate = ClipMate {
        type_: 10,
        ad_seq_num: vec![0; 25],
        cr4: Some(clipcr4_l3_clipcr4_clipcr4()),
        ..Default::default()
    };
    let mut chunk = Vec::new();
    chunk.extend_from_slice(b"@r1\n");
    chunk.extend_from_slice(b"AAAAAAAAAAAAAAAAAAAAAAAAA\n+\n");
    chunk.extend_from_slice(b"IIIIIIIIIIIIIIIIIIIIIIIII\n");
    chunk.extend_from_slice(b"@r2\n");
    chunk.extend_from_slice(b"TTTTTTTTTTTTTTTTTTTTTTTTT\n+\n");
    chunk.extend_from_slice(b"IIIIIIIIIIIIIIIIIIIIIIIII\n");
    let plus1 = b"@r1\nAAAAAAAAAAAAAAAAAAAAAAAAA\n".len();
    let plus2 = plus1 + b"+\nIIIIIIIIIIIIIIIIIIIIIIIII\n@r2\nTTTTTTTTTTTTTTTTTTTTTTTTT\n".len();
    let chunk_len = chunk.len() as u64;

    clipmate_clipchunk_l7_clipmate_clipchunk(&mut clip_mate, &mut chunk, chunk_len).unwrap();

    assert_eq!(chunk[plus1], 91);
    assert_eq!(chunk[plus2], 0);
    let cr4 = clip_mate.cr4.as_ref().unwrap();
    assert_eq!(cr4.store_clip[0] as usize, plus1);
    assert_eq!(cr4.store_clip[1] as usize, plus2);
}

#[test]
fn clipmate_clipchunk_noops_for_non_cr4_type() {
    let mut clip_mate = ClipMate {
        type_: 3,
        ..Default::default()
    };
    let mut chunk = b"@r\nACGT\n+\nIIII\n".to_vec();
    let original = chunk.clone();
    let chunk_len = chunk.len() as u64;

    clipmate_clipchunk_l7_clipmate_clipchunk(&mut clip_mate, &mut chunk, chunk_len).unwrap();

    assert_eq!(chunk, original);
}

#[test]
fn star_usage_matches_brief_and_full_output_shape() {
    let brief = star_l36_usage(0, "2.7.11b", "build host dir", b"--genomeDir test\n");
    assert!(brief.starts_with(
        "Usage: STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq\n"
    ));
    assert!(brief.contains(
        "Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2022\n\n"
    ));
    assert!(brief.contains("STAR version=2.7.11b\n"));
    assert!(brief.contains("STAR compilation time,server,dir=build host dir\n"));
    assert!(brief.ends_with("\nTo list all parameters, run STAR --help\n"));
    assert!(!brief.contains("--genomeDir test"));

    let full = star_l36_usage(1, "V", "C", b"--genomeDir test\n--runThreadN 4\n");
    assert!(full.contains("STAR version=V\n"));
    assert!(full.ends_with("--genomeDir test\n--runThreadN 4\n"));

    let other = star_l36_usage(2, "V", "C", b"ignored");
    assert!(other.contains("<https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>\n"));
    assert!(!other.contains("ignored"));
    assert!(!other.contains("To list all parameters"));
}

#[test]
fn star_main_returns_usage_for_empty_and_help_arguments() {
    let p = Parameters::default();
    let empty = star_l58_main(
        &["STAR".to_string()],
        p.clone(),
        b"--runThreadN 4\n",
        None,
        None,
        None,
        None,
        &std::collections::BTreeSet::new(),
        &[],
        &[],
    )
    .unwrap();
    assert_eq!(empty.exit_code, 0);
    assert!(empty.usage.contains("To list all parameters"));

    let help = star_l58_main(
        &["STAR".to_string(), "--help".to_string()],
        p,
        b"--runThreadN 4\n",
        None,
        None,
        None,
        None,
        &std::collections::BTreeSet::new(),
        &[],
        &[],
    )
    .unwrap();
    assert_eq!(help.exit_code, 0);
    assert!(help.usage.ends_with("--runThreadN 4\n"));
}

#[test]
fn star_main_reports_unknown_run_mode() {
    let result = star_l58_main(
        &[
            "STAR".to_string(),
            "--runMode".to_string(),
            "bad".to_string(),
        ],
        Parameters {
            run_mode_in: vec!["bad".to_string()],
            ..Default::default()
        },
        b"",
        None,
        None,
        None,
        None,
        &std::collections::BTreeSet::new(),
        &[],
        &[],
    )
    .unwrap();

    assert_eq!(result.exit_code, 1);
    assert!(
        result
            .log_main
            .contains("unknown value of input parameter runMode=bad")
    );
}

#[test]
fn star_main_genome_generate_writes_index_and_removes_tmp() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_main_genome_{}_{}",
        std::process::id(),
        unique
    ));
    let tmp = dir.join("tmp");
    std::fs::create_dir_all(&tmp).unwrap();
    let fasta = dir.join("tiny.fa");
    std::fs::write(&fasta, ">chr1\nACGT\n").unwrap();
    let p = Parameters {
        command_line: "STAR --runMode genomeGenerate".to_string(),
        command_line_full: "STAR --runMode genomeGenerate".to_string(),
        version_genome: "test-version".to_string(),
        run_mode_in: vec!["genomeGenerate".to_string()],
        out_file_tmp: tmp.to_str().unwrap().to_string(),
        limit_genome_generate_ram: 10_000,
        limit_sjdb_insert_nsj: 0,
        p_ge: ParametersGenome {
            g_dir: dir.to_str().unwrap().to_string(),
            g_type_string: "Full".to_string(),
            g_fasta_files: vec![fasta.to_str().unwrap().to_string()],
            g_chr_bin_nbits: 2,
            g_saindex_nbases: 0,
            g_sasparse_d: 1,
            g_suffix_length_max: 8,
            sjdb_file_chr_start_end: vec!["-".to_string()],
            sjdb_gtf_file: "-".to_string(),
            transform: ParametersGenomeTransform {
                type_string: "None".to_string(),
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };

    let result = star_l58_main(
        &[
            "STAR".to_string(),
            "--runMode".to_string(),
            "genomeGenerate".to_string(),
        ],
        p,
        b"",
        None,
        None,
        None,
        None,
        &std::collections::BTreeSet::new(),
        &[],
        &[],
    )
    .unwrap();

    assert_eq!(result.exit_code, 0);
    assert_eq!(result.genome_generate.len(), 1);
    assert!(result.log_main.contains("DONE: Genome generation, EXITING"));
    assert!(result.removed_tmp);
    assert!(!tmp.exists());
    assert_eq!(
        std::fs::read_to_string(dir.join("chrName.txt")).unwrap(),
        "chr1\n"
    );
    assert!(std::fs::metadata(dir.join("SA")).unwrap().len() > 0);

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn star_main_align_reads_with_injected_genome_finishes_and_closes_reads() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_main_align_{}_{}",
        std::process::id(),
        unique
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let tmp = dir.join("tmp");
    std::fs::create_dir_all(&tmp).unwrap();

    let p = Parameters {
        command_line: "STAR --genomeDir genome --readFilesIn r1.fq".to_string(),
        command_line_full: "STAR --genomeDir genome --readFilesIn r1.fq".to_string(),
        run_mode_in: vec!["alignReads".to_string()],
        run_thread_n: 1,
        read_nends: 1,
        out_sam_type: vec!["None".to_string()],
        out_sam_mode: "None".to_string(),
        out_file_tmp: tmp.to_str().unwrap().to_string(),
        out_tmp_keep: "None".to_string(),
        read_files_in: vec!["r1.fq".to_string()],
        read_in_open: vec![true],
        read_files_command_pid: vec![77],
        ..Default::default()
    };
    let genome = Genome {
        n_chr_real: 1,
        n_genome: 4,
        chr_name: vec!["chr1".to_string()],
        chr_length: vec![4],
        chr_start: vec![0, 4],
        ..Default::default()
    };
    let chunks = vec![ReadAlignChunk::default()];

    let result = star_l58_main(
        &[
            "STAR".to_string(),
            "--genomeDir".to_string(),
            "genome".to_string(),
        ],
        p,
        b"",
        Some(genome),
        Some(Transcriptome::default()),
        None,
        Some(chunks),
        &std::collections::BTreeSet::new(),
        &[],
        &[],
    )
    .unwrap();

    assert_eq!(result.exit_code, 0);
    assert!(result.log_stdout.contains("started STAR run"));
    assert!(result.log_stdout.contains("finished successfully"));
    assert!(result.log_progress.contains("ALL DONE!"));
    assert!(result.log_main.contains("ALL DONE!"));
    assert_eq!(result.killed_read_command_pids, vec![77]);
    assert!(result.removed_tmp);
    assert!(!tmp.exists());

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn clipmate_initialize_matches_adapter_sequence_modes() {
    let mut no_clip = ClipMate::default();
    clipmate_initialize_l5_clipmate_initialize(&mut no_clip, 0, "-", 3, 0.2);
    assert_eq!(no_clip.type_, -1);
    assert_eq!(no_clip.ad_seq, "");
    assert!(no_clip.ad_seq_num.is_empty());
    assert_eq!(no_clip.n_after_ad, 3);
    assert_eq!(no_clip.ad_mmp, 0.2);

    let mut adapter = ClipMate {
        type_: 1,
        ..Default::default()
    };
    clipmate_initialize_l5_clipmate_initialize(&mut adapter, 5, "ACGTN", 2, 0.1);
    assert_eq!(adapter.type_, 1);
    assert_eq!(adapter.n, 5);
    assert_eq!(adapter.ad_seq, "ACGTN");
    assert_eq!(adapter.ad_seq_num, vec![0, 1, 2, 3, 4]);

    let mut poly_a = ClipMate {
        type_: 0,
        ..Default::default()
    };
    clipmate_initialize_l5_clipmate_initialize(&mut poly_a, 1, "polyA", 0, 0.0);
    assert_eq!(poly_a.ad_seq_num.len(), DEF_READ_SEQ_LENGTH_MAX);
    assert!(poly_a.ad_seq_num.iter().all(|base| *base == 0));

    let mut cr4 = ClipMate {
        type_: 10,
        ..Default::default()
    };
    clipmate_initialize_l5_clipmate_initialize(&mut cr4, 1, "A", 0, 0.0);
    assert!(cr4.cr4.is_some());
}

#[test]
fn clipmate_clip_matches_fixed_base_adapter_and_cr4_modes() {
    let mut five_prime = ClipMate {
        type_: 0,
        n: 2,
        n_after_ad: 1,
        ..Default::default()
    };
    let mut seq = vec![0, 1, 2, 3, 0, 1];
    let mut l_read = 6;
    assert_eq!(
        clipmate_clip_l5_clipmate_clip(&mut five_prime, &mut l_read, &mut seq),
        3
    );
    assert_eq!(l_read, 3);
    assert_eq!(&seq[..3], &[3, 0, 1]);

    let mut too_short = ClipMate {
        type_: 1,
        n: 10,
        ..Default::default()
    };
    let mut short_seq = vec![0, 1, 2];
    let mut short_len = 3;
    assert_eq!(
        clipmate_clip_l5_clipmate_clip(&mut too_short, &mut short_len, &mut short_seq),
        3
    );
    assert_eq!(short_len, 0);

    let mut three_prime = ClipMate {
        type_: 1,
        ad_seq: "AC".to_string(),
        ad_seq_num: vec![0, 1],
        ad_mmp: 0.0,
        ..Default::default()
    };
    let mut seq3 = vec![2, 3, 0, 1, 2];
    let mut len3 = 5;
    assert_eq!(
        clipmate_clip_l5_clipmate_clip(&mut three_prime, &mut len3, &mut seq3),
        3
    );
    assert_eq!(len3, 2);
    assert_eq!(three_prime.clipped_ad_n, 3);

    let mut cr4_five = ClipMate {
        type_: 10,
        ad_seq: "A".to_string(),
        clipped_info: 2,
        ..Default::default()
    };
    let mut seq4 = vec![9, 8, 7, 6];
    let mut len4 = 4;
    assert_eq!(
        clipmate_clip_l5_clipmate_clip(&mut cr4_five, &mut len4, &mut seq4),
        2
    );
    assert_eq!(len4, 2);
    assert_eq!(&seq4[..2], &[7, 6]);

    let mut cr4_three = ClipMate {
        type_: 11,
        ad_seq: "A".to_string(),
        ..Default::default()
    };
    let mut seq5 = vec![2u8; 8];
    seq5.extend(std::iter::repeat_n(0u8, 22));
    let mut len5 = seq5.len() as u32;
    assert_eq!(
        clipmate_clip_l5_clipmate_clip(&mut cr4_three, &mut len5, &mut seq5),
        22
    );
    assert_eq!(len5, 8);
}

#[test]
fn parameters_clip_initialize_expands_defaults_and_builds_clip_mates() {
    let p = Parameters {
        read_nmates: 1,
        read_nends: 3,
        ..Default::default()
    };
    let mut clip = ParametersClip {
        adapter_type: vec!["CellRanger4".to_string()],
        in_: [
            ReadClipInput {
                n: vec![0],
                n_after_ad: vec![0],
                ad_seq: vec!["-".to_string()],
                ad_mmp: vec![0.0],
            },
            ReadClipInput {
                n: vec![0],
                n_after_ad: vec![0],
                ad_seq: vec!["-".to_string()],
                ad_mmp: vec![0.0],
            },
        ],
        ..Default::default()
    };
    parametersclip_initialize_l6_parametersclip_initialize(&mut clip, &p).unwrap();
    assert_eq!(clip.in_[0].ad_seq[0], "AAGCAGTGGTATCAACGCAGAGTACATGGG");
    assert_eq!(clip.in_[1].ad_seq, vec!["A"]);
    assert_eq!(clip.in_[0].n, vec![0]);
    assert_eq!(clip.in_[1].n_after_ad, vec![0]);

    let mut clip_mates = Vec::new();
    parametersclip_initialize_l84_parametersclip_initializeclipmates(&clip, &mut clip_mates);
    assert_eq!(clip_mates.len(), 3);
    assert_eq!(clip_mates[0][0].type_, 10);
    assert_eq!(clip_mates[0][1].type_, 11);
    assert!(clip_mates[0][0].cr4.is_some());
    assert_eq!(clip_mates[1][0].ad_seq, "");
    assert_eq!(clip_mates[2][0].type_, -1);
    assert_eq!(clip_mates[2][1].type_, -1);
}

#[test]
fn parameters_clip_initialize_rejects_invalid_or_unsupported_modes() {
    let p = Parameters {
        read_nmates: 2,
        read_nends: 2,
        ..Default::default()
    };
    let base_in = [
        ReadClipInput {
            n: vec![0],
            n_after_ad: vec![0],
            ad_seq: vec!["-".to_string()],
            ad_mmp: vec![0.0],
        },
        ReadClipInput {
            n: vec![0],
            n_after_ad: vec![0],
            ad_seq: vec!["-".to_string()],
            ad_mmp: vec![0.0],
        },
    ];

    let mut bad_type = ParametersClip {
        adapter_type: vec!["Bad".to_string()],
        in_: base_in.clone(),
        ..Default::default()
    };
    assert!(
        parametersclip_initialize_l6_parametersclip_initialize(&mut bad_type, &p)
            .unwrap_err()
            .contains("--clipAdapterType")
    );

    let mut hamming_unsupported_5p = ParametersClip {
        adapter_type: vec!["Hamming".to_string()],
        in_: [
            ReadClipInput {
                ad_seq: vec!["AC".to_string()],
                ..base_in[0].clone()
            },
            base_in[1].clone(),
        ],
        ..Default::default()
    };
    assert!(
        parametersclip_initialize_l6_parametersclip_initialize(&mut hamming_unsupported_5p, &p)
            .unwrap_err()
            .contains("--clip5pAdapterSeq")
    );

    let mut wrong_count = ParametersClip {
        adapter_type: vec!["Hamming".to_string()],
        in_: [
            base_in[0].clone(),
            ReadClipInput {
                ad_seq: vec!["A".to_string(), "C".to_string()],
                ad_mmp: vec![0.0],
                n: vec![0, 0],
                n_after_ad: vec![0, 0],
            },
        ],
        ..Default::default()
    };
    assert!(
        parametersclip_initialize_l6_parametersclip_initialize(&mut wrong_count, &p)
            .unwrap_err()
            .contains("AdapterMMp")
    );
}

#[test]
fn packed_array_primitives_match_bit_packing_contract() {
    let mut packed = packedarray_l3_packedarray_packedarray();
    assert!(!packed.array_allocated);
    assert!(packed.char_array.is_empty());

    packedarray_l8_packedarray_definebits(&mut packed, 3, 5);
    assert_eq!(packed.word_length, 3);
    assert_eq!(packed.word_comp_length, 61);
    assert_eq!(packed.bit_rec_mask, 7);
    assert_eq!(packed.length_byte, 9);

    packedarray_l31_packedarray_allocatearray(&mut packed);
    assert!(packed.array_allocated);
    assert_eq!(packed.char_array, vec![0; 9]);

    packedarray_l17_packedarray_writepacked(&mut packed, 0, 5);
    packedarray_l17_packedarray_writepacked(&mut packed, 1, 3);
    packedarray_l17_packedarray_writepacked(&mut packed, 2, 7);
    let word = u64::from_ne_bytes(packed.char_array[0..8].try_into().unwrap());
    assert_eq!(word & 0b111, 5);
    assert_eq!((word >> 3) & 0b111, 3);
    assert_eq!((word >> 6) & 0b111, 7);
    assert_eq!(packedarray_h18_packedarray_index(&packed, 0), 5);
    assert_eq!(packedarray_h18_packedarray_index(&packed, 1), 3);
    assert_eq!(packedarray_h18_packedarray_index(&packed, 2), 7);

    packedarray_l27_packedarray_pointarray(&mut packed, &[1, 2, 3, 4, 5, 6, 7, 8]);
    assert_eq!(packed.char_array, vec![1, 2, 3, 4, 5, 6, 7, 8]);
    packedarray_l37_packedarray_deallocatearray(&mut packed);
    assert!(!packed.array_allocated);
    assert!(packed.char_array.is_empty());

    let mut pointed = PackedArray::default();
    packedarray_l27_packedarray_pointarray(&mut pointed, &[9, 8]);
    packedarray_l37_packedarray_deallocatearray(&mut pointed);
    assert!(pointed.char_array.is_empty());
}

#[test]
fn solo_collapse_comparators_match_original_sort_keys() {
    assert_eq!(
        solofeature_collapseumiall_l540_funcomparesolo1(&[5, 2], &[1, 3]),
        -1
    );
    assert_eq!(
        solofeature_collapseumiall_l540_funcomparesolo1(&[5, 3], &[1, 3]),
        1
    );
    assert_eq!(
        solofeature_collapseumiall_l540_funcomparesolo1(&[5, 3], &[5, 3]),
        0
    );

    assert_eq!(
        solofeature_collapseumiall_l557_funcompare_uint32_1_2_0(&[9, 1, 2], &[8, 2, 0]),
        -1
    );
    assert_eq!(
        solofeature_collapseumiall_l557_funcompare_uint32_1_2_0(&[9, 2, 4], &[8, 2, 3]),
        1
    );
    assert_eq!(
        solofeature_collapseumiall_l557_funcompare_uint32_1_2_0(&[9, 2, 4], &[10, 2, 4]),
        -1
    );
}

#[test]
fn solo_read_flag_reset_clears_each_feature_flag() {
    let mut solo_read = SoloRead {
        read_feat: vec![
            SoloReadFeature {
                read_flag: star_rs::generated::structs::SoloReadFlagClass {
                    flag: 7,
                    ..Default::default()
                },
                ..Default::default()
            },
            SoloReadFeature {
                read_flag: star_rs::generated::structs::SoloReadFlagClass {
                    flag: 11,
                    ..Default::default()
                },
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    soloread_l18_soloread_readflagreset(&mut solo_read);
    assert_eq!(solo_read.read_feat[0].read_flag.flag, 0);
    assert_eq!(solo_read.read_feat[1].read_flag.flag, 0);
}

#[test]
fn solo_read_feature_and_solo_feature_constructors_follow_feature_type_switches() {
    let p = Parameters {
        out_file_tmp: "/tmp/star-solo".to_string(),
        run_thread_n: 3,
        p_solo: ParametersSolo {
            solo_type: 1,
            read_info_yes: vec![false, true, false, false, false, true, false, false],
            read_index_yes: vec![true, false, false, false, false, true, false, false],
            cb_wl_yes: true,
            cb_wl_size: 4,
            ..Default::default()
        },
        sj_all: [vec![10, 20, 30], vec![40]],
        ..Default::default()
    };

    let rf = soloreadfeature_l5_soloreadfeature_soloreadfeature(SOLO_FEATURE_TRANSCRIPT3P, &p, 7);
    assert_eq!(rf.feature_type, SOLO_FEATURE_TRANSCRIPT3P);
    assert!(rf.read_info_yes);
    assert!(!rf.read_index_yes);
    assert_eq!(rf.cb_read_count, vec![0; 4]);
    assert_eq!(rf.transcript_dist_count.len(), 10000);
    assert_eq!(
        rf.stream_reads_path,
        Some("/tmp/star-solo/soloTranscript3p_7".to_string())
    );

    let sf_gene = solofeature_l4_solofeature_solofeature(&p, SOLO_FEATURE_GENE, 11);
    assert_eq!(sf_gene.features_number, 11);
    assert_eq!(sf_gene.read_feat_all_len, 3);
    assert!(sf_gene.read_feat_sum.is_some());

    let sf_sj = solofeature_l4_solofeature_solofeature(&p, SOLO_FEATURE_SJ, 11);
    assert_eq!(sf_sj.features_number, 3);

    let p_off = Parameters {
        p_solo: ParametersSolo {
            solo_type: 0,
            cb_wl_yes: true,
            cb_wl_size: 5,
            ..Default::default()
        },
        ..Default::default()
    };
    let rf_off = soloreadfeature_l5_soloreadfeature_soloreadfeature(SOLO_FEATURE_GENE, &p_off, -1);
    assert_eq!(rf_off.cb_read_count.len(), 0);
}

#[test]
fn parameters_solo_complex_wl_strings_cycles_lengths_and_barcode_carry() {
    let mut p_solo = ParametersSolo {
        cb_wl_size: 6,
        cb_v: vec![
            SoloBarcode {
                min_len: 1,
                wl: vec![vec![], vec![0, 1], vec![6]],
                ..Default::default()
            },
            SoloBarcode {
                min_len: 1,
                wl: vec![vec![], vec![2, 3]],
                ..Default::default()
            },
        ],
        ..Default::default()
    };

    parameterssolo_l503_parameterssolo_complexwlstrings(&mut p_solo);

    assert_eq!(
        p_solo.cb_wl_str,
        vec!["A_G", "C_G", "CG_G", "A_T", "C_T", "CG_T"]
    );
    assert_eq!(p_solo.cb_v[0].i_len, 2);
    assert_eq!(p_solo.cb_v[0].i_cb, 1);
    assert_eq!(p_solo.cb_v[1].i_len, 1);
    assert_eq!(p_solo.cb_v[1].i_cb, 1);
}

#[test]
fn solo_feature_load_raw_matrix_sorts_counts_and_loads_side_files() {
    let p = Parameters {
        run_mode_in: vec![
            "soloCellFiltering".to_string(),
            "/raw/matrix".to_string(),
            "/filtered/out".to_string(),
        ],
        ..Default::default()
    };
    let mut p_solo = ParametersSolo {
        out_file_names: vec![
            "unused0".to_string(),
            "features.tsv".to_string(),
            "barcodes.tsv".to_string(),
            "matrix.mtx".to_string(),
        ],
        ..Default::default()
    };
    let mut sf = SoloFeature::default();
    let matrix = "%%MatrixMarket matrix coordinate integer general\n%\n4 3 4\n3 2 1.2\n1 1 2.6\n2 2 4.5\n1 3 1.0\n";

    solofeature_loadrawmatrix_l7_solofeature_loadrawmatrix(
        &mut sf,
        &p,
        &mut p_solo,
        matrix,
        "CB1\nCB2\nCB3\n",
        "geneA\tGene A\n",
    )
    .unwrap();

    assert_eq!(sf.features_number, 4);
    assert_eq!(sf.output_prefix, "/filtered/out");
    assert_eq!(sf.output_prefix_filtered, "/filtered/out");
    assert_eq!(sf.count_mat_stride, 3);
    assert_eq!(sf.n_cb, 2);
    assert_eq!(sf.ind_cb, vec![0, 1, 2]);
    assert_eq!(sf.count_cell_gene_umi_index, vec![0, 3, 9]);
    assert_eq!(sf.n_gene_per_cb, vec![1, 2, 1]);
    assert_eq!(sf.n_umi_per_cb, vec![3, 6, 1]);
    assert_eq!(
        sf.count_cell_gene_umi,
        vec![0, 3, 3, 1, 5, 5, 2, 1, 1, 0, 1, 1]
    );
    assert_eq!(p_solo.cb_wl_str, vec!["CB1", "CB2", "CB3"]);
    assert_eq!(sf.copied_features_tsv, "geneA\tGene A\n");
}

#[test]
fn solo_feature_load_raw_matrix_reports_missing_parameters_and_empty_matrix() {
    let p = Parameters {
        run_mode_in: vec!["soloCellFiltering".to_string()],
        ..Default::default()
    };
    let mut p_solo = ParametersSolo::default();
    let mut sf = SoloFeature::default();
    let err = solofeature_loadrawmatrix_l7_solofeature_loadrawmatrix(
        &mut sf,
        &p,
        &mut p_solo,
        "",
        "",
        "",
    )
    .unwrap_err();
    assert!(err.contains("--runMode soloCellFiltering"));

    let p = Parameters {
        run_mode_in: vec![
            "soloCellFiltering".to_string(),
            "in".to_string(),
            "out".to_string(),
        ],
        ..Default::default()
    };
    let mut p_solo = ParametersSolo {
        out_file_names: vec![
            "unused0".to_string(),
            "features.tsv".to_string(),
            "barcodes.tsv".to_string(),
            "matrix.mtx".to_string(),
        ],
        ..Default::default()
    };
    let err = solofeature_loadrawmatrix_l7_solofeature_loadrawmatrix(
        &mut sf,
        &p,
        &mut p_solo,
        "4 3 0\n",
        "",
        "",
    )
    .unwrap_err();
    assert!(err.contains("no counts detected"));
}

#[test]
fn solo_feature_output_results_writes_gene_barcodes_and_filtered_matrix() {
    let p = Parameters {
        out_file_name_prefix: "run/".to_string(),
        ..Default::default()
    };
    let p_solo = ParametersSolo {
        cb_wl_size: 3,
        cb_wl_str: vec!["CBA".to_string(), "CBB".to_string(), "CBC".to_string()],
        out_file_names: vec![
            String::new(),
            "features.tsv".to_string(),
            "barcodes.tsv".to_string(),
            "matrix.mtx".to_string(),
        ],
        out_format_features_gene_field3: "Gene Expression".to_string(),
        umi_dedup: UMIdedup {
            types: vec![1],
            ..Default::default()
        },
        ..Default::default()
    };
    let trans = Transcriptome {
        n_ge: 2,
        ge_id: vec!["g1".to_string(), "g2".to_string()],
        ge_name: vec!["Gene1".to_string(), String::new()],
        ..Default::default()
    };
    let mut sf = SoloFeature {
        feature_type: SOLO_FEATURE_GENE,
        features_number: 2,
        n_cb: 2,
        ind_cb: vec![0, 2],
        n_gene_per_cb: vec![2, 1],
        count_cell_gene_umi_index: vec![0, 4],
        count_cell_gene_umi: vec![0, 5, 1, 2, 1, 7],
        count_mat_stride: 2,
        filtered_cells: SoloFilteredCells {
            filt_vec_bool: vec![false, true],
            n_cells: 1,
            ..Default::default()
        },
        ..Default::default()
    };

    let unfiltered = solofeature_outputresults_l12_solofeature_outputresults(
        &mut sf, false, "out/", &p, &p_solo, &trans, "/cwd",
    )
    .unwrap();
    assert_eq!(unfiltered.created_directory, "out/");
    assert_eq!(
        unfiltered.files["out/features.tsv"],
        "g1\tGene1\tGene Expression\ng2\tg2\tGene Expression\n"
    );
    assert_eq!(unfiltered.files["out/barcodes.tsv"], "CBA\nCBB\nCBC\n");
    assert_eq!(
        unfiltered.files["out/matrix.mtx"],
        "%%MatrixMarket matrix coordinate integer general\n%\n2 3 3\n1 1 5\n2 1 2\n2 3 7\n"
    );

    let filtered = solofeature_outputresults_l12_solofeature_outputresults(
        &mut sf, true, "flt/", &p, &p_solo, &trans, "/cwd",
    )
    .unwrap();
    assert_eq!(filtered.files["flt/barcodes.tsv"], "CBC\n");
    assert_eq!(
        filtered.files["flt/matrix.mtx"],
        "%%MatrixMarket matrix coordinate integer general\n%\n2 1 1\n2 1 7\n"
    );
}

#[test]
fn solo_feature_output_results_merges_unique_and_multimapper_counts() {
    let mut multi_map = MultiMappers {
        types: vec![1],
        yes_multi: true,
        count_ind_i: [u32::MAX; 5],
        ..Default::default()
    };
    multi_map.count_ind_i[1] = 1;
    let p_solo = ParametersSolo {
        cb_wl_size: 1,
        cb_wl_str: vec!["CB1".to_string()],
        out_file_names: vec![
            String::new(),
            "features.tsv".to_string(),
            "barcodes.tsv".to_string(),
            "matrix.mtx".to_string(),
        ],
        out_format_features_gene_field3: "-".to_string(),
        umi_dedup: UMIdedup {
            types: vec![1],
            yes_n: 1,
            ..Default::default()
        },
        multi_map,
        ..Default::default()
    };
    let trans = Transcriptome {
        n_ge: 3,
        ge_id: vec!["g1".to_string(), "g2".to_string(), "g3".to_string()],
        ge_name: vec![String::new(), String::new(), String::new()],
        ..Default::default()
    };
    let mut sf = SoloFeature {
        feature_type: SOLO_FEATURE_GENE,
        features_number: 3,
        n_cb: 1,
        ind_cb: vec![0],
        n_gene_per_cb: vec![2],
        count_cell_gene_umi_index: vec![0, 4],
        count_cell_gene_umi: vec![0, 5, 2, 1],
        count_mat_stride: 2,
        count_mat_mult_i: vec![0, 6],
        count_mat_mult_m: vec![0.0, 1.5, 0.0, 1.0, 4.0, 0.0],
        count_mat_mult_s: 3,
        ..Default::default()
    };

    let out = solofeature_outputresults_l12_solofeature_outputresults(
        &mut sf,
        false,
        "multi/",
        &Parameters::default(),
        &p_solo,
        &trans,
        "/cwd",
    )
    .unwrap();

    assert_eq!(
        out.files["multi/UniqueAndMult-Uniform.mtx"],
        "%%MatrixMarket matrix coordinate real general\n%\n3 1 3\n1 1 6.5\n2 1 4\n3 1 1\n"
    );
    assert_eq!(sf.n_umi_per_cb_multi, vec![5]);
    assert_eq!(sf.n_gene_per_cb_multi, vec![1]);
}

#[test]
fn solo_feature_cell_filtering_top_cells_computes_stats_and_outputs_filtered_matrix() {
    let p = Parameters::default();
    let p_solo = ParametersSolo {
        cb_wl_size: 3,
        cb_wl_str: vec!["CB1".to_string(), "CB2".to_string(), "CB3".to_string()],
        out_file_names: vec![
            String::new(),
            "features.tsv".to_string(),
            "barcodes.tsv".to_string(),
            "matrix.mtx".to_string(),
        ],
        out_format_features_gene_field3: "-".to_string(),
        umi_dedup: UMIdedup {
            count_ind_main: 1,
            types: vec![1],
            ..Default::default()
        },
        cell_filter: SoloCellFilter {
            type_: vec!["TopCells".to_string()],
            top_cells: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    let trans = Transcriptome {
        n_ge: 3,
        ge_id: vec!["g1".to_string(), "g2".to_string(), "g3".to_string()],
        ge_name: vec![String::new(), String::new(), String::new()],
        ..Default::default()
    };
    let mut sf = SoloFeature {
        feature_type: SOLO_FEATURE_GENE,
        features_number: 3,
        n_cb: 3,
        ind_cb: vec![0, 1, 2],
        n_gene_per_cb: vec![2, 1, 1],
        n_read_per_cb_total: vec![10, 8, 2],
        n_read_per_cb_unique: vec![9, 7, 1],
        n_umi_per_cb: vec![12, 8, 3],
        count_cell_gene_umi_index: vec![0, 4, 6],
        count_mat_stride: 2,
        count_cell_gene_umi: vec![0, 5, 1, 7, 2, 8, 1, 3],
        output_prefix_filtered: "flt/".to_string(),
        ..Default::default()
    };

    let result = solofeature_cellfiltering_l5_solofeature_cellfiltering(
        &mut sf, &p_solo, None, &p, &trans, "/cwd",
    )
    .unwrap();

    assert_eq!(
        result.log_main,
        "cellFiltering: simple: nUMImax=0; nUMImin=8; nCellsSimple=2\n"
    );
    assert!(!result.empty_drops_requested);
    assert_eq!(sf.filtered_cells.filt_vec_bool, vec![true, true, false]);
    assert_eq!(sf.filtered_cells.n_cells_simple, 2);
    assert_eq!(sf.filtered_cells.n_cells, 2);
    assert_eq!(sf.filtered_cells.n_umi_in_cells, 20);
    assert_eq!(sf.filtered_cells.mean_umi_per_cell, 10);
    assert_eq!(sf.filtered_cells.median_umi_per_cell, 8);
    assert_eq!(sf.filtered_cells.n_read_in_cells_unique, 16);
    assert_eq!(sf.filtered_cells.mean_read_per_cell_unique, 8);
    assert_eq!(sf.filtered_cells.median_read_per_cell_unique, 9);
    assert_eq!(sf.filtered_cells.n_gene_detected, 3);
    assert_eq!(sf.filtered_cells.n_gene_in_cells, 3);
    assert_eq!(sf.filtered_cells.mean_gene_per_cell, 1);
    assert_eq!(sf.filtered_cells.median_gene_per_cell, 2);
    let out = result.output_results.unwrap();
    assert_eq!(out.files["flt/barcodes.tsv"], "CB1\nCB2\n");
    assert_eq!(
        out.files["flt/matrix.mtx"],
        "%%MatrixMarket matrix coordinate integer general\n%\n3 2 3\n1 1 5\n2 1 7\n3 2 8\n"
    );
}

#[test]
fn solo_feature_cell_filtering_knee_and_emptydrops_request_follow_simple_filter() {
    let p = Parameters::default();
    let p_solo = ParametersSolo {
        cb_wl_size: 3,
        cb_wl_str: vec!["CB1".to_string(), "CB2".to_string(), "CB3".to_string()],
        out_file_names: vec![
            String::new(),
            "features.tsv".to_string(),
            "barcodes.tsv".to_string(),
            "matrix.mtx".to_string(),
        ],
        out_format_features_gene_field3: "-".to_string(),
        umi_dedup: UMIdedup {
            count_ind_main: 1,
            types: vec![1],
            ..Default::default()
        },
        cell_filter: SoloCellFilter {
            type_: vec!["EmptyDrops_CR".to_string()],
            knee: star_rs::generated::structs::SoloCellFilterKnee {
                n_expected_cells: 2.0,
                max_percentile: 0.5,
                max_min_ratio: 2.0,
            },
            ed_cr: star_rs::generated::structs::SoloCellFilterEmptyDropsCr {
                ind_min: 10,
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    let trans = Transcriptome {
        n_ge: 1,
        ge_id: vec!["g1".to_string()],
        ge_name: vec![String::new()],
        ..Default::default()
    };
    let mut sf = SoloFeature {
        feature_type: SOLO_FEATURE_GENE,
        features_number: 1,
        n_cb: 3,
        ind_cb: vec![0, 1, 2],
        n_gene_per_cb: vec![1, 1, 1],
        n_umi_per_cb: vec![10, 4, 1],
        count_cell_gene_umi_index: vec![0, 2, 4],
        count_mat_stride: 2,
        count_cell_gene_umi: vec![0, 10, 0, 4, 0, 1],
        output_prefix_filtered: "flt/".to_string(),
        ..Default::default()
    };

    let result = solofeature_cellfiltering_l5_solofeature_cellfiltering(
        &mut sf, &p_solo, None, &p, &trans, "/cwd",
    )
    .unwrap();

    assert!(result.empty_drops_requested);
    assert!(result.empty_drops.is_some());
    assert!(
        result
            .log_main
            .starts_with("cellFiltering: simple: nUMImax=4; nUMImin=2; nCellsSimple=2\n")
    );
    assert!(
        result
            .log_main
            .contains("emptyDrops_CR filtering: total number of cells: nCB=3 is smaller")
    );
    assert_eq!(sf.filtered_cells.filt_vec_bool, vec![true, true, false]);
}

#[test]
fn solo_feature_emptydrops_cr_promotes_candidate_cells() {
    let p_solo = ParametersSolo {
        umi_dedup: UMIdedup {
            count_ind_main: 1,
            ..Default::default()
        },
        cell_filter: SoloCellFilter {
            ed_cr: star_rs::generated::structs::SoloCellFilterEmptyDropsCr {
                ind_min: 3,
                ind_max: 5,
                umi_min: 1,
                umi_min_frac_median: 0.0,
                cand_max_n: 2,
                fdr: 1.0,
                sim_n: 5,
            },
            ..Default::default()
        },
        ..Default::default()
    };
    let mut sf = SoloFeature {
        features_number: 4,
        n_cb: 5,
        n_gene_per_cb: vec![1, 1, 1, 1, 1],
        n_umi_per_cb: vec![20, 15, 5, 2, 1],
        n_umi_per_cb_sorted: vec![20, 15, 5, 2, 1],
        count_cell_gene_umi_index: vec![0, 2, 4, 6, 8],
        count_mat_stride: 2,
        count_cell_gene_umi: vec![0, 20, 1, 15, 2, 5, 0, 2, 1, 1],
        filtered_cells: SoloFilteredCells {
            filt_vec_bool: vec![true, true, false, false, false],
            n_cells_simple: 2,
            ..Default::default()
        },
        ..Default::default()
    };

    let result = solofeature_emptydrops_cr_l10_solofeature_emptydrops_cr(&mut sf, &p_solo);

    assert_eq!(result.feat_det_n, 3);
    assert_eq!(result.min_umi, 1);
    assert_eq!(result.candidate_first, 2);
    assert_eq!(result.candidate_last, 3);
    assert_eq!(result.extra_cells, 2);
    assert_eq!(
        sf.filtered_cells.filt_vec_bool,
        vec![true, true, true, true, false]
    );
    assert_eq!(result.p_values.len(), 2);
    assert!(result.log_main.contains("finished emptyDrops_CR filtering"));
}

#[test]
fn solo_feature_sum_threads_builds_undefined_whitelist_and_detected_cb_indexes() {
    let mut p_solo = ParametersSolo {
        cb_wl_yes: false,
        cb_l: 2,
        cb_type_type: 1,
        cb_match_wl: CBMatchWL {
            mm1_multi_pc: true,
            ..Default::default()
        },
        ..Default::default()
    };
    let p = Parameters {
        run_thread_n: 2,
        ..Default::default()
    };
    let mut sf = SoloFeature {
        feature_type: SOLO_FEATURE_GENE,
        read_feat_sum: Some(SoloReadFeature {
            cb_wl_yes: false,
            ..Default::default()
        }),
        ..Default::default()
    };
    let read_feat_all = vec![
        SoloReadFeature {
            cb_read_count_map: std::collections::BTreeMap::from([(1, 2)]),
            ..Default::default()
        },
        SoloReadFeature {
            cb_read_count_map: std::collections::BTreeMap::from([(1, 1), (6, 3)]),
            ..Default::default()
        },
    ];
    let mut read_bar_sum = SoloReadBarcode::default();

    solofeature_sumthreads_l8_solofeature_sumthreads(
        &mut sf,
        &p,
        &mut p_solo,
        &mut read_bar_sum,
        &read_feat_all,
        41,
    );

    assert_eq!(sf.n_reads_input, 42);
    assert_eq!(sf.read_feat_all_len, 2);
    assert_eq!(p_solo.cb_wl_size, 2);
    assert_eq!(p_solo.cb_wl, vec![1, 6]);
    assert_eq!(p_solo.cb_wl_str, vec!["AC", "CG"]);
    assert_eq!(sf.read_feat_sum.as_ref().unwrap().cb_read_count, vec![3, 3]);
    assert_eq!(read_bar_sum.cb_read_count_exact, vec![4, 4]);
    assert_eq!(sf.n_cb, 2);
    assert_eq!(sf.n_reads_mapped, 6);
    assert_eq!(sf.ind_cb, vec![0, 1]);
    assert_eq!(sf.ind_cb_wl, vec![0, 1]);
}

#[test]
fn solo_feature_sum_threads_restart_recounts_stream_cb_column() {
    let mut p_solo = ParametersSolo {
        cb_wl_yes: true,
        cb_wl_size: 3,
        ..Default::default()
    };
    let p = Parameters {
        run_restart_type: 1,
        ..Default::default()
    };
    let mut sf = SoloFeature {
        feature_type: SOLO_FEATURE_SJ,
        read_feat_sum: Some(SoloReadFeature {
            cb_wl_yes: true,
            cb_wl_size: 3,
            cb_read_count: vec![5, 0, 1],
            ..Default::default()
        }),
        ..Default::default()
    };
    let read_feat_all = vec![SoloReadFeature {
        cb_wl_yes: true,
        cb_wl_size: 3,
        cb_read_count: vec![0, 0, 0],
        stream_reads: "9 8 7 6 1\n4 3 2 1 2\n".to_string(),
        ..Default::default()
    }];
    let mut read_bar_sum = SoloReadBarcode::default();

    solofeature_sumthreads_l8_solofeature_sumthreads(
        &mut sf,
        &p,
        &mut p_solo,
        &mut read_bar_sum,
        &read_feat_all,
        0,
    );

    assert_eq!(
        sf.read_feat_sum.as_ref().unwrap().cb_read_count,
        vec![5, 1, 2]
    );
    assert_eq!(sf.n_cb, 3);
    assert_eq!(sf.n_reads_mapped, 8);
    assert_eq!(sf.ind_cb, vec![0, 1, 2]);
    assert_eq!(sf.ind_cb_wl, vec![0, 1, 2]);
}

#[test]
fn parameters_sam_attributes_expands_presets_and_quant_orders() {
    let mut p = Parameters {
        out_sam_attributes: vec!["Standard".to_string()],
        out_sam_attr_present: star_rs::generated::structs::SamAttrPresent {
            v_w: true,
            ..Default::default()
        },
        out_sam_attr_order: vec![99],
        out_sam_attr_order_quant: vec![98],
        out_sam_attr_rgline: vec!["-".to_string()],
        ..Default::default()
    };

    let log = parameters_samattributes_l4_parameters_samattributes(&mut p).unwrap();

    assert!(log.is_empty());
    assert_eq!(
        p.out_sam_attr_order,
        vec![ATTR_NH, ATTR_HI, ATTR_AS, ATTR_NM_LOWER]
    );
    assert_eq!(p.out_sam_attr_order_quant, vec![ATTR_NH, ATTR_HI]);
    assert!(p.out_sam_attr_present.nh);
    assert!(p.out_sam_attr_present.hi);
    assert!(p.out_sam_attr_present.as_);
    assert!(p.out_sam_attr_present.n_m);
    assert!(!p.out_sam_attr_present.v_w);
}

#[test]
fn parameters_sam_attributes_handles_rg_xs_wasp_and_bam_warnings() {
    let mut p = Parameters {
        out_sam_attributes: vec!["XS".to_string(), "CB".to_string()],
        out_sam_attr_rgline: vec!["ID:rg1".to_string()],
        out_bam_unsorted: true,
        out_sam_bool: true,
        wasp_yes: true,
        ..Default::default()
    };

    let log = parameters_samattributes_l4_parameters_samattributes(&mut p).unwrap();

    assert_eq!(
        p.out_sam_attr_order,
        vec![ATTR_XS, ATTR_CB, ATTR_RG, ATTR_VW]
    );
    assert_eq!(
        p.out_sam_attr_order_quant,
        vec![ATTR_NH, ATTR_HI, ATTR_CB, ATTR_RG, ATTR_VW]
    );
    assert_eq!(p.out_sam_strand_field_type, 1);
    assert!(p.out_sam_attr_present.xs);
    assert!(p.out_sam_attr_present.cb);
    assert!(p.out_sam_attr_present.rg);
    assert!(p.out_sam_attr_present.v_w);
    assert!(log.contains("outSAMattributes contains XS"));
    assert!(log.contains("outSAMattrRG defines a read group"));
    assert!(log.contains("waspOutputMode is set"));
    assert!(log.contains("CB tag. It will be output into BAM file"));
    assert!(log.contains("vW tag. It will be output into BAM file"));
}

#[test]
fn parameters_sam_attributes_reports_original_errors_and_enables_gene_quant() {
    let mut unknown = Parameters {
        out_sam_attributes: vec!["bad".to_string()],
        ..Default::default()
    };
    assert_eq!(
        parameters_samattributes_l4_parameters_samattributes(&mut unknown).unwrap_err(),
        "EXITING because of FATAL INPUT ERROR: unknown/unimplemented SAM atrribute (tag): bad\nSOLUTION: re-run STAR with --outSAMattributes that contains only implemented attributes\n"
    );

    let mut no_var = Parameters {
        out_sam_attributes: vec!["vA".to_string()],
        ..Default::default()
    };
    assert!(
        parameters_samattributes_l4_parameters_samattributes(&mut no_var)
            .unwrap_err()
            .contains("--varVCFfile is not set")
    );

    let mut rg_missing = Parameters {
        out_sam_attributes: vec!["RG".to_string()],
        out_sam_attr_rgline: vec!["-".to_string()],
        ..Default::default()
    };
    assert!(
        parameters_samattributes_l4_parameters_samattributes(&mut rg_missing)
            .unwrap_err()
            .contains("--outSAMattrRGline is not set")
    );

    let mut bam_required = Parameters {
        out_sam_attributes: vec!["CB".to_string()],
        ..Default::default()
    };
    assert!(
        parameters_samattributes_l4_parameters_samattributes(&mut bam_required)
            .unwrap_err()
            .contains("requires BAM output")
    );

    let mut gx = Parameters {
        out_sam_attributes: vec!["GX".to_string()],
        out_bam_coord: true,
        var_yes: true,
        wasp_yes: true,
        ..Default::default()
    };
    parameters_samattributes_l4_parameters_samattributes(&mut gx).unwrap();
    assert!(gx.quant_gene_yes);
    assert!(gx.quant_yes);
    assert_eq!(gx.p_solo.sam_attr_feature, SOLO_FEATURE_GENE);
}

#[test]
fn sjdb_load_from_files_reads_each_file_and_sets_priority() {
    let dir = std::env::temp_dir().join(format!("star-rs-sjdb-files-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let file1 = dir.join("sj1.tab");
    let file2 = dir.join("sj2.tab");
    std::fs::write(&file1, "chr1\t10\t20\t1\nchr2\t30\t40\t-\n").unwrap();
    std::fs::write(&file2, "chr3\t50\t60\t?\n").unwrap();

    let p = Parameters {
        p_ge: ParametersGenome {
            sjdb_file_chr_start_end: vec![
                file1.to_string_lossy().into_owned(),
                file2.to_string_lossy().into_owned(),
            ],
            ..Default::default()
        },
        ..Default::default()
    };
    let mut loci = SjdbClass {
        priority: vec![3],
        ..Default::default()
    };

    let log = sjdbloadfromfiles_l6_sjdbloadfromfiles(&p, &mut loci).unwrap();

    assert_eq!(loci.chr, vec!["chr1", "chr2", "chr3"]);
    assert_eq!(loci.start, vec![10, 30, 50]);
    assert_eq!(loci.end, vec![20, 40, 60]);
    assert_eq!(loci.str_, vec!['+', '-', '.']);
    assert_eq!(loci.priority, vec![3, 10, 10]);
    assert!(log.contains("total number of junctions:2"));
    assert!(log.contains("total number of junctions:3"));
}

#[test]
fn sjdb_load_from_files_obeys_dash_sentinel_and_reports_missing_files() {
    let mut loci = SjdbClass {
        chr: vec!["kept".to_string()],
        priority: vec![8],
        ..Default::default()
    };
    let p_skip = Parameters {
        p_ge: ParametersGenome {
            sjdb_file_chr_start_end: vec!["-".to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    assert_eq!(
        sjdbloadfromfiles_l6_sjdbloadfromfiles(&p_skip, &mut loci).unwrap(),
        ""
    );
    assert_eq!(loci.chr, vec!["kept"]);
    assert_eq!(loci.priority, vec![8]);

    let missing = "/tmp/star-rs-definitely-missing-sjdb-file.tab";
    let p_missing = Parameters {
        p_ge: ParametersGenome {
            sjdb_file_chr_start_end: vec![missing.to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    assert_eq!(
        sjdbloadfromfiles_l6_sjdbloadfromfiles(&p_missing, &mut SjdbClass::default()).unwrap_err(),
        format!(
            "FATAL INPUT error, could not open input file pGe.sjdbFileChrStartEnd={}\n",
            missing
        )
    );
}

#[test]
fn sjdb_load_from_stream_parses_rows_and_normalizes_strand() {
    let mut loci = SjdbClass::default();

    sjdbloadfromstream_l2_sjdbloadfromstream(
        "chr1\t10\t20\t1\nchr2 30 40 2\n\nchr3 50 60 +\nchr4 70 80 -\nchr5 90 100 ?\n",
        &mut loci,
    );

    assert_eq!(loci.chr, vec!["chr1", "chr2", "chr3", "chr4", "chr5"]);
    assert_eq!(loci.start, vec![10, 30, 50, 70, 90]);
    assert_eq!(loci.end, vec![20, 40, 60, 80, 100]);
    assert_eq!(loci.str_, vec!['+', '-', '+', '-', '.']);
}

#[test]
fn transcript_pe_overlap_se_to_pe_splits_and_copies_metadata() {
    let source = Transcript {
        exons: vec![[5, 100, 20, 0, 7], [25, 200, 10, 0, 9]],
        canon_sj: vec![1, 0],
        sj_annot: vec![1, 0],
        sj_str: vec![2, 0],
        shift_sj: vec![[3, 4], [0, 0]],
        intron_motifs: [1, 2, 3],
        sj_motif_strand: 6,
        n_exons: 2,
        l_read: 40,
        chr: 11,
        str_: 0,
        ro_str: 1,
        g_start: 100,
        g_length: 110,
        c_start: 99,
        n_gap: 4,
        l_gap: 5,
        n_del: 6,
        l_del: 66,
        n_ins: 7,
        l_ins: 8,
        n_unique: 9,
        n_anchor: 10,
        ..Default::default()
    };
    let mut out = Transcript {
        read_length: vec![30, 20],
        l_read: 51,
        exons: vec![[u32::MAX; 5]; 4],
        canon_sj: vec![99; 4],
        sj_annot: vec![9; 4],
        sj_str: vec![9; 4],
        shift_sj: vec![[9, 9]; 4],
        ..Default::default()
    };

    readalign_peoverlapmergemap_l136_transcript_peoverlapsetope(&mut out, &[0, 20], &source);

    assert_eq!(out.n_exons, 4);
    assert_eq!(out.exons[0], [5, 100, 20, 0, 7]);
    assert_eq!(out.exons[1], [25, 200, 5, 0, 9]);
    assert_eq!(out.exons[2], [31, 115, 5, 1, 7]);
    assert_eq!(out.exons[3], [36, 200, 10, 1, 9]);
    assert_eq!(out.canon_sj[0], 1);
    assert_eq!(out.canon_sj[1], -3);
    assert_eq!(out.canon_sj[2], 1);
    assert_eq!(out.canon_sj[3], -3);
    assert_eq!(out.sj_annot[..4], [1, 0, 1, 0]);
    assert_eq!(out.sj_str[..4], [2, 0, 2, 0]);
    assert_eq!(out.shift_sj[..4], [[3, 4], [0, 0], [3, 4], [0, 0]]);
    assert_eq!(out.intron_motifs, [1, 2, 3]);
    assert_eq!(out.sj_motif_strand, 6);
    assert_eq!(out.chr, 11);
    assert_eq!(out.str_, 0);
    assert_eq!(out.ro_start, 6);
    assert_eq!(out.r_length, 40);
    assert_eq!(out.mapped_length, 40);
    assert_eq!(out.r_start, 5);
    assert_eq!(out.l_del, 6);
    assert_eq!(out.l_ins, 8);
}

#[test]
fn read_align_pe_merge_mates_merges_best_overlap_and_rebuilds_reverse_buffers() {
    let mut ra = ReadAlign {
        l_read: 9,
        read_nmates: 2,
        read_length: vec![4, 4],
        read_length_original: vec![4, 4],
        read1: [
            vec![0, 1, 2, 3, 5, 2, 3, 0, 0],
            vec![3, 2, 1, 0, 5, 1, 0, 3, 3],
            vec![3, 3, 0, 1, 5, 0, 1, 2, 3],
        ],
        ..Default::default()
    };

    readalign_peoverlapmergemap_l79_readalign_pemergemates(&mut ra, 0.0, 2);

    assert_eq!(ra.pe_ov.n_ov, 2);
    assert_eq!(ra.pe_ov.mate_start, [0, 2]);
    assert_eq!(ra.l_read, 6);
    assert_eq!(ra.read_length, vec![6, 0]);
    assert_eq!(ra.read_length_original, vec![6, 0]);
    assert_eq!(ra.read_nmates, 1);
    assert_eq!(&ra.read1[0][..6], &[0, 1, 2, 3, 0, 0]);
    assert_eq!(&ra.read1[1][..6], &[3, 2, 1, 0, 3, 3]);
    assert_eq!(&ra.read1[2][..6], &[3, 3, 0, 1, 2, 3]);

    let mut ra = ReadAlign {
        l_read: 9,
        read_nmates: 2,
        read_length: vec![4, 4],
        read_length_original: vec![4, 4],
        read1: [vec![2, 3, 0, 0, 5, 0, 1, 2, 3], Vec::new(), Vec::new()],
        ..Default::default()
    };

    readalign_peoverlapmergemap_l79_readalign_pemergemates(&mut ra, 0.0, 2);

    assert_eq!(ra.pe_ov.n_ov, 2);
    assert_eq!(ra.pe_ov.mate_start, [2, 0]);
    assert_eq!(ra.l_read, 6);
    assert_eq!(&ra.read1[0][..6], &[0, 1, 2, 3, 0, 0]);

    let mut too_short = ReadAlign {
        l_read: 9,
        read_nmates: 2,
        read_length: vec![4, 4],
        read_length_original: vec![4, 4],
        read1: [
            vec![0, 1, 2, 3, 5, 2, 3, 0, 0],
            vec![3, 2, 1, 0, 5, 1, 0, 3, 3],
            vec![3, 3, 0, 1, 5, 0, 1, 2, 3],
        ],
        ..Default::default()
    };
    readalign_peoverlapmergemap_l79_readalign_pemergemates(&mut too_short, 0.0, 3);
    assert_eq!(too_short.pe_ov.n_ov, 0);
    assert_eq!(too_short.l_read, 9);
    assert_eq!(too_short.read_nmates, 2);
}

#[test]
fn read_align_pe_overlap_se_to_pe_converts_scores_and_selects_best_windows() {
    let tr_init = Transcript {
        read_length: vec![5, 5],
        l_read: 11,
        read_nmates: 2,
        exons: vec![[0; EX_SIZE]; 4],
        canon_sj: vec![0; 4],
        sj_annot: vec![0; 4],
        sj_str: vec![0; 4],
        shift_sj: vec![[0, 0]; 4],
        ..Default::default()
    };
    let mut ra = ReadAlign {
        l_read: 11,
        read1: [
            vec![0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0],
            Vec::new(),
            vec![0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0],
        ],
        pe_ov: star_rs::generated::structs::ReadAlignPeOverlap {
            mate_start: [0, 4],
            ..Default::default()
        },
        ..Default::default()
    };
    let bad = Transcript {
        n_exons: 1,
        l_read: 9,
        read_length: vec![5, 5],
        str_: 0,
        exons: vec![[0, 20, 5, 0, 0]],
        ..Default::default()
    };
    let good = Transcript {
        n_exons: 1,
        l_read: 9,
        read_length: vec![5, 5],
        str_: 0,
        exons: vec![[0, 10, 5, 0, 0]],
        ..Default::default()
    };
    let se_ra = ReadAlign {
        n_w: 2,
        n_win_tr: vec![2, 1],
        tr_all: vec![vec![bad, good.clone()], vec![good.clone()]],
        ..Default::default()
    };
    let mut genome = vec![1_u8; 32];
    for base in &mut genome[10..16] {
        *base = 0;
    }

    readalign_peoverlapmergemap_l266_readalign_peoverlapsetope(
        &mut ra, &se_ra, &tr_init, &genome, 0, -1, -2, -1, -2, -8, -4, -6, -8, 0.0,
    );

    assert_eq!(ra.n_w, 2);
    assert_eq!(ra.n_win_tr, vec![2, 1]);
    assert_eq!(ra.tr_all.len(), 2);
    assert_eq!(ra.tr_all[0].len(), 2);
    assert_eq!(ra.tr_all[0][0].exons[0][EX_G], 10);
    assert!(ra.tr_all[0][0].max_score > ra.tr_all[0][1].max_score);
    assert_eq!(ra.tr_best, ra.tr_all[0][0]);
    assert_eq!(ra.tr_best.n_exons, 2);
}

#[test]
fn read_align_pe_overlap_merge_map_requests_mapping_and_accepts_better_merged_alignment() {
    let p = Parameters {
        read_nmates: 2,
        pe_overlap_nbases_min: 2,
        pe_overlap_mmp: 0.0,
        p_ch: ParametersChimeric {
            segment_min: 1,
            multimap_nmax: 8,
            nonchim_score_drop_min: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    let tr_init = Transcript {
        read_length: vec![4, 4],
        l_read: 9,
        read_nmates: 2,
        exons: vec![[0; EX_SIZE]; 4],
        canon_sj: vec![0; 4],
        sj_annot: vec![0; 4],
        sj_str: vec![0; 4],
        shift_sj: vec![[0, 0]; 4],
        ..Default::default()
    };
    let mut ra = ReadAlign {
        l_read: 9,
        read_nmates: 2,
        read_length: vec![4, 4],
        read_length_original: vec![4, 4],
        tr_best: Transcript {
            max_score: 1,
            ..Default::default()
        },
        read1: [
            vec![0, 1, 2, 3, 5, 2, 3, 0, 0],
            vec![3, 2, 1, 0, 5, 1, 0, 3, 3],
            vec![3, 3, 0, 1, 5, 0, 1, 2, 3],
        ],
        ..Default::default()
    };
    let mut pe_merge_ra = ReadAlign::default();
    let empty_mapped = ReadAlign {
        n_w: 0,
        l_read: 6,
        ..Default::default()
    };
    let empty_genome = Genome {
        g: vec![0; 32],
        ..Default::default()
    };

    let pending = readalign_peoverlapmergemap_l4_readalign_peoverlapmergemap(
        &mut ra,
        &mut pe_merge_ra,
        Some(&empty_mapped),
        &p,
        &empty_genome,
        &tr_init,
        &[0; 32],
        0,
        -1,
        -2,
        -1,
        -2,
        -8,
        -4,
        -6,
        -8,
        0.0,
        Some(true),
    )
    .unwrap();
    assert!(pending.map_one_read_requested);
    assert_eq!(ra.pe_ov.n_ov, 2);
    assert!(!ra.pe_ov.yes);
    assert_eq!(pe_merge_ra.l_read, 6);
    ra.tr_best.max_score = -100;

    let mapped = ReadAlign {
        n_w: 1,
        n_win_tr: vec![1],
        tr_all: vec![vec![Transcript {
            n_exons: 1,
            l_read: 6,
            read_length: vec![6, 0],
            str_: 0,
            exons: vec![[0, 10, 6, 0, 0]],
            ..Default::default()
        }]],
        tr_best: Transcript {
            max_score: 6,
            ..Default::default()
        },
        read_length: vec![6, 0],
        ..pe_merge_ra.clone()
    };
    let mut genome = vec![1_u8; 32];
    for base in &mut genome[10..16] {
        *base = 0;
    }
    let done = readalign_peoverlapmergemap_l4_readalign_peoverlapmergemap(
        &mut ra,
        &mut pe_merge_ra,
        Some(&mapped),
        &p,
        &Genome {
            g: genome.clone(),
            ..Default::default()
        },
        &tr_init,
        &genome,
        0,
        -1,
        -2,
        -1,
        -2,
        -8,
        -4,
        -6,
        -8,
        0.0,
        None,
    )
    .unwrap();
    assert!(done.map_one_read_requested);
    assert!(ra.pe_ov.yes);
    assert_eq!(ra.n_w, 1);
    assert_eq!(ra.tr_best.n_exons, 2);
    assert!(done.chimeric_detection.is_some());

    let mut disabled = ReadAlign {
        pe_ov: ReadAlignPeOverlap {
            yes: true,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut merge = ReadAlign::default();
    let no_op = readalign_peoverlapmergemap_l4_readalign_peoverlapmergemap(
        &mut disabled,
        &mut merge,
        None,
        &Parameters::default(),
        &Genome::default(),
        &Transcript::default(),
        &[],
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0.0,
        None,
    )
    .unwrap();
    assert!(!no_op.map_one_read_requested);
    assert!(!disabled.pe_ov.yes);
}

#[test]
fn transcript_pe_overlap_se_to_pe_handles_reverse_strand_mate_coordinates() {
    let source = Transcript {
        exons: vec![[0, 100, 15, 0, 1], [20, 200, 10, 0, 2], [32, 300, 5, 0, 3]],
        canon_sj: vec![1, 2, 0],
        sj_annot: vec![1, 1, 0],
        sj_str: vec![1, 2, 0],
        shift_sj: vec![[1, 1], [2, 2], [0, 0]],
        n_exons: 3,
        l_read: 40,
        str_: 1,
        ro_str: 0,
        ..Default::default()
    };
    let mut out = Transcript {
        read_length: vec![30, 20],
        l_read: 51,
        ..Default::default()
    };

    readalign_peoverlapmergemap_l136_transcript_peoverlapsetope(&mut out, &[0, 20], &source);

    assert_eq!(out.n_exons, 4);
    assert_eq!(out.exons[0], [0, 100, 15, 1, 1]);
    assert_eq!(out.exons[1], [21, 110, 5, 0, 1]);
    assert_eq!(out.exons[2], [31, 200, 10, 0, 2]);
    assert_eq!(out.exons[3], [43, 300, 5, 0, 3]);
    assert_eq!(out.canon_sj[..4], [-3, 1, 2, -3]);
    assert_eq!(out.sj_annot[..4], [0, 1, 1, 0]);
    assert_eq!(out.sj_str[..4], [0, 1, 2, 0]);
    assert_eq!(out.r_length, 35);
    assert_eq!(out.r_start, 0);
    assert_eq!(out.ro_start, 0);
}

#[test]
fn read_align_pe_overlap_chimeric_se_to_pe_drops_shorter_mate_segment() {
    let tr_init = Transcript {
        read_length: vec![10, 10],
        l_read: 21,
        ..Default::default()
    };
    let se1 = Transcript {
        exons: vec![[0, 100, 5, 0, 0], [12, 200, 5, 0, 0]],
        canon_sj: vec![1, 0],
        sj_annot: vec![1, 0],
        sj_str: vec![2, 0],
        shift_sj: vec![[3, 4], [0, 0]],
        n_exons: 2,
        l_read: 21,
        str_: 0,
        ..Default::default()
    };
    let se2 = Transcript {
        exons: vec![[0, 300, 6, 0, 0], [12, 400, 4, 0, 0]],
        canon_sj: vec![2, 0],
        sj_annot: vec![1, 0],
        sj_str: vec![1, 0],
        shift_sj: vec![[5, 6], [0, 0]],
        n_exons: 2,
        l_read: 21,
        str_: 0,
        ..Default::default()
    };
    let mut out1 = Transcript::default();
    let mut out2 = Transcript::default();

    readalign_peoverlapmergemap_l308_readalign_peoverlapchimericsetope(
        &tr_init,
        &[0, 10],
        &[10, 10],
        &se1,
        &se2,
        &mut out1,
        &mut out2,
    );

    assert_eq!(out1.n_exons, 2);
    assert_eq!(out1.exons[0], [0, 100, 5, 0, 0]);
    assert_eq!(out1.exons[1], [13, 200, 5, 1, 0]);
    assert_eq!(out2.n_exons, 1);
    assert_eq!(out2.exons, vec![[0, 300, 6, 0, 0]]);
}

#[test]
fn output_read_cb_matches_gene_sj_transcript_and_read_info_formats() {
    let mut solo_bar = SoloReadBarcode {
        solo_type: 1,
        umi_b: 17,
        cb_match: 0,
        cb_match_string: "CB:QUAL".to_string(),
        ..Default::default()
    };
    let read_flag = SoloReadFlagClass {
        flag: 5,
        ..Default::default()
    };

    let mut out = String::new();
    assert_eq!(
        soloreadfeature_record_l206_outputreadcb(
            &mut out,
            123,
            SOLO_FEATURE_GENE,
            &mut solo_bar,
            &ReadSoloFeatures {
                gene: 42,
                ..Default::default()
            },
            &ReadAnnotations::default(),
            &read_flag,
        ),
        1
    );
    assert_eq!(out, "17 123 5 42 0 CB:QUAL\n");

    out.clear();
    assert_eq!(
        soloreadfeature_record_l206_outputreadcb(
            &mut out,
            u64::MAX,
            SOLO_FEATURE_GENE_FULL,
            &mut solo_bar,
            &ReadSoloFeatures {
                gene: 77,
                ..Default::default()
            },
            &ReadAnnotations::default(),
            &read_flag,
        ),
        1
    );
    assert_eq!(out, "17 77 0 CB:QUAL\n");

    out.clear();
    assert_eq!(
        soloreadfeature_record_l206_outputreadcb(
            &mut out,
            88,
            SOLO_FEATURE_GENE,
            &mut solo_bar,
            &ReadSoloFeatures {
                gene_mult: vec![101, 202],
                ..Default::default()
            },
            &ReadAnnotations::default(),
            &read_flag,
        ),
        2
    );
    assert_eq!(out, "17 88 5 101 0 CB:QUAL\n17 88 5 202 0 CB:QUAL\n");

    out.clear();
    assert_eq!(
        soloreadfeature_record_l206_outputreadcb(
            &mut out,
            u64::MAX,
            SOLO_FEATURE_SJ,
            &mut solo_bar,
            &ReadSoloFeatures {
                sj: vec![[100, 5], [200, 9]],
                ..Default::default()
            },
            &ReadAnnotations::default(),
            &read_flag,
        ),
        2
    );
    assert_eq!(out, "17 100 5 0 CB:QUAL\n17 200 9 0 CB:QUAL\n");

    out.clear();
    assert_eq!(
        soloreadfeature_record_l206_outputreadcb(
            &mut out,
            55,
            SOLO_FEATURE_TRANSCRIPT3P,
            &mut solo_bar,
            &ReadSoloFeatures::default(),
            &ReadAnnotations {
                transcript_concordant: vec![[7, 300], [8, 400]],
                ..Default::default()
            },
            &read_flag,
        ),
        1
    );
    assert_eq!(out, "CB:QUAL 17 2 7 300 8 400 55\n");

    out.clear();
    assert_eq!(
        soloreadfeature_record_l206_outputreadcb(
            &mut out,
            66,
            -1,
            &mut solo_bar,
            &ReadSoloFeatures::default(),
            &ReadAnnotations::default(),
            &read_flag,
        ),
        1
    );
    assert_eq!(out, "17 66 5 -1 0 CB:QUAL\n");
}

#[test]
fn output_read_cb_smartseq_replaces_umi_with_extended_alignment_locus() {
    let mut solo_bar = SoloReadBarcode {
        solo_type: 4,
        umi_b: 99,
        cb_match: 1,
        cb_match_string: "SMART".to_string(),
        ..Default::default()
    };
    let align = Transcript {
        c_start: 1000,
        l_read: 50,
        n_exons: 2,
        exons: vec![[10, 1000, 20, 0, 0], [35, 1100, 15, 0, 0]],
        ..Default::default()
    };
    let expected_umi = transcript_l53_transcript_chrstartlengthextended(&align);
    let mut out = String::new();

    soloreadfeature_record_l206_outputreadcb(
        &mut out,
        10,
        SOLO_FEATURE_GENE,
        &mut solo_bar,
        &ReadSoloFeatures {
            gene: 9,
            ind_annot_tr: 0,
            align_out: vec![align],
            ..Default::default()
        },
        &ReadAnnotations::default(),
        &SoloReadFlagClass {
            flag: 3,
            ..Default::default()
        },
    );

    assert_eq!(solo_bar.umi_b, expected_umi);
    assert_eq!(out, format!("{} 10 3 9 1 SMART\n", expected_umi));
}

#[test]
fn solo_read_feature_record_outputs_unique_gene_and_counts_cell_barcode() {
    let mut p = Parameters {
        p_solo: ParametersSolo {
            solo_type: 1,
            cb_wl_yes: true,
            cb_wl_size: 3,
            ..Default::default()
        },
        ..Default::default()
    };
    p.p_solo.read_stats_yes.resize(8, false);

    let mut rf = SoloReadFeature {
        feature_type: SOLO_FEATURE_GENE,
        read_index_yes: true,
        cb_read_count: vec![0; 3],
        ..Default::default()
    };
    let mut solo_bar = SoloReadBarcode {
        solo_type: 1,
        umi_b: 19,
        cb_match: 0,
        cb_match_string: "CB".to_string(),
        cb_match_ind: vec![2],
        ..Default::default()
    };
    let mut ann_features = vec![ReadAnnotFeature::default(); 8];
    ann_features[SOLO_FEATURE_GENE as usize].f_set = std::collections::BTreeSet::from([44]);

    soloreadfeature_record_l20_soloreadfeature_record(
        &mut rf,
        &p,
        &mut solo_bar,
        1,
        &[Transcript::default()],
        77,
        &ReadAnnotations {
            annot_features: ann_features,
            ..Default::default()
        },
    );

    assert_eq!(rf.stream_reads, "19 77 64 44 0 CB\n");
    assert_eq!(rf.cb_read_count, vec![0, 0, 1]);
}

#[test]
fn solo_read_feature_record_counts_no_cb_read_stats_before_returning() {
    let mut p = Parameters {
        p_solo: ParametersSolo {
            solo_type: 1,
            cb_wl_yes: true,
            ..Default::default()
        },
        p_ge: ParametersGenome {
            chr_set_mito: std::collections::BTreeSet::from([2]),
            ..Default::default()
        },
        ..Default::default()
    };
    p.p_solo.read_stats_yes.resize(8, false);
    p.p_solo.read_stats_yes[SOLO_FEATURE_GENE as usize] = true;

    let mut rf = SoloReadFeature {
        feature_type: SOLO_FEATURE_GENE,
        ..Default::default()
    };
    let mut solo_bar = SoloReadBarcode {
        cb_match: -1,
        ..Default::default()
    };
    let mut ann_features = vec![ReadAnnotFeature::default(); 8];
    ann_features[SOLO_FEATURE_GENE as usize].f_set = std::collections::BTreeSet::from([11]);
    ann_features[SOLO_FEATURE_GENE as usize].ov_type = 1;

    soloreadfeature_record_l20_soloreadfeature_record(
        &mut rf,
        &p,
        &mut solo_bar,
        1,
        &[Transcript {
            chr: 2,
            ..Default::default()
        }],
        3,
        &ReadAnnotations {
            annot_features: ann_features,
            ..Default::default()
        },
    );

    assert_eq!(rf.stream_reads, "");
    assert_eq!(
        rf.read_flag.flag_counts_no_cb[SOLO_READ_FLAG_CB_MATCH as usize],
        1
    );
    assert_eq!(
        rf.read_flag.flag_counts_no_cb[SOLO_READ_FLAG_GENOME_U as usize],
        1
    );
    assert_eq!(
        rf.read_flag.flag_counts_no_cb[SOLO_READ_FLAG_FEATURE_U as usize],
        1
    );
    assert_eq!(
        rf.read_flag.flag_counts_no_cb[SOLO_READ_FLAG_EXONIC as usize],
        1
    );
    assert_eq!(
        rf.read_flag.flag_counts_no_cb[SOLO_READ_FLAG_MITO as usize],
        1
    );
}

#[test]
fn solo_read_feature_record_outputs_splice_junctions_and_no_wl_counts() {
    let p = Parameters {
        p_solo: ParametersSolo {
            solo_type: 1,
            cb_wl_yes: false,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut rf = SoloReadFeature {
        feature_type: SOLO_FEATURE_SJ,
        read_index_yes: false,
        ..Default::default()
    };
    let mut solo_bar = SoloReadBarcode {
        umi_b: 33,
        cb_match: 0,
        cb_match_string: "CB".to_string(),
        cb_match_ind: vec![99],
        ..Default::default()
    };

    soloreadfeature_record_l20_soloreadfeature_record(
        &mut rf,
        &p,
        &mut solo_bar,
        1,
        &[Transcript {
            n_exons: 2,
            exons: vec![[0, 10, 5, 0, 0], [5, 90, 7, 0, 0]],
            canon_sj: vec![1],
            sj_annot: vec![0],
            ..Default::default()
        }],
        7,
        &ReadAnnotations::default(),
    );

    assert_eq!(rf.stream_reads, "33 15 75 0 CB\n");
    assert_eq!(rf.cb_read_count_map.get(&99), Some(&1));
}

#[test]
fn solo_read_feature_record_outputs_sorted_velocyto_records() {
    let p = Parameters {
        p_solo: ParametersSolo {
            solo_type: 1,
            cb_wl_yes: false,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut rf = SoloReadFeature {
        feature_type: SOLO_FEATURE_VELOCYTO,
        ..Default::default()
    };
    let mut solo_bar = SoloReadBarcode {
        cb_match: 0,
        cb_match_ind: vec![7],
        ..Default::default()
    };

    soloreadfeature_record_l20_soloreadfeature_record(
        &mut rf,
        &p,
        &mut solo_bar,
        1,
        &[Transcript::default()],
        55,
        &ReadAnnotations {
            tr_velocyto_type: vec![
                TrTypeStruct { tr: 4, type_: 3 },
                TrTypeStruct { tr: 1, type_: 8 },
            ],
            ..Default::default()
        },
    );

    assert_eq!(rf.stream_reads, "55 2 1 8 4 3\n");
    assert_eq!(rf.cb_read_count_map.get(&7), Some(&1));
}

#[test]
fn solo_input_feature_umi_matches_formatted_stream_extraction() {
    let sj_all = [vec![100, 200, 300], vec![10, 20, 30]];
    let mut iread = 0;
    let mut cbmatch = 0;
    let mut feature = 999;
    let mut umi = 0;
    let mut feat_vec = vec![99];
    let mut flags = star_rs::generated::structs::SoloReadFlagClass::default();

    let mut gene_tokens = "17 123 5 42 1".split_whitespace();
    assert!(soloinputfeatureumi_l5_soloinputfeatureumi(
        &mut gene_tokens,
        SOLO_FEATURE_GENE,
        true,
        &sj_all,
        &mut iread,
        &mut cbmatch,
        &mut feature,
        &mut umi,
        &mut feat_vec,
        &mut flags,
    ));
    assert_eq!(umi, 17);
    assert_eq!(iread, 123);
    assert_eq!(flags.flag, 5);
    assert_eq!(feature, 42);
    assert_eq!(cbmatch, 1);

    let mut sj_tokens = "91 200 20 -1".split_whitespace();
    assert!(soloinputfeatureumi_l5_soloinputfeatureumi(
        &mut sj_tokens,
        SOLO_FEATURE_SJ,
        false,
        &sj_all,
        &mut iread,
        &mut cbmatch,
        &mut feature,
        &mut umi,
        &mut feat_vec,
        &mut flags,
    ));
    assert_eq!(umi, 91);
    assert_eq!(feature, 1);
    assert_eq!(cbmatch, -1);

    let mut tr_tokens = "55 2 7 8 9 10 3".split_whitespace();
    assert!(soloinputfeatureumi_l5_soloinputfeatureumi(
        &mut tr_tokens,
        SOLO_FEATURE_TRANSCRIPT3P,
        false,
        &sj_all,
        &mut iread,
        &mut cbmatch,
        &mut feature,
        &mut umi,
        &mut feat_vec,
        &mut flags,
    ));
    assert_eq!(umi, 55);
    assert_eq!(feature, 0);
    assert_eq!(feat_vec, vec![7, 8, 9, 10]);
    assert_eq!(cbmatch, 3);

    let mut empty = "".split_whitespace();
    assert!(!soloinputfeatureumi_l5_soloinputfeatureumi(
        &mut empty,
        SOLO_FEATURE_GENE,
        false,
        &sj_all,
        &mut iread,
        &mut cbmatch,
        &mut feature,
        &mut umi,
        &mut feat_vec,
        &mut flags,
    ));
}

#[test]
fn solo_feature_clear_large_clears_owned_large_vectors_but_keeps_other_state() {
    let mut sf = SoloFeature {
        feature_type: SOLO_FEATURE_VELOCYTO,
        features_number: 9,
        cb_feature_umi_map: vec![()],
        count_cell_gene_umi: vec![1],
        count_cell_gene_umi_index: vec![1],
        count_mat_mult_i: vec![1],
        count_mat_mult_m: vec![1.0],
        ind_cb_wl: vec![1],
        n_gene_per_cb: vec![2],
        n_gene_per_cb_multi: vec![3],
        n_read_per_cb: vec![4],
        n_read_per_cb_total: vec![5],
        n_read_per_cb_unique: vec![6],
        n_umi_per_cb: vec![7],
        n_umi_per_cb_multi: vec![8],
        n_umi_per_cb_sorted: vec![9],
        sj_all: [vec![10], vec![11]],
        ..Default::default()
    };
    solofeature_l29_solofeature_clearlarge(&mut sf);
    assert_eq!(sf.feature_type, SOLO_FEATURE_VELOCYTO);
    assert_eq!(sf.features_number, 9);
    assert!(sf.cb_feature_umi_map.is_empty());
    assert!(sf.count_cell_gene_umi.is_empty());
    assert!(sf.count_cell_gene_umi_index.is_empty());
    assert!(sf.count_mat_mult_i.is_empty());
    assert!(sf.count_mat_mult_m.is_empty());
    assert!(sf.ind_cb_wl.is_empty());
    assert!(sf.n_gene_per_cb.is_empty());
    assert!(sf.n_gene_per_cb_multi.is_empty());
    assert!(sf.n_read_per_cb.is_empty());
    assert!(sf.n_read_per_cb_total.is_empty());
    assert!(sf.n_read_per_cb_unique.is_empty());
    assert!(sf.n_umi_per_cb.is_empty());
    assert!(sf.n_umi_per_cb_multi.is_empty());
    assert!(sf.n_umi_per_cb_sorted.is_empty());
    assert!(sf.sj_all[0].is_empty());
    assert!(sf.sj_all[1].is_empty());
}

#[test]
fn solo_read_barcode_counters_and_status_stats_match_original_methods() {
    let mut barcode = soloreadbarcode_l4_soloreadbarcode_soloreadbarcode(1, true, 3, 4);
    assert_eq!(barcode.cb_read_count_exact, vec![0, 0, 0]);
    assert_eq!(barcode.homo_polymer, [0, 85, 170, 255]);
    assert_eq!(barcode.qual_hist[42], 0);
    assert_eq!(barcode.stats.v, vec![0; SOLO_READ_BARCODE_N_STATS]);

    let mut incoming = soloreadbarcode_l4_soloreadbarcode_soloreadbarcode(1, true, 3, 4);
    incoming.cb_read_count_exact = vec![2, 3, 5];
    incoming.qual_hist[0] = 7;
    incoming.qual_hist[42] = 11;
    soloreadbarcode_l26_soloreadbarcode_addcounts(&mut barcode, &incoming);
    assert_eq!(barcode.cb_read_count_exact, vec![2, 3, 5]);
    assert_eq!(barcode.qual_hist[0], 7);
    assert_eq!(barcode.qual_hist[42], 11);

    incoming.stats.v = (1..=SOLO_READ_BARCODE_N_STATS as u64).collect();
    soloreadbarcode_l38_soloreadbarcode_addstats(&mut barcode, &incoming);
    assert_eq!(barcode.stats.v[0], 1);
    assert_eq!(barcode.stats.v[11], 12);
    let stats_out = soloreadbarcode_l44_soloreadbarcode_statsout(&barcode);
    assert!(stats_out.contains("noNoAdapter"));
    assert!(stats_out.contains("             1\n"));

    barcode.cb_match_ind = vec![1];
    soloreadbarcode_getcbandumi_l93_soloreadbarcode_addstats(&mut barcode, 0);
    soloreadbarcode_getcbandumi_l93_soloreadbarcode_addstats(&mut barcode, 1);
    soloreadbarcode_getcbandumi_l93_soloreadbarcode_addstats(&mut barcode, 4);
    soloreadbarcode_getcbandumi_l93_soloreadbarcode_addstats(&mut barcode, -24);
    assert_eq!(barcode.cb_read_count_exact, vec![2, 4, 5]);
    assert_eq!(barcode.stats.v[9], 10 + 1);
    assert_eq!(barcode.stats.v[10], 11 + 1);
    assert_eq!(barcode.stats.v[11], 12 + 1);
    assert_eq!(barcode.stats.v[5], 6 + 1);

    let mut no_wl = SoloReadBarcode {
        cb_wl_yes: false,
        stats: barcode.stats.clone(),
        ..Default::default()
    };
    soloreadbarcode_getcbandumi_l93_soloreadbarcode_addstats(&mut no_wl, -1);
    assert_eq!(no_wl.stats, barcode.stats);
}

#[test]
fn solo_read_barcode_convert_check_umi_matches_n_and_homopolymer_rules() {
    let mut barcode = soloreadbarcode_l4_soloreadbarcode_soloreadbarcode(1, true, 3, 4);

    barcode.umi_seq = "ACGT".to_string();
    assert!(soloreadbarcode_getcbandumi_l133_soloreadbarcode_convertcheckumi(&mut barcode));
    assert_eq!(barcode.umi_b, 27);

    barcode.umi_seq = "ANNT".to_string();
    assert!(!soloreadbarcode_getcbandumi_l133_soloreadbarcode_convertcheckumi(&mut barcode));
    assert_eq!(barcode.umi_check, -23);

    barcode.umi_seq = "CCCC".to_string();
    assert!(!soloreadbarcode_getcbandumi_l133_soloreadbarcode_convertcheckumi(&mut barcode));
    assert_eq!(barcode.umi_b, barcode.homo_polymer[1]);
    assert_eq!(barcode.umi_check, -24);
}

#[test]
fn solo_read_barcode_get_cb_and_umi_handles_simple_sam_and_complex_modes() {
    let mut simple_p = Parameters {
        p_solo: ParametersSolo {
            solo_type: SOLO_TYPE_CB_UMI_SIMPLE,
            barcode_read: 0,
            cb_wl_yes: true,
            cb_wl_size: 1,
            cb_wl: vec![27],
            cb_l: 4,
            cb_s: 1,
            umi_l: 2,
            umi_s: 5,
            b_l: 6,
            cb_type_type: 1,
            cb_match_wl: CBMatchWL {
                mm1: true,
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    let mut simple_bar = soloreadbarcode_l4_soloreadbarcode_soloreadbarcode(1, true, 1, 2);
    soloreadbarcode_getcbandumi_l147_soloreadbarcode_getcbandumi(
        &mut simple_bar,
        &mut simple_p,
        &["ACGTAC".to_string()],
        &["ABCDEF".to_string()],
        &[6],
        "",
        0,
        "read-simple",
    )
    .unwrap();
    assert_eq!(simple_bar.cb_seq, "ACGT");
    assert_eq!(simple_bar.umi_seq, "AC");
    assert_eq!(simple_bar.cb_qual, "ABCD");
    assert_eq!(simple_bar.umi_b, 1);
    assert_eq!(simple_bar.cb_match, 0);
    assert_eq!(simple_bar.cb_match_ind, vec![0]);
    assert_eq!(simple_bar.cb_match_string, "0");
    assert_eq!(simple_bar.cb_read_count_exact, vec![1]);
    assert_eq!(simple_bar.stats.v[9], 1);
    assert_eq!(simple_bar.qual_hist[b'A' as usize], 1);
    assert_eq!(simple_bar.qual_hist[b'E' as usize], 1);

    let mut sam_p = Parameters {
        read_files_type_n: 10,
        p_solo: ParametersSolo {
            solo_type: SOLO_TYPE_CB_SAM_TAG_OUT,
            cb_wl_yes: false,
            cb_l: 4,
            cb_s: 1,
            umi_l: 2,
            umi_s: 5,
            b_l: 6,
            sam_attr_barcode_seq: vec!["\tCR:Z:".to_string()],
            sam_attr_barcode_qual: vec!["\tCY:Z:".to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    let mut sam_bar = soloreadbarcode_l4_soloreadbarcode_soloreadbarcode(3, false, 0, 2);
    soloreadbarcode_getcbandumi_l147_soloreadbarcode_getcbandumi(
        &mut sam_bar,
        &mut sam_p,
        &[],
        &[],
        &[],
        "CR:Z:ACGTAC\tCY:Z:ABCDEF",
        0,
        "read-sam",
    )
    .unwrap();
    assert_eq!(sam_bar.b_seq, "ACGTAC");
    assert_eq!(sam_bar.b_qual, "ABCDEF");
    assert_eq!(sam_bar.cb_seq, "ACGT");
    assert_eq!(sam_bar.cb_seq_corrected, "ACGT");
    assert_eq!(sam_bar.cb_match, 0);
    assert_eq!(sam_bar.cb_match_ind, vec![27]);

    let mut complex_p = Parameters {
        p_solo: ParametersSolo {
            solo_type: SOLO_TYPE_CB_UMI_COMPLEX,
            barcode_read: 0,
            cb_wl_yes: true,
            cb_wl_size: 1,
            umi_l: 0,
            b_l: 6,
            umi_v: SoloBarcode {
                anchor_type: [0, 0],
                anchor_dist: [0, 1],
                ..Default::default()
            },
            cb_v: vec![SoloBarcode {
                anchor_type: [0, 0],
                anchor_dist: [2, 5],
                wl: vec![vec![], vec![], vec![], vec![], vec![27]],
                wl_factor: 1,
                wl_add: vec![0, 0, 0, 0, 0],
                min_len: 4,
                ..Default::default()
            }],
            cb_match_wl: CBMatchWL {
                mm1: true,
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    let mut complex_bar = soloreadbarcode_l4_soloreadbarcode_soloreadbarcode(2, true, 1, 2);
    soloreadbarcode_getcbandumi_l147_soloreadbarcode_getcbandumi(
        &mut complex_bar,
        &mut complex_p,
        &["TAACGT".to_string()],
        &["HHHHHH".to_string()],
        &[6],
        "",
        0,
        "read-complex",
    )
    .unwrap();
    assert_eq!(complex_bar.umi_seq, "TA");
    assert_eq!(complex_p.p_solo.umi_l, 2);
    assert_eq!(complex_bar.cb_seq, "ACGT");
    assert_eq!(complex_bar.cb_qual, "HHHH");
    assert_eq!(complex_bar.cb_match, 0);
    assert_eq!(complex_bar.cb_match_ind, vec![0]);
    assert_eq!(complex_bar.cb_match_string, "0");
    assert_eq!(complex_bar.cb_read_count_exact, vec![1]);
}

#[test]
fn solo_constructor_allocates_barcode_and_configured_features() {
    let trans = Transcriptome {
        n_ge: 11,
        ..Default::default()
    };
    let none = solo_l5_solo_solo(
        &Parameters {
            p_solo: ParametersSolo {
                solo_type: SOLO_TYPE_NONE,
                ..Default::default()
            },
            ..Default::default()
        },
        &trans,
    );
    assert!(none.read_bar_sum.is_none());
    assert!(none.solo_feat.is_empty());

    let sam_tag_out = solo_l5_solo_solo(
        &Parameters {
            p_solo: ParametersSolo {
                solo_type: SOLO_TYPE_CB_SAM_TAG_OUT,
                cb_wl_yes: true,
                cb_wl_size: 3,
                umi_l: 4,
                ..Default::default()
            },
            ..Default::default()
        },
        &trans,
    );
    assert_eq!(
        sam_tag_out
            .read_bar_sum
            .as_ref()
            .unwrap()
            .cb_read_count_exact,
        vec![0; 3]
    );
    assert!(sam_tag_out.solo_feat.is_empty());

    let p = Parameters {
        run_thread_n: 2,
        sj_all: [vec![10, 20, 30], vec![5, 6, 7]],
        p_solo: ParametersSolo {
            solo_type: SOLO_TYPE_CB_UMI_SIMPLE,
            cb_wl_yes: true,
            cb_wl_size: 2,
            umi_l: 3,
            n_features: 2,
            features: vec![SOLO_FEATURE_GENE as u32, SOLO_FEATURE_SJ as u32],
            ..Default::default()
        },
        ..Default::default()
    };
    let solo = solo_l5_solo_solo(&p, &trans);
    let read_bar = solo.read_bar_sum.as_ref().unwrap();
    assert_eq!(read_bar.solo_type, SOLO_TYPE_CB_UMI_SIMPLE);
    assert_eq!(read_bar.cb_read_count_exact, vec![0; 2]);
    assert_eq!(read_bar.homo_polymer[3], 63);
    assert_eq!(solo.solo_feat.len(), 2);
    assert_eq!(solo.solo_feat[0].feature_type, SOLO_FEATURE_GENE);
    assert_eq!(solo.solo_feat[0].features_number, 11);
    assert_eq!(solo.solo_feat[0].read_feat_all_len, 2);
    assert!(solo.solo_feat[0].read_feat_sum.is_some());
    assert_eq!(solo.solo_feat[1].feature_type, SOLO_FEATURE_SJ);
    assert_eq!(solo.solo_feat[1].features_number, 3);
}

#[test]
fn solo_cell_filtering_constructor_loads_matrix_filters_and_reports_exit() {
    let trans = Transcriptome {
        n_ge: 2,
        ge_id: vec!["g1".to_string(), "g2".to_string()],
        ge_name: vec![String::new(), String::new()],
        ..Default::default()
    };
    let p = Parameters {
        run_mode_in: vec![
            "soloCellFiltering".to_string(),
            "/raw".to_string(),
            "/filtered".to_string(),
        ],
        p_solo: ParametersSolo {
            out_file_names: vec![
                "unused0".to_string(),
                "features.tsv".to_string(),
                "barcodes.tsv".to_string(),
                "matrix.mtx".to_string(),
            ],
            umi_dedup: UMIdedup {
                count_ind_main: 2,
                types: vec![1],
                ..Default::default()
            },
            cell_filter: SoloCellFilter {
                type_: vec!["TopCells".to_string()],
                top_cells: 0,
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    let matrix =
        "%%MatrixMarket matrix coordinate integer general\n%\n2 2 3\n1 1 5\n2 1 2\n2 2 1\n";

    let result =
        solo_l23_solo_solo(&p, &trans, matrix, "CB1\nCB2\n", "g1\tG1\ng2\tG2\n", "/cwd").unwrap();

    assert!(result.exited);
    assert!(
        result
            .log_stdout
            .contains("..... starting SoloCellFiltering\n")
    );
    assert!(result.log_stdout.contains("..... finished successfully\n"));
    assert!(result.log_main.ends_with("ALL DONE!\n"));
    assert_eq!(result.solo.solo_feat.len(), 1);
    assert_eq!(result.solo.solo_feat[0].feature_type, -1);
    assert_eq!(result.solo.solo_feat[0].n_cb, 1);
    assert_eq!(
        result.solo.solo_feat[0].filtered_cells.filt_vec_bool,
        vec![true]
    );
    assert_eq!(result.solo.solo_feat[0].filtered_cells.n_cells, 1);
    assert_eq!(
        result
            .cell_filtering
            .as_ref()
            .unwrap()
            .output_results
            .as_ref()
            .unwrap()
            .files["/filteredmatrix.mtx"],
        "%%MatrixMarket matrix coordinate integer general\n%\n2 1 2\n1 1 5\n2 1 2\n"
    );
}

#[test]
fn read_align_reset_n_clears_new_read_counters() {
    let mut read_align = ReadAlign {
        read_nmates: 2,
        map_marker: 1,
        n_a: 2,
        n_p: 3,
        n_w: 4,
        n_tr: 5,
        n_um: [6, 7],
        stored_lmin: 8,
        uniq_lmax: 9,
        uniq_lmax_ind: 10,
        mult_lmax: 11,
        mult_lmax_n: 12,
        mult_nmin_l: 13,
        mult_nmin: 14,
        mult_nmax: 15,
        mult_nmax_l: 16,
        chim_n: 17,
        max_score_mate: vec![100, 200, 300],
        ..Default::default()
    };

    readalign_l113_readalign_resetn(&mut read_align);

    assert_eq!(read_align.map_marker, 0);
    assert_eq!(read_align.n_a, 0);
    assert_eq!(read_align.n_p, 0);
    assert_eq!(read_align.n_w, 0);
    assert_eq!(read_align.n_tr, 0);
    assert_eq!(read_align.n_um, [0, 0]);
    assert_eq!(read_align.stored_lmin, 0);
    assert_eq!(read_align.uniq_lmax, 0);
    assert_eq!(read_align.uniq_lmax_ind, 0);
    assert_eq!(read_align.mult_lmax, 0);
    assert_eq!(read_align.mult_lmax_n, 0);
    assert_eq!(read_align.mult_nmin_l, 0);
    assert_eq!(read_align.mult_nmin, 0);
    assert_eq!(read_align.mult_nmax, 0);
    assert_eq!(read_align.mult_nmax_l, 0);
    assert_eq!(read_align.chim_n, 0);
    assert_eq!(read_align.max_score_mate, vec![0, 0, 300]);
}

#[test]
fn read_align_constructor_allocates_core_buffers_and_resets_counters() {
    let read_clip = ReadClipInput {
        n: vec![1, 0],
        n_after_ad: vec![0, 0],
        ad_seq: vec!["-".to_string(), "-".to_string()],
        ad_mmp: vec![0.0, 0.0],
    };
    let p = Parameters {
        read_nmates: 2,
        read_nends: 2,
        run_rng_seed: 777,
        quant_tr_sam_yes: true,
        align_transcripts_per_read_nmax: 4,
        align_windows_per_read_nmax: 3,
        seed_per_window_nmax: 5,
        seed_per_read_nmax: 6,
        max_nsplit: 7,
        win_bin_n: 8,
        out_sam_mult_nmax: 2,
        p_ge: ParametersGenome {
            g_type_string: "SuperTranscriptome".to_string(),
            ..Default::default()
        },
        p_clip: ParametersClip {
            adapter_type: vec!["None".to_string()],
            read_nmates: 2,
            read_nends: 2,
            in_: [read_clip.clone(), read_clip],
        },
        ..Default::default()
    };
    let genome = Genome {
        genome_out: GenomeOut {
            conv_yes: true,
            ..Default::default()
        },
        ..Default::default()
    };

    let ra = readalign_l6_readalign_readalign(&p, &genome, None, 2);

    assert_eq!(ra.read_nmates, 2);
    assert_eq!(ra.rng_mult_order_seed, 2331);
    assert_eq!(ra.align_tr_all.len(), 4);
    assert_eq!(ra.win_bin[0], vec![u32::MAX; 8]);
    assert_eq!(ra.split_r[0].len(), 7);
    assert_eq!(ra.pc.len(), 6);
    assert_eq!(ra.wa.len(), 3);
    assert_eq!(ra.wa[0].len(), 5);
    assert_eq!(ra.tr_all.len(), 4);
    assert_eq!(ra.n_win_tr, vec![0, 0, 0]);
    assert_eq!(ra.tr_array.len(), 4);
    assert_eq!(ra.aligns_gen_out_al_mult.len(), 2);
    assert_eq!(ra.read0.len(), 2);
    assert_eq!(ra.read0[0].len(), DEF_READ_SEQ_LENGTH_MAX + 1);
    assert_eq!(ra.read_name_mates[0].len(), DEF_READ_NAME_LENGTH_MAX);
    assert_eq!(ra.out_bam_one_align.len(), 4);
    assert!(ra.chunk_out_chim_junction_opened);
    assert!(ra.spl_graph_present);
    assert_eq!(ra.clip_mates.len(), 2);
    assert_eq!(ra.clip_mates[0][0].n, 1);
    assert_eq!(ra.map_marker, 0);
    assert_eq!(ra.max_score_mate, vec![0, 0]);
}

#[test]
fn read_align_chunk_constructor_wires_outputs_and_thread_files() {
    let base = std::env::temp_dir().join(format!("star-rs-readalignchunk-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&base);
    std::fs::create_dir_all(&base).unwrap();
    let read_clip = ReadClipInput {
        n: vec![0, 0],
        n_after_ad: vec![0, 0],
        ad_seq: vec!["-".to_string(), "-".to_string()],
        ad_mmp: vec![0.0, 0.0],
    };
    let p = Parameters {
        out_file_tmp: base.to_string_lossy().to_string(),
        read_nmates: 2,
        read_nends: 2,
        run_rng_seed: 11,
        quant_yes: true,
        quant_ge_count_yes: true,
        quant_tr_sam_bam_yes: true,
        chunk_in_size_bytes_array: 5,
        chunk_out_bam_size_bytes: 90,
        out_sam_bool: true,
        out_bam_unsorted: true,
        out_bam_coord: true,
        out_bam_coord_nbins: 3,
        out_bam_sort_tmp_dir: format!("{}/sort", base.to_string_lossy()),
        out_sj: true,
        limit_out_sj_collapsed: 9,
        out_filter_by_sjout_stage: 1,
        out_filter_type: "BySJout".to_string(),
        out_reads_unmapped: "Fastx".to_string(),
        wasp_yes: true,
        pe_overlap_nbases_min: 1,
        p_ch: ParametersChimeric {
            segment_min: 1,
            out_sam_old: true,
            out_junctions: true,
            ..Default::default()
        },
        align_transcripts_per_read_nmax: 2,
        align_windows_per_read_nmax: 1,
        seed_per_window_nmax: 1,
        seed_per_read_nmax: 1,
        max_nsplit: 1,
        win_bin_n: 1,
        p_clip: ParametersClip {
            adapter_type: vec!["None".to_string()],
            read_nmates: 2,
            read_nends: 2,
            in_: [read_clip.clone(), read_clip],
        },
        ..Default::default()
    };
    let tr = Transcriptome {
        n_ge: 2,
        ..Default::default()
    };

    let chunk =
        readalignchunk_l5_readalignchunk_readalignchunk(&p, &Genome::default(), Some(&tr), 4)
            .unwrap();

    assert_eq!(chunk.i_thread, 4);
    assert_eq!(chunk.ra.i_read, 0);
    assert_eq!(chunk.chunk_in, vec![vec![b'\n'; 5], vec![b'\n'; 5]]);
    assert_eq!(chunk.read_in_stream_n, 2);
    assert_eq!(chunk.chunk_out_bam.len(), 90);
    assert!(chunk.chunk_out_bam_unsorted.is_some());
    assert_eq!(chunk.chunk_out_bam_coord.n_bins, 1);
    assert_eq!(chunk.chunk_out_bam_coord.bin_total_n.len(), 3);
    assert!(chunk.chunk_out_bam_quant.is_some());
    assert_eq!(chunk.chunk_out_sj.n_store, 9);
    assert_eq!(chunk.chunk_out_sj1.n_store, 9);
    assert_eq!(
        chunk.chunk_tr.as_ref().unwrap().quants.gene_counts.g_count[0].len(),
        2
    );
    assert!(
        chunk
            .chunk_out_chim_sam_path
            .as_ref()
            .unwrap()
            .ends_with("4")
    );
    assert!(
        chunk
            .chunk_out_chim_junction_path
            .as_ref()
            .unwrap()
            .ends_with("4")
    );
    assert_eq!(chunk.chunk_out_unmapped_reads_paths.len(), 2);
    assert_eq!(chunk.chunk_out_filter_by_sjout_files.len(), 2);
    assert!(chunk.wasp_ra_present);
    assert!(chunk.pe_merge_ra_present);
    assert!(chunk.log_main.contains("Chimeric.out.sam.thread4"));
}

#[test]
fn read_align_chunk_map_chunk_flushes_sam_collapses_sj_and_adds_stats() {
    let p = Parameters {
        out_sam_bool: true,
        out_sam_order: "Paired".to_string(),
        run_thread_n: 1,
        chunk_out_bam_size_bytes: 10,
        limit_out_sam_one_read_bytes: 4,
        out_sj: true,
        limit_out_sj_one_read: 2,
        out_filter_by_sjout_stage: 1,
        ..Default::default()
    };
    let mut chunk = ReadAlignChunk {
        no_reads_left: true,
        chunk_in: vec![b"@r\nAC\n+\nII\n".to_vec()],
        chunk_in_size_bytes_total: vec![11],
        chunk_out_bam: b"abcdefghi".to_vec(),
        chunk_out_sj: OutSJ {
            n: 2,
            n_store: 3,
            junctions: vec![
                JunctionRecord {
                    start: 20,
                    gap: 5,
                    count_unique: 1,
                    ..Default::default()
                },
                JunctionRecord {
                    start: 10,
                    gap: 2,
                    count_unique: 1,
                    ..Default::default()
                },
            ],
        },
        chunk_out_sj1: OutSJ {
            n: 1,
            n_store: 3,
            junctions: vec![JunctionRecord {
                start: 30,
                gap: 4,
                count_unique: 1,
                ..Default::default()
            }],
        },
        ra: ReadAlign {
            clip_mates: vec![vec![ClipMate {
                type_: 0,
                ..Default::default()
            }]],
            stats_ra: Stats {
                read_n: 99,
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    let mut stats_all = Stats {
        time_start_map: 1,
        time_last_report: 0,
        ..Default::default()
    };
    let mut calls = 0;

    let result = readalignchunk_mapchunk_l7_readalignchunk_mapchunk(
        &mut chunk,
        &p,
        &mut stats_all,
        70,
        |ra| {
            calls += 1;
            if calls == 1 {
                ra.out_bam_bytes = 7;
                ra.stats_ra.read_n = 1;
                ra.stats_ra.read_bases = 50;
                0
            } else {
                ra.out_bam_bytes = 0;
                -1
            }
        },
        None,
    )
    .unwrap();

    assert_eq!(calls, 2);
    assert_eq!(chunk.ra.i_read, 1);
    assert_eq!(result.reads_processed, 1);
    assert_eq!(result.direct_sam_output, b"abcdefg");
    assert_eq!(chunk.chunk_out_bam_total, 0);
    assert_eq!(chunk.chunk_out_sj.n, 2);
    assert_eq!(chunk.chunk_out_sj.junctions[0].start, 10);
    assert_eq!(chunk.chunk_out_sj.n_store, 6);
    assert_eq!(chunk.chunk_out_sj1.n, 1);
    assert_eq!(stats_all.read_n, 1);
    assert_eq!(stats_all.read_bases, 50);
    assert!(
        result
            .log_main
            .contains("Increased the size of chunkOutSJ to 6")
    );
    assert!(result.progress_report.is_some());
}

#[test]
fn read_align_chunk_map_chunk_uses_paired_keep_input_order_chunk_names() {
    let p = Parameters {
        out_file_tmp: "/tmp/star-mapchunk".to_string(),
        out_sam_bool: true,
        out_sam_order: "PairedKeepInputOrder".to_string(),
        run_thread_n: 2,
        chunk_out_bam_size_bytes: 100,
        limit_out_sam_one_read_bytes: 100,
        ..Default::default()
    };
    let mut chunk = ReadAlignChunk {
        i_thread: 3,
        no_reads_left: true,
        chunk_out_bam: b"xyz".to_vec(),
        ..Default::default()
    };
    let mut stats_all = Stats::default();

    let result = readalignchunk_mapchunk_l7_readalignchunk_mapchunk(
        &mut chunk,
        &p,
        &mut stats_all,
        0,
        |ra| {
            ra.out_bam_bytes = 0;
            -1
        },
        None,
    )
    .unwrap();

    assert_eq!(
        result.paired_keep_input_order_tmp_name.as_deref(),
        Some("/tmp/star-mapchunk/Aligned.tmp.sam.chunk3")
    );
    assert_eq!(
        result.paired_keep_input_order_final_name.as_deref(),
        Some("/tmp/star-mapchunk/Aligned.out.sam.chunk3")
    );
    assert!(result.direct_sam_output.is_empty());
    assert!(result.paired_keep_input_order_tmp.is_empty());
    assert_eq!(chunk.chunk_out_bam_total, 0);
}

#[test]
fn read_align_chunk_process_chunks_reads_fastq_chunks_and_finishes_outputs() {
    let mut p = Parameters {
        run_thread_n: 1,
        read_nends: 2,
        read_nmates: 2,
        chunk_in_size_bytes: 10_000,
        out_sam_bool: true,
        chunk_out_bam_size_bytes: 5,
        limit_out_sam_one_read_bytes: 3,
        out_bam_unsorted: true,
        out_bam_coord: true,
        quant_tr_sam_bam_yes: true,
        out_reads_unmapped: "Fastx".to_string(),
        p_ch: ParametersChimeric {
            segment_min: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut chunk = ReadAlignChunk {
        i_thread: 2,
        chunk_in: vec![Vec::new(), Vec::new()],
        chunk_in_size_bytes_total: vec![0, 0],
        chunk_out_bam: b"abcde".to_vec(),
        chunk_out_bam_quant: Some(Default::default()),
        chunk_out_chim_sam_path: Some("/tmp/Chimeric.out.sam.thread2".to_string()),
        chunk_out_chim_junction_path: Some("/tmp/Chimeric.out.junction.thread2".to_string()),
        chunk_out_unmapped_reads_paths: vec![
            "/tmp/Unmapped.out.mate0.thread2".to_string(),
            "/tmp/Unmapped.out.mate1.thread2".to_string(),
        ],
        ..Default::default()
    };
    let mut thread_chunks = ThreadControl::default();
    let mut stats = Stats::default();
    let input = vec![
        "@r1 1:Y:0:ACGT extra\nACGT\n+\nIIII\n".to_string(),
        "@r1 1:Y:0:ACGT extra\nTGCA\n+\nJJJJ\n".to_string(),
    ];
    let mut calls = 0;

    let result = readalignchunk_processchunks_l11_readalignchunk_processchunks(
        &mut chunk,
        &mut p,
        &mut thread_chunks,
        &mut stats,
        0,
        &input,
        |ra| {
            calls += 1;
            if calls == 1 {
                ra.out_bam_bytes = 3;
                ra.stats_ra.read_n = 1;
                0
            } else {
                ra.out_bam_bytes = 0;
                -1
            }
        },
        None,
    )
    .unwrap();

    assert_eq!(p.i_read_all, 1);
    assert_eq!(thread_chunks.chunk_in_n, 2);
    assert_eq!(result.chunks_read, 1);
    assert_eq!(result.map_chunks.len(), 2);
    assert_eq!(result.map_chunks[0].direct_sam_output, b"abc");
    assert_eq!(
        result.chunk_inputs[0][0],
        b"@r1 1 Y 0\nACGT\n+\nIIII\n\n".to_vec()
    );
    assert_eq!(
        result.chunk_inputs[0][1],
        b"@r1 1 Y 0\nTGCA\n+\nJJJJ\n\n".to_vec()
    );
    assert!(result.flushed_bam_unsorted);
    assert!(result.flushed_bam_coord);
    assert!(result.flushed_bam_quant);
    assert_eq!(
        result.chim_sam_cat_path.as_deref(),
        Some("/tmp/Chimeric.out.sam.thread2")
    );
    assert_eq!(result.unmapped_fastx_cat_paths.len(), 2);
    assert!(result.log_main.contains("Completed: thread #2"));
}

#[test]
fn read_align_chunk_process_chunks_converts_paired_sam_to_fastq_shape() {
    let mut p = Parameters {
        read_files_type_n: 10,
        read_nends: 2,
        read_nmates: 2,
        chunk_in_size_bytes: 10_000,
        out_sam_read_id: "Number".to_string(),
        read_files_index: 7,
        ..Default::default()
    };
    let mut chunk = ReadAlignChunk {
        chunk_in: vec![Vec::new(), Vec::new()],
        chunk_in_size_bytes_total: vec![0, 0],
        ..Default::default()
    };
    let mut thread_chunks = ThreadControl::default();
    let mut stats = Stats::default();
    let input = vec![
        concat!(
            "pairA 99 chr1 1 255 4M = 10 0 ACGT IIII NM:i:0\n",
            "pairA 147 chr1 10 255 4M = 1 0 AGTC JJJJ AS:i:4\n"
        )
        .to_string(),
        String::new(),
    ];

    let result = readalignchunk_processchunks_l11_readalignchunk_processchunks(
        &mut chunk,
        &mut p,
        &mut thread_chunks,
        &mut stats,
        0,
        &input,
        |_ra| -1,
        None,
    )
    .unwrap();

    assert_eq!(p.i_read_all, 1);
    assert_eq!(
        result.chunk_inputs[0][0],
        b"@1 1 N 7 NM:i:0\nACGT\n+\nIIII\n\n".to_vec()
    );
    assert_eq!(
        result.chunk_inputs[0][1],
        b"@1 1 N 7 AS:i:4\nGACT\n+\nJJJJ\n\n".to_vec()
    );
}

#[test]
fn read_align_chimeric_detection_old_accepts_unique_mate_bracketed_chimera() {
    let read_length = vec![10, 10];
    let left = Transcript {
        exons: vec![[0, 100, 6, 0, 0]],
        n_exons: 1,
        l_read: 21,
        read_length: read_length.clone(),
        r_length: 6,
        max_score: 8,
        chr: 0,
        str_: 0,
        ro_start: 0,
        ..Default::default()
    };
    let right = Transcript {
        exons: vec![[11, 200, 6, 1, 0]],
        n_exons: 1,
        l_read: 21,
        read_length: read_length.clone(),
        r_length: 6,
        max_score: 9,
        chr: 1,
        str_: 0,
        ro_start: 11,
        ..Default::default()
    };
    let p = Parameters {
        p_ch: ParametersChimeric {
            segment_min: 3,
            segment_read_gap_max: 0,
            junction_overhang_min: 2,
            score_min: 10,
            score_drop_max: 100,
            score_separation: 0,
            main_segment_mult_nmax: 10,
            ..Default::default()
        },
        align_mates_gap_max: 100,
        ..Default::default()
    };
    let map_gen = Genome {
        chr_start: vec![0, 0],
        g: vec![0; 256],
        ..Default::default()
    };
    let mut ra = ReadAlign {
        tr_best: left.clone(),
        n_tr: 1,
        n_w: 2,
        n_win_tr: vec![1, 1],
        tr_all: vec![vec![left], vec![right]],
        l_read: 21,
        read_length,
        tr_chim: vec![Transcript::default(); 2],
        ..Default::default()
    };

    assert!(
        readalign_chimericdetectionold_l7_readalign_chimericdetectionold(&mut ra, &p, &map_gen)
    );

    assert_eq!(ra.chim_n, 2);
    assert_eq!(ra.chim_motif, -1);
    assert_eq!(ra.chim_j0, 106);
    assert_eq!(ra.chim_j1, 199);
    assert_eq!(ra.tr_chim[0].chr, 0);
    assert_eq!(ra.tr_chim[1].chr, 1);
}

#[test]
fn read_align_chimeric_detection_old_rejects_nonunique_and_near_linear_pairs() {
    let read_length = vec![10, 10];
    let tr = |r_start, g_start, ifrag, score| Transcript {
        exons: vec![[r_start, g_start, 6, ifrag, 0]],
        n_exons: 1,
        l_read: 21,
        read_length: read_length.clone(),
        r_length: 6,
        max_score: score,
        chr: 0,
        str_: 0,
        ro_start: r_start,
        ..Default::default()
    };
    let p = Parameters {
        p_ch: ParametersChimeric {
            segment_min: 3,
            segment_read_gap_max: 0,
            junction_overhang_min: 2,
            score_min: 10,
            score_drop_max: 100,
            score_separation: 2,
            main_segment_mult_nmax: 10,
            ..Default::default()
        },
        align_mates_gap_max: 100,
        ..Default::default()
    };
    let map_gen = Genome {
        chr_start: vec![0],
        g: vec![0; 256],
        ..Default::default()
    };
    let left = tr(0, 100, 0, 8);
    let right = tr(11, 120, 1, 9);
    let alt_right = tr(11, 140, 1, 8);
    let mut nonunique = ReadAlign {
        tr_best: left.clone(),
        n_tr: 1,
        n_w: 2,
        n_win_tr: vec![1, 2],
        tr_all: vec![vec![left.clone()], vec![right.clone(), alt_right]],
        l_read: 21,
        read_length: read_length.clone(),
        tr_chim: vec![Transcript::default(); 2],
        ..Default::default()
    };
    assert!(
        !readalign_chimericdetectionold_l7_readalign_chimericdetectionold(
            &mut nonunique,
            &p,
            &map_gen
        )
    );

    let mut near_linear = ReadAlign {
        tr_best: left.clone(),
        n_tr: 1,
        n_w: 2,
        n_win_tr: vec![1, 1],
        tr_all: vec![vec![left], vec![right]],
        l_read: 21,
        read_length,
        tr_chim: vec![Transcript::default(); 2],
        ..Default::default()
    };
    assert!(
        !readalign_chimericdetectionold_l7_readalign_chimericdetectionold(
            &mut near_linear,
            &p,
            &map_gen
        )
    );
    assert_eq!(near_linear.chim_motif, -1);
}

#[test]
fn solo_feature_process_records_orchestrates_standard_counting_and_outputs() {
    let mut p = Parameters {
        out_file_name_prefix: "/tmp/star-solo/".to_string(),
        run_thread_n: 1,
        read_nmates: 1,
        read_quality_score_base: 33,
        ..Default::default()
    };
    let mut p_solo = ParametersSolo {
        solo_type: SOLO_TYPE_CB_UMI_SIMPLE,
        cb_wl_yes: true,
        cb_wl_size: 1,
        cb_wl_str: vec!["CELL1".to_string()],
        out_file_names: vec![
            "Solo.out/".to_string(),
            "features.tsv".to_string(),
            "barcodes.tsv".to_string(),
            "matrix.mtx".to_string(),
        ],
        out_format_features_gene_field3: "Gene Expression".to_string(),
        umi_dedup: UMIdedup {
            types: vec![0],
            ..Default::default()
        },
        cell_filter: SoloCellFilter {
            type_: vec!["None".to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    let stats_names: Vec<String> = (0..SOLO_READ_FEATURE_N_STATS)
        .map(|ii| format!("stat{ii}"))
        .collect();
    let mut stats_v = vec![0; SOLO_READ_FEATURE_N_STATS];
    stats_v[SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE] = 10;
    stats_v[SOLO_READ_FEATURE_STAT_YES_UMIS] = 5;
    let mut solo_feature = SoloFeature {
        feature_type: SOLO_FEATURE_GENE,
        features_number: 2,
        read_bar_sum: Some(SoloReadBarcode {
            qual_hist: {
                let mut q = vec![0; 256];
                q[70] = 10;
                q
            },
            stats: star_rs::generated::structs::SoloReadBarcodeStats {
                names: (0..12).map(|ii| format!("bc{ii}")).collect(),
                v: vec![0; 12],
            },
            ..Default::default()
        }),
        read_feat_sum: Some(SoloReadFeature {
            cb_read_count: vec![1],
            stats: star_rs::generated::structs::SoloReadFeatureStats {
                names: stats_names,
                v: stats_v,
            },
            ..Default::default()
        }),
        n_cb: 1,
        ind_cb: vec![0],
        n_gene_per_cb: vec![1],
        count_cell_gene_umi_index: vec![0],
        count_cell_gene_umi: vec![0, 5],
        count_mat_stride: 2,
        n_umi_per_cb: vec![5],
        n_read_per_cb_unique: vec![7],
        ..Default::default()
    };
    let mut read_bar_sum = solo_feature.read_bar_sum.clone().unwrap();
    let trans = Transcriptome {
        n_ge: 2,
        ge_id: vec!["g0".to_string(), "g1".to_string()],
        ge_name: vec!["G0".to_string(), "G1".to_string()],
        ..Default::default()
    };
    let g_stats = Stats {
        read_n: 20,
        mapped_reads_u: 12,
        mapped_reads_m: 3,
        ..Default::default()
    };
    let ra_chunks = vec![ReadAlign {
        qual_hist: {
            let mut h = vec![vec![0; 256]];
            h[0][70] = 20;
            h
        },
        ..Default::default()
    }];

    let result = solofeature_processrecords_l8_solofeature_processrecords(
        &mut solo_feature,
        &mut p,
        &mut p_solo,
        &mut read_bar_sum,
        &[],
        None,
        &trans,
        &g_stats,
        &ra_chunks,
        "/tmp",
        &[],
        "Apr 28 12:00:00",
        "Apr 28 12:00:01",
        "Apr 28 12:00:02",
        "Apr 28 12:00:03",
        "Apr 28 12:00:04",
        "Apr 28 12:00:05",
        "RAM",
        |sf| {
            sf.n_read_per_cb = vec![9];
            "counted\n".to_string()
        },
        |_sf| "quant\n".to_string(),
    )
    .unwrap();

    assert!(result.count_cb_gene_umi_called);
    assert!(!result.quant_transcript_called);
    assert_eq!(solo_feature.output_prefix, "/tmp/star-solo/Solo.out/Gene/");
    assert!(
        result
            .created_directories
            .contains(&"/tmp/star-solo/Solo.out/Gene/".to_string())
    );
    assert_eq!(
        result.files["/tmp/star-solo/Solo.out/Gene//raw/matrix.mtx"],
        "%%MatrixMarket matrix coordinate integer general\n%\n2 1 1\n1 1 5\n"
    );
    assert!(result.files["/tmp/star-solo/Solo.out/Gene/Features.stats"].contains("stat8"));
    assert!(
        result
            .files
            .contains_key("/tmp/star-solo/Solo.out/Gene/Summary.csv")
    );
    assert!(result.log_main.contains("Starting Solo post-map for Gene"));
    assert!(result.log_main.contains("counted"));
    assert!(result.log_main.ends_with("RAM"));
    assert!(solo_feature.count_cell_gene_umi.is_empty());
}

#[test]
fn solo_feature_process_records_loads_sj_rows_and_returns_after_quant_transcript() {
    let mut p = Parameters {
        out_file_name_prefix: "out/".to_string(),
        ..Default::default()
    };
    let mut p_solo = ParametersSolo {
        solo_type: SOLO_TYPE_CB_UMI_SIMPLE,
        cb_wl_yes: true,
        cb_wl_size: 0,
        out_file_names: vec![
            "Solo.out/".to_string(),
            "features.tsv".to_string(),
            "barcodes.tsv".to_string(),
            "matrix.mtx".to_string(),
        ],
        umi_dedup: UMIdedup {
            types: vec![0],
            ..Default::default()
        },
        cell_filter: SoloCellFilter {
            type_: vec!["None".to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    let mut solo_feature = SoloFeature {
        feature_type: SOLO_FEATURE_TRANSCRIPT3P,
        read_feat_sum: Some(SoloReadFeature::default()),
        ..Default::default()
    };
    let mut read_bar_sum = SoloReadBarcode::default();

    let result = solofeature_processrecords_l8_solofeature_processrecords(
        &mut solo_feature,
        &mut p,
        &mut p_solo,
        &mut read_bar_sum,
        &[],
        None,
        &Transcriptome::default(),
        &Stats {
            read_n: 1,
            ..Default::default()
        },
        &[],
        ".",
        &[(100, 20)],
        "t0",
        "t1",
        "t2",
        "t3",
        "t4",
        "t5",
        "mem",
        |_sf| "count\n".to_string(),
        |sf| {
            sf.n_reads_mapped = 42;
            "quant\n".to_string()
        },
    )
    .unwrap();

    assert!(result.returned_after_quant_transcript);
    assert!(result.quant_transcript_called);
    assert!(!result.count_cb_gene_umi_called);
    assert_eq!(solo_feature.n_reads_mapped, 42);
    assert!(
        !result
            .files
            .contains_key("out/Solo.out/Transcript3p/Features.stats")
    );

    solo_feature.feature_type = SOLO_FEATURE_SJ;
    solo_feature.read_bar_sum = Some(SoloReadBarcode {
        qual_hist: {
            let mut q = vec![0; 256];
            q[70] = 1;
            q
        },
        stats: star_rs::generated::structs::SoloReadBarcodeStats {
            names: (0..12).map(|ii| format!("bc{ii}")).collect(),
            v: vec![0; 12],
        },
        ..Default::default()
    });
    let mut stats_v = vec![0; SOLO_READ_FEATURE_N_STATS];
    stats_v[SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE] = 1;
    stats_v[SOLO_READ_FEATURE_STAT_YES_UMIS] = 1;
    solo_feature.read_feat_sum = Some(SoloReadFeature {
        stats: star_rs::generated::structs::SoloReadFeatureStats {
            names: (0..SOLO_READ_FEATURE_N_STATS)
                .map(|ii| format!("stat{ii}"))
                .collect(),
            v: stats_v,
        },
        ..Default::default()
    });
    read_bar_sum = solo_feature.read_bar_sum.clone().unwrap();
    let sj_result = solofeature_processrecords_l8_solofeature_processrecords(
        &mut solo_feature,
        &mut p,
        &mut p_solo,
        &mut read_bar_sum,
        &[],
        None,
        &Transcriptome::default(),
        &Stats {
            read_n: 1,
            ..Default::default()
        },
        &[ReadAlign {
            qual_hist: {
                let mut h = vec![vec![0; 256]];
                h[0][70] = 1;
                h
            },
            ..Default::default()
        }],
        ".",
        &[(100, 20)],
        "t0",
        "t1",
        "t2",
        "t3",
        "t4",
        "t5",
        "mem",
        |sf| {
            sf.count_mat_stride = 1;
            "count\n".to_string()
        },
        |_sf| "quant\n".to_string(),
    )
    .unwrap();
    assert_eq!(p.sj_all[0], vec![100]);
    assert!(
        sj_result
            .log_main
            .contains("Read splice junctions for Solo SJ feature: 1")
    );
}

#[test]
fn solo_process_and_output_collects_barcodes_and_dispatches_features() {
    let mut p = Parameters {
        out_file_name_prefix: "/tmp/star-solo/".to_string(),
        run_thread_n: 1,
        read_nmates: 1,
        read_quality_score_base: 33,
        ..Default::default()
    };
    let mut solo = Solo {
        p_solo: ParametersSolo {
            solo_type: SOLO_TYPE_CB_UMI_SIMPLE,
            cb_wl_yes: true,
            cb_wl_size: 1,
            cb_wl_str: vec!["CELL1".to_string()],
            out_file_names: vec![
                "Solo.out/".to_string(),
                "features.tsv".to_string(),
                "barcodes.tsv".to_string(),
                "matrix.mtx".to_string(),
            ],
            out_format_features_gene_field3: "Gene Expression".to_string(),
            n_features: 1,
            features: vec![SOLO_FEATURE_GENE as u32],
            umi_dedup: UMIdedup {
                types: vec![0],
                ..Default::default()
            },
            cell_filter: SoloCellFilter {
                type_: vec!["None".to_string()],
                ..Default::default()
            },
            cb_match_wl: CBMatchWL {
                mm1_multi_pc: true,
                ..Default::default()
            },
            ..Default::default()
        },
        read_bar_sum: Some(SoloReadBarcode {
            cb_wl_yes: true,
            cb_wl_size: 1,
            cb_read_count_exact: vec![2],
            qual_hist: vec![0; 256],
            stats: star_rs::generated::structs::SoloReadBarcodeStats {
                names: (0..12).map(|ii| format!("bc{ii}")).collect(),
                v: vec![0; 12],
            },
            ..Default::default()
        }),
        solo_feat: vec![SoloFeature {
            feature_type: SOLO_FEATURE_GENE,
            features_number: 1,
            read_bar_sum: Some(SoloReadBarcode {
                stats: star_rs::generated::structs::SoloReadBarcodeStats {
                    names: (0..12).map(|ii| format!("bc{ii}")).collect(),
                    v: vec![0; 12],
                },
                qual_hist: {
                    let mut q = vec![0; 256];
                    q[70] = 1;
                    q
                },
                ..Default::default()
            }),
            read_feat_sum: Some(SoloReadFeature {
                cb_read_count: vec![1],
                stats: star_rs::generated::structs::SoloReadFeatureStats {
                    names: (0..SOLO_READ_FEATURE_N_STATS)
                        .map(|ii| format!("stat{ii}"))
                        .collect(),
                    v: {
                        let mut v = vec![0; SOLO_READ_FEATURE_N_STATS];
                        v[SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE] = 1;
                        v[SOLO_READ_FEATURE_STAT_YES_UMIS] = 1;
                        v
                    },
                },
                ..Default::default()
            }),
            n_cb: 1,
            ind_cb: vec![0],
            n_gene_per_cb: vec![1],
            count_cell_gene_umi_index: vec![0],
            count_cell_gene_umi: vec![0, 3],
            count_mat_stride: 2,
            n_umi_per_cb: vec![3],
            n_read_per_cb_unique: vec![4],
            ..Default::default()
        }],
    };
    let mut ra_chunks = vec![ReadAlignChunk {
        ra: ReadAlign {
            solo_read: SoloRead {
                read_bar: Some(SoloReadBarcode {
                    cb_wl_yes: true,
                    cb_wl_size: 1,
                    cb_read_count_exact: vec![5],
                    qual_hist: {
                        let mut q = vec![0; 256];
                        q[70] = 2;
                        q
                    },
                    stats: star_rs::generated::structs::SoloReadBarcodeStats {
                        names: (0..12).map(|ii| format!("bc{ii}")).collect(),
                        v: vec![1; 12],
                    },
                    ..Default::default()
                }),
                read_feat: vec![SoloReadFeature {
                    feature_type: SOLO_FEATURE_GENE,
                    ..Default::default()
                }],
                ..Default::default()
            },
            qual_hist: {
                let mut h = vec![vec![0; 256]];
                h[0][70] = 5;
                h
            },
            ..Default::default()
        },
        ..Default::default()
    }];
    let trans = Transcriptome {
        n_ge: 1,
        ge_id: vec!["g0".to_string()],
        ge_name: vec!["G0".to_string()],
        ..Default::default()
    };
    let g_stats = Stats {
        read_n: 10,
        mapped_reads_u: 5,
        mapped_reads_m: 1,
        ..Default::default()
    };

    let result = solo_l48_solo_processandoutput(
        &mut solo,
        &mut p,
        &trans,
        &mut ra_chunks,
        &g_stats,
        "/tmp",
        &[],
        "t-start",
        "t-finish",
        "t-process",
        "t-raw",
        "t-filter",
        "t-redist",
        "t-collapse",
        "t-count",
        "mem",
        |_ii, sf| {
            sf.n_read_per_cb = vec![4];
            "feature-count\n".to_string()
        },
        |_ii, _sf| "quant\n".to_string(),
    )
    .unwrap();

    let read_bar = solo.read_bar_sum.as_ref().unwrap();
    assert_eq!(read_bar.cb_read_count_exact, vec![8]);
    assert!(ra_chunks[0].ra.solo_read.read_bar.is_none());
    assert!(result.files["/tmp/star-solo/Solo.out/Barcodes.stats"].contains("bc0"));
    assert_eq!(result.feature_results.len(), 1);
    assert!(result.feature_results[0].count_cb_gene_umi_called);
    assert!(
        result
            .log_stdout
            .contains("t-start ..... started Solo counting")
    );
    assert!(
        result
            .log_stdout
            .contains("t-finish ..... finished Solo counting")
    );
    assert!(result.log_main.contains("feature-count"));
}

#[test]
fn solo_process_and_output_returns_after_barcode_stats_for_cb_sam_tag_out() {
    let mut p = Parameters {
        out_file_name_prefix: "out/".to_string(),
        run_thread_n: 1,
        ..Default::default()
    };
    let mut solo = Solo {
        p_solo: ParametersSolo {
            solo_type: SOLO_TYPE_CB_SAM_TAG_OUT,
            cb_wl_yes: false,
            out_file_names: vec!["Solo.out/".to_string()],
            ..Default::default()
        },
        read_bar_sum: Some(SoloReadBarcode {
            stats: star_rs::generated::structs::SoloReadBarcodeStats {
                names: vec!["ok".to_string()],
                v: vec![7],
            },
            ..Default::default()
        }),
        ..Default::default()
    };

    let result = solo_l48_solo_processandoutput(
        &mut solo,
        &mut p,
        &Transcriptome::default(),
        &mut [],
        &Stats::default(),
        ".",
        &[],
        "start",
        "finish",
        "process",
        "raw",
        "filter",
        "redist",
        "collapse",
        "count",
        "mem",
        |_ii, _sf| "count\n".to_string(),
        |_ii, _sf| "quant\n".to_string(),
    )
    .unwrap();

    assert!(result.returned_after_barcode_output);
    assert!(result.feature_results.is_empty());
    assert_eq!(
        result.files["out/Solo.out/Barcodes.stats"],
        format!("{:>50}{:>15}\n", "ok", 7)
    );
}

#[test]
fn solo_feature_quant_transcript_builds_cluster_expression_matrix() {
    let mut counts = vec![0_u32; 1105];
    for ii in 0..1000 {
        counts[ii] = 10;
    }
    counts[1000] = 20;
    counts[1001] = 10;
    counts[1002] = 5;
    counts[1003] = 15;
    let mut solo_feature = SoloFeature {
        output_prefix: "/tmp/Transcript3p/".to_string(),
        read_feat_sum: Some(SoloReadFeature {
            transcript_dist_count: counts,
            ..Default::default()
        }),
        read_feat_all: vec![SoloReadFeature {
            stream_reads: ["0 7 2 0 10 1 10", "0 7 2 0 10 1 2000", "1 9 1 1 10"].join("\n"),
            ..Default::default()
        }],
        ..Default::default()
    };
    let p_solo = ParametersSolo {
        cluster_cb_file: "clusters.tsv".to_string(),
        cb_wl: vec![0, 1],
        out_file_names: vec![
            "Solo.out/".to_string(),
            "features.tsv".to_string(),
            "barcodes.tsv".to_string(),
            "matrix.mtx".to_string(),
        ],
        ..Default::default()
    };
    let trans = Transcriptome {
        n_tr: 2,
        n_ge: 2,
        tr_id: vec!["tr0".to_string(), "tr1".to_string()],
        tr_len: vec![20, 20],
        tr_gene: vec![0, 1],
        ge_name: vec!["gene0".to_string(), "gene1".to_string()],
        ..Default::default()
    };

    let result = solofeature_quanttranscript_l12_solofeature_quanttranscript(
        &mut solo_feature,
        &p_solo,
        &trans,
        1,
        "AAAAAAAAAAAAAAAA 1\nAAAAAAAAAAAAAAAC 2\nNAAA 9\nTTTTTTTTTTTTTTTT 3\n",
        "t-input",
        "t-cluster",
        "t-done",
    );

    assert!(!result.returned_no_cluster_file);
    assert!(
        result
            .log_main
            .contains("WARNING: cluster CB sequence contains non-ACGT base and is ignored: NAAA")
    );
    assert!(
        result
            .log_main
            .contains("t-input ... Transcript3p counting: finished input")
    );
    assert!(
        result
            .files
            .contains_key("/tmp/Transcript3p/transcriptEndDistanceDistribution.txt")
    );
    assert_eq!(
        result.files["/tmp/Transcript3p//features.tsv"],
        "tr0\t20\tgene0\ntr1\t20\tgene1\n"
    );
    let matrix = &result.files["/tmp/Transcript3p/matrix.mtx"];
    assert!(matrix.starts_with("%%MatrixMarket matrix coordinate real general\n%\n2 "));
    assert!(matrix.contains("\n1 1 "));
    assert!(matrix.contains("\n2 2 "));
    assert!(result.cluster_expression[&1][0] > 0.0);
    assert_eq!(result.cluster_expression[&1][1], 0.0);
    assert_eq!(result.cluster_expression[&2][0], 0.0);
    assert!(result.cluster_expression[&2][1] > 0.0);
}

#[test]
fn solo_feature_quant_transcript_returns_without_cluster_file() {
    let mut solo_feature = SoloFeature::default();
    let result = solofeature_quanttranscript_l12_solofeature_quanttranscript(
        &mut solo_feature,
        &ParametersSolo {
            cluster_cb_file: "-".to_string(),
            ..Default::default()
        },
        &Transcriptome::default(),
        0,
        "",
        "input",
        "cluster",
        "done",
    );
    assert!(result.returned_no_cluster_file);
    assert!(result.files.is_empty());
}

#[test]
fn read_align_store_aligns_inserts_sorted_updates_counters_and_skips_duplicates() {
    let p = Parameters {
        seed_multimap_nmax: 10,
        seed_per_read_nmax: 4,
        ..Default::default()
    };
    let mut read_align = ReadAlign {
        stored_lmin: 5,
        ..Default::default()
    };

    readalign_storealigns_l10_readalign_storealigns(&mut read_align, &p, 0, 10, 2, 4, [20, 22], 1)
        .unwrap();
    assert_eq!(read_align.n_p, 1);
    assert_eq!(read_align.n_a, 2);
    assert_eq!(read_align.n_um, [0, 2]);
    assert_eq!(read_align.pc[0], [10, 4, 0, 0, 2, 20, 22, 1]);
    assert_eq!(read_align.mult_nmin, 2);
    assert_eq!(read_align.mult_nmin_l, 5);
    assert_eq!(read_align.mult_lmax, 5);
    assert_eq!(read_align.mult_lmax_n, 2);
    assert_eq!(read_align.mult_nmax, 2);
    assert_eq!(read_align.mult_nmax_l, 5);

    readalign_storealigns_l10_readalign_storealigns(&mut read_align, &p, 0, 8, 1, 6, [30, 30], 0)
        .unwrap();
    assert_eq!(read_align.n_p, 2);
    assert_eq!(read_align.pc[0][PC_R_START], 8);
    assert_eq!(read_align.pc[1][PC_R_START], 10);
    assert_eq!(read_align.n_a, 3);
    assert_eq!(read_align.n_um, [1, 2]);
    assert_eq!(read_align.uniq_lmax, 6);
    assert_eq!(read_align.uniq_lmax_ind, 1);

    readalign_storealigns_l10_readalign_storealigns(&mut read_align, &p, 1, 15, 1, 6, [40, 41], 2)
        .unwrap();
    assert_eq!(read_align.pc[1][PC_R_START], 10);
    assert_eq!(read_align.pc[1][PC_DIR], 1);
    assert_eq!(read_align.pc[2][PC_LENGTH], 4);

    readalign_storealigns_l10_readalign_storealigns(&mut read_align, &p, 0, 8, 1, 6, [99, 99], 9)
        .unwrap();
    assert_eq!(read_align.n_p, 3);
    assert_eq!(read_align.n_a, 5);
    assert_eq!(read_align.n_um, [3, 2]);
    assert_eq!(read_align.pc[0][PC_SASTART], 30);

    readalign_storealigns_l10_readalign_storealigns(&mut read_align, &p, 0, 7, 25, 9, [0, 0], 0)
        .unwrap();
    assert_eq!(read_align.n_p, 3);
    assert_eq!(read_align.mult_nmin, 2);
    assert_eq!(read_align.n_a, 5);

    let mut too_many = ReadAlign::default();
    let p_limit = Parameters {
        seed_multimap_nmax: 10,
        seed_per_read_nmax: 0,
        ..Default::default()
    };
    assert!(
        readalign_storealigns_l10_readalign_storealigns(
            &mut too_many,
            &p_limit,
            0,
            1,
            1,
            1,
            [0, 0],
            0,
        )
        .unwrap_err()
        .contains("seedPerReadNmax")
    );
}

#[test]
fn read_align_copy_read_copies_only_read_level_state_and_lread_bytes() {
    let source = ReadAlign {
        l_read: 4,
        read_length: vec![4, 3, 99],
        read_length_original: vec![5, 6, 77],
        read_length_pair_original: 11,
        out_filter_mismatch_nmax_total: 8,
        read_name: "readA".to_string(),
        i_read_all: 123,
        read_filter: 2,
        read_files_index: 7,
        read1: [
            vec![1, 2, 3, 4, 99],
            vec![5, 6, 7, 8, 99],
            vec![9, 10, 11, 12, 99],
        ],
        n_a: 42,
        ..Default::default()
    };
    let mut target = ReadAlign {
        read_length: vec![0],
        read_length_original: vec![0],
        read1: [vec![100], vec![101], vec![102]],
        n_a: 5,
        ..Default::default()
    };

    readalign_waspmap_l115_readalign_copyread(&mut target, &source);

    assert_eq!(target.l_read, 4);
    assert_eq!(target.read_length, vec![4, 3]);
    assert_eq!(target.read_length_original, vec![5, 6]);
    assert_eq!(target.read_length_pair_original, 11);
    assert_eq!(target.out_filter_mismatch_nmax_total, 8);
    assert_eq!(target.read_name, "readA");
    assert_eq!(target.i_read_all, 123);
    assert_eq!(target.read_filter, 2);
    assert_eq!(target.read_files_index, 7);
    assert_eq!(
        target.read1,
        [vec![1, 2, 3, 4], vec![5, 6, 7, 8], vec![9, 10, 11, 12]]
    );
    assert_eq!(target.n_a, 5);
}

#[test]
fn read_align_wasp_map_classifies_early_exit_cases() {
    let mut ra = ReadAlign {
        tr_best: Transcript {
            var_allele: vec![1],
            ..Default::default()
        },
        ..Default::default()
    };
    let mut wasp_ra = ReadAlign::default();
    let disabled = readalign_waspmap_l3_readalign_waspmap(
        &mut ra,
        &mut wasp_ra,
        &Parameters::default(),
        &Genome::default(),
        None,
        None,
        &[],
    );
    assert_eq!(disabled.wasp_type, -1);
    assert_eq!(ra.wasp_type, -1);

    let p = Parameters {
        wasp_yes: true,
        ..Default::default()
    };
    let mut no_variant = ReadAlign {
        tr_best: Transcript::default(),
        ..Default::default()
    };
    let no_variant_result = readalign_waspmap_l3_readalign_waspmap(
        &mut no_variant,
        &mut wasp_ra,
        &p,
        &Genome::default(),
        None,
        None,
        &[],
    );
    assert_eq!(no_variant_result.wasp_type, -1);

    let mut multi = ReadAlign {
        n_tr: 2,
        tr_best: Transcript {
            var_allele: vec![1],
            ..Default::default()
        },
        ..Default::default()
    };
    let multi_result = readalign_waspmap_l3_readalign_waspmap(
        &mut multi,
        &mut wasp_ra,
        &p,
        &Genome::default(),
        None,
        None,
        &[],
    );
    assert_eq!(multi_result.wasp_type, 2);

    let mut n_base = ReadAlign {
        n_tr: 1,
        tr_best: Transcript {
            var_allele: vec![4],
            ..Default::default()
        },
        ..Default::default()
    };
    let n_base_result = readalign_waspmap_l3_readalign_waspmap(
        &mut n_base,
        &mut wasp_ra,
        &p,
        &Genome::default(),
        None,
        None,
        &[],
    );
    assert_eq!(n_base_result.wasp_type, 3);
}

#[test]
fn read_align_wasp_map_generates_remap_requests_and_accepts_same_alignment() {
    let p = Parameters {
        wasp_yes: true,
        ..Default::default()
    };
    let align = Transcript {
        n_exons: 1,
        exons: vec![[2, 100, 3, 0, 0]],
        var_ind: vec![0, 1],
        var_read_coord: vec![1, 3],
        var_allele: vec![1, 2],
        ..Default::default()
    };
    let mut ra = ReadAlign {
        n_tr: 1,
        tr_best: align.clone(),
        l_read: 5,
        read_length: vec![5, 0],
        read_length_original: vec![5, 0],
        read1: [
            vec![0, 1, 2, 3, 0],
            vec![3, 2, 1, 0, 3],
            vec![3, 0, 1, 2, 3],
        ],
        ..Default::default()
    };
    let mut wasp_ra = ReadAlign::default();
    let genome = Genome {
        var: Variation {
            yes: true,
            snp: SNP {
                nt: vec![[0, 1, 2], [0, 2, 3]],
                ..Default::default()
            },
        },
        ..Default::default()
    };
    let outcomes = vec![
        WaspMapOutcome {
            unmap_type: -1,
            n_tr: 1,
            align: align.clone(),
            ..Default::default()
        },
        WaspMapOutcome {
            unmap_type: -1,
            n_tr: 1,
            align: align.clone(),
            ..Default::default()
        },
        WaspMapOutcome {
            unmap_type: -1,
            n_tr: 1,
            align: align.clone(),
            ..Default::default()
        },
    ];

    let result = readalign_waspmap_l3_readalign_waspmap(
        &mut ra,
        &mut wasp_ra,
        &p,
        &genome,
        None,
        None,
        &outcomes,
    );

    assert_eq!(result.wasp_type, 1);
    assert_eq!(ra.wasp_type, 1);
    assert_eq!(
        result
            .requests
            .iter()
            .map(|r| r.alleles.clone())
            .collect::<Vec<_>>(),
        vec![vec![1, 1], vec![2, 1], vec![2, 2]]
    );
    assert_eq!(result.requests[0].read1[0], vec![0, 1, 2, 2, 0]);
    assert_eq!(result.requests[1].read1[0], vec![0, 2, 2, 2, 0]);
    assert_eq!(result.requests[2].read1[0], vec![0, 2, 2, 3, 0]);
}

#[test]
fn read_align_wasp_map_rejects_changed_remap_alignment() {
    let p = Parameters {
        wasp_yes: true,
        ..Default::default()
    };
    let align = Transcript {
        n_exons: 1,
        exons: vec![[2, 100, 3, 0, 0]],
        var_ind: vec![0],
        var_read_coord: vec![1],
        var_allele: vec![1],
        ..Default::default()
    };
    let changed = Transcript {
        n_exons: 1,
        exons: vec![[3, 100, 3, 0, 0]],
        ..Default::default()
    };
    let mut ra = ReadAlign {
        n_tr: 1,
        tr_best: align,
        l_read: 4,
        read_length: vec![4, 0],
        read_length_original: vec![4, 0],
        read1: [vec![0, 1, 2, 3], vec![3, 2, 1, 0], vec![0, 1, 2, 3]],
        ..Default::default()
    };
    let mut wasp_ra = ReadAlign::default();
    let genome = Genome {
        var: Variation {
            snp: SNP {
                nt: vec![[0, 1, 2]],
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    let outcomes = vec![WaspMapOutcome {
        unmap_type: -1,
        n_tr: 1,
        align: changed,
        ..Default::default()
    }];

    let result = readalign_waspmap_l3_readalign_waspmap(
        &mut ra,
        &mut wasp_ra,
        &p,
        &genome,
        None,
        None,
        &outcomes,
    );

    assert_eq!(result.wasp_type, 6);
    assert_eq!(ra.wasp_type, 6);
}

#[test]
fn transcriptome_gene_counts_add_align_updates_none_ambig_multi_and_stranded_counts() {
    let mut transcriptome = Transcriptome {
        ex_g: star_rs::generated::structs::TranscriptomeExG {
            n_ex: 3,
            s: vec![10, 30, 60],
            e: vec![50, 90, 120],
            e_max: vec![50, 90, 120],
            str_: vec![1, 2, 3],
            g: vec![0, 1, 2],
            t: vec![0, 0, 0],
        },
        quants: quantifications_l3_quantifications_quantifications(3),
        ..Default::default()
    };
    let align = Transcript {
        exons: vec![[0, 35, 10, 0, 0]],
        n_exons: 1,
        str_: 0,
        ..Default::default()
    };
    let mut gene1 = Vec::new();

    transcriptome_genecountsaddalign_l4_transcriptome_genecountsaddalign(
        &mut transcriptome,
        1,
        &[align],
        &mut gene1,
    );

    assert_eq!(gene1, vec![-2, 0, 1]);
    assert_eq!(transcriptome.quants.gene_counts.c_ambig, vec![1, 0, 0]);
    assert_eq!(transcriptome.quants.gene_counts.g_count[1][0], 1);
    assert_eq!(transcriptome.quants.gene_counts.g_count[2][1], 1);

    let no_gene = Transcript {
        exons: vec![[0, 200, 10, 0, 0]],
        n_exons: 1,
        str_: 0,
        ..Default::default()
    };
    transcriptome_genecountsaddalign_l4_transcriptome_genecountsaddalign(
        &mut transcriptome,
        1,
        &[no_gene],
        &mut gene1,
    );
    assert_eq!(gene1, vec![-1, -1, -1]);
    assert_eq!(transcriptome.quants.gene_counts.c_none, vec![1, 1, 1]);

    transcriptome_genecountsaddalign_l4_transcriptome_genecountsaddalign(
        &mut transcriptome,
        2,
        &[],
        &mut gene1,
    );
    assert_eq!(gene1, vec![-1, -1, -1]);
    assert_eq!(transcriptome.quants.gene_counts.c_multi, 1);
}

#[test]
fn align_blocks_overlap_exons_tracks_overlap_protrusion_and_splice_concordance() {
    let transcript = Transcript {
        exons: vec![[0, 110, 20, 0, 0], [0, 150, 10, 0, 0]],
        canon_sj: vec![1],
        n_exons: 2,
        ..Default::default()
    };
    assert_eq!(
        transcriptome_alignexonoverlap_l236_alignblocksoverlapexons(
            &transcript,
            2,
            &[0, 29, 50, 59],
            100,
        ),
        (30, true)
    );

    let shifted_junction = Transcript {
        exons: vec![[0, 110, 18, 0, 0], [0, 151, 9, 0, 0]],
        canon_sj: vec![1],
        n_exons: 2,
        ..Default::default()
    };
    assert_eq!(
        transcriptome_alignexonoverlap_l236_alignblocksoverlapexons(
            &shifted_junction,
            2,
            &[0, 29, 50, 59],
            100,
        ),
        (27, false)
    );

    let noncanonical_shift = Transcript {
        exons: vec![[0, 110, 18, 0, 0], [0, 151, 9, 0, 0]],
        canon_sj: vec![-1],
        n_exons: 2,
        ..Default::default()
    };
    assert_eq!(
        transcriptome_alignexonoverlap_l236_alignblocksoverlapexons(
            &noncanonical_shift,
            2,
            &[0, 29, 50, 59],
            100,
        ),
        (27, true)
    );

    let protruding = Transcript {
        exons: vec![[0, 95, 10, 0, 0]],
        n_exons: 1,
        ..Default::default()
    };
    assert_eq!(
        transcriptome_alignexonoverlap_l236_alignblocksoverlapexons(&protruding, 1, &[0, 29], 100,),
        (-1, true)
    );
}

#[test]
fn transcriptome_align_exon_overlap_prioritizes_overlap_types_and_counts_sense_only() {
    let transcriptome = Transcriptome {
        n_tr: 3,
        tr_s: vec![100, 120, 200],
        tr_e: vec![180, 190, 260],
        tr_e_max: vec![190, 190, 260],
        tr_ex_n: vec![2, 1, 1],
        tr_ex_i: vec![0, 2, 3],
        tr_str: vec![1, 2, 1],
        tr_gene: vec![10, 20, 30],
        ex_se: vec![0, 29, 50, 59, 0, 69, 0, 59],
        ..Default::default()
    };

    let exact = Transcript {
        exons: vec![[0, 110, 20, 0, 0], [0, 150, 10, 0, 0]],
        canon_sj: vec![1],
        n_exons: 2,
        str_: 0,
        ..Default::default()
    };
    let mut ann = star_rs::generated::structs::ReadAnnotFeature::default();
    transcriptome_alignexonoverlap_l10_transcriptome_alignexonoverlap(
        &transcriptome,
        1,
        &[exact],
        0,
        &mut ann,
    );
    assert_eq!(ann.ov_type, 1);
    assert_eq!(ann.f_set, std::collections::BTreeSet::from([10]));
    assert_eq!(ann.f_align[0], std::collections::BTreeSet::from([10]));

    let antisense = Transcript {
        exons: vec![[0, 125, 15, 0, 0]],
        n_exons: 1,
        str_: 0,
        ..Default::default()
    };
    let mut ann_as = star_rs::generated::structs::ReadAnnotFeature::default();
    transcriptome_alignexonoverlap_l10_transcriptome_alignexonoverlap(
        &transcriptome,
        1,
        &[antisense],
        0,
        &mut ann_as,
    );
    assert_eq!(ann_as.ov_type, 2);
    assert!(ann_as.f_set.is_empty());
    assert_eq!(ann_as.f_align.len(), 1);
    assert!(ann_as.f_align[0].is_empty());

    let partial = Transcript {
        exons: vec![[0, 110, 18, 0, 0], [0, 151, 9, 0, 0]],
        canon_sj: vec![1],
        n_exons: 2,
        str_: 0,
        ..Default::default()
    };
    let mut ann_partial = star_rs::generated::structs::ReadAnnotFeature::default();
    transcriptome_alignexonoverlap_l10_transcriptome_alignexonoverlap(
        &transcriptome,
        1,
        &[partial],
        0,
        &mut ann_partial,
    );
    assert_eq!(ann_partial.ov_type, 3);
    assert_eq!(ann_partial.f_set, std::collections::BTreeSet::from([10]));

    let outside = Transcript {
        exons: vec![[0, 500, 10, 0, 0]],
        n_exons: 1,
        ..Default::default()
    };
    let mut ann_outside = star_rs::generated::structs::ReadAnnotFeature::default();
    transcriptome_alignexonoverlap_l10_transcriptome_alignexonoverlap(
        &transcriptome,
        1,
        &[outside],
        0,
        &mut ann_outside,
    );
    assert_eq!(ann_outside.ov_type, 7);
    assert!(ann_outside.f_set.is_empty());
}

#[test]
fn transcriptome_quants_allocate_and_output_match_gene_count_table() {
    let mut transcriptome = Transcriptome {
        n_ge: 2,
        ge_id: vec!["geneA".to_string(), "geneB".to_string()],
        quants: quantifications_l3_quantifications_quantifications(1),
        ..Default::default()
    };

    transcriptome_l150_transcriptome_quantsallocate(&mut transcriptome, false);
    assert_eq!(transcriptome.quants.gene_counts.n_ge, 1);

    transcriptome_l150_transcriptome_quantsallocate(&mut transcriptome, true);
    assert_eq!(transcriptome.quants.gene_counts.n_ge, 2);
    assert_eq!(transcriptome.quants.gene_counts.n_type, 3);

    transcriptome.quants.gene_counts.c_multi = 7;
    transcriptome.quants.gene_counts.c_none = vec![1, 2, 3];
    transcriptome.quants.gene_counts.c_ambig = vec![4, 5, 6];
    transcriptome.quants.gene_counts.g_count = vec![vec![10, 11], vec![20, 21], vec![30, 31]];

    let dir = std::env::temp_dir().join(format!("star-rs-quants-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("ReadsPerGene.out.tab");
    let stats = Stats {
        unmapped_mismatch: 1,
        unmapped_short: 2,
        unmapped_other: 3,
        unmapped_multi: 4,
        ..Default::default()
    };

    transcriptome_l156_transcriptome_quantsoutput(&transcriptome, path.to_str().unwrap(), &stats)
        .unwrap();

    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        concat!(
            "N_unmapped\t10\t10\t10\n",
            "N_multimapping\t7\t7\t7\n",
            "N_noFeature\t1\t2\t3\n",
            "N_ambiguous\t4\t5\t6\n",
            "geneA\t10\t20\t30\n",
            "geneB\t11\t21\t31\n",
        )
    );
}

#[test]
fn transcriptome_constructor_loads_star_tables_and_gene_full_ranges() {
    let gene_info = concat!(
        "3\n",
        "geneB nameB typeB\n",
        "geneA nameA typeA\n",
        "geneC nameC typeC\n",
    );
    let transcript_info = concat!("2\n", "tr1 100 180 180 1 2 0 1\n", "tr2 20 40 40 2 1 2 0\n",);
    let exon_info = concat!("3\n", "100 120 0\n", "150 180 21\n", "20 40 0\n",);
    let exon_ge_tr_info = concat!(
        "4\n",
        "300 320 1 2 0\n",
        "100 120 1 1 0\n",
        "150 180 1 1 0\n",
        "20 40 2 0 1\n",
    );

    let (transcriptome, tr_info_dir, log) = transcriptome_l7_transcriptome_transcriptome(
        true,
        false,
        "annotations.gtf",
        "GenomeDir",
        "Pass2Dir",
        "TransformedDir",
        true,
        true,
        true,
        gene_info,
        Some(transcript_info),
        Some(exon_info),
        Some(exon_ge_tr_info),
    )
    .unwrap();

    assert_eq!(tr_info_dir, "Pass2Dir");
    assert_eq!(
        log,
        "Loaded transcript database, nTr=2\nLoaded exon database, nEx=3\n"
    );
    assert_eq!(transcriptome.n_ge, 3);
    assert_eq!(transcriptome.ge_id, ["geneB", "geneA", "geneC"]);
    assert_eq!(transcriptome.ge_name, ["nameB", "nameA", "nameC"]);
    assert_eq!(transcriptome.ge_biotype, ["typeB", "typeA", "typeC"]);

    assert_eq!(transcriptome.n_tr, 2);
    assert_eq!(transcriptome.tr_id, ["tr1", "tr2"]);
    assert_eq!(transcriptome.tr_s, [100, 20]);
    assert_eq!(transcriptome.tr_e, [180, 40]);
    assert_eq!(transcriptome.tr_e_max, [180, 40]);
    assert_eq!(transcriptome.tr_str, [1, 2]);
    assert_eq!(transcriptome.tr_ex_n, [2, 1]);
    assert_eq!(transcriptome.tr_ex_i, [0, 2]);
    assert_eq!(transcriptome.tr_gene, [1, 0]);
    assert_eq!(transcriptome.tr_len, [52, 21]);
    assert_eq!(transcriptome.ex_se, [100, 120, 150, 180, 20, 40]);
    assert_eq!(transcriptome.ex_len_cum, [0, 21, 0]);

    assert_eq!(transcriptome.ex_g.s, [300, 100, 150, 20]);
    assert_eq!(transcriptome.ex_g.e, [320, 120, 180, 40]);
    assert_eq!(transcriptome.ex_g.e_max, [320, 320, 320, 320]);
    assert_eq!(transcriptome.ex_g.str_, [1, 1, 1, 2]);
    assert_eq!(transcriptome.ex_g.g, [2, 1, 1, 0]);
    assert_eq!(transcriptome.ex_g.t, [0, 0, 0, 1]);

    assert_eq!(transcriptome.gene_full.s, [20, 100, 300]);
    assert_eq!(transcriptome.gene_full.e, [40, 180, 320]);
    assert_eq!(transcriptome.gene_full.e_max, [40, 180, 320]);
    assert_eq!(transcriptome.gene_full.str_, [2, 1, 1]);
    assert_eq!(transcriptome.gene_full.g, [0, 1, 2]);

    let (empty, tr_info_dir, log) = transcriptome_l7_transcriptome_transcriptome(
        false,
        true,
        "-",
        "GenomeDir",
        "Pass2Dir",
        "TransformedDir",
        true,
        true,
        true,
        gene_info,
        Some(transcript_info),
        Some(exon_info),
        Some(exon_ge_tr_info),
    )
    .unwrap();
    assert_eq!(empty, Transcriptome::default());
    assert!(tr_info_dir.is_empty());
    assert!(log.is_empty());
}

#[test]
fn transcriptome_quant_align_scans_candidate_transcripts_backwards() {
    let transcriptome = Transcriptome {
        n_tr: 2,
        tr_s: vec![100, 100],
        tr_e: vec![180, 180],
        tr_e_max: vec![180, 180],
        tr_ex_n: vec![2, 2],
        tr_ex_i: vec![0, 2],
        tr_str: vec![1, 2],
        ex_se: vec![0, 29, 50, 59, 0, 29, 50, 59],
        ex_len_cum: vec![0, 30, 0, 30],
        ..Default::default()
    };
    let mut genomic = Transcript {
        exons: vec![[3, 110, 20, 7, 0], [23, 150, 10, 7, 0]],
        n_exons: 2,
        canon_sj: vec![1, 0],
        l_read: 40,
        str_: 1,
        ..Default::default()
    };
    let mut transcript_aligns = Vec::new();

    assert_eq!(
        transcriptome_quantalign_l91_transcriptome_quantalign(
            &transcriptome,
            &mut genomic,
            &mut transcript_aligns,
        ),
        2
    );
    assert_eq!(transcript_aligns.len(), 2);
    assert_eq!(transcript_aligns[0].chr, 1);
    assert_eq!(transcript_aligns[0].str_, 0);
    assert_eq!(transcript_aligns[0].n_exons, 1);
    assert_eq!(transcript_aligns[0].exons[0][EX_G], 0);
    assert_eq!(transcript_aligns[1].chr, 0);
    assert_eq!(transcript_aligns[1].str_, 1);
    assert_eq!(transcript_aligns[1].n_exons, 1);
    assert_eq!(transcript_aligns[1].exons[0][EX_G], 10);
    assert_eq!(genomic.canon_sj[1], -999);
}

#[test]
fn transcriptome_quant_align_returns_zero_outside_transcript_starts() {
    let transcriptome = Transcriptome {
        n_tr: 1,
        tr_s: vec![100],
        tr_e: vec![180],
        tr_e_max: vec![180],
        tr_ex_n: vec![1],
        tr_ex_i: vec![0],
        tr_str: vec![1],
        ex_se: vec![0, 29],
        ex_len_cum: vec![0],
        ..Default::default()
    };
    let mut genomic = Transcript {
        exons: vec![[0, 90, 5, 0, 0]],
        n_exons: 1,
        canon_sj: vec![0],
        ..Default::default()
    };
    let mut transcript_aligns = vec![Transcript {
        chr: 77,
        ..Default::default()
    }];

    assert_eq!(
        transcriptome_quantalign_l91_transcriptome_quantalign(
            &transcriptome,
            &mut genomic,
            &mut transcript_aligns,
        ),
        0
    );
    assert_eq!(transcript_aligns[0].chr, 77);
}

#[test]
fn read_align_quant_transcriptome_filters_extends_and_records_bam_requests() {
    let transcriptome = Transcriptome {
        n_tr: 1,
        tr_s: vec![100],
        tr_e: vec![130],
        tr_e_max: vec![130],
        tr_ex_n: vec![1],
        tr_ex_i: vec![0],
        tr_str: vec![1],
        ex_se: vec![0, 30],
        ex_len_cum: vec![0],
        ..Default::default()
    };
    let mut genome = Genome {
        g: vec![0; 160],
        ..Default::default()
    };
    genome.g[100..108].copy_from_slice(&[0, 1, 2, 3, 0, 1, 2, 3]);

    let read_align = ReadAlign {
        l_read: 8,
        read_nmates: 1,
        read_length: vec![8, 0],
        out_filter_mismatch_nmax_total: 4,
        read1: [
            vec![0, 1, 2, 3, 0, 1, 2, 3],
            vec![3, 2, 1, 0, 3, 2, 1, 0],
            vec![3, 2, 1, 0, 3, 2, 1, 0],
        ],
        ..Default::default()
    };
    let p = Parameters {
        read_nmates: 1,
        quant_tr_sam_indel: false,
        quant_tr_sam_single_end: true,
        quant_tr_sam_soft_clip: false,
        quant_tr_sam_bam_yes: true,
        out_filter_mismatch_nover_lmax: 1.0,
        ..Default::default()
    };

    let good = Transcript {
        exons: vec![[1, 101, 6, 0, 0]],
        n_exons: 1,
        canon_sj: vec![0],
        l_read: 8,
        str_: 0,
        ..Default::default()
    };
    let indel_rejected = Transcript {
        n_del: 1,
        ..good.clone()
    };

    let out = readalign_quanttranscriptome_l7_readalign_quanttranscriptome(
        &read_align,
        &p,
        &genome,
        &transcriptome,
        &[indel_rejected, good],
        0.7,
    );

    assert_eq!(out.n_align_t, 1);
    assert_eq!(out.align_t.len(), 1);
    assert_eq!(out.align_t[0].exons[0][EX_R], 0);
    assert_eq!(out.align_t[0].exons[0][EX_G], 0);
    assert_eq!(out.align_t[0].exons[0][EX_L], 8);
    assert!(out.align_t[0].primary_flag);
    assert_eq!(out.bam_requests.len(), 1);
    assert_eq!(out.bam_requests[0].n_align_t, 1);
    assert_eq!(out.bam_requests[0].i_align_t, 0);
    assert!(out.bam_requests[0].transcript.primary_flag);

    let paired_single_mate = Transcript {
        exons: vec![[0, 100, 8, 0, 0]],
        n_exons: 1,
        canon_sj: vec![0],
        l_read: 8,
        ..Default::default()
    };
    let paired_ra = ReadAlign {
        l_read: 8,
        read_length: vec![4, 4],
        read1: read_align.read1.clone(),
        out_filter_mismatch_nmax_total: 4,
        ..Default::default()
    };
    let paired_p = Parameters {
        read_nmates: 2,
        quant_tr_sam_indel: true,
        quant_tr_sam_single_end: false,
        quant_tr_sam_soft_clip: true,
        out_filter_mismatch_nover_lmax: 1.0,
        ..Default::default()
    };
    let paired_out = readalign_quanttranscriptome_l7_readalign_quanttranscriptome(
        &paired_ra,
        &paired_p,
        &genome,
        &transcriptome,
        &[paired_single_mate],
        0.0,
    );
    assert_eq!(paired_out.n_align_t, 0);
    assert!(paired_out.align_t.is_empty());
    assert!(paired_out.bam_requests.is_empty());
}

#[test]
fn transcriptome_classify_align_records_gene_and_velocyto_annotations() {
    let transcriptome = Transcriptome {
        n_tr: 2,
        tr_s: vec![100, 100],
        tr_e: vec![180, 180],
        tr_e_max: vec![180, 180],
        tr_ex_n: vec![2, 2],
        tr_ex_i: vec![0, 2],
        tr_str: vec![1, 2],
        tr_gene: vec![7, 7],
        ex_se: vec![0, 29, 50, 59, 0, 29, 50, 59],
        ex_len_cum: vec![0, 30, 0, 30],
        ..Default::default()
    };
    let align = Transcript {
        exons: vec![[3, 110, 20, 0, 0], [23, 150, 10, 0, 0]],
        n_exons: 2,
        canon_sj: vec![1, 0],
        l_read: 40,
        str_: 1,
        ..Default::default()
    };
    let mut read_annot = ReadAnnotations {
        annot_features: vec![ReadAnnotFeature::default(); 8],
        ..Default::default()
    };
    let mut feature_yes = vec![false; 8];
    feature_yes[SOLO_FEATURE_VELOCYTO as usize] = true;

    transcriptome_classifyalign_l177_transcriptome_classifyalign(
        &transcriptome,
        &ParametersSolo {
            strand: -1,
            feature_yes,
            ..Default::default()
        },
        &[align],
        &mut read_annot,
    );

    let gene_ann = &read_annot.annot_features[SOLO_FEATURE_GENE as usize];
    assert_eq!(gene_ann.ov_type, 1);
    assert_eq!(gene_ann.f_set, std::collections::BTreeSet::from([7]));
    assert_eq!(gene_ann.f_align[0], std::collections::BTreeSet::from([7]));
    assert_eq!(read_annot.transcript_concordant, vec![[1, 10], [0, 0]]);
    assert_eq!(read_annot.gene_velocyto_simple[0], 7);
    assert_eq!(
        read_annot.gene_velocyto_simple[1],
        (1 << ALIGN_VS_TRANSCRIPT_EXON_INTRON_SPAN) | (1 << ALIGN_VS_TRANSCRIPT_CONCORDANT)
    );
    assert_eq!(read_annot.tr_velocyto_type.len(), 2);
    assert_eq!(read_annot.tr_velocyto_type[0].tr, 1);
    assert_eq!(
        read_annot.tr_velocyto_type[0].type_,
        1 << ALIGN_VS_TRANSCRIPT_CONCORDANT
    );

    let stranded_transcriptome = Transcriptome {
        tr_gene: vec![7, 8],
        ..transcriptome.clone()
    };
    let mut stranded = ReadAnnotations {
        annot_features: vec![ReadAnnotFeature::default(); 8],
        ..Default::default()
    };
    transcriptome_classifyalign_l177_transcriptome_classifyalign(
        &stranded_transcriptome,
        &ParametersSolo {
            strand: 1,
            feature_yes: vec![false; 8],
            ..Default::default()
        },
        &[Transcript {
            exons: vec![[0, 110, 5, 0, 0]],
            n_exons: 1,
            canon_sj: vec![0],
            str_: 0,
            ..Default::default()
        }],
        &mut stranded,
    );
    assert_eq!(stranded.transcript_concordant, vec![[1, 10]]);
    assert_eq!(
        stranded.annot_features[SOLO_FEATURE_GENE as usize].f_set,
        std::collections::BTreeSet::from([8])
    );
}

#[test]
fn quant_align_to_transcript_converts_genomic_blocks_to_transcript_coordinates() {
    let ex_se = [0, 29, 50, 59];
    let ex_len_cum = [0, 30];

    let mut genomic = Transcript {
        exons: vec![[3, 110, 20, 7, 0], [23, 150, 10, 7, 0]],
        n_exons: 2,
        canon_sj: vec![1, 0],
        l_read: 40,
        primary_flag: true,
        ..Default::default()
    };
    let mut transcript = Transcript {
        exons: vec![[u32::MAX; 5]; 3],
        canon_sj: vec![u32::MAX as i32; 3],
        sj_annot: vec![1, 1, 1],
        shift_sj: vec![[9, 9]; 3],
        sj_str: vec![2, 2, 2],
        primary_flag: true,
        ..Default::default()
    };

    assert_eq!(
        transcriptome_quantalign_l5_aligntotranscript(
            &mut genomic,
            100,
            1,
            &ex_se,
            &ex_len_cum,
            2,
            &mut transcript,
        ),
        1
    );
    assert_eq!(genomic.canon_sj[1], -999);
    assert!(!transcript.primary_flag);
    assert_eq!(transcript.n_exons, 1);
    assert_eq!(transcript.exons[0], [3, 10, 30, 7, u32::MAX]);
    assert_eq!(&transcript.sj_annot[..1], &[0]);
    assert_eq!(&transcript.shift_sj[..1], &[[0, 0]]);
    assert_eq!(&transcript.sj_str[..1], &[0]);
}

#[test]
fn quant_align_to_transcript_handles_negative_strand_mates_indels_and_rejections() {
    let ex_se = [0, 29, 50, 79];
    let ex_len_cum = [0, 30];

    let mut negative = Transcript {
        exons: vec![[0, 110, 10, 0, 0], [10, 120, 5, 0, 0], [20, 150, 12, 1, 0]],
        n_exons: 3,
        canon_sj: vec![-1, -3, 0],
        l_read: 40,
        ..Default::default()
    };
    let mut transcript = Transcript::default();
    assert_eq!(
        transcriptome_quantalign_l5_aligntotranscript(
            &mut negative,
            100,
            2,
            &ex_se,
            &ex_len_cum,
            2,
            &mut transcript,
        ),
        1
    );
    assert_eq!(transcript.n_exons, 3);
    assert_eq!(transcript.exons[0], [8, 18, 12, 1, 0]);
    assert_eq!(transcript.exons[1], [25, 35, 5, 0, 0]);
    assert_eq!(transcript.exons[2], [30, 40, 10, 0, 0]);
    assert_eq!(transcript.canon_sj[0], -3);
    assert_eq!(transcript.canon_sj[1], -1);

    let mut intron_start = Transcript {
        exons: vec![[0, 135, 5, 0, 0]],
        n_exons: 1,
        canon_sj: vec![0],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_quantalign_l5_aligntotranscript(
            &mut intron_start,
            100,
            1,
            &ex_se,
            &ex_len_cum,
            2,
            &mut Transcript::default(),
        ),
        0
    );

    let mut overhang = Transcript {
        exons: vec![[0, 120, 20, 0, 0]],
        n_exons: 1,
        canon_sj: vec![0],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_quantalign_l5_aligntotranscript(
            &mut overhang,
            100,
            1,
            &ex_se,
            &ex_len_cum,
            2,
            &mut Transcript::default(),
        ),
        0
    );

    let mut sj_mismatch = Transcript {
        exons: vec![[0, 110, 20, 0, 0], [20, 151, 10, 0, 0]],
        n_exons: 2,
        canon_sj: vec![1, 0],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_quantalign_l5_aligntotranscript(
            &mut sj_mismatch,
            100,
            1,
            &ex_se,
            &ex_len_cum,
            2,
            &mut Transcript::default(),
        ),
        0
    );
}

#[test]
fn align_to_transcript_classifies_blocks_and_sets_distances() {
    let ex_se = [0, 29, 50, 59];
    let ex_len_cum = [0, 30];

    let concordant = Transcript {
        exons: vec![[0, 110, 20, 0, 0], [0, 150, 10, 0, 0]],
        n_exons: 2,
        canon_sj: vec![1],
        ..Default::default()
    };
    let mut dist = [u32::MAX; 2];
    assert_eq!(
        transcriptome_classifyalign_l8_aligntotranscript(
            &concordant,
            100,
            2,
            &ex_se,
            &ex_len_cum,
            &mut dist,
        ),
        ALIGN_VS_TRANSCRIPT_CONCORDANT
    );
    assert_eq!(dist, [10, 0]);

    let span = Transcript {
        exons: vec![[0, 115, 35, 0, 0]],
        n_exons: 1,
        canon_sj: vec![],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_classifyalign_l8_aligntotranscript(
            &span,
            100,
            2,
            &ex_se,
            &ex_len_cum,
            &mut [0; 2],
        ),
        ALIGN_VS_TRANSCRIPT_EXON_INTRON_SPAN
    );

    let intronic = Transcript {
        exons: vec![[0, 130, 10, 0, 0]],
        n_exons: 1,
        canon_sj: vec![],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_classifyalign_l8_aligntotranscript(
            &intronic,
            100,
            2,
            &ex_se,
            &ex_len_cum,
            &mut [0; 2],
        ),
        ALIGN_VS_TRANSCRIPT_INTRON
    );

    let mixed = Transcript {
        exons: vec![[0, 105, 10, 0, 0], [0, 135, 8, 0, 0]],
        n_exons: 2,
        canon_sj: vec![-3],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_classifyalign_l8_aligntotranscript(
            &mixed,
            100,
            2,
            &ex_se,
            &ex_len_cum,
            &mut [0; 2],
        ),
        ALIGN_VS_TRANSCRIPT_EXON_INTRON
    );

    let sj_mismatch = Transcript {
        exons: vec![[0, 110, 18, 0, 0], [0, 151, 5, 0, 0]],
        n_exons: 2,
        canon_sj: vec![1],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_classifyalign_l8_aligntotranscript(
            &sj_mismatch,
            100,
            2,
            &ex_se,
            &ex_len_cum,
            &mut [0; 2],
        ),
        -1
    );

    let protrudes_left = Transcript {
        exons: vec![[0, 99, 5, 0, 0]],
        n_exons: 1,
        canon_sj: vec![],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_classifyalign_l8_aligntotranscript(
            &protrudes_left,
            100,
            2,
            &ex_se,
            &ex_len_cum,
            &mut [0; 2],
        ),
        -1
    );
}

#[test]
fn align_to_transcript_min_overlap_classifies_exonic_intronic_span_and_spliced_conflict() {
    let ex_se = [0, 19, 40, 59, 80, 99];

    let exonic = Transcript {
        exons: vec![[0, 105, 10, 0, 0]],
        n_exons: 1,
        canon_sj: vec![],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_classifyalign_l93_aligntotranscriptminoverlap(&exonic, 100, &ex_se, 3, 6),
        ALIGN_VS_TRANSCRIPT_CONCORDANT
    );

    let span = Transcript {
        exons: vec![[0, 115, 35, 0, 0]],
        n_exons: 1,
        canon_sj: vec![],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_classifyalign_l93_aligntotranscriptminoverlap(&span, 100, &ex_se, 3, 6),
        ALIGN_VS_TRANSCRIPT_EXON_INTRON_SPAN
    );

    let intronic = Transcript {
        exons: vec![[0, 125, 8, 0, 0]],
        n_exons: 1,
        canon_sj: vec![],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_classifyalign_l93_aligntotranscriptminoverlap(&intronic, 100, &ex_se, 3, 6),
        ALIGN_VS_TRANSCRIPT_INTRON
    );

    let mixed = Transcript {
        exons: vec![[0, 105, 10, 0, 0], [0, 125, 8, 0, 0]],
        n_exons: 2,
        canon_sj: vec![-3],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_classifyalign_l93_aligntotranscriptminoverlap(&mixed, 100, &ex_se, 3, 6),
        ALIGN_VS_TRANSCRIPT_EXON_INTRON
    );

    let spliced_conflict = Transcript {
        exons: vec![[0, 125, 8, 0, 0]],
        n_exons: 1,
        canon_sj: vec![],
        sj_yes: true,
        ..Default::default()
    };
    assert_eq!(
        transcriptome_classifyalign_l93_aligntotranscriptminoverlap(
            &spliced_conflict,
            100,
            &ex_se,
            3,
            6,
        ),
        -1
    );

    let large_intron = [0, 19, 1_000_030, 1_000_060];
    let swallowed = Transcript {
        exons: vec![[0, 125, 8, 0, 0]],
        n_exons: 1,
        canon_sj: vec![],
        ..Default::default()
    };
    assert_eq!(
        transcriptome_classifyalign_l93_aligntotranscriptminoverlap(
            &swallowed,
            100,
            &large_intron,
            2,
            6,
        ),
        -1
    );
}

#[test]
fn transcriptome_gene_full_align_overlap_scans_blocks_and_strand_filters() {
    let transcriptome = Transcriptome {
        n_ge: 4,
        gene_full: TranscriptomeGeneFull {
            s: vec![10, 30, 60, 150],
            e: vec![50, 90, 120, 180],
            e_max: vec![50, 90, 120, 180],
            str_: vec![1, 0, 1, 1],
            g: vec![101, 202, 303, 404],
        },
        ..Default::default()
    };
    let aligns = vec![
        Transcript {
            exons: vec![[0, 35, 10, 0, 0], [20, 45, 5, 0, 0]],
            n_exons: 2,
            str_: 1,
            ..Default::default()
        },
        Transcript {
            exons: vec![[0, 70, 10, 0, 0]],
            n_exons: 1,
            str_: 0,
            ..Default::default()
        },
    ];

    let mut ann = star_rs::generated::structs::ReadAnnotFeature::default();
    transcriptome_genefullalignoverlap_l5_transcriptome_genefullalignoverlap(
        &transcriptome,
        2,
        &aligns,
        -1,
        &mut ann,
    );

    assert_eq!(ann.f_set, std::collections::BTreeSet::from([101, 202, 303]));
    assert_eq!(ann.f_align[0], std::collections::BTreeSet::from([101, 202]));
    assert_eq!(ann.f_align[1], std::collections::BTreeSet::from([202, 303]));
    assert_eq!(ann.ov_type, 0);

    let mut strand_filtered = star_rs::generated::structs::ReadAnnotFeature::default();
    transcriptome_genefullalignoverlap_l5_transcriptome_genefullalignoverlap(
        &transcriptome,
        2,
        &aligns,
        0,
        &mut strand_filtered,
    );
    assert_eq!(
        strand_filtered.f_set,
        std::collections::BTreeSet::from([202, 303])
    );
    assert_eq!(
        strand_filtered.f_align[0],
        std::collections::BTreeSet::from([202])
    );
    assert_eq!(
        strand_filtered.f_align[1],
        std::collections::BTreeSet::from([303])
    );
}

#[test]
fn transcriptome_gene_full_exon_over_intron_prioritizes_exonic_then_intronic() {
    let mut concordant = star_rs::generated::structs::ReadAnnotFeature::default();
    concordant.f_set.insert(17);
    concordant.f_align = vec![std::collections::BTreeSet::from([17])];
    concordant.ov_type = 0;
    let mut ann = star_rs::generated::structs::ReadAnnotFeature::default();
    transcriptome_genefullalignoverlap_exonoverintron_l5_transcriptome_genefullalignoverlap_exonoverintron(
        &Transcriptome::default(),
        1,
        &[],
        -1,
        &mut ann,
        &concordant,
    );
    assert_eq!(ann.f_set, std::collections::BTreeSet::from([17]));
    assert_eq!(ann.ov_type, 1);

    let transcriptome = Transcriptome {
        n_ge: 3,
        gene_full: TranscriptomeGeneFull {
            s: vec![10, 30, 60],
            e: vec![50, 90, 120],
            e_max: vec![50, 90, 120],
            str_: vec![1, 0, 1],
            g: vec![101, 202, 303],
        },
        ..Default::default()
    };
    let aligns = vec![
        Transcript {
            exons: vec![[0, 35, 10, 0, 0], [20, 45, 5, 0, 0]],
            n_exons: 2,
            str_: 1,
            ..Default::default()
        },
        Transcript {
            exons: vec![[0, 70, 10, 0, 0]],
            n_exons: 1,
            str_: 0,
            ..Default::default()
        },
    ];
    let mut intronic = star_rs::generated::structs::ReadAnnotFeature::default();
    transcriptome_genefullalignoverlap_exonoverintron_l5_transcriptome_genefullalignoverlap_exonoverintron(
        &transcriptome,
        2,
        &aligns,
        -1,
        &mut intronic,
        &star_rs::generated::structs::ReadAnnotFeature::default(),
    );
    assert_eq!(intronic.ov_type, 5);
    assert_eq!(
        intronic.f_set,
        std::collections::BTreeSet::from([101, 202, 303])
    );
    assert_eq!(intronic.f_align.len(), 2);
    assert_eq!(
        intronic.f_align[0],
        std::collections::BTreeSet::from([101, 202])
    );
    assert_eq!(
        intronic.f_align[1],
        std::collections::BTreeSet::from([202, 303])
    );

    let mut strand_filtered = star_rs::generated::structs::ReadAnnotFeature::default();
    transcriptome_genefullalignoverlap_exonoverintron_l5_transcriptome_genefullalignoverlap_exonoverintron(
        &transcriptome,
        2,
        &aligns,
        0,
        &mut strand_filtered,
        &star_rs::generated::structs::ReadAnnotFeature::default(),
    );
    assert_eq!(
        strand_filtered.f_set,
        std::collections::BTreeSet::from([202, 303])
    );
    assert_eq!(strand_filtered.ov_type, 5);
}

#[test]
fn read_align_mapped_filter_matches_original_unmap_type_order() {
    let mut no_windows = ReadAlign::default();
    readalign_mappedfilter_l3_readalign_mappedfilter(&mut no_windows, 10, 0.5, 5, 0.25, 0.2, 10);
    assert_eq!(no_windows.unmap_type, 0);
    assert_eq!(no_windows.stats_ra.unmapped_other, 1);

    let mut too_short = ReadAlign {
        n_w: 1,
        l_read: 100,
        tr_best: Transcript {
            max_score: 40,
            n_match: 90,
            r_length: 90,
            ..Default::default()
        },
        ..Default::default()
    };
    readalign_mappedfilter_l3_readalign_mappedfilter(&mut too_short, 10, 0.5, 5, 0.25, 0.2, 10);
    assert_eq!(too_short.unmap_type, 1);
    assert_eq!(too_short.stats_ra.unmapped_short, 1);

    let mut too_many_mm = ReadAlign {
        n_w: 1,
        l_read: 100,
        out_filter_mismatch_nmax_total: 3,
        tr_best: Transcript {
            max_score: 80,
            n_match: 80,
            n_mm: 4,
            r_length: 80,
            ..Default::default()
        },
        ..Default::default()
    };
    readalign_mappedfilter_l3_readalign_mappedfilter(&mut too_many_mm, 10, 0.5, 5, 0.25, 0.2, 10);
    assert_eq!(too_many_mm.unmap_type, 2);
    assert_eq!(too_many_mm.stats_ra.unmapped_mismatch, 1);

    let mut too_multi = ReadAlign {
        n_w: 1,
        n_tr: 11,
        l_read: 100,
        out_filter_mismatch_nmax_total: 10,
        tr_best: Transcript {
            max_score: 80,
            n_match: 80,
            n_mm: 1,
            r_length: 80,
            ..Default::default()
        },
        ..Default::default()
    };
    readalign_mappedfilter_l3_readalign_mappedfilter(&mut too_multi, 10, 0.5, 5, 0.25, 0.2, 10);
    assert_eq!(too_multi.unmap_type, 3);
    assert_eq!(too_multi.stats_ra.unmapped_multi, 1);

    let mut mapped = too_multi.clone();
    mapped.n_tr = 2;
    readalign_mappedfilter_l3_readalign_mappedfilter(&mut mapped, 10, 0.5, 5, 0.25, 0.2, 10);
    assert_eq!(mapped.unmap_type, -1);
}

#[test]
fn read_align_out_reads_unmapped_formats_fastx_mates_like_original() {
    let read_name_mates = vec!["@r/1".to_string(), "@r/2".to_string()];
    let read_name_extra = vec!["extra1".to_string(), String::new()];
    let read0 = vec!["ACGT".to_string(), "TGCA".to_string()];
    let qual0 = vec!["IIII".to_string(), "HHHH".to_string()];
    let mut out = vec![String::new(), String::new()];

    readalign_outputalignments_l259_readalign_outreadsunmapped(
        "Fastx",
        2,
        2,
        &read_name_mates,
        7,
        &read_name_extra,
        &[true, false],
        &read0,
        2,
        &qual0,
        &mut out,
    );

    assert_eq!(out[0], "@r/1 0:7: extra1 10\nACGT\n+\nIIII\n");
    assert_eq!(out[1], "@r/2 1:7:  10\nTGCA\n+\nHHHH\n");

    readalign_outputalignments_l259_readalign_outreadsunmapped(
        "None",
        2,
        2,
        &read_name_mates,
        7,
        &read_name_extra,
        &[true, false],
        &read0,
        2,
        &qual0,
        &mut out,
    );
    assert_eq!(out[0], "@r/1 0:7: extra1 10\nACGT\n+\nIIII\n");
}

#[test]
fn read_align_out_filter_by_sjout_holds_unannotated_junction_reads() {
    let mut unmap_type = -1;
    let mut pass = true;
    let tr_mult = vec![Transcript {
        n_exons: 2,
        canon_sj: vec![1],
        sj_annot: vec![0],
        ..Default::default()
    }];
    let mut stats = Stats {
        read_n: 12,
        read_bases: 100,
        ..Default::default()
    };
    let mut streams = vec![String::new(), String::new()];
    let read_name_mates = vec!["@r/1".to_string(), "@r/2".to_string()];
    let read_name_extra = vec!["x1".to_string(), String::new()];
    let read0 = vec!["AAAA".to_string(), "CCCC".to_string()];
    let qual0 = vec!["IIII".to_string(), "HHHH".to_string()];
    let mut sj = OutSJ::default();

    readalign_outputalignments_l90_readalign_outfilterbysjout(
        &mut unmap_type,
        &mut pass,
        1,
        1,
        &tr_mult,
        &mut stats,
        &[4, 6],
        2,
        &mut streams,
        &read_name_mates,
        77,
        5,
        3,
        &read_name_extra,
        &read0,
        2,
        &qual0,
        false,
        "All",
        &mut sj,
    );

    assert!(!pass);
    assert_eq!(unmap_type, -3);
    assert_eq!(stats.read_n, 11);
    assert_eq!(stats.read_bases, 90);
    assert_eq!(streams[0], "@r/1 77 5 3 x1\nAAAA\n+\nIIII\n");
    assert_eq!(streams[1], "@r/2 77 5 3\nCCCC\n+\nHHHH\n");

    let mut early_unmap_type = 1;
    let mut early_pass = false;
    readalign_outputalignments_l90_readalign_outfilterbysjout(
        &mut early_unmap_type,
        &mut early_pass,
        1,
        1,
        &tr_mult,
        &mut stats,
        &[4, 6],
        2,
        &mut streams,
        &read_name_mates,
        77,
        5,
        3,
        &read_name_extra,
        &read0,
        2,
        &qual0,
        false,
        "All",
        &mut sj,
    );
    assert!(early_pass);
    assert_eq!(early_unmap_type, 1);
}

#[test]
fn read_align_output_splice_graph_sam_formats_mapped_and_unmapped_records() {
    let p = Parameters {
        out_sam_mapq_unique: 255,
        out_sam_flag_and: 0xffff,
        out_sam_flag_or: 0,
        out_sam_attr_ih_start: 1,
        out_sam_attr_order: vec![ATTR_NH, ATTR_HI, ATTR_AS, ATTR_NM_LOWER, ATTR_NM, ATTR_RG],
        out_sam_attr_rg: vec!["rgA".to_string()],
        out_sam_mode: "Full".to_string(),
        read_files_type_n: 10,
        ..Default::default()
    };
    let genome = Genome {
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![100],
        ..Default::default()
    };
    let tr = Transcript {
        cigar: vec![[0, 4], [3, 6], [0, 3]],
        chr: 0,
        g_start: 120,
        str_: 1,
        primary_flag: false,
        max_score: 42,
        n_mm: 2,
        l_ins: 1,
        l_del: 3,
        ..Default::default()
    };
    let read0 = vec!["ACGTAAA".to_string()];
    let qual0 = vec!["abcdefg".to_string()];
    let read_name_extra = vec!["XA:i:1".to_string()];

    let (sam, bytes) = readalign_outputsplicegraphsam_l5_readalign_outputsplicegraphsam(
        &tr,
        2,
        1,
        b'N',
        -1,
        "@read1",
        &read0,
        2,
        &qual0,
        &p,
        0,
        &read_name_extra,
        7,
        &genome,
    )
    .unwrap();

    assert_eq!(bytes as usize, sam.len());
    assert_eq!(
        sam,
        "read1\t272\tchr1\t21\t3\t4M6N3M\t*\t0\t0\tTTTACGT\tgfedcba\tNH:i:2\tHI:i:2\tAS:i:42\tnM:i:2\tNM:i:6\tRG:Z:rgA\tXA:i:1\n"
    );

    let (unmapped, unmapped_bytes) =
        readalign_outputsplicegraphsam_l5_readalign_outputsplicegraphsam(
            &tr,
            0,
            0,
            b'Y',
            4,
            "@read1",
            &read0,
            1,
            &qual0,
            &p,
            0,
            &read_name_extra,
            7,
            &genome,
        )
        .unwrap();
    assert_eq!(unmapped_bytes as usize, unmapped.len());
    assert_eq!(
        unmapped,
        "read1\t516\t*\t0\t0\t*\t*\t0\t0\tACGTAAA\t*\tNH:i:0\tHI:i:0\tAS:i:42\tnM:i:2\tuT:A:4\tRG:Z:rgA\tXA:i:1\n"
    );
}

#[test]
fn read_align_output_transcript_sam_formats_paired_mapped_records() {
    let mut p = Parameters {
        read_nmates: 2,
        out_sam_mapq_unique: 255,
        out_sam_flag_and: 0xffff,
        out_sam_attr_ih_start: 1,
        out_sam_attr_order: vec![
            ATTR_NH,
            ATTR_HI,
            ATTR_AS,
            ATTR_NM_LOWER,
            ATTR_JM,
            ATTR_JI,
            ATTR_XS,
            ATTR_NM,
            ATTR_MD,
            ATTR_MC,
            ATTR_RG,
            ATTR_HA,
        ],
        out_sam_attr_present: star_rs::generated::structs::SamAttrPresent {
            nm: true,
            md: true,
            mc: true,
            ..Default::default()
        },
        out_sam_attr_rg: vec!["rg1".to_string()],
        out_sam_mode: "Full".to_string(),
        read_files_type_n: 10,
        genome_num_to_nt: b"ACGTN".to_vec(),
        ..Default::default()
    };
    p.p_ge.transform.type_ = 2;

    let mut genome_bases = vec![4u8; 160];
    genome_bases[110..113].copy_from_slice(&[0, 1, 2]);
    genome_bases[120..122].copy_from_slice(&[3, 0]);
    genome_bases[140..144].copy_from_slice(&[0, 0, 1, 1]);
    let genome = Genome {
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![100],
        g: genome_bases,
        p_ge: p.p_ge.clone(),
        ..Default::default()
    };

    let tr = Transcript {
        exons: vec![[0, 110, 3, 0, 0], [3, 120, 2, 0, 0], [6, 140, 4, 1, 0]],
        canon_sj: vec![1, -3],
        sj_annot: vec![1, 0],
        n_exons: 3,
        chr: 0,
        str_: 0,
        ro_str: 0,
        primary_flag: true,
        max_score: 50,
        n_mm: 1,
        sj_motif_strand: 1,
        haplo_type: 2,
        ..Default::default()
    };
    let read0 = vec!["ACGTA".to_string(), "AACC".to_string()];
    let qual0 = vec!["IIIII".to_string(), "abcd".to_string()];
    let read_name_extra = vec!["XA:i:1".to_string(), String::new()];
    let mut read1 = [Vec::new(), Vec::new(), Vec::new()];
    read1[0] = vec![0, 1, 2, 3, 0, 4, 0, 0, 1, 1];
    let clip_mates = vec![
        vec![ClipMate::default(), ClipMate::default()],
        vec![ClipMate::default(), ClipMate::default()],
    ];
    let mut sam = String::new();

    let bytes = readalign_outputtranscriptsam_l5_readalign_outputtranscriptsam(
        &tr,
        2,
        0,
        u32::MAX,
        u32::MAX,
        -1,
        -1,
        None,
        &mut sam,
        b'N',
        "@read1",
        &read0,
        2,
        &qual0,
        &p,
        0,
        &read_name_extra,
        10,
        &[5, 4],
        &[5, 4],
        &clip_mates,
        &read1,
        &genome,
    )
    .unwrap();

    assert_eq!(bytes as usize, sam.len());
    assert_eq!(
        sam,
        "read1\t99\tchr1\t11\t3\t3M7N2M\t=\t41\t34\tACGTA\tIIIII\tNH:i:2\tHI:i:1\tAS:i:50\tnM:i:1\tjM:B:c,21\tjI:B:i,14,20\tXS:A:+\tNM:i:0\tMD:Z:5\tMC:Z:4M\tRG:Z:rg1\tha:i:2\tXA:i:1\nread1\t147\tchr1\t41\t3\t4M\t=\t11\t-34\tGGTT\tdcba\tNH:i:2\tHI:i:1\tAS:i:50\tnM:i:1\tjM:B:c,-1\tjI:B:i,-1\tXS:A:+\tNM:i:0\tMD:Z:4\tMC:Z:3M7N2M\tRG:Z:rg1\tha:i:2\n"
    );
}

#[test]
fn read_align_output_transcript_sam_formats_unmapped_mates() {
    let p = Parameters {
        read_nmates: 2,
        out_sam_attr_rg: vec!["rgA".to_string()],
        out_sam_mode: "Full".to_string(),
        read_files_type_n: 10,
        out_sam_unmapped_keep_pairs: true,
        ..Default::default()
    };
    let genome = Genome {
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![100],
        ..Default::default()
    };
    let tr = Transcript {
        exons: vec![[0, 120, 4, 0, 0]],
        chr: 0,
        str_: 1,
        primary_flag: false,
        max_score: 30,
        n_mm: 2,
        ..Default::default()
    };
    let clip_mates = vec![
        vec![ClipMate::default(), ClipMate::default()],
        vec![ClipMate::default(), ClipMate::default()],
    ];
    let read1 = [Vec::new(), Vec::new(), Vec::new()];
    let mut sam = String::new();

    let bytes = readalign_outputtranscriptsam_l5_readalign_outputtranscriptsam(
        &tr,
        0,
        0,
        0,
        0,
        -1,
        5,
        Some(&[true, false]),
        &mut sam,
        b'Y',
        "@read2",
        &["AAAA".to_string(), "CCCC".to_string()],
        1,
        &["IIII".to_string(), "HHHH".to_string()],
        &p,
        0,
        &["XA:i:1".to_string(), "XB:i:2".to_string()],
        4,
        &[4, 4],
        &[4, 4],
        &clip_mates,
        &read1,
        &genome,
    )
    .unwrap();

    assert_eq!(bytes as usize, sam.len());
    assert_eq!(
        sam,
        "read2\t933\t*\t0\t0\t*\tchr1\t21\t0\tCCCC\t*\tNH:i:0\tHI:i:0\tAS:i:30\tnM:i:2\tuT:A:5\tRG:Z:rgA\tXB:i:2\n"
    );
}

#[test]
fn read_align_write_sam_outputs_mapped_and_kept_unmapped_mate() {
    let p = Parameters {
        read_nmates: 2,
        out_sam_bool: true,
        out_sam_mapq_unique: 255,
        out_sam_flag_and: 0xffff,
        out_sam_attr_ih_start: 1,
        out_sam_attr_order: vec![ATTR_NH, ATTR_HI, ATTR_AS],
        out_sam_mode: "Full".to_string(),
        out_sam_mult_nmax: 10,
        out_sam_unmapped_keep_pairs: true,
        ..Default::default()
    };
    let genome = Genome {
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![100],
        ..Default::default()
    };
    let tr = Transcript {
        exons: vec![[0, 110, 4, 0, 0]],
        n_exons: 1,
        chr: 0,
        str_: 0,
        primary_flag: true,
        max_score: 20,
        ..Default::default()
    };
    let read_align = ReadAlign {
        l_read: 9,
        read_length: vec![4, 4],
        read_length_original: vec![4, 4],
        read_name: "@pair".to_string(),
        read_filter: b'N' as i32,
        clip_mates: vec![
            vec![ClipMate::default(), ClipMate::default()],
            vec![ClipMate::default(), ClipMate::default()],
        ],
        ..Default::default()
    };

    let result = readalign_outputalignments_l132_readalign_writesam(
        1,
        &[tr.clone()],
        &tr,
        &read_align,
        &p,
        &genome,
        -1,
        true,
        &["ACGT".to_string(), "TTTT".to_string()],
        2,
        &["IIII".to_string(), "HHHH".to_string()],
        &[String::new(), String::new()],
    )
    .unwrap();

    assert_eq!(result.unmap_type, 4);
    assert_eq!(result.mate_mapped, [true, false]);
    assert_eq!(result.bam_requests.len(), 0);
    assert_eq!(
        result.sam,
        "pair\t73\tchr1\t11\t255\t4M\t*\t0\t0\tACGT\tIIII\tNH:i:1\tHI:i:1\tAS:i:20\npair\t133\t*\t0\t0\t*\tchr1\t11\t0\tTTTT\tHHHH\tNH:i:0\tHI:i:0\tAS:i:20\tnM:i:0\tuT:A:4\n"
    );
    assert_eq!(result.out_bam_bytes as usize, result.sam.len());
}

#[test]
fn read_align_write_sam_applies_added_reference_filter_and_records_bam_requests() {
    let mut p = Parameters {
        read_nmates: 1,
        out_sam_bool: true,
        out_sam_mapq_unique: 255,
        out_sam_flag_and: 0xffff,
        out_sam_attr_ih_start: 1,
        out_sam_attr_order: vec![ATTR_NH],
        out_sam_mode: "Full".to_string(),
        out_sam_mult_nmax: 10,
        out_sam_filter_yes: true,
        out_sam_filter_keep_only_added_references: true,
        genome_insert_chr_ind_first: 2,
        ..Default::default()
    };
    let genome = Genome {
        n_chr_real: 3,
        chr_name: vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        chr_start: vec![100, 200, 300],
        ..Default::default()
    };
    let tr = Transcript {
        exons: vec![[0, 110, 4, 0, 0]],
        n_exons: 1,
        chr: 1,
        primary_flag: true,
        ..Default::default()
    };
    let read_align = ReadAlign {
        l_read: 4,
        read_length: vec![4],
        read_length_original: vec![4],
        read_name: "@r".to_string(),
        clip_mates: vec![vec![ClipMate::default(), ClipMate::default()]],
        ..Default::default()
    };

    let filtered = readalign_outputalignments_l132_readalign_writesam(
        1,
        &[tr.clone()],
        &tr,
        &read_align,
        &p,
        &genome,
        -1,
        true,
        &["ACGT".to_string()],
        1,
        &["IIII".to_string()],
        &[String::new()],
    )
    .unwrap();
    assert!(filtered.sam.is_empty());

    p.out_sam_filter_keep_only_added_references = false;
    p.out_sam_filter_keep_all_added_references = true;
    p.out_bam_unsorted = true;
    let added = Transcript {
        chr: 2,
        exons: vec![[0, 310, 4, 0, 0]],
        n_exons: 1,
        primary_flag: false,
        ..Default::default()
    };
    let written = readalign_outputalignments_l132_readalign_writesam(
        2,
        &[tr, added.clone()],
        &added,
        &read_align,
        &p,
        &genome,
        -1,
        true,
        &["ACGT".to_string()],
        1,
        &["IIII".to_string()],
        &[String::new()],
    )
    .unwrap();
    assert_eq!(written.bam_requests.len(), 1);
    assert_eq!(written.bam_requests[0].transcript.chr, 2);
    assert_eq!(written.bam_requests[0].n_tr_out, 1);
    assert_eq!(
        written.sam,
        "r\t0\tchr3\t11\t255\t4M\t*\t0\t0\tACGT\t*\tNH:i:1\n"
    );
}

#[test]
fn read_align_aligned_annotation_dispatches_enabled_quantifiers() {
    let mut transcriptome = Transcriptome {
        n_tr: 1,
        n_ge: 3,
        tr_s: vec![30],
        tr_e: vec![80],
        tr_e_max: vec![80],
        tr_ex_n: vec![1],
        tr_ex_i: vec![0],
        tr_str: vec![1],
        tr_gene: vec![0],
        ex_se: vec![0, 50],
        ex_len_cum: vec![0],
        ex_g: star_rs::generated::structs::TranscriptomeExG {
            n_ex: 1,
            s: vec![70],
            e: vec![90],
            e_max: vec![90],
            str_: vec![1],
            g: vec![2],
            t: vec![0],
        },
        gene_full: TranscriptomeGeneFull {
            s: vec![30, 70, 100],
            e: vec![80, 90, 120],
            e_max: vec![80, 90, 120],
            str_: vec![1, 1, 1],
            g: vec![0, 1, 2],
        },
        quants: quantifications_l3_quantifications_quantifications(3),
        ..Default::default()
    };
    let align = Transcript {
        exons: vec![[0, 35, 10, 0, 0]],
        n_exons: 1,
        str_: 0,
        canon_sj: vec![],
        l_read: 10,
        ..Default::default()
    };
    let converted = ReadAlignGenomeTransformResult {
        al_n: 1,
        al_mult: vec![Transcript {
            exons: vec![[0, 75, 5, 0, 0]],
            n_exons: 1,
            str_: 0,
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut read_annot = ReadAnnotations {
        annot_features: vec![ReadAnnotFeature::default(); 8],
        ..Default::default()
    };
    let mut feature_yes = vec![false; 8];
    feature_yes[SOLO_FEATURE_VELOCYTO as usize] = true;
    let p = Parameters {
        quant_ge_count_yes: true,
        quant_gene_full_yes: true,
        quant_gene_yes: true,
        quant_gene_full_exon_over_intron_yes: true,
        quant_gene_full_ex50p_as_yes: true,
        p_ge: ParametersGenome {
            transform: ParametersGenomeTransform {
                out_quant: true,
                ..Default::default()
            },
            ..Default::default()
        },
        p_solo: ParametersSolo {
            strand: -1,
            feature_yes,
            ..Default::default()
        },
        ..Default::default()
    };

    readalign_outputalignments_l298_readalign_alignedannotation(
        &mut transcriptome,
        &p,
        1,
        &[align],
        &converted,
        &mut read_annot,
    );

    assert_eq!(read_annot.gene_exon_overlap, vec![2, 2, -1]);
    assert_eq!(transcriptome.quants.gene_counts.g_count[0][2], 1);
    assert_eq!(transcriptome.quants.gene_counts.g_count[1][2], 1);
    assert_eq!(
        read_annot.annot_features[SOLO_FEATURE_GENE_FULL as usize].f_set,
        std::collections::BTreeSet::from([0])
    );
    assert_eq!(
        read_annot.annot_features[SOLO_FEATURE_GENE as usize].f_set,
        std::collections::BTreeSet::from([0])
    );
    assert_eq!(
        read_annot.annot_features[SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON as usize].f_set,
        std::collections::BTreeSet::from([0])
    );
    assert_eq!(
        read_annot.annot_features[SOLO_FEATURE_GENE_FULL_EX50P_AS as usize].f_set,
        std::collections::BTreeSet::from([0])
    );
    assert_eq!(read_annot.transcript_concordant, vec![[0, 36]]);
}

#[test]
fn read_align_output_alignments_orchestrates_mapped_sam_stats_and_sj() {
    let mut p = Parameters {
        read_nmates: 1,
        read_nends: 1,
        out_sam_bool: true,
        out_sam_mapq_unique: 255,
        out_sam_flag_and: 0xffff,
        out_sam_attr_ih_start: 1,
        out_sam_attr_order: vec![ATTR_NH, ATTR_HI, ATTR_AS],
        out_sam_mode: "Full".to_string(),
        out_sam_mult_nmax: 10,
        out_sj: true,
        ..Default::default()
    };
    let genome = Genome {
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![100],
        p_ge: ParametersGenome {
            g_type_string: "Full".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    let tr = Transcript {
        exons: vec![[0, 110, 4, 0, 0], [4, 130, 3, 0, 0]],
        n_exons: 2,
        canon_sj: vec![0],
        sj_annot: vec![1],
        chr: 0,
        str_: 0,
        primary_flag: true,
        max_score: 24,
        ..Default::default()
    };
    let mut read_align = ReadAlign {
        unmap_type: -1,
        n_tr: 1,
        tr_best: tr.clone(),
        l_read: 7,
        read_length: vec![7, 0],
        read_length_original: vec![7],
        read_name: "@mapped".to_string(),
        read_filter: b'N' as i32,
        clip_mates: vec![vec![ClipMate::default(), ClipMate::default()]],
        ..Default::default()
    };
    let mut transcriptome = Transcriptome::default();
    let mut solo_read = SoloRead::default();
    let mut chunk_out_sj = OutSJ::default();
    let mut chunk_out_sj1 = OutSJ::default();
    let mut held_fastq = vec![String::new()];
    let mut unmapped_fastx = vec![String::new()];
    let mut out_sam = String::new();

    let result = readalign_outputalignments_l5_readalign_outputalignments(
        &mut read_align,
        &mut p,
        &genome,
        &mut transcriptome,
        &mut solo_read,
        &[tr],
        &ReadAlignGenomeTransformResult::default(),
        &mut chunk_out_sj,
        &mut chunk_out_sj1,
        &mut held_fastq,
        &mut unmapped_fastx,
        &["@mapped".to_string()],
        &[String::new()],
        &["ACGTGTA".to_string()],
        &["IIIIIII".to_string()],
        2,
        &mut out_sam,
        0.0,
    )
    .unwrap();

    assert_eq!(read_align.unmap_type, -1);
    assert_eq!(read_align.stats_ra.mapped_reads_u, 1);
    assert_eq!(read_align.stats_ra.mapped_bases, 7);
    assert_eq!(read_align.stats_ra.splices_n[0], 1);
    assert_eq!(read_align.stats_ra.splices_nsjdb, 1);
    assert_eq!(chunk_out_sj.n, 1);
    assert!(chunk_out_sj1.junctions.is_empty());
    assert!(unmapped_fastx[0].is_empty());
    assert_eq!(result.out_bam_bytes as usize, out_sam.len());
    assert_eq!(
        out_sam,
        "mapped\t0\tchr1\t11\t255\t4M16N3M\t*\t0\t0\tACGTGTA\tIIIIIII\tNH:i:1\tHI:i:1\tAS:i:24\n"
    );
}

#[test]
fn read_align_output_alignments_holds_unannotated_sj_for_stage_two() {
    let mut p = Parameters {
        read_nmates: 1,
        read_nends: 1,
        out_filter_by_sjout_stage: 1,
        out_sj: true,
        out_reads_unmapped: "Fastx".to_string(),
        ..Default::default()
    };
    let genome = Genome {
        p_ge: ParametersGenome {
            g_type_string: "Full".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    let tr = Transcript {
        exons: vec![[0, 110, 4, 0, 0], [4, 130, 3, 0, 0]],
        n_exons: 2,
        canon_sj: vec![0],
        sj_annot: vec![0],
        chr: 0,
        ..Default::default()
    };
    let mut read_align = ReadAlign {
        unmap_type: -1,
        n_tr: 1,
        tr_best: tr.clone(),
        stats_ra: Stats {
            read_n: 10,
            read_bases: 7,
            ..Default::default()
        },
        read_length: vec![7, 0],
        read_length_original: vec![7],
        read_name: "@held".to_string(),
        i_read_all: 42,
        read_filter: b'N' as i32,
        ..Default::default()
    };
    let mut transcriptome = Transcriptome::default();
    let mut solo_read = SoloRead::default();
    let mut chunk_out_sj = OutSJ::default();
    let mut chunk_out_sj1 = OutSJ::default();
    let mut held_fastq = vec![String::new()];
    let mut unmapped_fastx = vec![String::new()];
    let mut out_sam = String::new();

    let result = readalign_outputalignments_l5_readalign_outputalignments(
        &mut read_align,
        &mut p,
        &genome,
        &mut transcriptome,
        &mut solo_read,
        &[tr],
        &ReadAlignGenomeTransformResult::default(),
        &mut chunk_out_sj,
        &mut chunk_out_sj1,
        &mut held_fastq,
        &mut unmapped_fastx,
        &["@held".to_string()],
        &["RX:Z:1".to_string()],
        &["ACGTGTA".to_string()],
        &["IIIIIII".to_string()],
        2,
        &mut out_sam,
        0.0,
    )
    .unwrap();

    assert_eq!(read_align.unmap_type, -3);
    assert!(!result.out_filter_by_sjout_pass);
    assert_eq!(read_align.stats_ra.read_n, 9);
    assert_eq!(read_align.stats_ra.read_bases, 0);
    assert_eq!(chunk_out_sj.n, 0);
    assert_eq!(chunk_out_sj1.n, 1);
    assert_eq!(held_fastq[0], "@held 42 78 0 RX:Z:1\nACGTGTA\n+\nIIIIIII\n");
    assert!(unmapped_fastx[0].is_empty());
    assert!(out_sam.is_empty());
}

#[test]
fn read_align_chimeric_detection_old_output_writes_sam_junction_and_bam_request() {
    let mut p = Parameters {
        read_nmates: 2,
        p_ch: ParametersChimeric {
            out_bam: true,
            out_sam_old: true,
            out_junctions: true,
            ..Default::default()
        },
        out_sam_mapq_unique: 255,
        out_sam_flag_and: 0xffff,
        out_sam_attr_ih_start: 1,
        out_sam_attr_order: vec![ATTR_NH, ATTR_HI, ATTR_AS, ATTR_NM_LOWER, ATTR_RG],
        out_sam_attr_present: star_rs::generated::structs::SamAttrPresent {
            rg: true,
            ..Default::default()
        },
        out_sam_attr_rg: vec!["rgC".to_string()],
        out_sam_mode: "Full".to_string(),
        ..Default::default()
    };
    p.genome_num_to_nt = b"ACGTN".to_vec();

    let mut genome_bases = vec![4u8; 180];
    genome_bases[110..114].copy_from_slice(&[0, 1, 2, 3]);
    genome_bases[130..134].copy_from_slice(&[0, 1, 2, 3]);
    let genome = Genome {
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![100],
        g: genome_bases,
        ..Default::default()
    };
    let read_align = ReadAlign {
        l_read: 9,
        read_nmates: 2,
        read_length: vec![4, 4],
        read_length_original: vec![4, 4],
        read_length_pair_original: 9,
        read_name: "@chim".to_string(),
        read_filter: b'N' as i32,
        read_files_index: 0,
        read1: [vec![0, 1, 2, 3], Vec::new(), vec![0, 1, 2, 3]],
        clip_mates: vec![
            vec![ClipMate::default(), ClipMate::default()],
            vec![ClipMate::default(), ClipMate::default()],
        ],
        ..Default::default()
    };
    let mut tr_chim = [
        Transcript {
            exons: vec![[0, 110, 4, 0, 0]],
            n_exons: 1,
            chr: 0,
            str_: 0,
            ro_str: 0,
            primary_flag: false,
            max_score: 0,
            n_mm: 7,
            ..Default::default()
        },
        Transcript {
            exons: vec![[0, 130, 4, 1, 0]],
            n_exons: 1,
            chr: 0,
            str_: 0,
            ro_str: 0,
            primary_flag: false,
            max_score: 0,
            n_mm: 7,
            ..Default::default()
        },
    ];

    let result = readalign_chimericdetectionoldoutput_l5_readalign_chimericdetectionoldoutput(
        true,
        &mut tr_chim,
        &read_align,
        &p,
        &genome,
        &["ACGT".to_string(), "ACGT".to_string()],
        2,
        &["IIII".to_string(), "HHHH".to_string()],
        &[String::new(), String::new()],
        2,
        113,
        130,
        1,
        2,
        3,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0.0,
    )
    .unwrap();

    assert_eq!(result.chim_n, 2);
    assert!(tr_chim[0].primary_flag);
    assert!(tr_chim[1].primary_flag);
    assert_eq!(tr_chim[0].max_score, 4);
    assert_eq!(tr_chim[1].max_score, 4);
    assert_eq!(result.bam_requests.len(), 1);
    assert!(result.bam_requests[0].is_best_chim_align);
    assert_eq!(
        result.chim_sam,
        "chim\t97\tchr1\t11\t3\t4M\tchr1\t31\t0\tACGT\tIIII\tNH:i:2\tHI:i:1\tAS:i:4\tnM:i:0\tRG:Z:rgC\nchim\t145\tchr1\t31\t3\t4M\tchr1\t11\t0\tACGT\tHHHH\tNH:i:2\tHI:i:2\tAS:i:4\tnM:i:0\tRG:Z:rgC\n"
    );
    assert_eq!(
        result.chim_junction,
        "chr1\t14\t+\tchr1\t31\t+\t1\t2\t3\tchim\t11\t4M\t31\t4M\trgC\n"
    );

    let skipped = readalign_chimericdetectionoldoutput_l5_readalign_chimericdetectionoldoutput(
        false,
        &mut tr_chim,
        &read_align,
        &p,
        &genome,
        &["ACGT".to_string(), "ACGT".to_string()],
        2,
        &["IIII".to_string(), "HHHH".to_string()],
        &[String::new(), String::new()],
        2,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0.0,
    )
    .unwrap();
    assert_eq!(skipped, Default::default());
}

#[test]
fn chimeric_align_bam_output_models_encompassing_and_supplementary_requests() {
    let p = Parameters {
        p_ch: ParametersChimeric {
            out_bam_hard_clip: true,
            ..Default::default()
        },
        ..Default::default()
    };
    let genome = Genome {
        n_chr_real: 2,
        chr_start: vec![100, 1000],
        ..Default::default()
    };
    let read_align = ReadAlign::default();
    let mate0 = Transcript {
        exons: vec![[0, 120, 8, 0, 0]],
        n_exons: 1,
        chr: 0,
        str_: 0,
        max_score: 30,
        ..Default::default()
    };
    let mate1 = Transcript {
        exons: vec![[8, 1040, 7, 1, 0]],
        n_exons: 1,
        chr: 1,
        str_: 1,
        max_score: 25,
        ..Default::default()
    };

    let encompassing = chimericalign_chimericbamoutput_l7_chimericalign_chimericbamoutput(
        &mate0,
        &mate1,
        &read_align,
        &genome,
        3,
        2,
        true,
        &p,
    );

    assert_eq!(encompassing.chim_represent, -1);
    assert_eq!(encompassing.chim_type, 2);
    assert_eq!(encompassing.bam_requests.len(), 2);
    assert_eq!(encompassing.bam_requests[0].align_type, -10);
    assert_eq!(encompassing.bam_requests[0].mate_chr, 1);
    assert_eq!(encompassing.bam_requests[0].mate_start, 40);
    assert_eq!(encompassing.bam_requests[0].mate_strand, 0);
    assert!(encompassing.bam_requests[0].transcript.primary_flag);
    assert_eq!(encompassing.bam_requests[1].mate_chr, 0);
    assert_eq!(encompassing.bam_requests[1].mate_start, 20);

    let split = Transcript {
        exons: vec![[0, 120, 5, 0, 0], [10, 180, 6, 1, 0]],
        n_exons: 2,
        chr: 0,
        str_: 0,
        max_score: 50,
        ..Default::default()
    };
    let supplementary = Transcript {
        exons: vec![[0, 1050, 9, 1, 0]],
        n_exons: 1,
        chr: 1,
        str_: 0,
        max_score: 20,
        ..Default::default()
    };

    let spanning = chimericalign_chimericbamoutput_l7_chimericalign_chimericbamoutput(
        &split,
        &supplementary,
        &read_align,
        &genome,
        4,
        7,
        false,
        &p,
    );

    assert_eq!(spanning.chim_represent, 0);
    assert_eq!(spanning.chim_type, 1);
    assert_eq!(spanning.representative_request_index, 0);
    assert_eq!(spanning.supplementary_request_index, 1);
    assert_eq!(spanning.bam_requests[0].align_type, -10);
    assert_eq!(spanning.bam_requests[0].n_tr_out, 7);
    assert_eq!(spanning.bam_requests[0].i_tr_out, 4);
    assert_eq!(spanning.bam_requests[1].align_type, -11);
    assert_eq!(spanning.bam_requests[1].mate_chr, 0);
    assert_eq!(spanning.bam_requests[1].mate_start, 20);
    assert_eq!(spanning.bam_requests[1].mate_strand, 0);
    assert!(!spanning.bam_requests[1].transcript.primary_flag);
}

#[test]
fn read_align_chimeric_detection_orchestrates_old_and_multimap_detectors() {
    let mut disabled = ReadAlign {
        chim_record: true,
        stats_ra: Stats {
            chimeric_all: 2,
            ..Default::default()
        },
        ..Default::default()
    };
    let no_request = readalign_chimericdetection_l16_readalign_chimericdetection(
        &mut disabled,
        &Parameters::default(),
        &Genome::default(),
        Some(true),
    )
    .unwrap();
    assert!(!disabled.chim_record);
    assert!(no_request.request.is_none());
    assert_eq!(disabled.stats_ra.chimeric_all, 2);

    let mut old = ReadAlign {
        n_w: 3,
        read_length: vec![50, 50],
        tr_best: Transcript {
            max_score: 80,
            ..Default::default()
        },
        ..Default::default()
    };
    let old_result = readalign_chimericdetection_l16_readalign_chimericdetection(
        &mut old,
        &Parameters {
            p_ch: ParametersChimeric {
                segment_min: 12,
                multimap_nmax: 0,
                ..Default::default()
            },
            ..Default::default()
        },
        &Genome::default(),
        Some(true),
    )
    .unwrap();
    assert!(old.chim_record);
    assert_eq!(old.stats_ra.chimeric_all, 1);
    assert!(old_result.old_output_requested);
    let old_request = old_result.request.unwrap();
    assert_eq!(old_request.detector, "chimericDetectionOld");
    assert_eq!(old_request.n_w, 3);
    assert_eq!(old_request.read_length, vec![50, 50]);

    let mut mult = ReadAlign {
        n_w: 4,
        read_length: vec![40, 40],
        tr_best: Transcript {
            max_score: 60,
            ..Default::default()
        },
        ..Default::default()
    };
    let mult_result = readalign_chimericdetection_l16_readalign_chimericdetection(
        &mut mult,
        &Parameters {
            p_ch: ParametersChimeric {
                segment_min: 12,
                multimap_nmax: 20,
                nonchim_score_drop_min: 15,
                ..Default::default()
            },
            ..Default::default()
        },
        &Genome::default(),
        Some(false),
    )
    .unwrap();
    assert!(!mult.chim_record);
    assert_eq!(mult.stats_ra.chimeric_all, 0);
    let mult_request = mult_result.request.unwrap();
    assert_eq!(mult_request.detector, "chimericDetectionMult");
    assert_eq!(mult_request.max_non_chim_align_score, 60);

    let read_length = vec![10, 10];
    let mut real_mult = ReadAlign {
        n_w: 2,
        n_win_tr: vec![1, 1],
        read_length: read_length.clone(),
        read1: [vec![0; 21], vec![0; 21], Vec::new()],
        tr_best: Transcript {
            max_score: 0,
            ..Default::default()
        },
        tr_all: vec![
            vec![Transcript {
                exons: vec![[0, 100, 6, 0, 0]],
                n_exons: 1,
                l_read: 21,
                read_length: read_length.clone(),
                read_length_original: read_length.clone(),
                read_length_pair_original: 21,
                read_nmates: 2,
                r_length: 6,
                max_score: 8,
                chr: 0,
                str_: 0,
                ro_start: 0,
                ..Default::default()
            }],
            vec![Transcript {
                exons: vec![[11, 200, 6, 1, 0]],
                n_exons: 1,
                l_read: 21,
                read_length: read_length.clone(),
                read_length_original: read_length.clone(),
                read_length_pair_original: 21,
                read_nmates: 2,
                r_length: 6,
                max_score: 9,
                chr: 1,
                str_: 0,
                ro_start: 11,
                ..Default::default()
            }],
        ],
        ..Default::default()
    };
    let real_mult_result = readalign_chimericdetection_l16_readalign_chimericdetection(
        &mut real_mult,
        &Parameters {
            p_ch: ParametersChimeric {
                segment_min: 3,
                junction_overhang_min: 2,
                score_min: 10,
                score_drop_max: 100,
                multimap_score_range: 0,
                multimap_nmax: 10,
                out_junctions: true,
                ..Default::default()
            },
            ..Default::default()
        },
        &Genome {
            chr_name: vec!["chr1".to_string(), "chr2".to_string()],
            chr_start: vec![0, 0],
            g: vec![0; 256],
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert!(real_mult.chim_record);
    assert_eq!(real_mult.stats_ra.chimeric_all, 1);
    let mult_output = real_mult_result.mult_output.unwrap();
    assert!(mult_output.chim_record);
    assert_eq!(mult_output.chim_n, 1);
    assert!(mult_output.chim_junction.contains("chr1\t107\t+"));

    let mut high_score = ReadAlign {
        read_length: vec![40, 40],
        tr_best: Transcript {
            max_score: 70,
            ..Default::default()
        },
        ..Default::default()
    };
    let skipped = readalign_chimericdetection_l16_readalign_chimericdetection(
        &mut high_score,
        &Parameters {
            p_ch: ParametersChimeric {
                segment_min: 12,
                multimap_nmax: 20,
                nonchim_score_drop_min: 15,
                ..Default::default()
            },
            ..Default::default()
        },
        &Genome::default(),
        Some(true),
    )
    .unwrap();
    assert!(skipped.request.is_none());
    assert!(!high_score.chim_record);
}

#[test]
fn read_align_chimeric_detection_pe_merged_models_old_and_new_paths() {
    let tr_init = Transcript {
        read_length: vec![10, 10],
        l_read: 21,
        ..Default::default()
    };
    let se1 = Transcript {
        exons: vec![[0, 100, 5, 0, 0], [12, 200, 5, 0, 0]],
        canon_sj: vec![1, 0],
        sj_annot: vec![1, 0],
        sj_str: vec![2, 0],
        shift_sj: vec![[3, 4], [0, 0]],
        n_exons: 2,
        l_read: 21,
        str_: 0,
        ..Default::default()
    };
    let se2 = Transcript {
        exons: vec![[0, 300, 6, 0, 0], [12, 400, 4, 0, 0]],
        canon_sj: vec![2, 0],
        sj_annot: vec![1, 0],
        sj_str: vec![1, 0],
        shift_sj: vec![[5, 6], [0, 0]],
        n_exons: 2,
        l_read: 21,
        str_: 0,
        ..Default::default()
    };
    let mut se_ra = ReadAlign {
        n_w: 2,
        read_length: vec![21, 0],
        tr_best: Transcript {
            max_score: 31,
            ..Default::default()
        },
        tr_chim: vec![se1, se2],
        ..Default::default()
    };
    let mut pe_ra = ReadAlign {
        tr_best: tr_init,
        pe_ov: ReadAlignPeOverlap {
            mate_start: [0, 10],
            ..Default::default()
        },
        read_length: vec![10, 10],
        read_length_original: vec![10, 10],
        ..Default::default()
    };

    let old = readalign_chimericdetectionpemerged_l5_readalign_chimericdetectionpemerged(
        &mut pe_ra,
        &mut se_ra,
        &Parameters {
            p_ch: ParametersChimeric {
                segment_min: 12,
                multimap_nmax: 0,
                ..Default::default()
            },
            ..Default::default()
        },
        &Genome::default(),
        Some(true),
    )
    .unwrap();
    assert!(old.mult_map_select_requested);
    assert!(old.mapped_filter_requested);
    assert!(old.old_output_requested);
    assert!(old.chim_record);
    assert_eq!(pe_ra.stats_ra.chimeric_all, 1);
    assert_eq!(
        old.request.as_ref().unwrap().detector,
        "chimericDetectionOld"
    );
    assert_eq!(old.pe_tr_chim.len(), 2);
    assert_eq!(pe_ra.tr_chim[0].n_exons, 2);

    let mut no_chim = pe_ra.clone();
    no_chim.stats_ra.chimeric_all = 0;
    let skipped_old = readalign_chimericdetectionpemerged_l5_readalign_chimericdetectionpemerged(
        &mut no_chim,
        &mut se_ra,
        &Parameters {
            p_ch: ParametersChimeric {
                segment_min: 12,
                multimap_nmax: 0,
                ..Default::default()
            },
            ..Default::default()
        },
        &Genome::default(),
        Some(false),
    )
    .unwrap();
    assert!(!skipped_old.chim_record);
    assert!(!skipped_old.old_output_requested);
    assert_eq!(no_chim.stats_ra.chimeric_all, 0);

    let mut mult = ReadAlign {
        read_length: vec![10, 10],
        tr_best: Transcript {
            max_score: 12,
            ..Default::default()
        },
        ..Default::default()
    };
    let mult_result = readalign_chimericdetectionpemerged_l5_readalign_chimericdetectionpemerged(
        &mut mult,
        &mut se_ra,
        &Parameters {
            p_ch: ParametersChimeric {
                segment_min: 12,
                multimap_nmax: 8,
                nonchim_score_drop_min: 5,
                ..Default::default()
            },
            ..Default::default()
        },
        &Genome::default(),
        Some(true),
    )
    .unwrap();
    assert!(mult_result.chim_record);
    assert_eq!(mult.stats_ra.chimeric_all, 1);
    assert_eq!(
        mult_result.request.as_ref().unwrap().detector,
        "chimericDetectionMult"
    );
    assert_eq!(mult_result.request.unwrap().max_non_chim_align_score, 31);

    let read_length = vec![10, 10];
    let mut se_real = ReadAlign {
        n_w: 2,
        n_win_tr: vec![1, 1],
        read_length: read_length.clone(),
        read1: [vec![0; 21], vec![0; 21], Vec::new()],
        tr_best: Transcript {
            max_score: 0,
            ..Default::default()
        },
        tr_all: vec![
            vec![Transcript {
                exons: vec![[0, 100, 6, 0, 0]],
                n_exons: 1,
                l_read: 21,
                read_length: read_length.clone(),
                read_length_original: read_length.clone(),
                read_length_pair_original: 21,
                read_nmates: 2,
                r_length: 6,
                max_score: 8,
                chr: 0,
                str_: 0,
                ro_start: 0,
                ..Default::default()
            }],
            vec![Transcript {
                exons: vec![[11, 200, 6, 1, 0]],
                n_exons: 1,
                l_read: 21,
                read_length: read_length.clone(),
                read_length_original: read_length.clone(),
                read_length_pair_original: 21,
                read_nmates: 2,
                r_length: 6,
                max_score: 9,
                chr: 1,
                str_: 0,
                ro_start: 11,
                ..Default::default()
            }],
        ],
        ..Default::default()
    };
    let mut pe_real = ReadAlign {
        read_length: read_length.clone(),
        read_length_original: read_length.clone(),
        tr_best: Transcript {
            max_score: 0,
            ..Default::default()
        },
        pe_ov: ReadAlignPeOverlap {
            mate_start: [0, 10],
            ..Default::default()
        },
        ..Default::default()
    };
    let real_pe_mult = readalign_chimericdetectionpemerged_l5_readalign_chimericdetectionpemerged(
        &mut pe_real,
        &mut se_real,
        &Parameters {
            p_ch: ParametersChimeric {
                segment_min: 3,
                junction_overhang_min: 2,
                score_min: 10,
                score_drop_max: 100,
                multimap_score_range: 0,
                multimap_nmax: 10,
                out_junctions: true,
                ..Default::default()
            },
            ..Default::default()
        },
        &Genome {
            chr_name: vec!["chr1".to_string(), "chr2".to_string()],
            chr_start: vec![0, 0],
            g: vec![0; 256],
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert!(pe_real.chim_record);
    assert_eq!(pe_real.stats_ra.chimeric_all, 1);
    let pe_mult_output = real_pe_mult.mult_output.unwrap();
    assert!(pe_mult_output.chim_record);
    assert_eq!(pe_mult_output.chim_n, 1);
    assert!(pe_mult_output.chim_junction.contains("chr1\t107\t+"));
}

#[test]
fn read_align_transform_genome_converts_marks_primary_and_removes_diploid_duplicates() {
    let mut tr_mult = vec![
        Transcript {
            exons: vec![[0, 10, 5, 0, 0]],
            n_exons: 1,
            chr: 0,
            str_: 1,
            max_score: 30,
            ..Default::default()
        },
        Transcript {
            exons: vec![[0, 110, 5, 0, 0]],
            n_exons: 1,
            chr: 3,
            str_: 1,
            max_score: 25,
            ..Default::default()
        },
    ];
    let read_align = ReadAlign {
        n_tr: 2,
        tr_best: tr_mult[0].clone(),
        ..Default::default()
    };
    let genome = Genome {
        n_chr_real: 4,
        genome_out: GenomeOut {
            conv_yes: true,
            conv_blocks: vec![[0, 100, 0], [100, 100, 0], [u64::MAX, 0, 0]],
            ..Default::default()
        },
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 4,
            transform: ParametersGenomeTransform {
                type_: 2,
                ..Default::default()
            },
            ..Default::default()
        },
        chr_bin: vec![0; 16],
        ..Default::default()
    };

    let transformed = readalign_transformgenome_l5_readalign_transformgenome(
        &read_align,
        &genome,
        &mut tr_mult,
        10,
        30,
        false,
        false,
        "",
        &[],
    );

    assert_eq!(tr_mult[0].haplo_type, 1);
    assert_eq!(tr_mult[1].haplo_type, 2);
    assert_eq!(transformed.al_n, 1);
    assert_eq!(transformed.al_mult[0].exons, vec![[0, 10, 5, 0, 0]]);
    assert_eq!(transformed.al_mult[0].haplo_type, 0);
    assert!(transformed.al_mult[0].primary_flag);
    assert_eq!(transformed.al_best.exons, vec![[0, 10, 5, 0, 0]]);

    let no_conversion = readalign_transformgenome_l5_readalign_transformgenome(
        &ReadAlign {
            n_tr: 1,
            tr_best: tr_mult[0].clone(),
            ..Default::default()
        },
        &Genome::default(),
        &mut tr_mult,
        10,
        30,
        false,
        false,
        "",
        &[],
    );
    assert_eq!(no_conversion.al_best, tr_mult[0]);
    assert_eq!(no_conversion.al_n, 0);
}

#[test]
fn read_align_splice_graph_write_sam_appends_all_transcripts_and_keeps_conv_noop() {
    let p = Parameters {
        out_sam_mapq_unique: 255,
        out_sam_flag_and: 0xffff,
        out_sam_attr_ih_start: 1,
        out_sam_attr_order: vec![ATTR_NH, ATTR_HI, ATTR_AS],
        out_sam_mode: "Full".to_string(),
        ..Default::default()
    };
    let genome = Genome {
        chr_name: vec!["chrA".to_string()],
        chr_start: vec![100],
        ..Default::default()
    };
    let tr_mult = vec![
        Transcript {
            cigar: vec![[0, 4]],
            chr: 0,
            g_start: 110,
            primary_flag: true,
            max_score: 20,
            ..Default::default()
        },
        Transcript {
            cigar: vec![[0, 3], [4, 1]],
            chr: 0,
            g_start: 130,
            primary_flag: false,
            max_score: 10,
            ..Default::default()
        },
    ];
    let read0 = vec!["ACGT".to_string()];
    let qual0 = vec!["abcd".to_string()];
    let read_name_extra = vec![String::new()];
    let mut out = String::new();

    let bytes = readalign_outputalignments_l277_readalign_splicegraphwritesam(
        &tr_mult,
        2,
        &mut out,
        b'N',
        -1,
        "@sg",
        &read0,
        2,
        &qual0,
        &p,
        0,
        &read_name_extra,
        4,
        &genome,
    )
    .unwrap();

    assert_eq!(bytes as usize, out.len());
    assert_eq!(
        out,
        "sg\t0\tchrA\t11\t3\t4M\t*\t0\t0\tACGT\tabcd\tNH:i:2\tHI:i:1\tAS:i:20\nsg\t256\tchrA\t31\t3\t3M1S\t*\t0\t0\tACGT\tabcd\tNH:i:2\tHI:i:2\tAS:i:10\n"
    );

    let mut converted_genome = genome.clone();
    converted_genome.genome_out.conv_yes = true;
    out.clear();
    let bytes = readalign_outputalignments_l277_readalign_splicegraphwritesam(
        &tr_mult,
        2,
        &mut out,
        b'N',
        -1,
        "@sg",
        &read0,
        2,
        &qual0,
        &p,
        0,
        &read_name_extra,
        4,
        &converted_genome,
    )
    .unwrap();
    assert_eq!(bytes, 0);
    assert_eq!(out, "");
}

#[test]
fn read_align_assign_align_to_window_records_replaces_and_prunes_like_original() {
    let mut ra = ReadAlign {
        l_read: 100,
        n_w: 1,
        win_bin: [vec![0], vec![UINT_WIN_BIN_MAX]],
        wa: vec![vec![[0; WA_SIZE]; 4]],
        n_wa: vec![0],
        n_wap: vec![0],
        wal_rec: vec![0],
        w_last_anchor: vec![0],
        ..Default::default()
    };

    readalign_assignaligntowindow_l6_readalign_assignaligntowindow(
        &mut ra, 100, 10, 0, 1, 0, 20, false, 0, 10, 3,
    )
    .unwrap();
    readalign_assignaligntowindow_l6_readalign_assignaligntowindow(
        &mut ra, 90, 5, 0, 1, 0, 5, false, 0, 10, 3,
    )
    .unwrap();
    assert_eq!(ra.n_wa[0], 2);
    assert_eq!(ra.wa[0][0][WA_R_START], 5);
    assert_eq!(ra.wa[0][1][WA_R_START], 20);
    assert_eq!(ra.n_wap[0], 2);

    readalign_assignaligntowindow_l6_readalign_assignaligntowindow(
        &mut ra, 102, 15, 0, 2, 0, 22, false, 0, 10, 3,
    )
    .unwrap();
    assert_eq!(ra.n_wa[0], 2);
    assert_eq!(ra.wa[0][1][WA_R_START], 22);
    assert_eq!(ra.wa[0][1][WA_LENGTH], 15);
    assert_eq!(ra.wa[0][1][WA_N_REP], 2);

    readalign_assignaligntowindow_l6_readalign_assignaligntowindow(
        &mut ra, 130, 3, 0, 1, 0, 60, false, 0, 10, 3,
    )
    .unwrap();
    assert_eq!(ra.n_wa[0], 3);
    readalign_assignaligntowindow_l6_readalign_assignaligntowindow(
        &mut ra, 150, 7, 0, 1, 0, 80, false, 0, 10, 3,
    )
    .unwrap();
    assert_eq!(ra.wal_rec[0], 3);
    assert_eq!(ra.n_wa[0], 3);
    assert_eq!(ra.wa[0][0][WA_LENGTH], 5);
    assert_eq!(ra.wa[0][1][WA_LENGTH], 15);
    assert_eq!(ra.wa[0][2][WA_LENGTH], 7);
    assert_eq!(ra.n_wap[0], 4);
}

#[test]
fn read_align_assign_align_to_window_handles_anchor_overflow_marker() {
    let mut ra = ReadAlign {
        l_read: 30,
        n_w: 1,
        win_bin: [vec![0], vec![UINT_WIN_BIN_MAX]],
        wa: vec![vec![[0; WA_SIZE]; 3]],
        n_wa: vec![2],
        n_wap: vec![2],
        wal_rec: vec![0],
        w_last_anchor: vec![1],
        ..Default::default()
    };
    ra.wa[0][0][WA_LENGTH] = 5;
    ra.wa[0][0][WA_ANCHOR] = 1;
    ra.wa[0][1][WA_LENGTH] = 6;
    ra.wa[0][1][WA_ANCHOR] = 1;

    readalign_assignaligntowindow_l6_readalign_assignaligntowindow(
        &mut ra, 20, 7, 0, 1, 0, 3, true, 0, 10, 2,
    )
    .unwrap();
    assert_eq!(ra.map_marker, MARKER_TOO_MANY_ANCHORS_PER_WINDOW);
    assert_eq!(ra.n_w, 0);
}

#[test]
fn read_align_create_extend_windows_merges_neighbors_and_limits_window_count() {
    let p = Parameters {
        win_bin_nbits: 0,
        win_anchor_dist_nbins: 1,
        win_bin_chr_nbits: 2,
        win_bin_n: 8,
        align_windows_per_read_nmax: 4,
        ..Default::default()
    };
    let genome = Genome {
        chr_bin: vec![0, 0, 0, 0, 1, 1, 1, 1],
        ..Default::default()
    };
    let mut ra = ReadAlign {
        n_w: 0,
        win_bin: [
            [UINT_WIN_BIN_MAX; 8].to_vec(),
            [UINT_WIN_BIN_MAX; 8].to_vec(),
        ],
        wc: vec![[0; WC_SIZE]; 4],
        ..Default::default()
    };

    assert_eq!(
        readalign_createextendwindowswithalign_l7_readalign_createextendwindowswithalign(
            &mut ra, &genome, &p, 2, 0,
        ),
        0
    );
    assert_eq!(ra.n_w, 1);
    assert_eq!(ra.win_bin[0][2], 0);
    assert_eq!(ra.wc[0], [0, 0, 2, 2]);

    assert_eq!(
        readalign_createextendwindowswithalign_l7_readalign_createextendwindowswithalign(
            &mut ra, &genome, &p, 3, 0,
        ),
        0
    );
    assert_eq!(ra.n_w, 1);
    assert_eq!(ra.win_bin[0][3], 0);
    assert_eq!(ra.wc[0], [0, 0, 2, 3]);

    assert_eq!(
        readalign_createextendwindowswithalign_l7_readalign_createextendwindowswithalign(
            &mut ra, &genome, &p, 0, 0,
        ),
        0
    );
    assert_eq!(ra.n_w, 2);
    assert_eq!(ra.wc[1], [0, 0, 0, 0]);

    assert_eq!(
        readalign_createextendwindowswithalign_l7_readalign_createextendwindowswithalign(
            &mut ra, &genome, &p, 1, 0,
        ),
        0
    );
    assert_eq!(ra.n_w, 2);
    assert_eq!(ra.win_bin[0][0..=3], [1, 1, 1, 1]);
    assert_eq!(ra.wc[1], [0, 0, 0, 3]);
    assert_eq!(ra.wc[0][WC_G_START], 1);
    assert_eq!(ra.wc[0][WC_G_END], 0);

    let p_limited = Parameters {
        win_bin_nbits: 0,
        win_anchor_dist_nbins: 0,
        win_bin_chr_nbits: 2,
        win_bin_n: 8,
        align_windows_per_read_nmax: 2,
        ..Default::default()
    };
    let mut limited = ReadAlign {
        n_w: 1,
        win_bin: [
            [UINT_WIN_BIN_MAX; 8].to_vec(),
            [UINT_WIN_BIN_MAX; 8].to_vec(),
        ],
        wc: vec![[0; WC_SIZE]; 2],
        ..Default::default()
    };
    assert_eq!(
        readalign_createextendwindowswithalign_l7_readalign_createextendwindowswithalign(
            &mut limited,
            &genome,
            &p_limited,
            5,
            0,
        ),
        EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS
    );
    assert_eq!(limited.n_w, 1);
}

#[test]
fn read_align_stitch_pieces_records_single_anchor_window() {
    let p = Parameters {
        win_bin_nbits: 0,
        win_anchor_dist_nbins: 1,
        win_anchor_multimap_nmax: 1,
        win_flank_nbins: 0,
        win_bin_chr_nbits: 8,
        win_bin_n: 64,
        align_windows_per_read_nmax: 4,
        seed_per_window_nmax: 4,
        align_transcripts_per_window_nmax: 4,
        align_transcripts_per_read_nmax: 16,
        out_filter_intron_motifs: "None".to_string(),
        out_filter_multimap_score_range: 0,
        align_soft_clip_at_reference_ends_yes: true,
        ..Default::default()
    };
    let genome = Genome {
        sa: vec![10],
        n_genome: 256,
        gstrand_bit: 31,
        gstrand_mask: 0x7fff_ffff,
        sj_gstart: 1_000,
        chr_bin: vec![0; 64],
        g: vec![0; 256],
        chr_start: vec![0],
        chr_length: vec![256],
        ..Default::default()
    };
    let mut pc = vec![[0u32; PC_SIZE]; 1];
    pc[0][PC_R_START] = 0;
    pc[0][PC_LENGTH] = 4;
    pc[0][PC_DIR] = 0;
    pc[0][PC_NREP] = 1;
    pc[0][PC_SASTART] = 0;
    pc[0][PC_SAEND] = 0;
    pc[0][PC_IFRAG] = 0;
    let mut ra = ReadAlign {
        n_p: 1,
        l_read: 4,
        pc,
        out_filter_mismatch_nmax_total: 10,
        read_length: vec![4],
        max_score_mate: vec![i32::MIN],
        tr_init: Box::new(Transcript::default()),
        ..Default::default()
    };
    let reads = vec![vec![0, 1, 2, 3], vec![3, 2, 1, 0]];

    readalign_stitchpieces_l12_readalign_stitchpieces(&mut ra, &reads, 4, &genome, &p).unwrap();

    assert_eq!(ra.map_marker, 0);
    assert_eq!(ra.n_wall, 1);
    assert_eq!(ra.n_w, 1);
    assert_eq!(ra.n_tr, 1);
    assert_eq!(ra.n_wa[0], 1);
    assert_eq!(ra.n_win_tr[0], 1);
    assert_eq!(ra.tr_all[0][0].chr, 0);
    assert_eq!(ra.tr_all[0][0].str_, 0);
    assert_eq!(ra.tr_all[0][0].ro_str, 0);
    assert_eq!(ra.tr_all[0][0].exons[0][EX_R], 0);
    assert_eq!(ra.tr_all[0][0].exons[0][EX_G], 10);
    assert_eq!(ra.tr_all[0][0].exons[0][EX_L], 4);
    assert_eq!(ra.tr_best.max_score, 4);
}

#[test]
fn read_align_stitch_pieces_converts_reverse_strand_coordinates() {
    let p = Parameters {
        win_bin_nbits: 0,
        win_anchor_dist_nbins: 1,
        win_anchor_multimap_nmax: 1,
        win_flank_nbins: 0,
        win_bin_chr_nbits: 8,
        win_bin_n: 220,
        align_windows_per_read_nmax: 4,
        seed_per_window_nmax: 4,
        align_transcripts_per_window_nmax: 4,
        align_transcripts_per_read_nmax: 16,
        out_filter_intron_motifs: "None".to_string(),
        out_filter_multimap_score_range: 0,
        align_soft_clip_at_reference_ends_yes: true,
        ..Default::default()
    };
    let genome = Genome {
        sa: vec![(1u32 << 31) | 20],
        n_genome: 200,
        gstrand_bit: 31,
        gstrand_mask: 0x7fff_ffff,
        sj_gstart: 1_000,
        chr_bin: vec![0; 220],
        g: vec![0; 220],
        chr_start: vec![0],
        chr_length: vec![220],
        ..Default::default()
    };
    let mut pc = vec![[0u32; PC_SIZE]; 1];
    pc[0][PC_R_START] = 0;
    pc[0][PC_LENGTH] = 4;
    pc[0][PC_DIR] = 0;
    pc[0][PC_NREP] = 1;
    pc[0][PC_SASTART] = 0;
    pc[0][PC_SAEND] = 0;
    pc[0][PC_IFRAG] = 0;
    let mut ra = ReadAlign {
        n_p: 1,
        l_read: 4,
        pc,
        out_filter_mismatch_nmax_total: 10,
        read_length: vec![4],
        max_score_mate: vec![i32::MIN],
        tr_init: Box::new(Transcript::default()),
        ..Default::default()
    };
    let reads = vec![vec![0, 1, 2, 3], vec![3, 2, 1, 0]];

    readalign_stitchpieces_l12_readalign_stitchpieces(&mut ra, &reads, 4, &genome, &p).unwrap();

    assert_eq!(ra.map_marker, 0);
    assert_eq!(ra.n_w, 1);
    assert_eq!(ra.wc[0][WC_STR], 1);
    assert_eq!(ra.wa[0][0][WA_G_START], 176);
    assert_eq!(ra.wa[0][0][WA_R_START], 0);
    assert_eq!(ra.tr_all[0][0].str_, 1);
    assert_eq!(ra.tr_all[0][0].ro_str, 1);
    assert_eq!(ra.tr_all[0][0].exons[0][EX_G], 176);
    assert_eq!(ra.tr_best.max_score, 4);
}

#[test]
fn read_align_one_read_loads_fastq_and_finishes_injected_standard_mapping() {
    let mut input = std::io::Cursor::new(b"@SRR000001.1 77 A 3 extra\nACGT\n+\nIIII\n".to_vec());
    let mut streams: Vec<&mut dyn std::io::BufRead> = vec![&mut input];
    let mut p = Parameters {
        read_nends: 1,
        read_nmates: 1,
        out_filter_multimap_nmax: 10,
        out_filter_mismatch_nmax: 10,
        out_filter_mismatch_nover_read_lmax: 1.0,
        out_filter_score_min: 0,
        out_filter_score_min_over_lread: 0.0,
        out_filter_match_nmin: 0,
        out_filter_match_nmin_over_lread: 0.0,
        out_filter_mismatch_nover_lmax: 1.0,
        out_filter_multimap_score_range: 0,
        out_sam_mult_nmax: 10,
        out_sam_primary_flag: "OneBestScore".to_string(),
        read_name_separator_char: vec![' '],
        ..Default::default()
    };
    let genome = Genome {
        chr_start: vec![0],
        chr_length: vec![100],
        chr_name: vec!["chr1".to_string()],
        g: vec![0; 100],
        ..Default::default()
    };
    let mut tr = Transcript {
        exons: vec![[0, 5, 4, 0, u32::MAX]],
        n_exons: 1,
        r_start: 0,
        r_length: 4,
        g_start: 5,
        g_length: 4,
        chr: 0,
        str_: 0,
        ro_str: 0,
        max_score: 4,
        n_match: 4,
        ..Default::default()
    };
    tr.read_length = vec![4, 0];
    let mapped = ReadAlign {
        n_w: 1,
        tr_best: tr.clone(),
        tr_all: vec![vec![tr]],
        n_wap: vec![1],
        n_win_tr: vec![1],
        ..Default::default()
    };
    let mut ra = ReadAlign {
        read_nmates: 1,
        read_length: vec![0, 0],
        read_length_original: vec![0, 0],
        read0: vec![Vec::new()],
        qual0: vec![Vec::new()],
        read_name_mates: vec![Vec::new()],
        read_name_extra: vec![String::new()],
        clip_mates: vec![vec![ClipMate::default(); 2]],
        qual_hist: vec![vec![0; 256]],
        max_score_mate: vec![0],
        tr_init: Box::new(Transcript::default()),
        ..Default::default()
    };
    let mut pe_merge_ra = ReadAlign::default();
    let mut wasp_ra = ReadAlign::default();
    let mut transcriptome = Transcriptome::default();
    let mut out_sj = OutSJ::default();
    let mut out_sj1 = OutSJ::default();
    let mut filter_files = Vec::<String>::new();
    let mut unmapped = vec![String::new()];
    let mut sam = String::new();

    let result = readalign_oneread_l8_readalign_oneread(
        &mut ra,
        &mut streams,
        &mut p,
        &genome,
        &mut transcriptome,
        None,
        Some(&mapped),
        None,
        &mut pe_merge_ra,
        &mut wasp_ra,
        None,
        &[],
        &mut out_sj,
        &mut out_sj1,
        &mut filter_files,
        &mut unmapped,
        &mut sam,
        0.0,
        None,
    )
    .unwrap();

    assert_eq!(result.status, 0);
    assert!(!result.map_one_read_requested);
    assert_eq!(ra.read_name, "@SRR000001.1");
    assert_eq!(ra.i_read_all, 77);
    assert_eq!(ra.read_filter, b'A' as i32);
    assert_eq!(ra.read_files_index, 3);
    assert_eq!(ra.read_name_extra[0], "extra");
    assert_eq!(ra.l_read, 4);
    assert_eq!(ra.read1[0][..4], [0, 1, 2, 3]);
    assert_eq!(ra.read1[1][..4], [3, 2, 1, 0]);
    assert_eq!(ra.read1[2][..4], [0, 1, 2, 3]);
    assert_eq!(ra.qual_hist[0][b'I' as usize], 4);
    assert_eq!(ra.out_filter_mismatch_nmax_total, 4);
    assert_eq!(ra.stats_ra.read_n, 1);
    assert_eq!(ra.stats_ra.mapped_reads_u, 1);
    assert_eq!(ra.unmap_type, -1);
    assert_eq!(ra.n_tr, 1);
    assert_eq!(result.tr_mult[0].primary_flag, true);
    assert!(result.output_alignments.is_some());
}

#[test]
fn read_align_one_read_reports_inconsistent_mate_eof() {
    let mut input1 = std::io::Cursor::new(b"@r1\nACGT\n+\nIIII\n".to_vec());
    let mut input2 = std::io::Cursor::new(Vec::<u8>::new());
    let mut streams: Vec<&mut dyn std::io::BufRead> = vec![&mut input1, &mut input2];
    let mut p = Parameters {
        read_nends: 2,
        read_nmates: 2,
        read_name_separator_char: vec![' '],
        ..Default::default()
    };
    let mut ra = ReadAlign {
        read_nmates: 2,
        read_length: vec![0, 0],
        read_length_original: vec![0, 0],
        clip_mates: vec![vec![ClipMate::default(); 2], vec![ClipMate::default(); 2]],
        qual_hist: vec![vec![0; 256], vec![0; 256]],
        tr_init: Box::new(Transcript::default()),
        ..Default::default()
    };
    let mut transcriptome = Transcriptome::default();
    let mut pe_merge_ra = ReadAlign::default();
    let mut wasp_ra = ReadAlign::default();
    let mut out_sj = OutSJ::default();
    let mut out_sj1 = OutSJ::default();
    let mut sam = String::new();
    let err = readalign_oneread_l8_readalign_oneread(
        &mut ra,
        &mut streams,
        &mut p,
        &Genome::default(),
        &mut transcriptome,
        None,
        None,
        None,
        &mut pe_merge_ra,
        &mut wasp_ra,
        None,
        &[],
        &mut out_sj,
        &mut out_sj1,
        &mut [],
        &mut [String::new(), String::new()],
        &mut sam,
        0.0,
        None,
    )
    .unwrap_err();

    assert!(err.contains("read files are not consistent"));
}

#[test]
fn variation_snp_block_helpers_match_original_window_logic() {
    let snp = SNP {
        n: 3,
        loci: vec![10, 20, 30],
        nt: vec![[0, 1, 0], [2, 2, 3], [1, 2, 3]],
        ..Default::default()
    };
    let mut out = vec![Vec::new(), Vec::new()];
    variation_l124_snp_snponblocks(&snp, 9, 15, 100, &mut out);
    assert_eq!(out[0], vec![[101, 1]]);
    assert_eq!(out[1], vec![[111, 3]]);

    let variation = Variation {
        yes: true,
        snp: SNP {
            n: 4,
            loci: vec![90, 105, 205, 210],
            nt: vec![[0, 1, 0], [0, 0, 2], [1, 2, 1], [3, 0, 1]],
            ..Default::default()
        },
    };
    assert_eq!(
        variation_l139_variation_sjdbsnp(&variation, 100, 200, 10),
        vec![vec![[0, 1], [14, 2], [19, 0]], vec![[19, 1]]]
    );

    let no_variation = variation_l8_variation_variation(false);
    assert_eq!(
        variation_l139_variation_sjdbsnp(&no_variation, 100, 200, 10),
        vec![Vec::<[i32; 2]>::new()]
    );
}

#[test]
fn transcript_variation_adjust_records_overlapping_snp_alleles() {
    let genome = Genome {
        chr_start: vec![100],
        var: Variation {
            yes: true,
            snp: SNP {
                n: 4,
                loci: vec![100, 105, 111, 130],
                nt: vec![[0, 1, 2], [2, 0, 3], [1, 2, 3], [3, 0, 1]],
                ..Default::default()
            },
        },
        ..Default::default()
    };
    let mut transcript = Transcript {
        n_exons: 2,
        chr: 0,
        exons: vec![
            {
                let mut e = [0u32; EX_SIZE];
                e[EX_R] = 0;
                e[EX_G] = 103;
                e[EX_L] = 5;
                e
            },
            {
                let mut e = [0u32; EX_SIZE];
                e[EX_R] = 10;
                e[EX_G] = 110;
                e[EX_L] = 4;
                e
            },
        ],
        ..Default::default()
    };
    let read = vec![0, 1, 3, 3, 4, 0, 0, 0, 0, 0, 0, 2, 4, 0];

    assert_eq!(
        transcript_variationadjust_l4_transcript_variationadjust(&mut transcript, &genome, &read),
        0
    );
    assert_eq!(transcript.var_ind, vec![1, 2]);
    assert_eq!(transcript.var_gen_coord, vec![5, 11]);
    assert_eq!(transcript.var_read_coord, vec![2, 11]);
    assert_eq!(transcript.var_allele, vec![2, 1]);

    let no_var_genome = Genome {
        var: Variation::default(),
        ..Default::default()
    };
    assert_eq!(
        transcript_variationadjust_l4_transcript_variationadjust(
            &mut transcript,
            &no_var_genome,
            &read
        ),
        0
    );
    assert_eq!(transcript.var_ind, vec![1, 2]);
}

#[test]
fn scan_vcf_records_supported_snv_genotypes_in_genome_coordinates() {
    let mut chr_index = std::collections::BTreeMap::new();
    chr_index.insert("chr1".to_string(), 0);
    chr_index.insert("chr2".to_string(), 1);

    let vcf = "\
##fileformat=VCFv4.2
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
chr1 11 . A C . . . GT 0/1
chr1 12 . A C . . . GT 0/0
chr1 13 . G T . . . GT 1/1
chr2 4 . C A,G . . . GT 1/2
chr2 5 . C CC . . . GT 0/1
chr3 9 . A C . . . GT 0/1
chr1 14 . N C . . . GT 0/1
chr1 15 . A A . . . GT 1/1
";
    let mut snp = SNP::default();
    let n_homoz = variation_l23_scanvcf(vcf, &mut snp, &[100, 200], &chr_index, true);
    assert_eq!(n_homoz, 1);
    assert_eq!(snp.n, 2);
    assert_eq!(snp.loci_v, vec![110, 203]);
    assert_eq!(snp.nt, vec![[0, 0, 1], [1, 0, 2]]);

    let mut snp_with_homoz = SNP::default();
    let n_homoz = variation_l23_scanvcf(vcf, &mut snp_with_homoz, &[100, 200], &chr_index, false);
    assert_eq!(n_homoz, 0);
    assert_eq!(snp_with_homoz.n, 3);
    assert_eq!(snp_with_homoz.loci_v, vec![110, 112, 203]);
    assert_eq!(snp_with_homoz.nt[1], [2, 3, 3]);
}

#[test]
fn variation_load_vcf_copies_clears_and_sorts_scanned_snps() {
    let mut chr_index = std::collections::BTreeMap::new();
    chr_index.insert("chr1".to_string(), 0);
    chr_index.insert("chr2".to_string(), 1);

    let vcf = "\
chr2 3 . C T . . . GT 0/1
chr1 9 . A G . . . GT 0/1
chr1 4 . T C . . . GT 0/1
";
    let mut variation = variation_l8_variation_variation(true);
    assert_eq!(
        variation_l81_variation_loadvcf(&mut variation, vcf, &[100, 200], &chr_index, false),
        Ok(0)
    );
    assert_eq!(variation.snp.n, 3);
    assert!(variation.snp.loci_v.is_empty());
    assert_eq!(variation.snp.loci, vec![103, 108, 202]);
    assert_eq!(variation.snp.nt, vec![[3, 3, 1], [0, 0, 2], [1, 1, 3]]);

    let empty = variation_l81_variation_loadvcf(
        &mut variation_l8_variation_variation(true),
        "chr1 1 . A C . . . GT 0/0\n",
        &[100],
        &chr_index,
        false,
    );
    assert!(empty.is_err());
}

#[test]
fn collapse_umi_with_one_mm_low_half_marks_colors_edges_and_directional_duplicates() {
    let mut umi_arr = vec![
        0b0000,
        10,
        u32::MAX,
        0,
        0b0001,
        3,
        u32::MAX,
        0,
        0b0010,
        4,
        u32::MAX,
        0,
        0b0100,
        20,
        u32::MAX,
        0,
        0b0101,
        2,
        u32::MAX,
        0,
    ];
    let mut n_u1 = 5;
    let mut n_u2 = 5;
    let mut n_c = 0;
    let mut edges = Vec::new();

    solofeature_collapseumi_graph_l80_collapseumiwith1mmlowhalf(
        &mut umi_arr,
        4,
        0b0011,
        5,
        &mut n_u1,
        &mut n_u2,
        &mut n_c,
        &mut edges,
    );

    assert_eq!(n_c, 2);
    assert_eq!(n_u1, 0);
    assert!(edges.is_empty());
    assert_eq!(umi_arr[2], 0);
    assert_eq!(umi_arr[6], 0);
    assert_eq!(umi_arr[10], 0);
    assert_eq!(umi_arr[14], 1);
    assert_eq!(umi_arr[18], 1);
    assert_eq!(umi_arr[5] & (1 << 31), 1 << 31);
    assert_eq!(umi_arr[9] & (1 << 31), 1 << 31);
    assert_eq!(umi_arr[17] & (1 << 31), 1 << 31);
    assert_eq!(n_u2, 2);
}

#[test]
fn graph_number_of_connected_components_colors_only_linked_nodes_like_original() {
    let edges = vec![[0, 1], [1, 2], [3, 4]];
    let mut node_color = Vec::new();

    let n = solofeature_collapseumi_graph_l142_graphnumberofconnectedcomponents(
        6,
        &edges,
        &mut node_color,
    );

    assert_eq!(n, 3);
    assert_eq!(node_color, vec![0, 0, 0, 3, 3, u32::MAX]);

    assert_eq!(
        solofeature_collapseumi_graph_l142_graphnumberofconnectedcomponents(
            4,
            &[],
            &mut node_color
        ),
        4
    );
    assert_eq!(node_color, vec![u32::MAX; 4]);
}

#[test]
fn solo_feature_umi_array_correct_graph_collapses_and_records_best_umi() {
    let p_solo = ParametersSolo {
        umi_l: 2,
        umi_mask_low: 0b0011,
        ..Default::default()
    };
    let mut umi_arr = vec![0b0010, 4, 123, 0, 0b0000, 10, 123, 0, 0b0001, 3, 123, 0];
    let mut umi_corr = std::collections::BTreeMap::new();

    let n_u = solofeature_collapseumi_graph_l16_solofeature_umiarraycorrect_graph(
        &p_solo,
        3,
        &mut umi_arr,
        4,
        true,
        true,
        &mut umi_corr,
    );

    assert_eq!(n_u, 1);
    assert_eq!(
        umi_arr,
        vec![0, 10, 0, 0, 1, 1 << 31 | 3, 0, 0, 2, 1 << 31 | 4, 0, 0]
    );
    assert_eq!(umi_corr.get(&0), Some(&0));
    assert_eq!(umi_corr.get(&1), Some(&0));
    assert_eq!(umi_corr.get(&2), Some(&0));
}

#[test]
fn solo_feature_umi_array_correct_graph_obeys_n_umi_flag_and_leaves_swapped_without_info() {
    let p_solo = ParametersSolo {
        umi_l: 2,
        umi_mask_low: 0b0011,
        ..Default::default()
    };
    let mut umi_arr = vec![0b0001, 2, 0, 0, 0b0100, 7, 0, 0];
    let mut umi_corr = std::collections::BTreeMap::new();

    let n_u = solofeature_collapseumi_graph_l16_solofeature_umiarraycorrect_graph(
        &p_solo,
        2,
        &mut umi_arr,
        4,
        false,
        false,
        &mut umi_corr,
    );

    assert_eq!(n_u, 0);
    assert!(umi_corr.is_empty());
    assert_eq!(umi_arr[0], 1);
    assert_eq!(umi_arr[4], 4);
}

#[test]
fn solo_feature_umi_array_correct_cr_replaces_one_mm_with_higher_count_umi() {
    let mut umi_arr = vec![0b0001, 2, 99, 0b1111, 9, 99, 0b0000, 5, 99];
    let mut umi_corr = std::collections::BTreeMap::new();

    let n_u = solofeature_collapseumiall_l580_solofeature_umiarraycorrect_cr(
        3,
        &mut umi_arr,
        3,
        true,
        true,
        &mut umi_corr,
    );

    assert_eq!(n_u, 2);
    assert_eq!(umi_arr, vec![1, 2, 0, 0, 5, 0, 15, 9, 15]);
    assert_eq!(umi_corr.get(&1), Some(&0));
}

#[test]
fn solo_feature_umi_array_correct_directional_obeys_count_threshold() {
    let mut umi_arr = vec![0b0001, 4, 99, 0b0010, 6, 99, 0b0000, 10, 99];
    let mut umi_corr = std::collections::BTreeMap::new();

    let n_u = solofeature_collapseumiall_l617_solofeature_umiarraycorrect_directional(
        3,
        &mut umi_arr,
        3,
        true,
        true,
        &mut umi_corr,
        1,
    );

    assert_eq!(n_u, 2);
    assert_eq!(umi_arr, vec![0, 10, 0, 2, 6, 2, 1, 4, 0]);
    assert_eq!(umi_corr.get(&1), Some(&0));
    assert!(!umi_corr.contains_key(&2));
}

#[test]
fn genome_parameters_write_emits_original_key_order_and_spacing() {
    let dir = std::env::temp_dir().join(format!("star-rs-genome-par-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("genomeParameters.txt");
    let _ = std::fs::remove_file(&path);

    let p = Parameters {
        command_line_full: "STAR --runMode genomeGenerate".to_string(),
        version_genome: "20201".to_string(),
        ..Default::default()
    };
    let genome = Genome {
        gstrand_bit: 33,
        sjdb_overhang: 149,
        p_ge: ParametersGenome {
            g_type_string: "Full".to_string(),
            g_fasta_files: vec!["chr1.fa".to_string(), "chr2.fa".to_string()],
            g_saindex_nbases: 14,
            g_chr_bin_nbits: 18,
            g_sasparse_d: 3,
            transform: ParametersGenomeTransform {
                type_string: "None".to_string(),
                vcf_file: "variants.vcf".to_string(),
                ..Default::default()
            },
            sjdb_file_chr_start_end: vec!["sj1.tab".to_string(), "sj2.tab".to_string()],
            sjdb_gtf_file: "genes.gtf".to_string(),
            sjdb_gtf_chr_prefix: "chr".to_string(),
            sjdb_gtf_feature_exon: "exon".to_string(),
            sjdb_gtf_tag_exon_parent_transcript: "transcript_id".to_string(),
            sjdb_gtf_tag_exon_parent_gene: "gene_id".to_string(),
            sjdb_insert_save: "Basic".to_string(),
            g_file_sizes: vec![11, 22, 33],
            ..Default::default()
        },
        ..Default::default()
    };

    genomeparameterswrite_l4_genomeparameterswrite(path.to_str().unwrap(), &p, "ERR", &genome)
        .unwrap();

    let actual = std::fs::read_to_string(&path).unwrap();
    assert_eq!(
        actual,
        concat!(
            "### STAR --runMode genomeGenerate\n",
            "### GstrandBit 33\n",
            "versionGenome\t20201\n",
            "genomeType\tFull\n",
            "genomeFastaFiles\tchr1.fa chr2.fa \n",
            "genomeSAindexNbases\t14\n",
            "genomeChrBinNbits\t18\n",
            "genomeSAsparseD\t3\n",
            "genomeTransformType\tNone\n",
            "genomeTransformVCF\tvariants.vcf\n",
            "sjdbOverhang\t149\n",
            "sjdbFileChrStartEnd\tsj1.tab sj2.tab \n",
            "sjdbGTFfile\tgenes.gtf\n",
            "sjdbGTFchrPrefix\tchr\n",
            "sjdbGTFfeatureExon\texon\n",
            "sjdbGTFtagExonParentTranscript\ttranscript_id\n",
            "sjdbGTFtagExonParentGene\tgene_id\n",
            "sjdbInsertSave\tBasic\n",
            "genomeFileSizes\t11 22 33\n",
        )
    );
}

#[test]
fn genome_constructor_sets_genome_bin_and_sjdb_lengths() {
    let genome = genome_l15_genome_genome(ParametersGenome {
        g_chr_bin_nbits: 4,
        sjdb_overhang: 10,
        ..Default::default()
    });

    assert_eq!(genome.genome_chr_bin_nbases, 16);
    assert_eq!(genome.sjdb_overhang, 10);
    assert_eq!(genome.sjdb_length, 21);

    let genome_without_sjdb = genome_l15_genome_genome(ParametersGenome {
        g_chr_bin_nbits: 2,
        sjdb_overhang: 0,
        ..Default::default()
    });
    assert_eq!(genome_without_sjdb.genome_chr_bin_nbases, 4);
    assert_eq!(genome_without_sjdb.sjdb_length, 0);
}

#[test]
fn genome_chr_bin_fill_maps_bins_to_chromosomes_with_terminal_start() {
    let mut genome = Genome {
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 4,
            ..Default::default()
        },
        chr_name: vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        chr_start: vec![0, 20, 40, 70],
        ..Default::default()
    };

    genome_l209_genome_chrbinfill(&mut genome);

    assert_eq!(genome.n_chr_real, 3);
    assert_eq!(genome.genome_chr_bin_nbases, 16);
    assert_eq!(genome.chr_bin_n, 5);
    assert_eq!(genome.chr_bin, vec![0, 0, 1, 2, 2]);
}

#[test]
fn genome_free_memory_clears_index_arrays_only_without_shared_memory() {
    let mut genome = Genome {
        p_ge: ParametersGenome {
            g_load: "NoSharedMemory".to_string(),
            ..Default::default()
        },
        g: vec![1, 2, 3],
        sa: vec![4, 5],
        sa_pass2: vec![6, 7],
        sai: vec![8, 9],
        chr_bin: vec![10],
        ..Default::default()
    };

    assert!(genome_l33_genome_freememory(&mut genome));
    assert!(genome.g.is_empty());
    assert!(genome.sa.is_empty());
    assert!(genome.sa_pass2.is_empty());
    assert!(genome.sai.is_empty());
    assert_eq!(genome.chr_bin, vec![10]);

    let mut shared = Genome {
        p_ge: ParametersGenome {
            g_load: "LoadAndKeep".to_string(),
            ..Default::default()
        },
        g: vec![1],
        sa: vec![2],
        sa_pass2: vec![3],
        sai: vec![4],
        ..Default::default()
    };
    assert!(!genome_l33_genome_freememory(&mut shared));
    assert_eq!(shared.g, vec![1]);
    assert_eq!(shared.sa, vec![2]);
    assert_eq!(shared.sa_pass2, vec![3]);
    assert_eq!(shared.sai, vec![4]);
}

#[test]
fn solo_read_constructor_allocates_barcode_and_requested_features() {
    let p = Parameters {
        out_file_tmp: "/tmp/star-solo".to_string(),
        p_solo: ParametersSolo {
            solo_type: 1,
            cb_wl_yes: true,
            cb_wl_size: 3,
            umi_l: 4,
            features: vec![5, 1],
            n_features: 2,
            read_info_yes: vec![false, false, true, false, false, true],
            read_index_yes: vec![false; 6],
            ..Default::default()
        },
        ..Default::default()
    };

    let solo_read = soloread_l3_soloread_soloread(&p, 7);

    assert_eq!(solo_read.i_chunk, 7);
    assert_eq!(
        solo_read.read_bar.as_ref().unwrap().cb_read_count_exact,
        vec![0; 3]
    );
    assert_eq!(solo_read.read_feat.len(), 2);
    assert_eq!(solo_read.read_feat[0].feature_type, 5);
    assert!(solo_read.read_feat[0].read_info_yes);
    assert_eq!(
        solo_read.read_feat[0].stream_reads_path.as_deref(),
        Some("/tmp/star-solo/soloGene_7")
    );
    assert_eq!(solo_read.read_feat[1].feature_type, 1);
    assert_eq!(solo_read.read_feat[1].transcript_dist_count.len(), 10000);

    let no_solo = soloread_l3_soloread_soloread(
        &Parameters {
            p_solo: ParametersSolo {
                solo_type: 0,
                features: vec![5],
                n_features: 1,
                ..Default::default()
            },
            ..Default::default()
        },
        0,
    );
    assert!(no_solo.read_bar.is_some());
    assert!(no_solo.read_feat.is_empty());

    let sam_tag_out = soloread_l3_soloread_soloread(
        &Parameters {
            p_solo: ParametersSolo {
                solo_type: 3,
                features: vec![5],
                n_features: 1,
                ..Default::default()
            },
            ..Default::default()
        },
        0,
    );
    assert!(sam_tag_out.read_feat.is_empty());
}

#[test]
fn solo_read_record_resets_flags_and_dispatches_configured_features() {
    let mut p = Parameters {
        p_solo: ParametersSolo {
            solo_type: 1,
            cb_wl_yes: false,
            features: vec![SOLO_FEATURE_GENE as u32, SOLO_FEATURE_SJ as u32],
            n_features: 2,
            read_stats_yes: vec![false; 8],
            ..Default::default()
        },
        ..Default::default()
    };
    p.p_solo.read_stats_yes[SOLO_FEATURE_GENE as usize] = true;

    let mut solo_read = SoloRead {
        read_bar: Some(SoloReadBarcode {
            umi_b: 12,
            cb_match: 0,
            cb_match_string: "CB".to_string(),
            cb_match_ind: vec![5],
            ..Default::default()
        }),
        read_feat: vec![
            SoloReadFeature {
                feature_type: SOLO_FEATURE_GENE,
                read_index_yes: true,
                read_flag: SoloReadFlagClass {
                    flag: 999,
                    ..Default::default()
                },
                ..Default::default()
            },
            SoloReadFeature {
                feature_type: SOLO_FEATURE_SJ,
                read_index_yes: false,
                read_flag: SoloReadFlagClass {
                    flag: 888,
                    ..Default::default()
                },
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    let mut ann_features = vec![ReadAnnotFeature::default(); 8];
    ann_features[SOLO_FEATURE_GENE as usize].f_set = std::collections::BTreeSet::from([9]);

    soloread_record_l3_soloread_record(
        &mut solo_read,
        &p,
        1,
        &[Transcript {
            n_exons: 2,
            exons: vec![[0, 10, 5, 0, 0], [5, 21, 6, 0, 0]],
            canon_sj: vec![1],
            sj_annot: vec![0],
            ..Default::default()
        }],
        44,
        &ReadAnnotations {
            annot_features: ann_features,
            ..Default::default()
        },
    );

    assert_eq!(solo_read.read_feat[0].stream_reads, "12 44 80 9 0 CB\n");
    assert_eq!(solo_read.read_feat[1].stream_reads, "12 15 6 0 CB\n");
    assert_eq!(solo_read.read_feat[0].cb_read_count_map.get(&5), Some(&1));
    assert_eq!(solo_read.read_feat[1].cb_read_count_map.get(&5), Some(&1));

    let mut tag_only = solo_read.clone();
    let tag_only_before = tag_only.clone();
    soloread_record_l3_soloread_record(
        &mut tag_only,
        &Parameters {
            p_solo: ParametersSolo {
                solo_type: 3,
                n_features: 2,
                ..Default::default()
            },
            ..Default::default()
        },
        1,
        &[],
        1,
        &ReadAnnotations::default(),
    );
    assert_eq!(tag_only, tag_only_before);
}

#[test]
fn solo_read_feature_input_records_counts_exact_rescued_and_rejected_records() {
    let mut rf = SoloReadFeature {
        feature_type: SOLO_FEATURE_GENE,
        read_info_yes: true,
        read_index_yes: true,
        stream_reads: [
            format!("101 0 {} 7 0 0", 1_u32 << SOLO_READ_FLAG_FEATURE_U),
            format!("102 1 {} 8 1 1", 1_u32 << SOLO_READ_FLAG_FEATURE_M),
            format!("103 2 {} 9 2 0 I 2 !", 1_u32 << SOLO_READ_FLAG_FEATURE_U),
            format!("104 3 {} 10 2 0 ! 2 !", 1_u32 << SOLO_READ_FLAG_FEATURE_U),
            "105 4 0 4294967295 0 0".to_string(),
        ]
        .join("\n"),
        stats: star_rs::generated::structs::SoloReadFeatureStats {
            v: vec![0; SOLO_READ_FEATURE_N_STATS],
            ..Default::default()
        },
        ..Default::default()
    };
    let p_solo = ParametersSolo {
        cb_wl_yes: true,
        cb_match_wl: CBMatchWL {
            one_exact: true,
            ..Default::default()
        },
        read_stats_yes: {
            let mut v = vec![false; SOLO_FEATURE_GENE as usize + 1];
            v[SOLO_FEATURE_GENE as usize] = true;
            v
        },
        qs_base: 33,
        qs_max: 40,
        cb_min_p: 0.8,
        ..Default::default()
    };
    let mut cb_p = vec![Vec::<u32>::new(); 3];
    let cb_read_count_total = vec![5, 0, 10];
    let mut read_info = vec![SoloFeatureReadInfo::default(); 6];
    let mut read_flag_counts = SoloReadFlagClass::default();
    let mut n_read_per_cb_unique = vec![0; 3];
    let mut n_read_per_cb_multi = vec![0; 3];

    soloreadfeature_inputrecords_l8_soloreadfeature_inputrecords(
        &mut rf,
        &p_solo,
        &[Vec::new(), Vec::new()],
        &mut cb_p,
        3,
        &cb_read_count_total,
        &mut read_info,
        &mut read_flag_counts,
        &mut n_read_per_cb_unique,
        &mut n_read_per_cb_multi,
    );

    assert_eq!(cb_p[0], vec![7, 101, 0]);
    assert_eq!(cb_p[1], Vec::<u32>::new());
    assert_eq!(cb_p[2], vec![9, 103, 2]);
    assert_eq!(read_info[4], SoloFeatureReadInfo { cb: 0, umi: 105 });
    assert_eq!(rf.stats.v[SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_EXACT], 1);
    assert_eq!(
        rf.stats.v[SOLO_READ_FEATURE_STAT_NO_MM_TO_WL_WITHOUT_EXACT],
        1
    );
    assert_eq!(rf.stats.v[SOLO_READ_FEATURE_STAT_NO_TOO_MANY_WL_MATCHES], 1);
    assert_eq!(n_read_per_cb_unique, vec![1, 0, 1]);
    assert_eq!(n_read_per_cb_multi, vec![0, 0, 0]);

    let cb0_counts = read_flag_counts.flag_counts.get(&0).unwrap();
    assert_eq!(cb0_counts[SOLO_READ_FLAG_CB_MATCH as usize], 2);
    assert_eq!(cb0_counts[SOLO_READ_FLAG_CB_PERFECT as usize], 2);
    assert_eq!(cb0_counts[SOLO_READ_FLAG_COUNTED_U as usize], 1);

    let cb2_counts = read_flag_counts.flag_counts.get(&2).unwrap();
    assert_eq!(cb2_counts[SOLO_READ_FLAG_CB_MM_MULTIPLE as usize], 1);
    assert_eq!(cb2_counts[SOLO_READ_FLAG_COUNTED_U as usize], 1);
    assert_eq!(
        read_flag_counts.flag_counts_no_cb[SOLO_READ_FLAG_CB_MATCH as usize],
        2
    );
}

#[test]
fn solo_feature_count_cbgene_umi_reads_records_and_collapses_exact_umis() {
    let mut umi_dedup = UMIdedup {
        types_in: vec!["Exact".to_string()],
        types: vec![1],
        type_main: 1,
        yes_n: 1,
        yes_b: [false, true, false, false, false, false],
        count_ind_i: [u32::MAX, 1, u32::MAX, u32::MAX, u32::MAX, u32::MAX],
        count_ind_main: 1,
    };
    let p_solo = ParametersSolo {
        read_info_yes: {
            let mut v = vec![false; SOLO_FEATURE_GENE as usize + 1];
            v[SOLO_FEATURE_GENE as usize] = true;
            v
        },
        read_index_yes: {
            let mut v = vec![false; SOLO_FEATURE_GENE as usize + 1];
            v[SOLO_FEATURE_GENE as usize] = true;
            v
        },
        read_stats_yes: vec![false; SOLO_FEATURE_GENE as usize + 1],
        cb_wl_yes: true,
        cb_wl_size: 2,
        cb_match_wl: CBMatchWL {
            one_exact: true,
            ..Default::default()
        },
        umi_dedup: {
            umi_dedup.count_ind_i[1] = 1;
            umi_dedup
        },
        umi_filtering: SoloUmiFiltering::default(),
        qs_base: 33,
        qs_max: 40,
        cb_min_p: 0.8,
        ..Default::default()
    };
    let read_flag_unique = 1_u32 << SOLO_READ_FLAG_FEATURE_U;
    let read_feat = SoloReadFeature {
        feature_type: SOLO_FEATURE_GENE,
        read_info_yes: true,
        read_index_yes: true,
        stream_reads: [
            format!("101 0 {} 7 0 0", read_flag_unique),
            format!("101 1 {} 7 0 0", read_flag_unique),
            format!("102 2 {} 8 0 0", read_flag_unique),
            format!("201 3 {} 9 0 1", read_flag_unique),
            format!("202 4 {} 9 0 1", read_flag_unique),
        ]
        .join("\n"),
        stats: star_rs::generated::structs::SoloReadFeatureStats {
            v: vec![0; SOLO_READ_FEATURE_N_STATS],
            ..Default::default()
        },
        ..Default::default()
    };
    let mut solo_feature = SoloFeature {
        feature_type: SOLO_FEATURE_GENE,
        features_number: 10,
        read_bar_sum: Some(SoloReadBarcode {
            cb_read_count_exact: vec![3, 2],
            ..Default::default()
        }),
        read_feat_sum: Some(SoloReadFeature {
            cb_read_count: vec![3, 2],
            stats: star_rs::generated::structs::SoloReadFeatureStats {
                v: vec![0; SOLO_READ_FEATURE_N_STATS],
                ..Default::default()
            },
            ..Default::default()
        }),
        read_feat_all: vec![read_feat],
        n_reads_mapped: 5,
        n_reads_input: 5,
        n_cb: 2,
        ind_cb: vec![0, 1],
        ..Default::default()
    };

    let log = solofeature_countcbgeneumi_l7_solofeature_countcbgeneumi(
        &mut solo_feature,
        &Parameters {
            run_thread_n: 1,
            ..Default::default()
        },
        &p_solo,
        0,
    )
    .unwrap();

    assert!(log.contains("Finished reading reads from Solo files nCB=2, nReadPerCBmax=3"));
    assert_eq!(solo_feature.rgu_stride, 3);
    assert_eq!(solo_feature.n_read_per_cb, vec![3, 2]);
    assert_eq!(solo_feature.n_read_per_cb_unique, vec![3, 2]);
    assert_eq!(solo_feature.n_read_per_cb_total, vec![3, 2]);
    assert_eq!(solo_feature.n_gene_per_cb, vec![2, 1]);
    assert_eq!(solo_feature.n_umi_per_cb, vec![2, 2]);
    assert_eq!(solo_feature.count_mat_stride, 2);
    assert_eq!(solo_feature.count_cell_gene_umi_index, vec![0, 4, 6]);
    assert_eq!(&solo_feature.count_cell_gene_umi[..6], &[7, 1, 8, 1, 9, 2]);
    assert_eq!(
        solo_feature.read_info[0],
        SoloFeatureReadInfo { cb: 0, umi: 101 }
    );
    assert_eq!(
        solo_feature.read_info[1],
        SoloFeatureReadInfo { cb: 0, umi: 101 }
    );
    assert_eq!(
        solo_feature.read_info[4],
        SoloFeatureReadInfo { cb: 1, umi: 202 }
    );
    let stats = &solo_feature.read_feat_sum.as_ref().unwrap().stats.v;
    assert_eq!(stats[SOLO_READ_FEATURE_STAT_YES_UMIS], 4);
    assert_eq!(stats[SOLO_READ_FEATURE_STAT_YES_CELL_BARCODES], 2);
    assert_eq!(stats[SOLO_READ_FEATURE_STAT_YES_WL_MATCH], 5);
    assert_eq!(
        stats[SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE],
        5
    );
}

#[test]
fn solo_feature_redistribute_reads_by_cb_partitions_raw_records() {
    let mut solo_feature = SoloFeature {
        feature_type: SOLO_FEATURE_GENE,
        n_cb: 3,
        ind_cb: vec![0, 2, 4],
        ind_cb_wl: vec![0, 0, 1, 0, 2],
        read_feat_sum: Some(SoloReadFeature {
            cb_read_count: vec![2, 0, 4, 0, 6],
            ..Default::default()
        }),
        read_feat_all: vec![
            SoloReadFeature {
                stream_reads: "100 7 9 0\n101 8 9 2\n".to_string(),
                ..Default::default()
            },
            SoloReadFeature {
                stream_reads: "102 9 9 4\n103 10 9 0\n".to_string(),
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    let p_solo = ParametersSolo {
        redistr_reads_nfiles: 2,
        ..Default::default()
    };

    let log = solofeature_redistributereadsbycb_l8_solofeature_redistributereadsbycb(
        &mut solo_feature,
        &p_solo,
        2,
    );

    assert_eq!(
        log,
        "     Redistributing reads into 2files; nReadRec=12;   nReadRecBin=6\n"
    );
    assert_eq!(solo_feature.redistr_files_cb_first, vec![0, 2, 3]);
    assert_eq!(solo_feature.redistr_files_cb_index, vec![0, 0, 1]);
    assert_eq!(solo_feature.redistr_files_nreads, vec![6, 6]);
    assert_eq!(
        solo_feature.redistr_files_streams,
        vec![
            "100 7 9 0\n101 8 9 2\n103 10 9 0\n".to_string(),
            "102 9 9 4\n".to_string()
        ]
    );

    solo_feature.feature_type = SOLO_FEATURE_SJ;
    solo_feature.redistr_files_cb_first.clear();
    solo_feature.redistr_files_cb_index.clear();
    solo_feature.redistr_files_nreads.clear();
    solo_feature.redistr_files_streams.clear();
    solo_feature.read_feat_all = vec![SoloReadFeature {
        stream_reads: "200 1 2 3 4\n201 1 2 3 0\n".to_string(),
        ..Default::default()
    }];
    solofeature_redistributereadsbycb_l8_solofeature_redistributereadsbycb(
        &mut solo_feature,
        &p_solo,
        1,
    );
    assert_eq!(
        solo_feature.redistr_files_streams,
        vec!["201 1 2 3 0\n".to_string(), "200 1 2 3 4\n".to_string()]
    );
}

#[test]
fn sjdb_prepare_detects_motifs_collapses_duplicates_and_builds_sj_sequences() {
    let mut genome = Genome {
        g: vec![4; 32],
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![0, 32],
        chr_bin: vec![0; 32],
        sjdb_overhang: 2,
        sjdb_length: 5,
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    genome.g[2] = 0;
    genome.g[3] = 1;
    genome.g[4] = 2;
    genome.g[5] = 3;
    genome.g[8] = 0;
    genome.g[9] = 2;
    genome.g[10] = 1;
    genome.g[11] = 2;
    let p = Parameters {
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    let sjdb_loci = SjdbClass {
        chr: vec!["chr1".to_string(), "chr1".to_string(), "chr1".to_string()],
        start: vec![5, 5, 5],
        end: vec![10, 10, 10],
        str_: vec!['+', '+', '.'],
        priority: vec![10, 20, 20],
        ..Default::default()
    };

    let out = sjdbprepare_l5_sjdbprepare(&sjdb_loci, &p, 32, "/tmp", &mut genome).unwrap();

    assert_eq!(genome.sjdb_n, 1);
    assert_eq!(genome.sjdb_start, vec![4]);
    assert_eq!(genome.sjdb_end, vec![9]);
    assert_eq!(genome.sjdb_motif, vec![1]);
    assert_eq!(genome.sjdb_shift_left, vec![0]);
    assert_eq!(genome.sjdb_shift_right, vec![0]);
    assert_eq!(genome.sjdb_strand, vec![1]);
    assert_eq!(genome.sj_dstart, vec![2]);
    assert_eq!(genome.sj_astart, vec![10]);
    assert_eq!(out.gsj, vec![0, 1, 1, 2, GENOME_SPACING_CHAR]);
    assert_eq!(out.sjdb_info_txt, "1\t2\n4\t9\t1\t0\t0\t1\n");
    assert_eq!(out.sjdb_list_out_tab, "chr1\t5\t10\t+\n");
    assert_eq!(out.log_main, "");
}

#[test]
fn sjdb_insert_junctions_prepares_builds_and_recomputes_win_bins() {
    let mut sa = packedarray_l3_packedarray_packedarray();
    packedarray_l8_packedarray_definebits(&mut sa, 33, 1);
    packedarray_l31_packedarray_allocatearray(&mut sa);
    packedarray_l17_packedarray_writepacked(&mut sa, 0, 0);

    let mut genome = Genome {
        g: vec![4; 40],
        sa: vec![0],
        sa_packed: sa,
        sa_insert: PackedArray::default(),
        sai_packed: PackedArray::default(),
        n_genome: 32,
        n_sa: 1,
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![0, 32],
        chr_bin: vec![0; 40],
        gstrand_bit: 32,
        gstrand_mask: u32::MAX,
        sjdb_overhang: 2,
        sjdb_length: 5,
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    genome.g[2] = 0;
    genome.g[3] = 1;
    genome.g[4] = 2;
    genome.g[5] = 3;
    genome.g[8] = 0;
    genome.g[9] = 2;
    genome.g[10] = 1;
    genome.g[11] = 2;

    let genome1 = Genome {
        n_genome: 32,
        n_sa: 1,
        sjdb_n: 0,
        ..Default::default()
    };
    let mut p = Parameters {
        run_mode_in: vec!["genomeGenerate".to_string()],
        limit_sjdb_insert_nsj: 10,
        win_bin_nbits: 2,
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut loci = SjdbClass {
        chr: vec!["chr1".to_string()],
        start: vec![5],
        end: vec![10],
        str_: vec!['+'],
        priority: vec![10],
        ..Default::default()
    };

    let out = sjdbinsertjunctions_l11_sjdbinsertjunctions(&mut p, &mut genome, &genome1, &mut loci)
        .unwrap();

    assert_eq!(genome.sjdb_n, 1);
    assert_eq!(genome.sjdb_start, vec![4]);
    assert_eq!(genome.sjdb_end, vec![9]);
    assert_eq!(genome.n_genome, 37);
    assert_eq!(&genome.g[32..37], &[0, 1, 1, 2, GENOME_SPACING_CHAR]);
    assert_eq!(p.win_bin_n, 10);
    assert_eq!(out.sjdb_prepare.sjdb_info_txt, "1\t2\n4\t9\t1\t0\t0\t1\n");
    assert!(
        out.sjdb_build_index
            .log_main
            .contains("inserting junctions into the genome indices")
    );
    assert!(out.log_main.contains("Finished preparing junctions"));
}

#[test]
fn sjdb_insert_junctions_reports_limit_after_prepare() {
    let mut genome = Genome {
        g: vec![4; 32],
        n_genome: 32,
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![0, 32],
        chr_bin: vec![0; 32],
        sjdb_overhang: 2,
        sjdb_length: 5,
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    genome.g[2] = 0;
    genome.g[3] = 1;
    genome.g[4] = 2;
    genome.g[5] = 3;
    genome.g[8] = 0;
    genome.g[9] = 2;
    let genome1 = Genome::default();
    let mut p = Parameters {
        run_mode_in: vec!["genomeGenerate".to_string()],
        limit_sjdb_insert_nsj: 0,
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut loci = SjdbClass {
        chr: vec!["chr1".to_string()],
        start: vec![5],
        end: vec![10],
        str_: vec!['+'],
        priority: vec![10],
        ..Default::default()
    };

    let err = sjdbinsertjunctions_l11_sjdbinsertjunctions(&mut p, &mut genome, &genome1, &mut loci)
        .unwrap_err();
    assert!(err.contains("limitSjdbInsertNsj=0"));
    assert_eq!(genome.sjdb_n, 1);
}

#[test]
fn solo_feature_count_velocyto_intersects_umis_and_classifies_counts() {
    let mut solo_feature = SoloFeature {
        feature_type: SOLO_FEATURE_VELOCYTO,
        n_cb: 2,
        n_reads_mapped: 20,
        ind_cb_wl: vec![0, 1],
        read_feat_sum: Some(SoloReadFeature {
            cb_read_count: vec![4, 3],
            stats: star_rs::generated::structs::SoloReadFeatureStats {
                v: vec![0; SOLO_READ_FEATURE_N_STATS],
                ..Default::default()
            },
            ..Default::default()
        }),
        read_feat_all: vec![SoloReadFeature {
            stream_reads: [
                format!("0 1 0 {}", 1_u8 << ALIGN_VS_TRANSCRIPT_CONCORDANT as u32),
                format!("0 1 0 {}", 1_u8 << ALIGN_VS_TRANSCRIPT_CONCORDANT as u32),
                format!("1 1 1 {}", 1_u8 << ALIGN_VS_TRANSCRIPT_INTRON as u32),
                format!(
                    "2 2 0 {} 2 {}",
                    1_u8 << ALIGN_VS_TRANSCRIPT_CONCORDANT as u32,
                    1_u8 << ALIGN_VS_TRANSCRIPT_INTRON as u32
                ),
                format!(
                    "3 2 3 {} 4 {}",
                    1_u8 << ALIGN_VS_TRANSCRIPT_CONCORDANT as u32,
                    1_u8 << ALIGN_VS_TRANSCRIPT_INTRON as u32
                ),
                format!("4 1 0 {}", 1_u8 << ALIGN_VS_TRANSCRIPT_CONCORDANT as u32),
            ]
            .join("\n"),
            ..Default::default()
        }],
        ..Default::default()
    };
    let gene_read_info = vec![
        SoloFeatureReadInfo { cb: 0, umi: 10 },
        SoloFeatureReadInfo { cb: 0, umi: 11 },
        SoloFeatureReadInfo { cb: 0, umi: 12 },
        SoloFeatureReadInfo { cb: 1, umi: 20 },
        SoloFeatureReadInfo {
            cb: -1,
            umi: u64::MAX,
        },
    ];
    let transcriptome = Transcriptome {
        tr_gene: vec![0, 1, 2, 3, 3],
        ..Default::default()
    };

    let log = solofeature_countvelocyto_l12_solofeature_countvelocyto(
        &mut solo_feature,
        &gene_read_info,
        &transcriptome,
        1,
        "T0",
        "T1",
        "T2",
        "MEM\n",
    );

    assert_eq!(
        log,
        "T0 ... Velocyto counting: allocated arrays\nT1 ... Velocyto counting: finished input\nT2 ... Velocyto counting: finished collapsing UMIs\nRAM for solo feature Velocyto\nMEM\n"
    );
    assert_eq!(solo_feature.n_read_per_cb, vec![4, 1]);
    assert_eq!(solo_feature.n_read_per_cb_total, vec![4, 1]);
    assert_eq!(solo_feature.n_read_per_cb_unique, vec![4, 1]);
    assert_eq!(solo_feature.n_umi_per_cb, vec![2, 1]);
    assert_eq!(solo_feature.n_gene_per_cb, vec![2, 1]);
    assert_eq!(solo_feature.count_mat_stride, 4);
    assert_eq!(solo_feature.count_cell_gene_umi_index, vec![0, 8, 12]);
    assert_eq!(
        &solo_feature.count_cell_gene_umi[..12],
        &[0, 1, 0, 0, 1, 0, 1, 0, 3, 0, 0, 1]
    );
    let stats = &solo_feature.read_feat_sum.as_ref().unwrap().stats.v;
    assert_eq!(stats[SOLO_READ_FEATURE_STAT_YES_UMIS], 3);
    assert_eq!(stats[SOLO_READ_FEATURE_STAT_YES_CELL_BARCODES], 2);
}

#[test]
fn solo_feature_count_smartseq_collapses_feature_umis_per_cell() {
    let p = Parameters {
        run_thread_n: 2,
        p_solo: ParametersSolo {
            redistr_reads_nfiles: 2,
            umi_dedup: UMIdedup {
                yes_n: 2,
                yes_b: [true, true, false, false, false, false],
                count_ind_i: [1, 2, u32::MAX, u32::MAX, u32::MAX, u32::MAX],
                count_ind_main: 2,
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    let mut read_feat_sum = soloreadfeature_l5_soloreadfeature_soloreadfeature(
        SOLO_FEATURE_GENE,
        &Parameters::default(),
        -1,
    );
    read_feat_sum.cb_read_count = vec![4, 2];
    read_feat_sum.stats = star_rs::generated::structs::SoloReadFeatureStats {
        names: vec![],
        v: vec![0; SOLO_READ_FEATURE_N_STATS],
    };
    let read_feat_all = vec![
        SoloReadFeature {
            stream_reads: concat!("10 5 0 0\n", "10 5 0 0\n", "11 5 0 0\n", "20 3 0 0\n",)
                .to_string(),
            stats: star_rs::generated::structs::SoloReadFeatureStats {
                names: vec![],
                v: {
                    let mut v = vec![0; SOLO_READ_FEATURE_N_STATS];
                    v[SOLO_READ_FEATURE_STAT_NO_NO_FEATURE] = 2;
                    v
                },
            },
            ..Default::default()
        },
        SoloReadFeature {
            stream_reads: concat!("30 5 0 1\n", "30 5 0 1\n").to_string(),
            stats: star_rs::generated::structs::SoloReadFeatureStats {
                names: vec![],
                v: {
                    let mut v = vec![0; SOLO_READ_FEATURE_N_STATS];
                    v[SOLO_READ_FEATURE_STAT_NO_UNMAPPED] = 1;
                    v
                },
            },
            ..Default::default()
        },
    ];
    let mut solo_feature = SoloFeature {
        feature_type: SOLO_FEATURE_GENE,
        n_cb: 2,
        ind_cb: vec![0, 1],
        ind_cb_wl: vec![0, 1],
        read_feat_sum: Some(read_feat_sum),
        read_feat_all,
        ..Default::default()
    };

    let log = solofeature_countsmartseq_l9_solofeature_countsmartseq(
        &mut solo_feature,
        &p,
        "T1",
        "T2",
        "T3",
    );

    assert!(log.contains("Redistributing reads into 2files; nReadRec=6"));
    assert!(log.contains("T1 ... Finished redistribution"));
    assert!(log.contains("T2 ... Finished reading / collapsing"));
    assert!(log.contains("T3 ... Finished SmartSeq counting"));
    assert_eq!(solo_feature.redistr_files_cb_first, vec![0, 1, 2]);
    assert_eq!(solo_feature.n_read_per_cb, vec![4, 2]);
    assert_eq!(solo_feature.n_read_per_cb_total, vec![4, 2]);
    assert_eq!(solo_feature.n_read_per_cb_unique, vec![4, 2]);
    assert_eq!(solo_feature.count_mat_stride, 3);
    assert_eq!(solo_feature.n_gene_per_cb, vec![2, 1]);
    assert_eq!(solo_feature.n_umi_per_cb, vec![3, 1]);
    assert_eq!(solo_feature.count_cell_gene_umi_index, vec![0, 6, 9]);
    assert_eq!(
        solo_feature.count_cell_gene_umi,
        vec![3, 1, 1, 5, 3, 2, 5, 2, 1]
    );
    let stats = &solo_feature.read_feat_sum.as_ref().unwrap().stats.v;
    assert_eq!(stats[SOLO_READ_FEATURE_STAT_NO_NO_FEATURE], 2);
    assert_eq!(stats[SOLO_READ_FEATURE_STAT_NO_UNMAPPED], 1);
    assert_eq!(stats[SOLO_READ_FEATURE_STAT_YES_WL_MATCH], 6);
    assert_eq!(stats[SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_EXACT], 6);
    assert_eq!(
        stats[SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE],
        6
    );
    assert_eq!(stats[SOLO_READ_FEATURE_STAT_YES_UMIS], 4);
    assert_eq!(stats[SOLO_READ_FEATURE_STAT_YES_CELL_BARCODES], 2);
}

#[test]
fn sam_headers_builds_main_and_transcriptome_bam_headers() {
    let mut p = Parameters {
        command_line_full: "STAR --genomeDir real --readFilesIn reads.fq".to_string(),
        command_line: "STAR --runThreadN 8".to_string(),
        out_sam_mode: "Full".to_string(),
        out_sam_type: vec!["SAM".to_string()],
        out_sam_bool: true,
        out_bam_unsorted: true,
        quant_tr_sam_bam_yes: true,
        out_sam_header_hd: vec!["-".to_string()],
        out_sam_header_pg: vec![
            "@PG".to_string(),
            "ID:pre".to_string(),
            "PN:preprocessor".to_string(),
        ],
        out_sam_header_comment_file: "comments.sam".to_string(),
        out_sam_attr_rgline_split: vec!["ID:rg1\tSM:sample1".to_string()],
        sam_header_extra: "@CO\textra header line\n".to_string(),
        ..Default::default()
    };
    let mut genome = Genome {
        n_chr_real: 2,
        chr_name: vec!["chr1".to_string(), "chr2".to_string()],
        chr_length: vec![248_956_422, 242_193_529],
        ..Default::default()
    };
    let transcriptome = Transcriptome {
        tr_id: vec!["ENST00000335137".to_string(), "ENST00000448914".to_string()],
        tr_ex_i: vec![0, 2],
        tr_ex_n: vec![2, 1],
        ex_len_cum: vec![0, 50, 0],
        ex_se: vec![100, 149, 200, 249, 300, 319],
        ..Default::default()
    };

    samheaders_l5_samheaders(
        &mut p,
        &mut genome,
        &transcriptome,
        "\n@SQ\tSN:phiX\tLN:5386\n",
        "  \n@CO\tlibrary: real RNA-seq\n\t\n@CO\tinstrument: NovaSeq\n",
    );

    let expected_tail = concat!(
        "@SQ\tSN:chr1\tLN:248956422\n",
        "@SQ\tSN:chr2\tLN:242193529\n",
        "@SQ\tSN:phiX\tLN:5386\n",
        "@PG\tID:pre\tPN:preprocessor\n",
        "@PG\tID:STAR\tPN:STAR\tVN:2.7.11b\tCL:STAR --genomeDir real --readFilesIn reads.fq\n",
        "@CO\tlibrary: real RNA-seq\n",
        "@CO\tinstrument: NovaSeq\n",
        "@RG\tID:rg1\tSM:sample1\n",
        "@CO\tuser command line: STAR --runThreadN 8\n",
        "@CO\textra header line\n",
    );
    let expected_header = format!("@HD\tVN:1.4\n{}", expected_tail);

    assert_eq!(p.sam_header_hd, "@HD\tVN:1.4");
    assert_eq!(p.sam_header, expected_header);
    assert_eq!(p.out_sam_contents, p.sam_header);
    assert_eq!(
        p.sam_header_sorted_coord,
        format!("@HD\tVN:1.4\tSO:coordinate\n{}", expected_tail)
    );
    assert_eq!(
        genome.chr_name_all,
        vec!["chr1".to_string(), "chr2".to_string(), "phiX".to_string()]
    );
    assert_eq!(genome.chr_length_all, vec![248_956_422, 242_193_529, 5386]);
    assert_eq!(
        p.out_bam_unsorted_header,
        bamfunctions_l77_outbamwriteheader(
            &p.sam_header,
            &genome.chr_name_all,
            &genome.chr_length_all
        )
    );

    let expected_tr_header = concat!(
        "@SQ\tSN:ENST00000335137\tLN:100\n",
        "@SQ\tSN:ENST00000448914\tLN:20\n",
        "@RG\tID:rg1\tSM:sample1\n",
    );
    assert_eq!(
        p.out_quant_bam_header,
        bamfunctions_l77_outbamwriteheader(expected_tr_header, &transcriptome.tr_id, &[100, 20],)
    );
}

#[test]
fn gtf_transcript_gene_sj_writes_tables_and_collapses_duplicate_junctions() {
    let genome = Genome {
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 10,
            sjdb_gtf_file: "real_chr1.gtf".to_string(),
            ..Default::default()
        },
        chr_bin: vec![0, 0, 0],
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![1000],
        ..Default::default()
    };
    let mut gtf = GTF {
        gtf_yes: true,
        exon_loci: vec![
            [1, 1400, 1449, 1],
            [0, 1200, 1249, 0],
            [0, 1000, 1099, 0],
            [2, 1000, 1099, 2],
            [2, 1200, 1249, 2],
        ],
        transcript_strand: vec![1, 2, 1],
        transcript_id: vec![
            "ENST_A".to_string(),
            "ENST_B".to_string(),
            "ENST_C".to_string(),
        ],
        gene_id: vec![
            "GENE_A".to_string(),
            "GENE_B".to_string(),
            "GENE_C".to_string(),
        ],
        gene_attr: vec![
            ["GeneA".to_string(), "protein_coding".to_string()],
            ["GeneB".to_string(), "lncRNA".to_string()],
            ["GeneC".to_string(), "protein_coding".to_string()],
        ],
        ..Default::default()
    };
    let mut sjdb = SjdbClass {
        chr: vec!["chrM".to_string()],
        start: vec![1],
        end: vec![10],
        str_: vec!['+'],
        priority: vec![10],
        ..Default::default()
    };
    let mut log = String::new();

    let out = gtf_transcriptgenesj_l23_gtf_transcriptgenesj(
        &mut gtf, &genome, &mut sjdb, "/tmp", &mut log,
    );

    assert_eq!(out.n_junctions_added, 1);
    assert_eq!(gtf.exon_n, 5);
    assert_eq!(
        out.gene_info_tab,
        "3\nGENE_A\tGeneA\tprotein_coding\nGENE_B\tGeneB\tlncRNA\nGENE_C\tGeneC\tprotein_coding\n"
    );
    assert_eq!(
        out.exon_ge_tr_info_tab,
        "5\n1000\t1099\t1\t0\t0\n1000\t1099\t1\t2\t2\n1200\t1249\t1\t0\t0\n1200\t1249\t1\t2\t2\n1400\t1449\t2\t1\t1\n"
    );
    assert_eq!(
        out.exon_info_tab,
        "5\n0\t99\t0\n200\t249\t100\n0\t99\t0\n200\t249\t100\n0\t49\t0\n"
    );
    assert_eq!(
        out.transcript_info_tab,
        "3\nENST_A\t1000\t1249\t1249\t1\t2\t0\t0\nENST_C\t1000\t1249\t1249\t1\t2\t2\t2\nENST_B\t1400\t1449\t1249\t2\t1\t4\t1\n"
    );
    assert_eq!(out.sjdb_list_from_gtf_out_tab, "chr1\t101\t200\t+\t1,3\n");
    assert_eq!(sjdb.chr, vec!["chrM".to_string(), "chr1".to_string()]);
    assert_eq!(sjdb.start, vec![1, 101]);
    assert_eq!(sjdb.end, vec![10, 200]);
    assert_eq!(sjdb.str_, vec!['+', '+']);
    assert_eq!(sjdb.priority, vec![10, 20]);
    assert_eq!(sjdb.gene[1].iter().copied().collect::<Vec<_>>(), vec![1, 3]);
    assert!(log.contains("3 transcripts"));
    assert!(log.contains("1 collapsed junctions"));
}

#[test]
fn gtf_supertranscript_condenses_exons_and_builds_transcriptome_outputs() {
    let mut genome = Genome {
        n_chr_real: 1,
        genome_chr_bin_nbases: 4,
        n_genome: 8,
        chr_name: vec!["chr1".to_string()],
        chr_start: vec![0, 8],
        chr_length: vec![8],
        g: vec![0, 1, 2, 3, 0, 1, 2, 3, 5, 5, 4, 0, 3, 2, 1, 0],
        ..Default::default()
    };
    let mut p = Parameters {
        p_ge: ParametersGenome {
            g_type_string: "Transcriptome".to_string(),
            g_dir: "/tmp/star_rs_gtf_st".to_string(),
            ..Default::default()
        },
        sjdb_insert_yes: true,
        ..Default::default()
    };
    let mut gtf = GTF {
        gtf_yes: true,
        exon_loci: vec![[0, 1, 3, 0], [1, 2, 4, 1], [0, 6, 7, 0]],
        transcript_strand: vec![1, 1],
        transcript_id: vec!["ENST_A".to_string(), "ENST_B".to_string()],
        gene_id: vec!["GENE_A".to_string(), "GENE_B".to_string()],
        gene_attr: vec![
            ["GeneA".to_string(), "protein_coding".to_string()],
            ["GeneB".to_string(), "lncRNA".to_string()],
        ],
        ..Default::default()
    };

    let out = gtf_supertranscript_l9_gtf_supertranscript(&mut gtf, &mut genome, &mut p);

    assert_eq!(out.full_genome_chr_name_txt, "chr1\n");
    assert_eq!(out.full_genome_chr_start_txt, "0\n8\n");
    assert_eq!(out.full_genome_chr_length_txt, "8\n");
    assert_eq!(out.full_genome_sequence, vec![0, 1, 2, 3, 0, 1, 2, 3]);
    assert_eq!(gtf.super_trome.seq_concat, vec![1, 2, 3, 0, 2, 3]);
    assert_eq!(gtf.super_trome.seq, vec![vec![1, 2, 3, 0, 2, 3]]);
    assert_eq!(gtf.transcript_start_end, vec![[0, 5], [1, 3]]);
    assert_eq!(gtf.super_trome.tr_index, vec![0, 0]);
    assert_eq!(gtf.super_trome.tr_start_end, vec![[0, 5], [1, 3]]);
    assert_eq!(gtf.transcript_seq, vec![vec![1, 2, 3, 2, 3], vec![2, 3, 0]]);
    assert_eq!(
        out.transcript_sequences_fasta,
        ">ENST_A\nCGTGT\n>ENST_B\nGTA\n"
    );
    assert_eq!(out.super_transcript_sequences_fasta, ">st0\nCGTAGT\n");
    assert_eq!(out.super_transcript_sj_tsv, "0\t2\t4\n");
    assert_eq!(
        out.conversion_to_full_genome_tsv,
        "2\t8\n0\t4\t1\n4\t2\t6\n"
    );
    assert!(!gtf.gtf_yes);
    assert!(!p.sjdb_insert_yes);
    assert_eq!(
        genome.chr_name,
        vec!["ENST_A".to_string(), "ENST_B".to_string()]
    );
    assert_eq!(genome.chr_length, vec![5, 3]);
    assert_eq!(genome.chr_start, vec![0, 8, 16]);
    assert_eq!(&genome.g[0..5], &[1, 2, 3, 2, 3]);
    assert_eq!(&genome.g[8..11], &[2, 3, 0]);
    assert!(out.log_main.contains("condensed) genome length = 6"));
    assert!(out.log_main.contains("Number of superTranscripts = 1"));
}

#[test]
fn gtf_constructor_parses_exons_attributes_and_warnings() {
    let mut genome = Genome {
        sjdb_overhang: 100,
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_length: vec![100],
        chr_start: vec![1000, 1100],
        p_ge: ParametersGenome {
            sjdb_gtf_file: "genes.gtf".to_string(),
            sjdb_gtf_chr_prefix: "chr".to_string(),
            sjdb_gtf_feature_exon: "exon".to_string(),
            sjdb_gtf_tag_exon_parent_transcript: "transcript_id".to_string(),
            sjdb_gtf_tag_exon_parent_gene: "gene_id".to_string(),
            sjdb_gtf_tag_exon_parent_gene_name: vec!["gene_name".to_string()],
            sjdb_gtf_tag_exon_parent_gene_type: vec![
                "gene_biotype".to_string(),
                "gene_type".to_string(),
            ],
            ..Default::default()
        },
        ..Default::default()
    };
    let p = Parameters::default();
    let gtf_text = concat!(
        "# comment\n",
        "1\tsrc\texon\t10\t19\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Name1\"; gene_biotype \"coding\";\n",
        "2\tsrc\texon\t1\t5\t.\t+\t.\tgene_id \"GX\"; transcript_id \"TX\";\n",
        "1\tsrc\texon\t95\t101\t.\t+\t.\tgene_id \"Gbad\"; transcript_id \"Tbad\";\n",
        "1\tsrc\texon\t20\t25\t.\t-\t.\tgene_name \"NoIds\";\n",
        "1\tsrc\tCDS\t30\t35\t.\t+\t.\tgene_id \"Gcds\"; transcript_id \"Tcds\";\n",
    );

    let (gtf, log) = gtf_l7_gtf_gtf(&mut genome, &p, "/tmp", Some(gtf_text)).unwrap();

    assert!(gtf.gtf_yes);
    assert_eq!(gtf.exon_n, 2);
    assert_eq!(
        gtf.transcript_id,
        vec!["T1".to_string(), "tr_chr1_20_25_1".to_string()]
    );
    assert_eq!(gtf.transcript_strand, vec![1, 2]);
    assert_eq!(
        gtf.gene_id,
        vec!["G1".to_string(), "MissingGeneID".to_string()]
    );
    assert_eq!(
        gtf.gene_attr,
        vec![
            ["Name1".to_string(), "coding".to_string()],
            ["NoIds".to_string(), "MissingGeneType".to_string()],
        ]
    );
    assert_eq!(gtf.exon_loci, vec![[0, 1009, 1018, 0], [1, 1019, 1024, 1]]);
    assert_eq!(genome.chr_name_index.get("chr1"), Some(&0));
    assert!(log.contains("processing annotations GTF"));
    assert!(log.contains("chromosome 'chr2' not found"));
    assert!(log.contains("exon end = 101 is larger"));
    assert!(log.contains("no transcript_id"));
    assert!(log.contains("no gene_id"));
}
