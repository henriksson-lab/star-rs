#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Parameters::Parameters` at STAR/source/Parameters.cpp:19. Args: "]
pub fn parameters_l19_parameters_parameters() -> crate::parameter_info::ParametersScanState {
    let par_specs = [
        (-1, -1, "versionGenome"),
        (-1, 2, "parametersFiles"),
        (-1, -1, "sysShell"),
        (-1, -1, "runMode"),
        (-1, -1, "runThreadN"),
        (-1, -1, "runDirPerm"),
        (-1, -1, "runRNGseed"),
        (-1, -1, "genomeType"),
        (-1, -1, "genomeDir"),
        (-1, -1, "genomeLoad"),
        (-1, -1, "genomeFastaFiles"),
        (-1, -1, "genomeChainFiles"),
        (-1, -1, "genomeSAindexNbases"),
        (-1, -1, "genomeChrBinNbits"),
        (-1, -1, "genomeSAsparseD"),
        (-1, -1, "genomeSuffixLengthMax"),
        (-1, -1, "genomeFileSizes"),
        (-1, -1, "genomeTransformType"),
        (-1, -1, "genomeTransformVCF"),
        (-1, -1, "genomeTransformOutput"),
        (-1, -1, "genomeChrSetMitochondrial"),
        (-1, -1, "readFilesType"),
        (-1, -1, "readFilesIn"),
        (-1, -1, "readFilesPrefix"),
        (-1, -1, "readFilesCommand"),
        (-1, -1, "readMatesLengthsIn"),
        (-1, -1, "readMapNumber"),
        (-1, -1, "readNameSeparator"),
        (-1, -1, "readQualityScoreBase"),
        (-1, -1, "readFilesManifest"),
        (-1, -1, "readFilesSAMattrKeep"),
        (-1, -1, "inputBAMfile"),
        (-1, -1, "bamRemoveDuplicatesType"),
        (-1, -1, "bamRemoveDuplicatesMate2basesN"),
        (-1, -1, "limitGenomeGenerateRAM"),
        (-1, -1, "limitIObufferSize"),
        (-1, -1, "limitOutSAMoneReadBytes"),
        (-1, -1, "limitOutSJcollapsed"),
        (-1, -1, "limitOutSJoneRead"),
        (-1, -1, "limitBAMsortRAM"),
        (-1, -1, "limitSjdbInsertNsj"),
        (-1, -1, "limitNreadsSoft"),
        (-1, 2, "outFileNamePrefix"),
        (-1, 2, "outTmpDir"),
        (-1, 2, "outTmpKeep"),
        (-1, 2, "outStd"),
        (-1, -1, "outReadsUnmapped"),
        (-1, -1, "outQSconversionAdd"),
        (-1, -1, "outMultimapperOrder"),
        (-1, -1, "outSAMtype"),
        (-1, -1, "outSAMmode"),
        (-1, -1, "outSAMstrandField"),
        (-1, -1, "outSAMattributes"),
        (-1, -1, "outSAMunmapped"),
        (-1, -1, "outSAMorder"),
        (-1, -1, "outSAMprimaryFlag"),
        (-1, -1, "outSAMreadID"),
        (-1, -1, "outSAMmapqUnique"),
        (-1, -1, "outSAMflagOR"),
        (-1, -1, "outSAMflagAND"),
        (-1, -1, "outSAMattrRGline"),
        (-1, -1, "outSAMheaderHD"),
        (-1, -1, "outSAMheaderPG"),
        (-1, -1, "outSAMheaderCommentFile"),
        (-1, -1, "outBAMcompression"),
        (-1, -1, "outBAMsortingThreadN"),
        (-1, -1, "outBAMsortingBinsN"),
        (-1, -1, "outSAMfilter"),
        (-1, -1, "outSAMmultNmax"),
        (-1, -1, "outSAMattrIHstart"),
        (-1, -1, "outSAMtlen"),
        (-1, -1, "outSJtype"),
        (-1, -1, "outSJfilterReads"),
        (-1, -1, "outSJfilterCountUniqueMin"),
        (-1, -1, "outSJfilterCountTotalMin"),
        (-1, -1, "outSJfilterOverhangMin"),
        (-1, -1, "outSJfilterDistToOtherSJmin"),
        (-1, -1, "outSJfilterIntronMaxVsReadN"),
        (-1, -1, "outWigType"),
        (-1, -1, "outWigStrand"),
        (-1, -1, "outWigReferencesPrefix"),
        (-1, -1, "outWigNorm"),
        (-1, -1, "outFilterType"),
        (-1, -1, "outFilterMultimapNmax"),
        (-1, -1, "outFilterMultimapScoreRange"),
        (-1, -1, "outFilterScoreMin"),
        (-1, -1, "outFilterScoreMinOverLread"),
        (-1, -1, "outFilterMatchNmin"),
        (-1, -1, "outFilterMatchNminOverLread"),
        (-1, -1, "outFilterMismatchNmax"),
        (-1, -1, "outFilterMismatchNoverLmax"),
        (-1, -1, "outFilterMismatchNoverReadLmax"),
        (-1, -1, "outFilterIntronMotifs"),
        (-1, -1, "outFilterIntronStrands"),
        (-1, -1, "clipAdapterType"),
        (-1, -1, "clip5pNbases"),
        (-1, -1, "clip3pNbases"),
        (-1, -1, "clip5pAfterAdapterNbases"),
        (-1, -1, "clip3pAfterAdapterNbases"),
        (-1, -1, "clip5pAdapterSeq"),
        (-1, -1, "clip3pAdapterSeq"),
        (-1, -1, "clip5pAdapterMMp"),
        (-1, -1, "clip3pAdapterMMp"),
        (-1, -1, "winBinNbits"),
        (-1, -1, "winAnchorDistNbins"),
        (-1, -1, "winFlankNbins"),
        (-1, -1, "winAnchorMultimapNmax"),
        (-1, -1, "winReadCoverageRelativeMin"),
        (-1, -1, "winReadCoverageBasesMin"),
        (-1, -1, "scoreGap"),
        (-1, -1, "scoreGapNoncan"),
        (-1, -1, "scoreGapGCAG"),
        (-1, -1, "scoreGapATAC"),
        (-1, -1, "scoreStitchSJshift"),
        (-1, -1, "scoreGenomicLengthLog2scale"),
        (-1, -1, "scoreDelBase"),
        (-1, -1, "scoreDelOpen"),
        (-1, -1, "scoreInsOpen"),
        (-1, -1, "scoreInsBase"),
        (-1, -1, "seedSearchLmax"),
        (-1, -1, "seedSearchStartLmax"),
        (-1, -1, "seedSearchStartLmaxOverLread"),
        (-1, -1, "seedPerReadNmax"),
        (-1, -1, "seedPerWindowNmax"),
        (-1, -1, "seedNoneLociPerWindow"),
        (-1, -1, "seedMultimapNmax"),
        (-1, -1, "seedSplitMin"),
        (-1, -1, "seedMapMin"),
        (-1, -1, "alignIntronMin"),
        (-1, -1, "alignIntronMax"),
        (-1, -1, "alignMatesGapMax"),
        (-1, -1, "alignTranscriptsPerReadNmax"),
        (-1, -1, "alignSJoverhangMin"),
        (-1, -1, "alignSJDBoverhangMin"),
        (-1, -1, "alignSJstitchMismatchNmax"),
        (-1, -1, "alignSplicedMateMapLmin"),
        (-1, -1, "alignSplicedMateMapLminOverLmate"),
        (-1, -1, "alignWindowsPerReadNmax"),
        (-1, -1, "alignTranscriptsPerWindowNmax"),
        (-1, -1, "alignEndsType"),
        (-1, -1, "alignSoftClipAtReferenceEnds"),
        (-1, -1, "alignEndsProtrude"),
        (-1, -1, "alignInsertionFlush"),
        (-1, -1, "peOverlapNbasesMin"),
        (-1, -1, "peOverlapMMp"),
        (-1, -1, "chimSegmentMin"),
        (-1, -1, "chimScoreMin"),
        (-1, -1, "chimScoreDropMax"),
        (-1, -1, "chimScoreSeparation"),
        (-1, -1, "chimScoreJunctionNonGTAG"),
        (-1, -1, "chimMainSegmentMultNmax"),
        (-1, -1, "chimJunctionOverhangMin"),
        (-1, -1, "chimOutType"),
        (-1, -1, "chimFilter"),
        (-1, -1, "chimSegmentReadGapMax"),
        (-1, -1, "chimMultimapNmax"),
        (-1, -1, "chimMultimapScoreRange"),
        (-1, -1, "chimNonchimScoreDropMin"),
        (-1, -1, "chimOutJunctionFormat"),
        (-1, -1, "sjdbFileChrStartEnd"),
        (-1, -1, "sjdbGTFfile"),
        (-1, -1, "sjdbGTFchrPrefix"),
        (-1, -1, "sjdbGTFfeatureExon"),
        (-1, -1, "sjdbGTFtagExonParentTranscript"),
        (-1, -1, "sjdbGTFtagExonParentGene"),
        (-1, -1, "sjdbGTFtagExonParentGeneName"),
        (-1, -1, "sjdbGTFtagExonParentGeneType"),
        (-1, -1, "sjdbOverhang"),
        (-1, -1, "sjdbScore"),
        (-1, -1, "sjdbInsertSave"),
        (-1, -1, "varVCFfile"),
        (-1, -1, "waspOutputMode"),
        (-1, -1, "quantMode"),
        (-1, -1, "quantTranscriptomeBAMcompression"),
        (-1, -1, "quantTranscriptomeSAMoutput"),
        (-1, -1, "twopass1readsN"),
        (-1, -1, "twopassMode"),
        (-1, -1, "soloType"),
        (-1, -1, "soloCBstart"),
        (-1, -1, "soloUMIstart"),
        (-1, -1, "soloCBlen"),
        (-1, -1, "soloUMIlen"),
        (-1, -1, "soloBarcodeReadLength"),
        (-1, -1, "soloBarcodeMate"),
        (-1, -1, "soloCBwhitelist"),
        (-1, -1, "soloStrand"),
        (-1, -1, "soloOutFileNames"),
        (-1, -1, "soloFeatures"),
        (-1, -1, "soloUMIdedup"),
        (-1, -1, "soloAdapterSequence"),
        (-1, -1, "soloAdapterMismatchesNmax"),
        (-1, -1, "soloCBmatchWLtype"),
        (-1, -1, "soloCBposition"),
        (-1, -1, "soloUMIposition"),
        (-1, -1, "soloCellFilter"),
        (-1, -1, "soloUMIfiltering"),
        (-1, -1, "soloMultiMappers"),
        (-1, -1, "soloClusterCBfile"),
        (-1, -1, "soloOutFormatFeaturesGeneField3"),
        (-1, -1, "soloInputSAMattrBarcodeSeq"),
        (-1, -1, "soloInputSAMattrBarcodeQual"),
        (-1, -1, "soloCellReadStats"),
        (-1, -1, "soloCBtype"),
    ];

    crate::parameter_info::ParametersScanState {
        par_array: par_specs
            .iter()
            .map(|(input_level, input_level_allowed, name)| {
                crate::parameter_info::ParameterScanEntry {
                    name_string: (*name).to_string(),
                    input_level_allowed: *input_level_allowed,
                    input_level: *input_level,
                    value_line: String::new(),
                }
            })
            .collect(),
        parameter_input_name: vec![
            "Default".to_string(),
            "Command-Line-Initial".to_string(),
            "Command-Line".to_string(),
            "genomeParameters.txt".to_string(),
        ],
        log_main: String::new(),
    }
}

#[doc = "Original `Parameters::inputParameters` at STAR/source/Parameters.cpp:310. Args: argInN: int, argIn: char"]
pub fn parameters_l310_parameters_inputparameters(
    p: &mut crate::parameters_chimeric::Parameters,
    arg_in: &[String],
    parameters_default: &str,
    parameter_files: &[(String, String)],
    manifest_contents: Option<&str>,
    whitelist_contents: &[(String, String)],
) -> Result<crate::parameter_info::ParametersScanState, String> {
    const PAR_NAME_PRINT_WIDTH: usize = 30;

    p.run_restart_type = 0;
    let mut state = parameters_l19_parameters_parameters();

    parameters_l1197_parameters_scanalllines(parameters_default, 0, -1, &mut state)?;
    for par in state.par_array.iter() {
        if par.input_level < 0 {
            return Err(format!(
                "BUG: DEFAULT parameter value not defined: {}",
                par.name_string
            ));
        }
    }

    p.command_line.clear();
    let mut command_line_file = String::new();
    if arg_in.len() > 1 {
        p.command_line.push_str(&arg_in[0]);
        for one_arg_in in arg_in.iter().skip(1) {
            let mut one_arg = one_arg_in.clone();
            if one_arg == "--version" {
                p.command_line_full = "STAR_VERSION".to_string();
                return Ok(state);
            }
            if let Some(found) = one_arg.find('=') {
                if one_arg.starts_with("--") {
                    let key = &one_arg[2..found];
                    let mut val = one_arg[found + 1..].to_string();
                    if val.contains(' ') || val.contains('\t') {
                        val = format!("\"{}\"", val);
                    }
                    command_line_file.push('\n');
                    command_line_file.push_str(key);
                    command_line_file.push(' ');
                    command_line_file.push_str(&val);
                } else {
                    if one_arg.contains(' ') || one_arg.contains('\t') {
                        one_arg = format!("\"{}\"", one_arg);
                    }
                    command_line_file.push(' ');
                    command_line_file.push_str(&one_arg);
                }
            } else if one_arg.starts_with("--") {
                command_line_file.push('\n');
                command_line_file.push_str(&one_arg[2..]);
            } else {
                if one_arg.contains(' ') || one_arg.contains('\t') {
                    one_arg = format!("\"{}\"", one_arg);
                }
                command_line_file.push(' ');
                command_line_file.push_str(&one_arg);
            }
            p.command_line.push(' ');
            p.command_line.push_str(&one_arg);
        }
        parameters_l1197_parameters_scanalllines(&command_line_file, 1, 2, &mut state)?;
    }

    state.log_main.push_str("##### Command Line:\n");
    state.log_main.push_str(&p.command_line);
    state.log_main.push('\n');
    state
        .log_main
        .push_str("##### Initial USER parameters from Command Line:\n");
    for par in state.par_array.iter() {
        if par.input_level == 1 {
            state.log_main.push_str(&format!(
                "{:<width$}    {}\n",
                par.name_string,
                par.value_line,
                width = PAR_NAME_PRINT_WIDTH
            ));
        }
    }

    for (file_name, file_contents) in parameter_files.iter() {
        state.parameter_input_name.push(file_name.clone());
        state.log_main.push_str(&format!(
            "##### USER parameters from user-defined parameters file {}:\n",
            file_name
        ));
        parameters_l1197_parameters_scanalllines(
            file_contents,
            state.parameter_input_name.len() as i32 - 1,
            -1,
            &mut state,
        )?;
    }

    if arg_in.len() > 1 {
        state
            .log_main
            .push_str("###### All USER parameters from Command Line:\n");
        parameters_l1197_parameters_scanalllines(&command_line_file, 2, -1, &mut state)?;
    }

    state
        .log_main
        .push_str("##### Finished reading parameters from all sources\n\n");
    state
        .log_main
        .push_str("##### Final user re-defined parameters-----------------:\n");

    let mut cl_full = String::new();
    if let Some(first) = arg_in.first() {
        cl_full.push_str(first);
    }
    for par in state.par_array.iter() {
        if par.input_level > 0 {
            state.log_main.push_str(&format!(
                "{:<width$}    {}\n",
                par.name_string,
                par.value_line,
                width = PAR_NAME_PRINT_WIDTH
            ));
            if par.name_string != "parametersFiles" {
                cl_full.push_str("   --");
                cl_full.push_str(&par.name_string);
                cl_full.push(' ');
                cl_full.push_str(&par.value_line);
            }
        }
    }
    p.command_line_full = cl_full.clone();
    state
        .log_main
        .push_str("\n-------------------------------\n##### Final effective command line:\n");
    state.log_main.push_str(&cl_full);
    state
        .log_main
        .push_str("\n----------------------------------------\n\n");

    let words = |value: &str| -> Vec<String> {
        let mut out = Vec::new();
        let mut current = String::new();
        let mut in_quotes = false;
        for ch in value.chars() {
            if ch == '"' {
                in_quotes = !in_quotes;
            } else if ch.is_whitespace() && !in_quotes {
                if !current.is_empty() {
                    out.push(std::mem::take(&mut current));
                }
            } else {
                current.push(ch);
            }
        }
        if !current.is_empty() {
            out.push(current);
        }
        out
    };
    let parse_i32 = |value: &str| -> Result<i32, String> {
        value
            .parse::<i32>()
            .map_err(|_| format!("invalid i32 parameter value: {}", value))
    };
    let parse_u32 = |value: &str| -> Result<u32, String> {
        value
            .parse::<u32>()
            .map_err(|_| format!("invalid u32 parameter value: {}", value))
    };
    let parse_u32_or_max = |value: &str| -> Result<u32, String> {
        if value == "-1" {
            Ok(u32::MAX)
        } else {
            value
                .parse::<u32>()
                .map_err(|_| format!("invalid u32 parameter value: {}", value))
        }
    };
    let parse_u64 = |value: &str| -> Result<u64, String> {
        value
            .parse::<u64>()
            .map_err(|_| format!("invalid u64 parameter value: {}", value))
    };
    let parse_u64_or_max = |value: &str| -> Result<u64, String> {
        if value == "-1" {
            Ok(u64::MAX)
        } else {
            value
                .parse::<u64>()
                .map_err(|_| format!("invalid u64 parameter value: {}", value))
        }
    };
    let parse_u64_words = |value: &str| -> Result<Vec<u64>, String> {
        value
            .split_whitespace()
            .map(|word| {
                word.parse::<u64>()
                    .map_err(|_| format!("invalid u64 parameter value: {}", word))
            })
            .collect()
    };
    let parse_u8 = |value: &str| -> Result<u8, String> {
        value
            .parse::<u8>()
            .map_err(|_| format!("invalid u8 parameter value: {}", value))
    };
    let parse_f64 = |value: &str| -> Result<f64, String> {
        value
            .parse::<f64>()
            .map_err(|_| format!("invalid f64 parameter value: {}", value))
    };
    let parse_i8 = |value: &str| -> Result<i8, String> {
        value
            .parse::<i8>()
            .map_err(|_| format!("invalid i8 parameter value: {}", value))
    };

    let mut out_tmp_dir = "-".to_string();
    let mut out_wig_type: Vec<String> = Vec::new();
    let mut out_wig_strand: Vec<String> = Vec::new();
    let mut out_wig_norm: Vec<String> = Vec::new();
    let mut quant_mode: Vec<String> = Vec::new();
    let mut quant_transcriptome_bam_compression = 1_i32;
    let mut quant_transcriptome_sam_output: Vec<String> = Vec::new();
    let mut out_sj_type: Vec<String> = Vec::new();
    let mut read_name_separator: Vec<String> = Vec::new();
    let mut out_sam_unmapped_mode: Vec<String> = Vec::new();
    let mut out_sam_filter_mode: Vec<String> = Vec::new();
    let mut out_multimapper_order = "Old_2.4".to_string();
    let mut out_filter_type = String::new();
    let mut twopass_mode = "None".to_string();
    let mut twopass1reads_n_defined = false;

    for par in state.par_array.iter() {
        if par.input_level < 0 {
            continue;
        }
        let value_raw = par.value_line.trim();
        let value_unquoted;
        let value =
            if value_raw.len() >= 2 && value_raw.starts_with('"') && value_raw.ends_with('"') {
                value_unquoted = value_raw[1..value_raw.len() - 1].to_string();
                value_unquoted.as_str()
            } else {
                value_raw
            };
        if par.name_string == "twopass1readsN" && par.input_level > 0 {
            twopass1reads_n_defined = true;
        }
        match par.name_string.as_str() {
            "versionGenome" => p.version_genome = value.to_string(),
            "runMode" => p.run_mode_in = words(value),
            "runThreadN" => p.run_thread_n = parse_i32(value)?,
            "runRNGseed" => p.run_rng_seed = parse_i32(value)?,
            "sysShell" => p.sys_shell = value.to_string(),
            "runDirPerm" => {
                p.run_dir_perm = if value == "User_RWX" {
                    0o700
                } else if value == "All_RWX" {
                    0o777
                } else {
                    return Err(format!(
                        "EXITING because of FATAL INPUT ERROR: unrecognized option in --runDirPerm={}\nSOLUTION: use one of the allowed values of --runDirPerm : 'User_RWX' or 'All_RWX' \n",
                        value
                    ));
                };
            }
            "genomeDir" => {
                p.genome_dir = value.to_string();
                p.p_ge.g_dir = value.to_string();
            }
            "genomeLoad" => p.p_ge.g_load = value.to_string(),
            "genomeFastaFiles" => p.p_ge.g_fasta_files = words(value),
            "genomeChainFiles" => p.p_ge.g_chain_files = words(value),
            "genomeFileSizes" => p.p_ge.g_file_sizes = parse_u64_words(value)?,
            "genomeType" => p.p_ge.g_type_string = value.to_string(),
            "genomeChrBinNbits" => p.p_ge.g_chr_bin_nbits = parse_u32(value)?,
            "genomeSAindexNbases" => p.p_ge.g_saindex_nbases = parse_u32(value)?,
            "genomeSAsparseD" => p.p_ge.g_sasparse_d = parse_u32(value)?,
            "genomeSuffixLengthMax" => p.p_ge.g_suffix_length_max = parse_u32_or_max(value)?,
            "genomeTransformType" => p.p_ge.transform.type_string = value.to_string(),
            "genomeTransformVCF" => p.p_ge.transform.vcf_file = value.to_string(),
            "genomeTransformOutput" => p.p_ge.transform.output = words(value),
            "genomeChrSetMitochondrial" => p.p_ge.chr_set_mito_strings = words(value),
            "sjdbFileChrStartEnd" => p.p_ge.sjdb_file_chr_start_end = words(value),
            "sjdbGTFfile" => p.p_ge.sjdb_gtf_file = value.to_string(),
            "sjdbGTFchrPrefix" => p.p_ge.sjdb_gtf_chr_prefix = value.to_string(),
            "sjdbGTFfeatureExon" => p.p_ge.sjdb_gtf_feature_exon = value.to_string(),
            "sjdbGTFtagExonParentTranscript" => {
                p.p_ge.sjdb_gtf_tag_exon_parent_transcript = value.to_string()
            }
            "sjdbGTFtagExonParentGene" => p.p_ge.sjdb_gtf_tag_exon_parent_gene = value.to_string(),
            "sjdbGTFtagExonParentGeneName" => {
                p.p_ge.sjdb_gtf_tag_exon_parent_gene_name = words(value)
            }
            "sjdbGTFtagExonParentGeneType" => {
                p.p_ge.sjdb_gtf_tag_exon_parent_gene_type = words(value)
            }
            "sjdbOverhang" => p.p_ge.sjdb_overhang = parse_u32(value)?,
            "sjdbScore" => p.p_ge.sjdb_score = parse_i32(value)?,
            "sjdbInsertSave" => p.p_ge.sjdb_insert_save = value.to_string(),
            "readFilesIn" => p.read_files_in = words(value),
            "readFilesType" => p.read_files_type = words(value),
            "readFilesPrefix" => p.read_files_prefix = value.to_string(),
            "readFilesManifest" => p.read_files_manifest = words(value),
            "readFilesSAMattrKeep" => p.read_files_sam_attr_keep_in = words(value),
            "readMatesLengthsIn" => p.read_mates_lengths_in = value.to_string(),
            "readFilesCommand" => {
                p.read_files_command = words(value);
                p.read_files_command_string = value.to_string();
            }
            "inputBAMfile" => p.input_bam_file = value.to_string(),
            "bamRemoveDuplicatesType" => p.bam_remove_duplicates_type = value.to_string(),
            "bamRemoveDuplicatesMate2basesN" => {
                p.bam_remove_duplicates_mate2bases_n = parse_u32(value)?
            }
            "readMapNumber" => p.read_map_number = parse_u64_or_max(value)?,
            "readQualityScoreBase" => {
                p.read_quality_score_base = parse_u8(value)?;
                p.p_solo.qs_base = parse_i8(value)?;
            }
            "readNameSeparator" => read_name_separator = words(value),
            "clipAdapterType" => p.p_clip.adapter_type = words(value),
            "clip5pNbases" => {
                p.p_clip.in_[0].n = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<u32>()
                            .map_err(|_| format!("invalid u32 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "clip3pNbases" => {
                p.p_clip.in_[1].n = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<u32>()
                            .map_err(|_| format!("invalid u32 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "clip5pAfterAdapterNbases" => {
                p.p_clip.in_[0].n_after_ad = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<u32>()
                            .map_err(|_| format!("invalid u32 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "clip3pAfterAdapterNbases" => {
                p.p_clip.in_[1].n_after_ad = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<u32>()
                            .map_err(|_| format!("invalid u32 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "clip5pAdapterSeq" => p.p_clip.in_[0].ad_seq = words(value),
            "clip3pAdapterSeq" => p.p_clip.in_[1].ad_seq = words(value),
            "clip5pAdapterMMp" => {
                p.p_clip.in_[0].ad_mmp = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<f64>()
                            .map_err(|_| format!("invalid f64 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "clip3pAdapterMMp" => {
                p.p_clip.in_[1].ad_mmp = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<f64>()
                            .map_err(|_| format!("invalid f64 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "outFileNamePrefix" => p.out_file_name_prefix = value.to_string(),
            "outStd" => p.out_std = value.to_string(),
            "outReadsUnmapped" => p.out_reads_unmapped = value.to_string(),
            "outTmpDir" => out_tmp_dir = value.to_string(),
            "outTmpKeep" => p.out_tmp_keep = value.to_string(),
            "outSAMtype" => p.out_sam_type = words(value),
            "outSAMmode" => p.out_sam_mode = value.to_string(),
            "outSAMattributes" => p.out_sam_attributes = words(value),
            "outMultimapperOrder" => out_multimapper_order = value.to_string(),
            "outBAMcompression" => p.out_bam_compression = parse_i32(value)?,
            "outBAMsortingThreadN" => p.out_bam_sorting_thread_n = parse_u32(value)?,
            "outBAMsortingBinsN" => p.out_bam_sorting_bins_n = parse_u32(value)?,
            "outSAMattrRGline" => p.out_sam_attr_rgline = words(value),
            "outSAMheaderHD" => p.out_sam_header_hd = words(value),
            "outSAMheaderPG" => p.out_sam_header_pg = words(value),
            "outSAMheaderCommentFile" => p.out_sam_header_comment_file = value.to_string(),
            "outSAMfilter" => out_sam_filter_mode = words(value),
            "outSAMattrIHstart" => p.out_sam_attr_ih_start = parse_u32(value)?,
            "outSAMtlen" => p.out_sam_tlen = parse_i32(value)?,
            "outSAMmapqUnique" => p.out_sam_mapq_unique = parse_i32(value)?,
            "outSAMmultNmax" => p.out_sam_mult_nmax = parse_u64_or_max(value)?,
            "outSAMprimaryFlag" => p.out_sam_primary_flag = value.to_string(),
            "outSAMflagOR" => p.out_sam_flag_or = parse_u32(value)? as u16,
            "outSAMflagAND" => p.out_sam_flag_and = parse_u32(value)? as u16,
            "outSAMorder" => p.out_sam_order = value.to_string(),
            "outSAMstrandField" => {
                p.out_sam_strand_field_type = if value == "None" {
                    0
                } else if value == "intronMotif" {
                    1
                } else {
                    return Err(format!(
                        "EXITING because of fatal INPUT error: unrecognized option in outSAMstrandField={}\nSOLUTION: use one of the allowed values of --outSAMstrandField : None or intronMotif \n",
                        value
                    ));
                };
            }
            "outSAMunmapped" => out_sam_unmapped_mode = words(value),
            "outSAMreadID" => {
                p.out_sam_read_id = value.to_string();
                p.out_sam_read_id_number = value == "Number";
            }
            "outQSconversionAdd" => p.out_qs_conversion_add = parse_i32(value)?,
            "outSJtype" => out_sj_type = words(value),
            "outSJfilterReads" => p.out_sj_filter_reads = value.to_string(),
            "outSJfilterCountUniqueMin" => {
                p.out_sjfilter_count_unique_min = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<i32>()
                            .map_err(|_| format!("invalid i32 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "outSJfilterCountTotalMin" => {
                p.out_sjfilter_count_total_min = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<i32>()
                            .map_err(|_| format!("invalid i32 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "outSJfilterOverhangMin" => {
                p.out_sjfilter_overhang_min = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<i32>()
                            .map_err(|_| format!("invalid i32 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "outSJfilterDistToOtherSJmin" => {
                p.out_sjfilter_dist_to_other_sj_min = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<i32>()
                            .map_err(|_| format!("invalid i32 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "outSJfilterIntronMaxVsReadN" => {
                p.out_sjfilter_intron_max_vs_read_n = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<u32>()
                            .map_err(|_| format!("invalid u32 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "outFilterType" => {
                p.out_filter_type = value.to_string();
                out_filter_type = value.to_string();
            }
            "outFilterIntronMotifs" => p.out_filter_intron_motifs = value.to_string(),
            "outFilterIntronStrands" => p.out_filter_intron_strands = value.to_string(),
            "outFilterMultimapNmax" => p.out_filter_multimap_nmax = parse_u64(value)?,
            "outFilterMultimapScoreRange" => p.out_filter_multimap_score_range = parse_i32(value)?,
            "outFilterMismatchNmax" => p.out_filter_mismatch_nmax = parse_u32(value)?,
            "outFilterMismatchNoverLmax" => p.out_filter_mismatch_nover_lmax = parse_f64(value)?,
            "outFilterMismatchNoverReadLmax" => {
                p.out_filter_mismatch_nover_read_lmax = parse_f64(value)?
            }
            "outFilterScoreMin" => p.out_filter_score_min = parse_i32(value)?,
            "outFilterScoreMinOverLread" => p.out_filter_score_min_over_lread = parse_f64(value)?,
            "outFilterMatchNmin" => p.out_filter_match_nmin = parse_u32(value)?,
            "outFilterMatchNminOverLread" => p.out_filter_match_nmin_over_lread = parse_f64(value)?,
            "limitGenomeGenerateRAM" => p.limit_genome_generate_ram = parse_u64(value)?,
            "limitBAMsortRAM" => p.limit_bam_sort_ram = parse_u64(value)?,
            "limitOutSAMoneReadBytes" => p.limit_out_sam_one_read_bytes = parse_u64(value)?,
            "limitOutSJoneRead" => p.limit_out_sj_one_read = parse_u64(value)?,
            "limitOutSJcollapsed" => p.limit_out_sj_collapsed = parse_u64(value)?,
            "limitSjdbInsertNsj" => p.limit_sjdb_insert_nsj = parse_u32(value)?,
            "limitNreadsSoft" => p.limit_nreads_soft = parse_u64_or_max(value)?,
            "limitIObufferSize" => {
                let values = parse_u64_words(value)?;
                if values.len() != 2 {
                    return Err("EXITING because of FATAL input ERROR: --limitIObufferSize requires 2 numbers since 2.7.9a.\nSOLUTION: specify 2 numbers in --limitIObufferSize : size of input and output buffers in bytes.\n".to_string());
                }
                p.chunk_in_size_bytes_array = values[0];
                p.chunk_out_bam_size_bytes = values[1];
            }
            "seedSearchLmax" => p.seed_search_lmax = parse_u32(value)?,
            "seedSearchStartLmax" => p.seed_search_start_lmax = parse_u32(value)?,
            "seedSearchStartLmaxOverLread" => {
                p.seed_search_start_lmax_over_lread = parse_f64(value)?
            }
            "seedSplitMin" => p.seed_split_min = parse_u32(value)?,
            "seedMapMin" => p.seed_map_min = parse_u32(value)?,
            "seedMultimapNmax" => p.seed_multimap_nmax = parse_u32(value)?,
            "seedPerReadNmax" => p.seed_per_read_nmax = parse_u32(value)?,
            "seedPerWindowNmax" => p.seed_per_window_nmax = parse_u32(value)?,
            "seedNoneLociPerWindow" => p.seed_none_loci_per_window = parse_u32(value)?,
            "scoreGap" => p.score_gap = parse_i32(value)?,
            "scoreGapNoncan" => p.score_gap_noncan = parse_i32(value)?,
            "scoreGapGCAG" => p.score_gap_gcag = parse_i32(value)?,
            "scoreGapATAC" => p.score_gap_atac = parse_i32(value)?,
            "scoreDelBase" => p.score_del_base = parse_i32(value)?,
            "scoreDelOpen" => p.score_del_open = parse_i32(value)?,
            "scoreInsBase" => p.score_ins_base = parse_i32(value)?,
            "scoreInsOpen" => p.score_ins_open = parse_i32(value)?,
            "scoreStitchSJshift" => p.score_stitch_sj_shift = parse_i32(value)?,
            "scoreGenomicLengthLog2scale" => p.score_genomic_length_log2scale = parse_f64(value)?,
            "winBinNbits" => p.win_bin_nbits = parse_u32(value)?,
            "winReadCoverageRelativeMin" => p.win_read_coverage_relative_min = parse_f64(value)?,
            "winReadCoverageBasesMin" => p.win_read_coverage_bases_min = parse_u32(value)?,
            "winAnchorDistNbins" => p.win_anchor_dist_nbins = parse_u32(value)?,
            "winFlankNbins" => p.win_flank_nbins = parse_u32(value)?,
            "winAnchorMultimapNmax" => p.win_anchor_multimap_nmax = parse_u32(value)?,
            "alignIntronMin" => p.align_intron_min = parse_u32(value)?,
            "alignIntronMax" => p.align_intron_max = parse_u32(value)?,
            "alignMatesGapMax" => p.align_mates_gap_max = parse_u32(value)?,
            "alignSJoverhangMin" => p.align_sj_overhang_min = parse_u32(value)?,
            "alignSJDBoverhangMin" => p.align_sjdb_overhang_min = parse_u32(value)?,
            "alignSJstitchMismatchNmax" => {
                p.align_sj_stitch_mismatch_nmax = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<i32>()
                            .map_err(|_| format!("invalid i32 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "alignSplicedMateMapLmin" => p.align_spliced_mate_map_lmin = parse_u32(value)?,
            "alignSplicedMateMapLminOverLmate" => {
                p.align_spliced_mate_map_lmin_over_lmate = parse_f64(value)?
            }
            "alignTranscriptsPerReadNmax" => p.align_transcripts_per_read_nmax = parse_u32(value)?,
            "alignTranscriptsPerWindowNmax" => {
                p.align_transcripts_per_window_nmax = parse_u32(value)?
            }
            "alignWindowsPerReadNmax" => p.align_windows_per_read_nmax = parse_u32(value)?,
            "alignEndsType" => p.align_ends_type.in_ = value.to_string(),
            "alignEndsProtrude" => p.align_ends_protrude.in_ = words(value),
            "alignInsertionFlush" => p.align_insertion_flush.in_ = value.to_string(),
            "alignSoftClipAtReferenceEnds" => {
                p.align_soft_clip_at_reference_ends_yes = value == "Yes"
            }
            "peOverlapNbasesMin" => p.pe_overlap_nbases_min = parse_u32(value)?,
            "peOverlapMMp" => p.pe_overlap_mmp = parse_f64(value)?,
            "chimSegmentMin" => p.p_ch.segment_min = parse_u64(value)?,
            "chimJunctionOverhangMin" => p.p_ch.junction_overhang_min = parse_u64(value)?,
            "chimSegmentReadGapMax" => p.p_ch.segment_read_gap_max = parse_u64(value)?,
            "chimScoreMin" => p.p_ch.score_min = parse_i32(value)?,
            "chimScoreDropMax" => p.p_ch.score_drop_max = parse_i32(value)?,
            "chimScoreSeparation" => p.p_ch.score_separation = parse_i32(value)?,
            "chimScoreJunctionNonGTAG" => p.p_ch.score_junction_non_gtag = parse_i32(value)?,
            "chimMainSegmentMultNmax" => p.p_ch.main_segment_mult_nmax = parse_u32(value)?,
            "chimOutType" => p.p_ch.out_type = words(value),
            "chimOutJunctionFormat" => {
                p.p_ch.out_junction_format = value
                    .split_whitespace()
                    .map(|word| {
                        word.parse::<i32>()
                            .map_err(|_| format!("invalid i32 parameter value: {}", word))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            "chimFilter" => p.p_ch.filter_string_in = words(value),
            "chimMultimapNmax" => p.p_ch.multimap_nmax = parse_u64(value)?,
            "chimMultimapScoreRange" => p.p_ch.multimap_score_range = parse_u64(value)?,
            "chimNonchimScoreDropMin" => p.p_ch.nonchim_score_drop_min = parse_u64(value)?,
            "varVCFfile" => {
                p.var_yes = value != "-";
            }
            "quantMode" => quant_mode = words(value),
            "quantTranscriptomeBAMcompression" => {
                quant_transcriptome_bam_compression = parse_i32(value)?
            }
            "quantTranscriptomeSAMoutput" => quant_transcriptome_sam_output = words(value),
            "twopassMode" => twopass_mode = value.to_string(),
            "twopass1readsN" => p.two_pass_pass1reads_n = parse_u64_or_max(value)?,
            "soloType" => {
                p.p_solo.type_str = value.to_string();
                p.p_solo.yes = value != "None";
            }
            "soloCBstart" => p.p_solo.cb_s = parse_u32(value)?,
            "soloUMIstart" => p.p_solo.umi_s = parse_u32(value)?,
            "soloCBlen" => p.p_solo.cb_l = parse_u32(value)?,
            "soloUMIlen" => p.p_solo.umi_l = parse_u32(value)?,
            "soloBarcodeReadLength" => p.p_solo.b_l = parse_u32(value)?,
            "soloBarcodeMate" => p.p_solo.barcode_read = parse_u32(value)?,
            "soloCBwhitelist" => p.p_solo.solo_cb_whitelist = words(value),
            "soloCBposition" => p.p_solo.cb_position_str = words(value),
            "soloUMIposition" => p.p_solo.umi_position_str = value.to_string(),
            "soloAdapterSequence" => p.p_solo.adapter_seq = value.to_string(),
            "soloAdapterMismatchesNmax" => p.p_solo.adapter_mismatches_nmax = parse_u32(value)?,
            "soloCBmatchWLtype" => p.p_solo.cb_match_wl.type_ = value.to_string(),
            "soloCBtype" => {
                p.p_solo.cb_type_type = match value {
                    "Sequence" => 1,
                    "String" => 2,
                    _ => 0,
                };
            }
            "soloInputSAMattrBarcodeSeq" => p.p_solo.sam_attr_barcode_seq = words(value),
            "soloInputSAMattrBarcodeQual" => p.p_solo.sam_attr_barcode_qual = words(value),
            "soloStrand" => {
                p.p_solo.strand = match value {
                    "Unstranded" => 0,
                    "Forward" => 1,
                    "Reverse" => -1,
                    _ => -1,
                };
            }
            "soloFeatures" => {
                p.p_solo.features.clear();
                for feature in words(value) {
                    p.p_solo.features.push(match feature.as_str() {
                        "SJ" => SOLO_FEATURE_SJ as u32,
                        "Transcript3p" => SOLO_FEATURE_TRANSCRIPT3P as u32,
                        "GeneFull" => SOLO_FEATURE_GENE_FULL as u32,
                        "GeneFull_ExonOverIntron" => {
                            SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON as u32
                        }
                        "GeneFull_Ex50pAS" => SOLO_FEATURE_GENE_FULL_EX50P_AS as u32,
                        "Gene" => SOLO_FEATURE_GENE as u32,
                        "VelocytoSimple" => SOLO_FEATURE_VELOCYTO_SIMPLE as u32,
                        "Velocyto" => SOLO_FEATURE_VELOCYTO as u32,
                        other => {
                            return Err(format!(
                                "EXITING because of fatal PARAMETERS error: unrecognized option(s) in --soloFeatures {}\nSOLUTION: use allowed option: Gene",
                                other
                            ));
                        }
                    });
                }
            }
            "soloUMIdedup" => p.p_solo.umi_dedup.types_in = words(value),
            "soloUMIfiltering" => p.p_solo.umi_filtering.type_ = words(value),
            "soloMultiMappers" => p.p_solo.multi_map.types_in = words(value),
            "soloOutFileNames" => p.p_solo.out_file_names = words(value),
            "soloOutFormatFeaturesGeneField3" => {
                p.p_solo.out_format_features_gene_field3 = value.to_string()
            }
            "soloCellFilter" => p.p_solo.cell_filter.type_ = words(value),
            "soloCellReadStats" => p.p_solo.read_stats_type = value.to_string(),
            "soloClusterCBfile" => p.p_solo.cluster_cb_file = value.to_string(),
            "outWigType" => out_wig_type = words(value),
            "outWigStrand" => out_wig_strand = words(value),
            "outWigReferencesPrefix" => p.out_wig_references_prefix = value.to_string(),
            "outWigNorm" => out_wig_norm = words(value),
            "waspOutputMode" => p.wasp_output_mode = value.to_string(),
            _ => {}
        }
    }

    p.max_nsplit = 10;
    p.i_read_all = 0;
    parametersgenome_l5_parametersgenome_initialize(&mut p.p_ge)?;

    if p.limit_genome_generate_ram == 0 {
        return Err("EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM=0\nSOLUTION: please specify a >0 value for limitGenomeGenerateRAM\n".to_string());
    } else if p.limit_genome_generate_ram > 1_000_000_000_000 {
        state.log_main.push_str(&format!(
            "WARNING: specified limitGenomeGenerateRAM={} bytes appears to be too large, if you do not have enough memory the code will crash!\n",
            p.limit_genome_generate_ram
        ));
    }

    if p.out_file_name_prefix.is_empty() {
        p.out_file_name_prefix = "./".to_string();
    }
    p.out_file_tmp = if out_tmp_dir == "-" {
        format!("{}_STARtmp/", p.out_file_name_prefix)
    } else {
        format!("{}/", out_tmp_dir)
    };
    match p.out_std.as_str() {
        "Log" | "SAM" | "BAM_Unsorted" | "BAM_SortedByCoordinate" | "BAM_Quant" => {}
        other => {
            return Err(format!(
                "EXITING because of FATAL PARAMETER error: outStd={} is not a valid value of the parameter\nSOLUTION: provide a valid value fot outStd: Log / SAM / BAM_Unsorted / BAM_SortedByCoordinate",
                other
            ));
        }
    }

    p.out_wig_flags.yes = false;
    if out_wig_type.first().map(|s| s.as_str()) == Some("bedGraph") {
        p.out_wig_flags.yes = true;
        p.out_wig_flags.format = 0;
    } else if out_wig_type.first().map(|s| s.as_str()) == Some("wiggle") {
        p.out_wig_flags.yes = true;
        p.out_wig_flags.format = 1;
    } else if out_wig_type.first().map(|s| s.as_str()) != Some("None") {
        return Err(format!(
            "EXITING because of FATAL INPUT ERROR: unrecognized option in --outWigType={}\nSOLUTION: use one of the allowed values of --outWigType : 'None' or 'bedGraph' \n",
            out_wig_type.first().map(String::as_str).unwrap_or("")
        ));
    }
    p.out_wig_flags.strand = match out_wig_strand.first().map(|s| s.as_str()) {
        Some("Stranded") => true,
        Some("Unstranded") => false,
        Some(other) => {
            return Err(format!(
                "EXITING because of FATAL INPUT ERROR: unrecognized option in --outWigStrand={}\nSOLUTION: use one of the allowed values of --outWigStrand : 'Stranded' or 'Unstranded' \n",
                other
            ));
        }
        None => false,
    };
    p.out_wig_flags.type_ = if out_wig_type.len() == 1 {
        0
    } else if out_wig_type.get(1).map(|s| s.as_str()) == Some("read1_5p") {
        1
    } else if out_wig_type.get(1).map(|s| s.as_str()) == Some("read2") {
        2
    } else {
        return Err(format!(
            "EXITING because of FATAL INPUT ERROR: unrecognized second option in --outWigType={}\nSOLUTION: use one of the allowed values of --outWigType : 'read1_5p' \n",
            out_wig_type.get(1).map(String::as_str).unwrap_or("")
        ));
    };
    p.out_wig_flags.norm = match out_wig_norm.first().map(|s| s.as_str()) {
        Some("None") => 0,
        Some("RPM") => 1,
        Some(other) => {
            return Err(format!(
                "EXITING because of fatal parameter ERROR: unrecognized option in --outWigNorm={}\nSOLUTION: use one of the allowed values of --outWigNorm : 'None' or 'RPM' \n",
                other
            ));
        }
        None => 0,
    };

    p.bam_remove_duplicates_yes = false;
    p.bam_remove_duplicates_mark_multi = false;
    match p.bam_remove_duplicates_type.as_str() {
        "UniqueIdentical" => {
            p.bam_remove_duplicates_yes = true;
            p.bam_remove_duplicates_mark_multi = true;
        }
        "UniqueIdenticalNotMulti" => {
            p.bam_remove_duplicates_yes = true;
        }
        "-" => {}
        other => {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: unrecognized option in of --bamRemoveDuplicatesType={}\nSOLUTION: use allowed option: - or UniqueIdentical or UniqueIdenticalNotMulti",
                other
            ));
        }
    }

    p.out_sam_bool = false;
    p.out_bam_unsorted = false;
    p.out_bam_coord = false;
    if p.out_sam_mode != "None" {
        match p.out_sam_type.first().map(|s| s.as_str()) {
            Some("BAM") => {
                if p.out_sam_type.len() < 2 {
                    return Err("EXITING because of fatal PARAMETER error: missing BAM option\nSOLUTION: re-run STAR with one of the allowed values of --outSAMtype BAM Unsorted OR SortedByCoordinate OR both\n".to_string());
                }
                for type1 in p.out_sam_type.iter().skip(1) {
                    if type1 == "Unsorted" {
                        p.out_bam_unsorted = true;
                    } else if type1 == "SortedByCoordinate" {
                        p.out_bam_coord = true;
                    } else {
                        return Err(format!(
                            "EXITING because of fatal input ERROR: unknown value for outSAMtype: {}\nSOLUTION: re-run STAR with one of the allowed values of --outSAMtype BAM Unsorted or SortedByCoordinate or both\n",
                            type1
                        ));
                    }
                }
            }
            Some("SAM") => {
                if p.out_sam_type.len() > 1 {
                    return Err(format!(
                        "EXITING because of fatal PARAMETER error: --outSAMtype SAM can cannot be combined with {} or any other options\nSOLUTION: re-run STAR with with '--outSAMtype SAM' only, or with --outSAMtype BAM Unsorted|SortedByCoordinate\n",
                        p.out_sam_type[1]
                    ));
                }
                p.out_sam_bool = true;
            }
            Some("None") => {}
            Some(other) => {
                return Err(format!(
                    "EXITING because of fatal input ERROR: unknown value for the first word of outSAMtype: {}\nSOLUTION: re-run STAR with one of the allowed values of outSAMtype: BAM or SAM \n",
                    other
                ));
            }
            None => {}
        }
    }
    if p.out_bam_coord {
        p.out_bam_sorting_thread_nactual = if p.out_bam_sorting_thread_n == 0 {
            std::cmp::min(6, p.run_thread_n.max(0) as u32)
        } else {
            p.out_bam_sorting_thread_n
        };
        p.out_bam_coord_nbins = std::cmp::max(
            p.out_bam_sorting_thread_nactual * 3,
            p.out_bam_sorting_bins_n,
        );
        p.out_bam_sorting_bin_start = vec![0; p.out_bam_coord_nbins as usize];
        if !p.out_bam_sorting_bin_start.is_empty() {
            p.out_bam_sorting_bin_start[0] = 1;
        }
        p.out_bam_sort_tmp_dir = format!("{}/BAMsort/", p.out_file_tmp);
        if p.limit_bam_sort_ram == 0 && p.p_ge.g_load != "NoSharedMemory" {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: limitBAMsortRAM=0 (default) cannot be used with --genomeLoad={}, or any other shared memory options\nSOLUTION: please use default --genomeLoad NoSharedMemory, \n        OR specify --limitBAMsortRAM the amount of RAM (bytes) that can be allocated for BAM sorting in addition to shared memory allocated for the genome.\n        --limitBAMsortRAM typically has to be > 10000000000 (i.e 10GB).\n",
                p.p_ge.g_load
            ));
        }
        if p.limit_bam_sort_ram == 0 {
            state.log_main.push_str(
                "WARNING: --limitBAMsortRAM=0, will use genome size as RAM limit for BAM sorting\n",
            );
        }
    }
    if !p.out_sam_bool && p.out_sam_order == "PairedKeepInputOrder" {
        return Err("EXITING: fatal input ERROR: --outSAMorder=PairedKeepInputOrder is presently only compatible with SAM output, i.e. default --outSMAtype SAM\nSOLUTION: re-run STAR without --outSAMorder=PairedKeepInputOrder, or with --outSAMorder=PairedKeepInputOrder --outSMAtype SAM .\n".to_string());
    }

    if p.out_sam_mapq_unique < 0 || p.out_sam_mapq_unique > 255 {
        return Err(format!(
            "EXITING because of FATAL input ERROR: out of range value for outSAMmapqUnique={}\nSOLUTION: specify outSAMmapqUnique within the range of 0 to 255\n",
            p.out_sam_mapq_unique
        ));
    }

    p.out_sam_filter_keep_only_added_references = false;
    p.out_sam_filter_keep_all_added_references = false;
    p.out_sam_filter_yes = true;
    match out_sam_filter_mode.first().map(|s| s.as_str()) {
        Some("KeepOnlyAddedReferences") => p.out_sam_filter_keep_only_added_references = true,
        Some("KeepAllAddedReferences") => p.out_sam_filter_keep_all_added_references = true,
        Some("None") | None => p.out_sam_filter_yes = false,
        Some(other) => {
            return Err(format!(
                "EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --outSAMfilter: {}\nSOLUTION: specify one of the allowed values: KeepOnlyAddedReferences or None\n",
                other
            ));
        }
    }
    if (p.out_sam_filter_keep_only_added_references || p.out_sam_filter_keep_all_added_references)
        && p.p_ge
            .g_fasta_files
            .first()
            .map(|s| s.as_str())
            .unwrap_or("-")
            == "-"
    {
        return Err("EXITING because of FATAL INPUT ERROR: --outSAMfilter KeepOnlyAddedReferences OR KeepAllAddedReferences options can only be used if references are added on-the-fly with --genomeFastaFiles\nSOLUTION: use default --outSAMfilter None, OR add references with --genomeFataFiles\n".to_string());
    }

    p.out_multimapper_order_random = false;
    match out_multimapper_order.as_str() {
        "Old_2.4" => {}
        "Random" => p.out_multimapper_order_random = true,
        "SortedByCoordinate" => {}
        other => {
            return Err(format!(
                "EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --outMultimapperOrder: {}\nSOLUTION: specify one of the allowed values: Old_2.4 or SortedByCoordinate or Random\n",
                other
            ));
        }
    }
    p.out_sam_mult_nmax_is_limited = p.out_sam_mult_nmax != u64::MAX;

    p.quant_yes = false;
    p.quant_ge_count_yes = false;
    p.quant_tr_sam_yes = false;
    p.quant_tr_sam_bam_yes = false;
    if quant_mode.first().map(|s| s.as_str()) != Some("-") {
        p.quant_yes = true;
        for mode in quant_mode.iter() {
            if mode == "TranscriptomeSAM" {
                p.quant_tr_sam_yes = true;
            } else if mode == "GeneCounts" {
                p.quant_ge_count_yes = true;
            } else if mode != "-" {
                return Err(format!(
                    "EXITING because of fatal INPUT error: unrecognized option in --quantMode={}\nSOLUTION: use one of the allowed values of --quantMode : TranscriptomeSAM or GeneCounts or - .\n",
                    mode
                ));
            }
        }
    }
    p.quant_tr_sam_bam_yes = p.quant_tr_sam_yes && quant_transcriptome_bam_compression != -2;
    p.quant_tr_sam_indel = true;
    p.quant_tr_sam_single_end = true;
    p.quant_tr_sam_soft_clip = true;
    match quant_transcriptome_sam_output.first().map(|s| s.as_str()) {
        Some("BanSingleEnd_BanIndels_ExtendSoftclip") | None => {
            p.quant_tr_sam_indel = false;
            p.quant_tr_sam_single_end = false;
            p.quant_tr_sam_soft_clip = false;
        }
        Some("BanSingleEnd") => {
            p.quant_tr_sam_single_end = false;
        }
        Some("BanSingleEnd_ExtendSoftclip") => {
            p.quant_tr_sam_single_end = false;
            p.quant_tr_sam_soft_clip = false;
        }
        Some(other) => {
            return Err(format!(
                "EXITING because of fatal INPUT error: unrecognized option in --quantTranscriptomeSAMoutput={}\nSOLUTION: use one of the allowed values of --quantTranscriptomeSAMoutput : BanSingleEnd_BanIndels_ExtendSoftclip, BanSingleEnd, BanSingleEnd_ExtendSoftclip .\n",
                other
            ));
        }
    }
    p.quant_gene_full_yes = false;
    p.quant_gene_yes = false;

    if twopass1reads_n_defined && twopass_mode == "None" {
        return Err("EXITING because of fatal PARAMETERS error: --twopass1readsN is defined, but --twoPassMode is not defined\nSOLUTION: to activate the 2-pass mode, use --twopassMode Basic".to_string());
    }
    p.two_pass_yes = false;
    p.two_pass_pass2 = false;
    if twopass_mode != "None" {
        let run_mode = p
            .run_mode_in
            .first()
            .map(String::as_str)
            .unwrap_or("alignReads");
        if run_mode != "alignReads" {
            return Err("EXITING because of fatal PARAMETERS error: 2-pass mapping option  can only be used with --runMode alignReads\nSOLUTION: remove --twopassMode option".to_string());
        }
        if twopass_mode != "Basic" {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: unrecognized value of --twopassMode={}\nSOLUTION: for the 2-pass mode, use allowed values --twopassMode: Basic",
                twopass_mode
            ));
        }
        if p.two_pass_pass1reads_n == 0 {
            return Err("EXITING because of fatal PARAMETERS error: --twopass1readsN = 0 in the 2-pass mode\nSOLUTION: for the 2-pass mode, specify --twopass1readsN > 0. Use a very large number or -1 to map all reads in the 1st pass.\n".to_string());
        }
        if p.p_ge.g_load != "NoSharedMemory" {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: 2-pass method is not compatible with --genomeLoad {}\nSOLUTION: re-run STAR with --genomeLoad NoSharedMemory ; this is the only option compatible with --twopassMode Basic .\n",
                p.p_ge.g_load
            ));
        }
        p.two_pass_yes = true;
        p.two_pass_dir = format!("{}_STARpass1/", p.out_file_name_prefix);
    }

    let sam_attr_log = parameters_samattributes_l4_parameters_samattributes(p)?;
    state.log_main.push_str(&sam_attr_log);

    let read_files_log =
        parameters_readfilesinit_l8_parameters_readfilesinit(p, manifest_contents)?;
    state.log_main.push_str(&read_files_log);

    p.read_in_open = vec![true; p.read_files_in.len()];
    p.read_files_command_pid = vec![0; p.read_files_in.len()];

    let parameters_for_clip = p.clone();
    parametersclip_initialize_l6_parametersclip_initialize(&mut p.p_clip, &parameters_for_clip)?;

    let solo_log = parameterssolo_l10_parameterssolo_initialize(p, whitelist_contents, 0)?;
    state.log_main.push_str(&solo_log);

    p.align_ends_type.ext = [[false; 2]; 2];
    match p.align_ends_type.in_.as_str() {
        "EndToEnd" => p.align_ends_type.ext = [[true; 2]; 2],
        "Extend5pOfRead1" => p.align_ends_type.ext[0][0] = true,
        "Extend5pOfReads12" => {
            p.align_ends_type.ext[0][0] = true;
            p.align_ends_type.ext[1][0] = true;
        }
        "Extend3pOfRead1" => p.align_ends_type.ext[0][1] = true,
        "Local" => {}
        other => {
            return Err(format!(
                "EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --alignEndsType: {}\nSOLUTION: re-run STAR with --alignEndsType Local OR EndToEnd OR Extend5pOfRead1 OR Extend3pOfRead1\n",
                other
            ));
        }
    }

    p.genome_num_to_nt = b"ACGTN".to_vec();
    p.sjdb_insert_yes = p.p_ge.sjdb_file_chr_start_end.first().map(|s| s.as_str()) != Some("-")
        || p.p_ge.sjdb_gtf_file != "-"
        || p.two_pass_yes;
    if p.sjdb_insert_yes {
        p.sjdb_insert_out_dir = format!("{}_STARgenome/", p.out_file_name_prefix);
    }

    p.read_name_separator_char.clear();
    for separator in read_name_separator.iter() {
        if separator == "space" {
            p.read_name_separator_char.push(' ');
        } else if separator == "none" {
        } else if separator.len() == 1 {
            p.read_name_separator_char
                .push(separator.chars().next().unwrap());
        } else {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: unrecognized value of --readNameSeparator={}\nSOLUTION: use allowed values: space OR single characters",
                separator
            ));
        }
    }

    p.out_sam_unmapped_within = false;
    p.out_sam_unmapped_keep_pairs = false;
    match out_sam_unmapped_mode.as_slice() {
        [one] if one == "None" => {}
        [one] if one == "Within" => p.out_sam_unmapped_within = true,
        [one, two] if one == "Within" && two == "KeepPairs" => {
            p.out_sam_unmapped_within = true;
            if p.read_nmates == 2 {
                p.out_sam_unmapped_keep_pairs = true;
            }
        }
        _ => {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: unrecognized option for --outSAMunmapped= {}\nSOLUTION: use allowed options: None OR Within OR Within KeepPairs",
                out_sam_unmapped_mode.join(" ")
            ));
        }
    }

    p.align_ends_protrude.n_bases_max = p
        .align_ends_protrude
        .in_
        .first()
        .map(|v| v.parse::<i32>())
        .transpose()
        .map_err(|_| "invalid i32 parameter value for alignEndsProtrude".to_string())?
        .unwrap_or(0);
    p.align_ends_protrude.concordant_pair = false;
    if p.align_ends_protrude.n_bases_max > 0 {
        match p.align_ends_protrude.in_.get(1).map(|s| s.as_str()) {
            Some("ConcordantPair") => p.align_ends_protrude.concordant_pair = true,
            Some("DiscordantPair") => p.align_ends_protrude.concordant_pair = false,
            Some(other) => {
                return Err(format!(
                    "EXITING because of fatal PARAMETERS error: unrecognized option in of --alignEndsProtrude={}\nSOLUTION: use allowed option: ConcordantPair or DiscordantPair",
                    other
                ));
            }
            None => {}
        }
    }

    p.align_insertion_flush.flush_right = match p.align_insertion_flush.in_.as_str() {
        "None" => false,
        "Right" => true,
        other => {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: unrecognized option in of --alignInsertionFlush={}\nSOLUTION: use allowed option: None or Right",
                other
            ));
        }
    };
    p.align_soft_clip_at_reference_ends_yes = match p.align_soft_clip_at_reference_ends_yes {
        value => value,
    };

    if p.pe_overlap_nbases_min > 0 {
        /* peOverlap.yes is represented by pe_overlap_nbases_min > 0 */
    }

    if p.chunk_in_size_bytes_array > 0 {
        p.chunk_in_size_bytes_array /= p.read_nends.max(1) as u64;
        p.chunk_in_size_bytes = p.chunk_in_size_bytes_array.saturating_sub(
            2 * (DEF_READ_SEQ_LENGTH_MAX as u64 + 1) + 2 * DEF_READ_NAME_LENGTH_MAX as u64,
        );
    }

    p.out_sj = match out_sj_type.first().map(|s| s.as_str()) {
        Some("None") => false,
        Some("Standard") => true,
        Some(other) => {
            return Err(format!(
                "EXITING because of FATAL input ERROR: unrecognized option in --outSJtype   {}\nSOLUTION: use one of the allowed options: --outSJtype   Standard    OR    None\n",
                other
            ));
        }
        None => false,
    };
    if out_filter_type == "Normal" {
        p.out_filter_by_sjout_stage = 0;
    } else if out_filter_type == "BySJout" {
        if !p.out_sj {
            return Err("EXITING because of FATAL input ERROR: --outFilterType BySJout requires --outSJtype Standard\nSOLUTION: --outFilterType Normal    OR   --outFilterType BySJout --outSJtype Standard\n".to_string());
        }
        if p.out_sam_order == "PairedKeepInputOrder" {
            return Err("EXITING: fatal input ERROR: --outFilterType=BySJout is not presently compatible with --outSAMorder=PairedKeepInputOrder\nSOLUTION: re-run STAR without --outSAMorder=PairedKeepInputOrder, or without --outFilterType=BySJout\n".to_string());
        }
        p.out_filter_by_sjout_stage = 1;
    } else if !out_filter_type.is_empty() {
        return Err(format!(
            "EXITING because of FATAL input ERROR: unknown value of parameter outFilterType: {}\nSOLUTION: specify one of the allowed values: Normal | BySJout\n",
            out_filter_type
        ));
    }

    state
        .log_main
        .push_str("Finished loading and checking parameters\n");

    Ok(state)
}

#[doc = "Original `Parameters::scanAllLines` at STAR/source/Parameters.cpp:1197. Args: streamIn: istream, inputLevel: int, inputLevelRequested: int"]
pub fn parameters_l1197_parameters_scanalllines(
    stream_in: &str,
    input_level: i32,
    input_level_requested: i32,
    state: &mut crate::parameter_info::ParametersScanState,
) -> Result<(), String> {
    for line_in in stream_in.lines() {
        let mut line_in = line_in.to_string();
        parameters_l1205_parameters_scanoneline(
            &mut line_in,
            input_level,
            input_level_requested,
            state,
        )?;
    }
    Ok(())
}

#[doc = "Original `Parameters::scanOneLine` at STAR/source/Parameters.cpp:1205. Args: lineIn: string, inputLevel: int, inputLevelRequested: int"]
pub fn parameters_l1205_parameters_scanoneline(
    line_in: &mut str,
    input_level: i32,
    input_level_requested: i32,
    state: &mut crate::parameter_info::ParametersScanState,
) -> Result<i32, String> {
    const PAR_NAME_PRINT_WIDTH: usize = 30;

    if line_in.is_empty() {
        return Ok(0);
    }

    if input_level == 0 && (line_in.starts_with(' ') || line_in.starts_with('\t')) {
        return Ok(0);
    }

    let mut words = line_in.split_whitespace();
    let par_in = words.next().unwrap_or("");
    if par_in.is_empty() || par_in.starts_with("//") || par_in.starts_with('#') {
        return Ok(0);
    }

    let mut i_par = state.par_array.len();
    for (ii, par) in state.par_array.iter().enumerate() {
        if par_in == par.name_string {
            if input_level_requested < 0 || input_level_requested == par.input_level_allowed {
                i_par = ii;
                break;
            } else {
                return Ok(1);
            }
        }
    }

    let par_v = words.next().unwrap_or("");
    if par_v.is_empty() {
        let input_name = state
            .parameter_input_name
            .get(input_level as usize)
            .map(String::as_str)
            .unwrap_or("");
        return Err(format!(
            "EXITING: FATAL INPUT ERROR: empty value for parameter \"{}\" in input \"{}\"\nSOLUTION: use non-empty value for this parameter\n",
            par_in, input_name
        ));
    }

    let value_start = line_in
        .find(par_in)
        .map(|pos| pos + par_in.len())
        .unwrap_or(0);
    let par_v_all = line_in[value_start..].to_string();

    if i_par == state.par_array.len() {
        let input_name = state
            .parameter_input_name
            .get(input_level as usize)
            .map(String::as_str)
            .unwrap_or("");
        return Err(format!(
            "EXITING: FATAL INPUT ERROR: unrecognized parameter name \"{}\" in input \"{}\"\nSOLUTION: use correct parameter name (check the manual)\n",
            par_in, input_name
        ));
    }

    if input_level == 0 && state.par_array[i_par].input_level > 0 {
        state.log_main.push_str(&format!(
            "{:<width$}{} ... is RE-DEFINED on Command Line as: {}\n",
            state.par_array[i_par].name_string,
            par_v_all,
            state.par_array[i_par].value_line,
            width = PAR_NAME_PRINT_WIDTH
        ));
    } else if state.par_array[i_par].input_level_allowed > 0
        && state.par_array[i_par].input_level_allowed < input_level
    {
        let input_name = state
            .parameter_input_name
            .get(input_level as usize)
            .map(String::as_str)
            .unwrap_or("");
        let allowed_name = state
            .parameter_input_name
            .get(state.par_array[i_par].input_level_allowed as usize)
            .map(String::as_str)
            .unwrap_or("");
        return Err(format!(
            "EXITING: FATAL INPUT ERROR: parameter \"{}\" cannot be defined at the input level \"{}\"\nSOLUTION: define parameter \"{}\" in \"{}\"\n",
            par_in, input_name, par_in, allowed_name
        ));
    } else if state.par_array[i_par].input_level == input_level {
        let input_name = state
            .parameter_input_name
            .get(input_level as usize)
            .map(String::as_str)
            .unwrap_or("");
        return Err(format!(
            "EXITING: FATAL INPUT ERROR: duplicate parameter \"{}\" in input \"{}\"\nSOLUTION: keep only one definition of input parameters in each input source\n",
            par_in, input_name
        ));
    } else {
        state.par_array[i_par].value_line = par_v_all.trim_start().to_string();
        state.par_array[i_par].input_level = input_level;
        state.log_main.push_str(&format!(
            "{:<width$}{}",
            state.par_array[i_par].name_string,
            state.par_array[i_par].value_line,
            width = PAR_NAME_PRINT_WIDTH
        ));
        if state.par_array[i_par].input_level > 0 {
            state.log_main.push_str("     ~RE-DEFINED");
        }
        state.log_main.push('\n');
    }

    Ok(1)
}
