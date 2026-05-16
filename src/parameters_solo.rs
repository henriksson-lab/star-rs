#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `MultiMappers` at STAR/source/ParametersSolo.h:45."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MultiMappers {
    pub types_in: Vec<String>,
    pub types: Vec<i32>,
    pub type_main: i32,
    pub yes_multi: bool,
    pub yes_n: u32,
    pub yes_b: [bool; 5],
    pub count_ind_i: [u32; 5],
    pub count_ind_main: u32,
}

#[doc = "Original class `ParametersSolo` at STAR/source/ParametersSolo.h:14."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParametersSolo {
    pub yes: bool,
    pub type_str: String,
    pub solo_type: i32,
    pub barcode_read: u32,
    pub barcode_read_separate: bool,
    pub barcode_start: u32,
    pub barcode_end: u32,
    pub sam_attr_yes: bool,
    pub read_info_yes: Vec<bool>,
    pub read_index_yes: Vec<bool>,
    pub read_stats_yes: Vec<bool>,
    pub read_stats_type: String,
    pub cb_wl_yes: bool,
    pub cb_wl_size: u32,
    pub cb_wl: Vec<u64>,
    pub cb_wl_str: Vec<String>,
    pub solo_cb_whitelist: Vec<String>,
    pub cb_position_str: Vec<String>,
    pub umi_position_str: String,
    pub cb_v: Vec<SoloBarcode>,
    pub umi_v: SoloBarcode,
    pub adapter_yes: bool,
    pub adapter_seq: String,
    pub adapter_mismatches_nmax: u32,
    pub cb_l: u32,
    pub cb_s: u32,
    pub umi_s: u32,
    pub b_l: u32,
    pub cbumi_l: u32,
    pub cb_type_type: i32,
    pub cb_type_str_map: std::collections::BTreeMap<String, u32>,
    pub sam_attr_barcode_seq: Vec<String>,
    pub sam_attr_barcode_qual: Vec<String>,
    pub features: Vec<u32>,
    pub feature_yes: Vec<bool>,
    pub feature_ind: Vec<i32>,
    pub feature_first: i32,
    pub n_features: u32,
    pub out_file_names: Vec<String>,
    pub out_format_features_gene_field3: String,
    pub umi_l: u32,
    pub umi_mask_low: u32,
    pub umi_dedup: UMIdedup,
    pub umi_filtering: SoloUmiFiltering,
    pub multi_map: MultiMappers,
    pub cb_match_wl: CBMatchWL,
    pub sam_attr_feature: i32,
    pub qs_base: i8,
    pub qs_max: i8,
    pub cb_min_p: f64,
    pub redistr_reads_nfiles: u32,
    pub cluster_cb_file: String,
    pub cell_filter: SoloCellFilter,
    pub strand: i32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloCellFilter {
    pub type_: Vec<String>,
    pub top_cells: u32,
    pub knee: SoloCellFilterKnee,
    pub ed_cr: SoloCellFilterEmptyDropsCr,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloCellFilterKnee {
    pub n_expected_cells: f64,
    pub max_percentile: f64,
    pub max_min_ratio: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloCellFilterEmptyDropsCr {
    pub ind_min: u32,
    pub ind_max: u32,
    pub umi_min: u32,
    pub umi_min_frac_median: f64,
    pub cand_max_n: u32,
    pub fdr: f64,
    pub sim_n: u32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CBMatchWL {
    pub type_: String,
    pub mm1: bool,
    pub mm1_multi: bool,
    pub one_exact: bool,
    pub mm1_multi_pc: bool,
    pub mm1_multi_nbase: bool,
    pub edit_dist_2: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloUmiFiltering {
    pub type_: Vec<String>,
    pub yes: bool,
    pub multi_gene_umi: bool,
    pub multi_gene_umi_all: bool,
    pub multi_gene_umi_cr: bool,
}

#[doc = "Original class `UMIdedup` at STAR/source/ParametersSolo.h:16."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct UMIdedup {
    pub types_in: Vec<String>,
    pub types: Vec<i32>,
    pub type_main: i32,
    pub yes_n: u32,
    pub yes_b: [bool; 6],
    pub count_ind_i: [u32; 6],
    pub count_ind_main: u32,
}

#[doc = "Original `ParametersSolo::initialize` at STAR/source/ParametersSolo.cpp:10. Args: pPin: Parameters"]
pub fn parameterssolo_l10_parameterssolo_initialize(
    p: &mut crate::parameters_chimeric::Parameters,
    whitelist_contents: &[(String, String)],
    raw_time: libc::time_t,
) -> Result<String, String> {
    let mut log_main = String::new();
    let mut p_solo = std::mem::take(&mut p.p_solo);

    log_main.push_str(&parameterssolo_l533_parameterssolo_cellfiltering(
        &mut p_solo,
    )?);
    if p_solo.cell_filter.type_.is_empty() {
        p_solo.cell_filter.type_.push("None".to_string());
    }

    if p.out_sam_mode == "soloCellFiltering" {
        p_solo.yes = true;
        p_solo.umi_dedup.types_in = vec!["NoDedup".to_string()];
        let p_solo_for_dedup = p_solo.clone();
        parameterssolo_l585_umidedup_initialize(&mut p_solo.umi_dedup, &p_solo_for_dedup)?;
        p.p_solo = p_solo;
        return Ok(log_main);
    }

    p_solo.redistr_reads_nfiles = 3_u32.saturating_mul(p.run_thread_n.max(0) as u32);
    p_solo.barcode_start = 0;
    p_solo.barcode_end = 0;
    p_solo.yes = true;

    let solo_feature_n = SOLO_FEATURE_VELOCYTO as usize + 1;
    if p_solo.feature_yes.len() < solo_feature_n {
        p_solo.feature_yes.resize(solo_feature_n, false);
    }
    if p_solo.feature_ind.len() < solo_feature_n {
        p_solo.feature_ind.resize(solo_feature_n, -1);
    }
    if p_solo.read_info_yes.len() < solo_feature_n {
        p_solo.read_info_yes.resize(solo_feature_n, false);
    }
    if p_solo.read_index_yes.len() < solo_feature_n {
        p_solo.read_index_yes.resize(solo_feature_n, false);
    }
    if p_solo.read_stats_yes.len() < solo_feature_n {
        p_solo.read_stats_yes.resize(solo_feature_n, false);
    }

    if p_solo.type_str == "None" || p_solo.type_str == "SmartSeq" {
        let solo_tags_requested = p.out_sam_attr_present.cr
            || p.out_sam_attr_present.cy
            || p.out_sam_attr_present.ur
            || p.out_sam_attr_present.uy
            || p.out_sam_attr_present.cb
            || p.out_sam_attr_present.ub
            || p.out_sam_attr_present.s_s
            || p.out_sam_attr_present.s_q
            || p.out_sam_attr_present.s_m
            || p.out_sam_attr_present.s_f;
        if solo_tags_requested {
            return Err(format!(
                "EXITING because of FATAL INPUT ERROR: --outSAMattributes contains CR/CY/UR/UY/CB/UB tags which are not allowed for --soloType {}\nSOLUTION: re-run STAR without these attribures\n",
                p_solo.type_str
            ));
        }
    }

    if p_solo.type_str == "None" {
        p_solo.solo_type = SOLO_TYPE_NONE;
        p_solo.yes = false;
        p_solo.sam_attr_yes = false;
        p.p_solo = p_solo;
        return Ok(log_main);
    } else if p_solo.type_str == "CB_UMI_Simple" || p_solo.type_str == "Droplet" {
        p_solo.solo_type = SOLO_TYPE_CB_UMI_SIMPLE;
        if p_solo.umi_l > 16 {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: UMI length is too long: --soloUMIlen={}\nSOLUTION: UMI length cannot be longer than 16",
                p_solo.umi_l
            ));
        }
        if p_solo.cb_l > 31 {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: CB length is too long: --soloCBlen={}\nSOLUTION: CB length cannot be longer than 31",
                p_solo.cb_l
            ));
        }
        p_solo.cbumi_l = p_solo.cb_l + p_solo.umi_l;
        if p_solo.b_l == 1 {
            p_solo.b_l = p_solo.cbumi_l;
        }
        p_solo.barcode_start = std::cmp::min(p_solo.cb_s, p_solo.umi_s).saturating_sub(1);
        p_solo.barcode_end =
            std::cmp::max(p_solo.cb_s + p_solo.cb_l, p_solo.umi_s + p_solo.umi_l).saturating_sub(2);
    } else if p_solo.type_str == "CB_UMI_Complex" {
        p_solo.solo_type = SOLO_TYPE_CB_UMI_COMPLEX;
        p_solo.b_l = 0;
        p_solo.cbumi_l = 0;
    } else if p_solo.type_str == "CB_samTagOut" {
        p_solo.solo_type = SOLO_TYPE_CB_SAM_TAG_OUT;
        p_solo.cbumi_l = p_solo.cb_l + p_solo.umi_l;
        if p_solo.b_l == 1 {
            p_solo.b_l = p_solo.cbumi_l;
        }
        p_solo.barcode_start = std::cmp::min(p_solo.cb_s, p_solo.umi_s).saturating_sub(1);
        p_solo.barcode_end =
            std::cmp::max(p_solo.cb_s + p_solo.cb_l, p_solo.umi_s + p_solo.umi_l).saturating_sub(2);
    } else if p_solo.type_str == "SmartSeq" {
        p_solo.solo_type = SOLO_TYPE_SMART_SEQ;
    } else {
        return Err(format!(
            "EXITING because of fatal PARAMETERS error: unrecognized option in --soloType={}\nSOLUTION: use allowed option: None OR CB_UMI_Simple OR CB_UMI_Complex\nObsolete option Droplet should be replaced with CB_UMI_Simple",
            p_solo.type_str
        ));
    }

    if p_solo.cb_type_type == 0 {
        p_solo.cb_type_type = 1;
    }

    let barcode_read_in = p_solo.barcode_read;
    p_solo.barcode_read = u32::MAX;
    p_solo.barcode_read_separate = false;
    if p.read_files_type_n != 10 {
        if p_solo.solo_type != SOLO_TYPE_SMART_SEQ {
            if barcode_read_in == 0 || barcode_read_in == u32::MAX {
                if p.read_nends < 2 {
                    return Err("EXITING because of fatal PARAMETERS error: --soloType (except SmartSeq) with --soloBarcodeMate 0 (default) require 2 reads or 3 reads, where the last read is the barcode read.\nSOLUTION: if barcode is in a separate mate, specify it as the last file in --readFilesIn. If barcode sequence is a part of one of the mates, specify that mate with --soloBarcodeMate 1 (or 2 or 3)".to_string());
                }
                p.read_nmates = p.read_nends - 1;
                p_solo.barcode_read = p.read_nends - 1;
            } else if barcode_read_in > p.read_nends {
                return Err(format!(
                    "EXITING because of fatal PARAMETERS error: --soloBarcodeMate {}is larger than number of mates {}\nSOLUTION: specify --soloBarcodeMate <= than the number of mates.",
                    barcode_read_in, p.read_nends
                ));
            } else {
                if p_solo.solo_type != SOLO_TYPE_CB_UMI_SIMPLE {
                    return Err(format!(
                        "EXITING because of fatal PARAMETERS error: --soloBarcodeMate {}>0 for is not allowed for --soloType {}\nSOLUTION: specify --soloBarcodeMate 0   or   --soloType CB_UMI_Simple",
                        barcode_read_in, p_solo.type_str
                    ));
                }
                p_solo.barcode_read = barcode_read_in - 1;
                p_solo.barcode_read_separate = true;
                p_solo.b_l = 0;
                let br = p_solo.barcode_read as usize;
                let clip_5 = p.p_clip.in_[0].n.get(br).copied().unwrap_or(0);
                let clip_3 = p.p_clip.in_[1].n.get(br).copied().unwrap_or(0);
                if clip_5 == 0 && clip_3 == 0 {
                    return Err(format!(
                        "EXITING because of fatal PARAMETERS error: --soloBarcodeMate {} specifies that barcode sequence is a part of the mate {}, which requires clipping the barcode off this mate.\nSOLUTION: clip the barcode sequence from 5' or/and 3' with --clip5pNbases   or/and --clip3pNbases . The values for mate1 and mate2 have to specified, specify 0 for no clipping.",
                        br + 1,
                        br + 1
                    ));
                }
            }
        }
    } else {
        if p_solo.type_str == "SmartSeq" {
            return Err("EXITING because of fatal PARAMETERS error: --readFilesType SAM SE/PE cannot be used with --soloType SmartSeq\nSOLUTION: for Smart-seq input from BAM files, use --soloType CB_UMI_Simple , create whitelist of SmartSeq file names, and specify the SAM tag that records these file names in --soloInputSAMattrBarcodeSeq".to_string());
        }
        if p_solo.sam_attr_barcode_seq.first().map(String::as_str) == Some("-") {
            return Err("EXITING because of fatal PARAMETERS error: --readFilesType SAM SE/PE requires --soloInputSAMattrBarcodeSeq.\nSOLUTION: specify input SAM attributes for barcode sequence in --soloInputSAMattrBarcodeSeq, and (optionally) quality with --soloInputSAMattrBarcodeQual".to_string());
        }
        if p_solo.sam_attr_barcode_qual.first().map(String::as_str) == Some("-") {
            log_main.push_str("WARNING: since --readFilesType SAM SE/PE --soloInputSAMattrBarcodeQual - : qualities for barcode read will be replaced with 'H'\n");
            p_solo.sam_attr_barcode_qual.clear();
        }
        for tag in p_solo.sam_attr_barcode_seq.iter_mut() {
            if tag.len() != 2 {
                return Err("EXITING because of fatal PARAMETERS error: --soloInputSAMattrBarcodeSeq attributes have to be two-letter strings.\nSOLUTION: specify correct two-letter strings in --soloInputSAMattrBarcodeSeq".to_string());
            }
            tag.insert(0, '\t');
        }
        for tag in p_solo.sam_attr_barcode_qual.iter_mut() {
            if tag.len() != 2 {
                return Err("EXITING because of fatal PARAMETERS error: --soloInputSAMattrBarcodeQual attributes have to be two-letter strings.\nSOLUTION: specify correct two-letter strings in --soloInputSAMattrBarcodeQual".to_string());
            }
            tag.insert(0, '\t');
        }
    }

    p_solo.strand = match p_solo.strand {
        -1 | 0 | 1 => p_solo.strand,
        _ => -1,
    };

    p_solo.feature_ind.fill(-1);
    p_solo.feature_yes.fill(false);
    p_solo.feature_first = -1;
    let feature_names = [
        "SJ",
        "Transcript3p",
        "GeneFull",
        "GeneFull_ExonOverIntron",
        "GeneFull_Ex50pAS",
        "Gene",
        "VelocytoSimple",
        "Velocyto",
    ];
    let feature_in: Vec<String> = if p_solo.features.is_empty() {
        vec!["Gene".to_string()]
    } else {
        p_solo
            .features
            .iter()
            .map(|f| feature_names.get(*f as usize).unwrap_or(&"").to_string())
            .collect()
    };
    p_solo.features.clear();
    for fin in feature_in.iter() {
        let mut found = false;
        for (ii, name) in feature_names.iter().enumerate() {
            if fin == name {
                p_solo.feature_yes[ii] = true;
                p_solo.features.push(ii as u32);
                if p_solo.feature_first == -1 {
                    p_solo.feature_first = ii as i32;
                }
                found = true;
                break;
            }
        }
        if !found {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: unrecognized option(s) in --soloFeatures {}\nSOLUTION: use allowed option: Gene",
                feature_in.join(" ")
            ));
        }
    }
    p_solo.n_features = p_solo.features.len() as u32;
    p_solo.features.sort_unstable();
    for (ii, feature) in p_solo.features.iter().enumerate() {
        p_solo.feature_ind[*feature as usize] = ii as i32;
    }

    if p_solo.feature_yes[SOLO_FEATURE_VELOCYTO as usize] && p_solo.solo_type == SOLO_TYPE_SMART_SEQ
    {
        return Err("EXITING because of fatal PARAMETERS error: --soloFeatures Velocyto is presently not compatible with --soloType SmartSeq .\nSOLUTION: re-run without --soloFeatures Velocyto .".to_string());
    }

    if p_solo.feature_yes[SOLO_FEATURE_GENE as usize] {
        p.quant_gene_yes = true;
        p.quant_yes = true;
    }
    if p_solo.feature_yes[SOLO_FEATURE_GENE_FULL as usize] {
        p.quant_gene_full_yes = true;
        p.quant_yes = true;
        if !p_solo.feature_yes[SOLO_FEATURE_GENE as usize] {
            p.quant_gene_yes = false;
        }
    }
    if p_solo.feature_yes[SOLO_FEATURE_GENE_FULL_EX50P_AS as usize] {
        p.quant_gene_full_ex50p_as_yes = true;
        p.quant_yes = true;
        if !p_solo.feature_yes[SOLO_FEATURE_GENE as usize] {
            p.quant_gene_yes = false;
        }
    }
    if p_solo.feature_yes[SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON as usize] {
        p.quant_gene_full_exon_over_intron_yes = true;
        p.quant_gene_yes = true;
        p.quant_yes = true;
    }

    parameterssolo_l675_parameterssolo_init_cbmatchwl(&mut p_solo)?;
    let p_solo_for_dedup = p_solo.clone();
    parameterssolo_l585_umidedup_initialize(&mut p_solo.umi_dedup, &p_solo_for_dedup)?;

    p_solo.qs_base = 33;
    p_solo.qs_max = 33;
    p_solo.cb_min_p = 0.975;
    p_solo.umi_mask_low = if p_solo.umi_l >= 32 {
        u32::MAX
    } else {
        ((1_u64 << p_solo.umi_l) - 1) as u32
    };

    if p_solo.solo_type == SOLO_TYPE_CB_UMI_SIMPLE || p_solo.solo_type == SOLO_TYPE_CB_SAM_TAG_OUT {
        if p_solo.solo_cb_whitelist.len() > 1 {
            return Err("EXITING because of FATAL ERROR in INPUT parameters: --soloCBwhitelist contains more than one file which is not allowed with --soloType CB_UMI_Simple \nSOLUTION: in --soloCBwhitelist specify only one whitelist file \n".to_string());
        } else if p_solo.solo_cb_whitelist.first().map(String::as_str) == Some("-") {
            return Err("EXITING because of FATAL ERROR in INPUT parameters: --soloCBwhitelist is not defined\nSOLUTION: in --soloCBwhitelist specify path to and name of the whitelist file, or None for CB demultiplexing without whitelist \n".to_string());
        } else if p_solo.solo_cb_whitelist.first().map(String::as_str) == Some("None") {
            p_solo.cb_wl_yes = false;
        } else if let Some(wl_name) = p_solo.solo_cb_whitelist.first() {
            p_solo.cb_wl_yes = true;
            let wl_contents = whitelist_contents
                .iter()
                .find(|(name, _)| name == wl_name)
                .map(|(_, contents)| contents.as_str())
                .ok_or_else(|| {
                    format!(
                        "EXITING because of FATAL ERROR: could not open CB whitelist file {}",
                        wl_name
                    )
                })?;
            p_solo.cb_wl.clear();
            for seq1 in wl_contents.split_whitespace() {
                if seq1.len() != p_solo.cb_l as usize {
                    return Err(format!(
                        "EXITING because of FATAL ERROR in input CB whitelist file: {} the total length of barcode sequence is {} not equal to expected {}\nSOLUTION: make sure that the barcode read is the second in --readFilesIn and check that is has the correct formatting\n",
                        wl_name,
                        seq1.len(),
                        p_solo.b_l
                    ));
                }
                let mut cb1 = 0_u64;
                if sequencefuns_l249_convertnuclstrtoint64(seq1, &mut cb1) == -1 {
                    p_solo.cb_wl.push(cb1);
                } else {
                    log_main.push_str(&format!(
                        "WARNING: CB whitelist sequence contains non-ACGT base and is ignored: {}\n",
                        seq1
                    ));
                }
            }
            if p_solo.cb_wl.is_empty() {
                return Err(format!(
                    "EXITING because of FATAL ERROR: CB whitelist file {} is empty. \nSOLUTION: provide non-empty whitelist.\n",
                    wl_name
                ));
            }
        }
        p_solo.cb_wl.sort_unstable();
        p_solo.cb_wl.dedup();
        p_solo.cb_wl_size = p_solo.cb_wl.len() as u32;
        log_main.push_str(&format!(
            "Number of CBs in the whitelist = {}\n",
            p_solo.cb_wl_size
        ));
        p_solo
            .cb_wl_str
            .resize(p_solo.cb_wl_size as usize, String::new());
        for ii in 0..p_solo.cb_wl_size as usize {
            p_solo.cb_wl_str[ii] =
                sequencefuns_l267_convertnuclint64tostring(p_solo.cb_wl[ii], p_solo.cb_l);
        }
    } else if p_solo.solo_type == SOLO_TYPE_SMART_SEQ {
        p_solo.cb_wl_str = p.out_sam_attr_rg.clone();
        p_solo.cb_wl_size = p_solo.cb_wl_str.len() as u32;
        p_solo.cb_wl_yes = true;
    } else if p_solo.solo_type == SOLO_TYPE_CB_UMI_COMPLEX {
        p_solo.cb_wl_yes = true;
        p_solo.adapter_yes = p_solo.adapter_seq != "-";
        if p_solo.cb_position_str.len() != p_solo.solo_cb_whitelist.len() {
            return Err(format!(
                "EXITING because of fatal PARAMETER error: number of barcodes in --soloCBposition : {} is not equal to the number of WhiteLists in --soloCBwhitelist : {}\nSOLUTION: make sure that the number of CB whitelists and CB positions are the same\n",
                p_solo.cb_position_str.len(),
                p_solo.solo_cb_whitelist.len()
            ));
        }
        p_solo
            .cb_v
            .resize(p_solo.cb_position_str.len(), Default::default());
        p_solo.cb_wl_size = 1;
        for icb in 0..p_solo.cb_v.len() {
            p_solo.cb_v[icb].adapter_length = p_solo.adapter_seq.len() as i32;
            let wl_name = &p_solo.solo_cb_whitelist[icb];
            let wl_contents = whitelist_contents
                .iter()
                .find(|(name, _)| name == wl_name)
                .map(|(_, contents)| contents.as_str())
                .ok_or_else(|| {
                    format!(
                        "EXITING because of FATAL ERROR: could not open CB whitelist file {}",
                        wl_name
                    )
                })?;
            for seq1 in wl_contents.split_whitespace() {
                let mut cb1 = 0_u64;
                if sequencefuns_l249_convertnuclstrtoint64(seq1, &mut cb1) != -1 {
                    log_main.push_str(&format!(
                        "WARNING: CB whitelist sequence contains non-ACGT base and is ignored: {}\n",
                        seq1
                    ));
                    continue;
                }
                let len1 = seq1.len();
                if len1 >= p_solo.cb_v[icb].wl.len() {
                    p_solo.cb_v[icb].wl.resize(len1 + 1, Vec::new());
                }
                p_solo.cb_v[icb].wl[len1].push(cb1);
            }
            solobarcode_l9_solobarcode_sortwhitelist(
                &mut p_solo.cb_v[icb],
                p_solo.cb_match_wl.edit_dist_2,
            );
            p_solo.cb_v[icb].wl_factor = p_solo.cb_wl_size as u64;
            p_solo.cb_wl_size = p_solo
                .cb_wl_size
                .saturating_mul(p_solo.cb_v[icb].total_size);
        }
        parameterssolo_l503_parameterssolo_complexwlstrings(&mut p_solo);
    }

    log_main.push_str(&format!(
        "{} ... Finished reading, sorting and deduplicating CB whitelist sequences.\n",
        timefunctions_l14_timemonthdaytime(raw_time)
    ));

    p_solo.sam_attr_yes = false;
    if (p.out_sam_attr_present.cb || p.out_sam_attr_present.ub)
        && p_solo.solo_type != SOLO_TYPE_CB_SAM_TAG_OUT
    {
        p_solo.sam_attr_yes = true;
        if !p.out_bam_coord {
            return Err("EXITING because of fatal PARAMETERS error: CB and/or UB attributes in --outSAMattributes can only be output in the sorted BAM file.\nSOLUTION: re-run STAR with --outSAMtype BAM SortedByCoordinate ...\n".to_string());
        }
    } else if p.out_sam_attr_present.ub && p_solo.solo_type == SOLO_TYPE_CB_SAM_TAG_OUT {
        return Err("EXITING because of fatal PARAMETERS error: UB attribute (corrected UMI) in --outSAMattributes cannot be used with --soloType CB_samTagOut \nSOLUTION: instead, use UR (uncorrected UMI) in --outSAMattributes\n".to_string());
    }

    p_solo.read_info_yes.fill(false);
    if p_solo.feature_yes[SOLO_FEATURE_VELOCYTO_SIMPLE as usize]
        || p_solo.feature_yes[SOLO_FEATURE_VELOCYTO as usize]
    {
        p_solo.read_info_yes[SOLO_FEATURE_GENE as usize] = true;
    }
    p_solo.sam_attr_feature = p_solo.feature_first;
    if p_solo.sam_attr_yes {
        match p_solo.feature_first {
            SOLO_FEATURE_GENE
            | SOLO_FEATURE_GENE_FULL
            | SOLO_FEATURE_GENE_FULL_EX50P_AS
            | SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON => {}
            _ => {
                return Err("EXITING because of fatal PARAMETERS error: CB and/or UB attributes in --outSAMattributes require --soloFeatures Gene OR/AND GeneFull OR/AND GeneFull_Ex50pAS.\nSOLUTION: re-run STAR adding Gene AND/OR GeneFull OR/AND GeneFull_Ex50pAS OR/AND GeneFull_ExonOverIntron to --soloFeatures\n".to_string());
            }
        }
        p_solo.read_info_yes[p_solo.sam_attr_feature as usize] = true;
    }
    p_solo.read_index_yes = p_solo.read_info_yes.clone();

    p_solo.read_stats_yes.fill(false);
    if p_solo.read_stats_type == "Standard" {
        p_solo.read_stats_yes.fill(true);
        p_solo.read_stats_yes[SOLO_FEATURE_VELOCYTO_SIMPLE as usize] = false;
        p_solo.read_stats_yes[SOLO_FEATURE_VELOCYTO as usize] = false;
        p_solo.read_stats_yes[SOLO_FEATURE_SJ as usize] = false;
        for ff in 0..p_solo.read_index_yes.len() {
            p_solo.read_index_yes[ff] |= p_solo.read_stats_yes[ff];
        }
    } else if !p_solo.read_stats_type.is_empty() && p_solo.read_stats_type != "None" {
        return Err(format!(
            "EXITING because of fatal PARAMETERS error: unrecognized option in --soloCellReadStats{}\nSOLUTION: use allowed options: None OR Standard \n",
            p_solo.read_stats_type
        ));
    }

    if p_solo.umi_filtering.type_.is_empty() {
        p_solo.umi_filtering.type_.push("-".to_string());
    }
    if p_solo.umi_filtering.type_[0] == "MultiGeneUMI" {
        p_solo.umi_filtering.multi_gene_umi = true;
        p_solo.umi_filtering.yes = true;
    } else if p_solo.umi_filtering.type_[0] == "MultiGeneUMI_All" {
        p_solo.umi_filtering.multi_gene_umi_all = true;
        p_solo.umi_filtering.yes = true;
    } else if p_solo.umi_filtering.type_[0] == "MultiGeneUMI_CR" {
        p_solo.umi_filtering.multi_gene_umi_cr = true;
        if p_solo.umi_dedup.types_in.len() > 1
            || p_solo.umi_dedup.types_in.first().map(String::as_str) != Some("1MM_CR")
        {
            return Err("EXITING because of fatal PARAMETERS error: --soloUMIfiltering MultiGeneUMI_CR only works with --soloUMIdedup 1MM_CR\nSOLUTION: rerun with --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \n".to_string());
        }
    } else if p_solo.umi_filtering.type_[0] != "-" {
        return Err(format!(
            "EXITING because of fatal PARAMETERS error: unrecognized option in --soloUMIfiltering={}\nSOLUTION: use allowed options: - or MultiGeneUMI or MultiGeneUMI_CR \n",
            p_solo.umi_filtering.type_[0]
        ));
    }

    let p_solo_for_multi = p_solo.clone();
    parameterssolo_l624_multimappers_initialize(&mut p_solo.multi_map, &p_solo_for_multi)?;
    if p_solo.multi_map.yes_multi {
        if p_solo.solo_type == SOLO_TYPE_CB_SAM_TAG_OUT || p_solo.solo_type == SOLO_TYPE_SMART_SEQ {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: multimapping options do not work for --soloType {}\nSOLUTION: use default option --soloMultiMappers Unique\n",
                p_solo.type_str
            ));
        }
        p_solo.read_index_yes[SOLO_FEATURE_GENE as usize] = true;
        p_solo.read_index_yes[SOLO_FEATURE_GENE_FULL as usize] = true;
        p_solo.read_index_yes[SOLO_FEATURE_GENE_FULL_EX50P_AS as usize] = true;
        p_solo.read_index_yes[SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON as usize] = true;
    }

    p.p_solo = p_solo;
    Ok(log_main)
}

#[doc = "Original `ParametersSolo::umiSwapHalves` at STAR/source/ParametersSolo.cpp:496. Args: umi: uint32"]
pub fn parameterssolo_l496_parameterssolo_umiswaphalves(
    p_solo: &crate::parameters_solo::ParametersSolo,
    umi: &mut u32,
) {
    let high = *umi >> p_solo.umi_l;
    *umi &= p_solo.umi_mask_low;
    *umi <<= p_solo.umi_l;
    *umi |= high;
}

#[doc = "Original `ParametersSolo::complexWLstrings` at STAR/source/ParametersSolo.cpp:503. Args: "]
pub fn parameterssolo_l503_parameterssolo_complexwlstrings(
    p_solo: &mut crate::parameters_solo::ParametersSolo,
) {
    p_solo
        .cb_wl_str
        .resize(p_solo.cb_wl_size as usize, String::new());

    for cb in p_solo.cb_v.iter_mut() {
        cb.i_cb = 0;
        cb.i_len = cb.min_len;
    }

    for ii in 0..p_solo.cb_wl_size as usize {
        for ii in 0..p_solo.cb_v.len() {
            if p_solo.cb_v[ii].i_cb
                == p_solo.cb_v[ii].wl[p_solo.cb_v[ii].i_len as usize].len() as u32
            {
                p_solo.cb_v[ii].i_len += 1;
                p_solo.cb_v[ii].i_cb = 0;
            }
            if p_solo.cb_v[ii].i_len == p_solo.cb_v[ii].wl.len() as u32 {
                p_solo.cb_v[ii + 1].i_cb += 1;
                p_solo.cb_v[ii].i_len = p_solo.cb_v[ii].min_len;
            }
        }

        for cb in p_solo.cb_v.iter() {
            p_solo.cb_wl_str[ii].push_str(&sequencefuns_l267_convertnuclint64tostring(
                cb.wl[cb.i_len as usize][cb.i_cb as usize],
                cb.i_len,
            ));
            p_solo.cb_wl_str[ii].push('_');
        }
        p_solo.cb_wl_str[ii].pop();

        p_solo.cb_v[0].i_cb += 1;
    }
}

#[doc = "Original `ParametersSolo::cellFiltering` at STAR/source/ParametersSolo.cpp:533. Args: "]
pub fn parameterssolo_l533_parameterssolo_cellfiltering(
    p_solo: &mut crate::parameters_solo::ParametersSolo,
) -> Result<String, String> {
    let mut log_main = String::new();
    let mut pars1 = String::new();
    for s in p_solo.cell_filter.type_.iter().skip(1) {
        pars1.push(' ');
        pars1.push_str(s);
    }

    if p_solo.cell_filter.type_[0] == "CellRanger2.2" {
        if p_solo.cell_filter.type_.len() == 1 {
            log_main.push_str("ParametersSolo: using hardcoded filtering parameters for --soloCellFilterType CellRanger2.2\n");
            pars1 = "3000 0.99 10".to_string();
        } else if p_solo.cell_filter.type_.len() < 4 {
            return Err("EXITING because of fatal PARAMETERS error: --soloCellFilterType CellRanger2.2 requires exactly 3 numerical parameters\nSOLUTION: re-run with --soloCellFilterType CellRanger2.2 <nExpectedCells> <maxPercentile> <maxMinRatio>\n".to_string());
        }

        log_main
            .push_str("ParametersSolo: --soloCellFilterType CellRanger2.2 filtering parameters: ");
        log_main.push_str(&pars1);
        log_main.push('\n');
        let pars: Vec<&str> = pars1.split_whitespace().collect();
        p_solo.cell_filter.knee.n_expected_cells = pars[0].parse().unwrap_or(0.0);
        p_solo.cell_filter.knee.max_percentile = pars[1].parse().unwrap_or(0.0);
        p_solo.cell_filter.knee.max_min_ratio = pars[2].parse().unwrap_or(0.0);
    } else if p_solo.cell_filter.type_[0] == "EmptyDrops_CR" {
        if p_solo.cell_filter.type_.len() == 1 {
            log_main
                .push_str("ParametersSolo: using hardcoded filtering parameters for --soloCellFilterType EmptyDrops_CR\n");
            pars1 = "3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000".to_string();
        } else if p_solo.cell_filter.type_.len() < 11 {
            return Err("EXITING because of fatal PARAMETERS error: --soloCellFilterType EmptyDrops_CR requires exactly 10 numerical parameters\nSOLUTION: re-run with --soloCellFilterType EmptyDrops_CR <nExpectedCells> <maxPercentile> <maxMinRatio> <indMin> <indMax> <umiMin> <umiMinFracMedian> <candMaxN> <FDR> <simN>\n".to_string());
        }

        log_main
            .push_str("ParametersSolo: --soloCellFilterType EmptyDrops_CR filtering parameters: ");
        log_main.push_str(&pars1);
        log_main.push('\n');
        let pars: Vec<&str> = pars1.split_whitespace().collect();
        p_solo.cell_filter.knee.n_expected_cells = pars[0].parse().unwrap_or(0.0);
        p_solo.cell_filter.knee.max_percentile = pars[1].parse().unwrap_or(0.0);
        p_solo.cell_filter.knee.max_min_ratio = pars[2].parse().unwrap_or(0.0);
        p_solo.cell_filter.ed_cr.ind_min = pars[3].parse().unwrap_or(0);
        p_solo.cell_filter.ed_cr.ind_max = pars[4].parse().unwrap_or(0);
        p_solo.cell_filter.ed_cr.umi_min = pars[5].parse().unwrap_or(0);
        p_solo.cell_filter.ed_cr.umi_min_frac_median = pars[6].parse().unwrap_or(0.0);
        p_solo.cell_filter.ed_cr.cand_max_n = pars[7].parse().unwrap_or(0);
        p_solo.cell_filter.ed_cr.fdr = pars[8].parse().unwrap_or(0.0);
        p_solo.cell_filter.ed_cr.sim_n = pars[9].parse().unwrap_or(0);
    } else if p_solo.cell_filter.type_[0] == "TopCells" {
        if p_solo.cell_filter.type_.len() < 2 {
            return Err("EXITING because of fatal PARAMETERS error: number of cells not specified for --soloCellFilterType TopCells\nSOLUTION: --soloCellFilterType TopCells <NumberOfCells>\n".to_string());
        }
        p_solo.cell_filter.top_cells = p_solo.cell_filter.type_[1].parse().unwrap_or(0);
    } else if p_solo.cell_filter.type_[0] == "None" {
    } else {
        return Err(format!(
            "EXITING because of fatal PARAMETERS error: unrecognized option in --soloCellFilterType={}\nSOLUTION: use allowed options: CellRanger2.2 or None\n",
            p_solo.cell_filter.type_[0]
        ));
    }

    Ok(log_main)
}

#[doc = "Original `UMIdedup::initialize` at STAR/source/ParametersSolo.cpp:585. Args: pS: ParametersSolo"]
pub fn parameterssolo_l585_umidedup_initialize(
    umi_dedup: &mut crate::parameters_solo::UMIdedup,
    p_solo: &crate::parameters_solo::ParametersSolo,
) -> Result<(), String> {
    let type_names = [
        "NoDedup",
        "Exact",
        "1MM_All",
        "1MM_Directional",
        "1MM_CR",
        "1MM_Directional_UMItools",
    ];

    umi_dedup.yes_n = 0;
    umi_dedup.count_ind_i = [u32::MAX; 6];
    umi_dedup.yes_b = [false; 6];
    umi_dedup.types.clear();

    for (iin, type_in) in umi_dedup.types_in.iter().enumerate() {
        let mut itype = type_names.len();
        for (ii, type_name) in type_names.iter().enumerate() {
            if type_in == type_name {
                itype = ii;
                break;
            }
        }

        if itype == type_names.len() {
            let mut tall = String::new();
            for type_name in type_names.iter() {
                tall.push(' ');
                tall.push_str(type_name);
            }
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: unrecognzied option --soloUMIdedup = {}\nSOLUTION: use allowed values: {}\n",
                type_in, tall
            ));
        }

        umi_dedup.types.push(itype as i32);
        umi_dedup.yes_b[itype] = true;
        umi_dedup.yes_n += 1;
        umi_dedup.count_ind_i[itype] = iin as u32 + 1;

        if p_solo.solo_type == 4 && (umi_dedup.yes_b[2] || umi_dedup.yes_b[3] || umi_dedup.yes_b[4])
        {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: --soloUMIdedup = {} is not allowed for --soloType SmartSeq\nSOLUTION: use allowed options: Exact and/or NoDedup\n",
                type_in
            ));
        }
    }

    umi_dedup.type_main = umi_dedup.types[0];
    umi_dedup.count_ind_main = 1;
    Ok(())
}

#[doc = "Original `MultiMappers::initialize` at STAR/source/ParametersSolo.cpp:624. Args: pS: ParametersSolo"]
pub fn parameterssolo_l624_multimappers_initialize(
    multi_mappers: &mut crate::parameters_solo::MultiMappers,
    p_solo: &crate::parameters_solo::ParametersSolo,
) -> Result<(), String> {
    let type_names = ["Unique", "Uniform", "Rescue", "PropUnique", "EM"];

    multi_mappers.yes_n = 0;
    multi_mappers.count_ind_i = [u32::MAX; 5];
    multi_mappers.yes_b = [false; 5];
    multi_mappers.types.clear();

    for type_in in multi_mappers.types_in.iter() {
        let mut itype = type_names.len();
        for (ii, type_name) in type_names.iter().enumerate() {
            if type_in == type_name {
                itype = ii;
                break;
            }
        }

        if itype == type_names.len() {
            let mut tall = String::new();
            for type_name in type_names.iter() {
                tall.push(' ');
                tall.push_str(type_name);
            }
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: unrecognzied option --soloMultiMappers = {}\nSOLUTION: use allowed values: {}\n",
                type_in, tall
            ));
        }

        if itype == 0 {
            continue;
        }

        multi_mappers.types.push(itype as i32);
        multi_mappers.yes_b[itype] = true;
        multi_mappers.yes_n += 1;
    }

    if multi_mappers.yes_n == 0 {
        multi_mappers.yes_multi = false;
        return Ok(());
    }

    let mut ind1 = 1u32;
    for &itype in multi_mappers.types.iter() {
        multi_mappers.count_ind_i[itype as usize] = ind1;
        ind1 += p_solo.umi_dedup.yes_n;
    }

    multi_mappers.type_main = multi_mappers.types[0];
    multi_mappers.count_ind_main = 1;
    multi_mappers.yes_multi = multi_mappers.yes_b[1]
        || multi_mappers.yes_b[2]
        || multi_mappers.yes_b[3]
        || multi_mappers.yes_b[4];
    Ok(())
}

#[doc = "Original `ParametersSolo::init_CBmatchWL` at STAR/source/ParametersSolo.cpp:675. Args: "]
pub fn parameterssolo_l675_parameterssolo_init_cbmatchwl(
    p_solo: &mut crate::parameters_solo::ParametersSolo,
) -> Result<(), String> {
    let incomp1 = (p_solo.type_str == "CB_UMI_Complex"
        && p_solo.cb_match_wl.type_ != "Exact"
        && p_solo.cb_match_wl.type_ != "1MM"
        && p_solo.cb_match_wl.type_ != "EditDist_2")
        || (p_solo.type_str == "CB_samTagOut"
            && p_solo.cb_match_wl.type_ != "Exact"
            && p_solo.cb_match_wl.type_ != "1MM")
        || (p_solo.type_str != "CB_UMI_Complex" && p_solo.cb_match_wl.type_ == "EditDist_2");

    if incomp1 {
        return Err(format!(
            "EXITING because of fatal PARAMETERS error: --soloCBmatchWLtype {} does not work with --soloType {}\nSOLUTION: use allowed option: use --soloCBmatchWLtype Exact (exact matches only) OR 1MM (one match with 1 mismatched base)\n",
            p_solo.cb_match_wl.type_, p_solo.type_str
        ));
    }

    p_solo.cb_match_wl.mm1 = false;
    p_solo.cb_match_wl.mm1_multi = false;
    p_solo.cb_match_wl.mm1_multi_pc = false;
    p_solo.cb_match_wl.mm1_multi_nbase = false;
    p_solo.cb_match_wl.one_exact = false;
    p_solo.cb_match_wl.edit_dist_2 = false;

    match p_solo.cb_match_wl.type_.as_str() {
        "Exact" => {
            p_solo.cb_match_wl.one_exact = true;
        }
        "1MM" => {
            p_solo.cb_match_wl.mm1 = true;
            p_solo.cb_match_wl.one_exact = true;
        }
        "1MM_multi" => {
            p_solo.cb_match_wl.mm1 = true;
            p_solo.cb_match_wl.mm1_multi = true;
            p_solo.cb_match_wl.one_exact = true;
        }
        "1MM_multi_pseudocounts" => {
            p_solo.cb_match_wl.mm1 = true;
            p_solo.cb_match_wl.mm1_multi = true;
            p_solo.cb_match_wl.mm1_multi_pc = true;
        }
        "1MM_multi_Nbase_pseudocounts" => {
            p_solo.cb_match_wl.mm1 = true;
            p_solo.cb_match_wl.mm1_multi = true;
            p_solo.cb_match_wl.mm1_multi_pc = true;
            p_solo.cb_match_wl.mm1_multi_nbase = true;
        }
        "EditDist_2" => {
            p_solo.cb_match_wl.edit_dist_2 = true;
        }
        _ => {
            return Err(format!(
                "EXITING because of fatal PARAMETERS error: unrecognized option in --soloCBmatchWLtype {}\nSOLUTION: use allowed options: Exact or 1MM or 1MM_multi or 1MM_multi_pseudocounts 1MM_multi_Nbase_pseudocounts\n",
                p_solo.cb_match_wl.type_
            ));
        }
    }

    Ok(())
}
