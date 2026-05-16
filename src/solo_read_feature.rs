#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SoloReadFeature` at STAR/source/SoloReadFeature.h:15."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloReadFeature {
    pub feature_type: i32,
    pub read_info_yes: bool,
    pub read_index_yes: bool,
    pub stream_reads_path: Option<String>,
    pub stream_reads: String,
    pub cb_wl_yes: bool,
    pub cb_wl_size: u32,
    pub cb_read_count: Vec<u32>,
    pub cb_read_count_map: std::collections::BTreeMap<u64, u32>,
    pub transcript_dist_count: Vec<u32>,
    pub read_flag: SoloReadFlagClass,
    pub stats: SoloReadFeatureStats,
}

#[doc = "Original `SoloReadFeature::SoloReadFeature` at STAR/source/SoloReadFeature.cpp:5. Args: feTy: int32, Pin: Parameters, iChunk: int"]
pub fn soloreadfeature_l5_soloreadfeature_soloreadfeature(
    feature_type: i32,
    p: &crate::parameters_chimeric::Parameters,
    i_chunk: i32,
) -> crate::solo_read_feature::SoloReadFeature {
    let p_solo = &p.p_solo;
    let mut rf = crate::solo_read_feature::SoloReadFeature {
        feature_type,
        cb_wl_yes: p_solo.cb_wl_yes,
        cb_wl_size: p_solo.cb_wl_size,
        ..Default::default()
    };
    if p_solo.solo_type == 0 {
        return rf;
    }

    rf.read_info_yes = p_solo
        .read_info_yes
        .get(feature_type as usize)
        .copied()
        .unwrap_or(false);
    rf.read_index_yes = p_solo
        .read_index_yes
        .get(feature_type as usize)
        .copied()
        .unwrap_or(false);

    if p_solo.cb_wl_yes {
        rf.cb_read_count.resize(p_solo.cb_wl_size as usize, 0);
    }
    if i_chunk >= 0 {
        let names = [
            "SJ",
            "Transcript3p",
            "GeneFull",
            "GeneFull_ExonOverIntron",
            "GeneFull_Ex50pAS",
            "Gene",
            "VelocytoSimple",
            "Velocyto",
        ];
        let name = names.get(feature_type as usize).copied().unwrap_or("");
        rf.stream_reads_path = Some(format!("{}/solo{}_{}", p.out_file_tmp, name, i_chunk));
    }
    if feature_type == SOLO_FEATURE_TRANSCRIPT3P {
        rf.transcript_dist_count.resize(10000, 0);
    }
    rf
}

#[doc = "Original `SoloReadFeature::addCounts` at STAR/source/SoloReadFeature.cpp:29. Args: rfIn: SoloReadFeature"]
pub fn soloreadfeature_l29_soloreadfeature_addcounts(
    rf: &mut crate::solo_read_feature::SoloReadFeature,
    rf_in: &crate::solo_read_feature::SoloReadFeature,
) {
    if rf.cb_wl_yes {
        for ii in 0..rf.cb_wl_size as usize {
            rf.cb_read_count[ii] += rf_in.cb_read_count[ii];
        }
    } else {
        for (cb, count) in rf_in.cb_read_count_map.iter() {
            *rf.cb_read_count_map.entry(*cb).or_insert(0) += *count;
        }
    }

    if !rf.transcript_dist_count.is_empty() {
        for ii in 0..rf.transcript_dist_count.len() {
            rf.transcript_dist_count[ii] += rf_in.transcript_dist_count[ii];
        }
    }
}

#[doc = "Original `SoloReadFeature::addStats` at STAR/source/SoloReadFeature.cpp:47. Args: rfIn: SoloReadFeature"]
pub fn soloreadfeature_l47_soloreadfeature_addstats(
    rf: &mut crate::solo_read_feature::SoloReadFeature,
    rf_in: &crate::solo_read_feature::SoloReadFeature,
) {
    for ii in 0..rf.stats.v.len() {
        rf.stats.v[ii] += rf_in.stats.v[ii];
    }
    for ii in 0..SOLO_READ_FLAG_N_BITS {
        rf.read_flag.flag_counts_no_cb[ii] += rf_in.read_flag.flag_counts_no_cb[ii];
    }
}

#[doc = "Original `SoloReadFeature::statsOut` at STAR/source/SoloReadFeature.cpp:56. Args: streamOut: ofstream"]
pub fn soloreadfeature_l56_soloreadfeature_statsout(
    rf: &crate::solo_read_feature::SoloReadFeature,
) -> String {
    use std::fmt::Write;

    let mut out = String::new();
    for ii in 0..rf.stats.v.len() {
        writeln!(out, "{:>50}{:>15}", rf.stats.names[ii], rf.stats.v[ii]).unwrap();
    }
    out
}
