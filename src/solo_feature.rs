#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SoloFeature` at STAR/source/SoloFeature.h:18."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloFeature {
    pub feature_type: i32,
    pub features_number: i64,
    pub read_bar_sum: Option<SoloReadBarcode>,
    pub read_feat_sum: Option<SoloReadFeature>,
    pub read_feat_all: Vec<SoloReadFeature>,
    pub read_feat_all_len: usize,
    pub n_reads_mapped: u64,
    pub n_reads_input: u64,
    pub n_cb: u32,
    pub rgu_stride: u32,
    pub r_cbp: Vec<Vec<u32>>,
    pub cb_feature_umi_map: Vec<()>,
    pub count_cell_gene_umi: Vec<u32>,
    pub count_cell_gene_umi_index: Vec<u32>,
    pub count_mat_mult_i: Vec<u32>,
    pub count_mat_mult_m: Vec<f64>,
    pub count_mat_mult_s: u32,
    pub count_mat_stride: u32,
    pub output_prefix: String,
    pub output_prefix_filtered: String,
    pub copied_features_tsv: String,
    pub ind_cb_wl: Vec<u32>,
    pub ind_cb: Vec<u32>,
    pub n_gene_per_cb: Vec<u32>,
    pub n_gene_per_cb_multi: Vec<u32>,
    pub n_read_per_cb: Vec<u32>,
    pub n_read_per_cb_total: Vec<u32>,
    pub n_read_per_cb_unique: Vec<u32>,
    pub n_umi_per_cb: Vec<u32>,
    pub n_umi_per_cb_multi: Vec<u32>,
    pub n_umi_per_cb_sorted: Vec<u32>,
    pub read_info: Vec<SoloFeatureReadInfo>,
    pub read_flag_counts: SoloReadFlagClass,
    pub sj_all: [Vec<u64>; 2],
    pub filtered_cells: SoloFilteredCells,
    pub redistr_files_cb_first: Vec<u32>,
    pub redistr_files_cb_index: Vec<u32>,
    pub redistr_files_nreads: Vec<u64>,
    pub redistr_files_streams: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloFeatureReadInfo {
    pub cb: i64,
    pub umi: u64,
}

#[doc = "Original `SoloFeature::SoloFeature` at STAR/source/SoloFeature.cpp:4. Args: Pin: Parameters, RAchunk: ReadAlignChunk, inTrans: Transcriptome, feTy: int32, readBarSumIn: SoloReadBarcode, soloFeatAll: SoloFeature"]
pub fn solofeature_l4_solofeature_solofeature(
    p: &crate::parameters_chimeric::Parameters,
    feature_type: i32,
    trans_n_ge: u32,
) -> crate::solo_feature::SoloFeature {
    let read_feat_sum = if feature_type >= 0 {
        Some(soloreadfeature_l5_soloreadfeature_soloreadfeature(
            feature_type,
            p,
            -1,
        ))
    } else {
        None
    };
    let features_number = match feature_type {
        SOLO_FEATURE_GENE
        | SOLO_FEATURE_GENE_FULL
        | SOLO_FEATURE_GENE_FULL_EX50P_AS
        | SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON
        | SOLO_FEATURE_VELOCYTO => trans_n_ge as i64,
        SOLO_FEATURE_SJ => p.sj_all[0].len() as i64,
        _ => -1,
    };
    crate::solo_feature::SoloFeature {
        feature_type,
        features_number,
        read_feat_sum,
        read_feat_all_len: if feature_type >= 0 {
            p.run_thread_n as usize
        } else {
            0
        },
        ..Default::default()
    }
}

#[doc = "Original `SoloFeature::clearLarge` at STAR/source/SoloFeature.cpp:29. Args: "]
pub fn solofeature_l29_solofeature_clearlarge(
    solo_feature: &mut crate::solo_feature::SoloFeature,
) {
    solo_feature.cb_feature_umi_map.clear();
    solo_feature.cb_feature_umi_map.shrink_to_fit();
    solo_feature.count_cell_gene_umi.clear();
    solo_feature.count_cell_gene_umi.shrink_to_fit();
    solo_feature.count_cell_gene_umi_index.clear();
    solo_feature.count_cell_gene_umi_index.shrink_to_fit();
    solo_feature.count_mat_mult_i.clear();
    solo_feature.count_mat_mult_i.shrink_to_fit();
    solo_feature.count_mat_mult_m.clear();
    solo_feature.count_mat_mult_m.shrink_to_fit();
    solo_feature.ind_cb_wl.clear();
    solo_feature.ind_cb_wl.shrink_to_fit();
    solo_feature.ind_cb.clear();
    solo_feature.ind_cb.shrink_to_fit();
    solo_feature.n_gene_per_cb.clear();
    solo_feature.n_gene_per_cb.shrink_to_fit();
    solo_feature.n_gene_per_cb_multi.clear();
    solo_feature.n_gene_per_cb_multi.shrink_to_fit();
    solo_feature.n_read_per_cb.clear();
    solo_feature.n_read_per_cb.shrink_to_fit();
    solo_feature.n_read_per_cb_total.clear();
    solo_feature.n_read_per_cb_total.shrink_to_fit();
    solo_feature.n_read_per_cb_unique.clear();
    solo_feature.n_read_per_cb_unique.shrink_to_fit();
    solo_feature.n_umi_per_cb.clear();
    solo_feature.n_umi_per_cb.shrink_to_fit();
    solo_feature.n_umi_per_cb_multi.clear();
    solo_feature.n_umi_per_cb_multi.shrink_to_fit();
    solo_feature.n_umi_per_cb_sorted.clear();
    solo_feature.n_umi_per_cb_sorted.shrink_to_fit();
    solo_feature.sj_all[0].clear();
    solo_feature.sj_all[0].shrink_to_fit();
    solo_feature.sj_all[1].clear();
    solo_feature.sj_all[1].shrink_to_fit();
}
