#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SoloFilteredCells` at STAR/source/SoloFilteredCells.h:4."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloFilteredCells {
    pub filt_vec_bool: Vec<bool>,
    pub n_cells: u32,
    pub n_read_in_cells: u64,
    pub median_read_per_cell_unique: u64,
    pub mean_read_per_cell_unique: u64,
    pub n_umi_in_cells: u64,
    pub median_umi_per_cell: u64,
    pub mean_umi_per_cell: u64,
    pub n_gene_in_cells: u64,
    pub median_gene_per_cell: u64,
    pub mean_gene_per_cell: u64,
    pub n_gene_detected: u64,
    pub n_cells_simple: u64,
    pub n_read_in_cells_unique: u64,
    pub n_read_per_cell_unique: Vec<u32>,
    pub n_gene_per_cell: Vec<u32>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloOutputResults {
    pub created_directory: String,
    pub files: std::collections::BTreeMap<String, String>,
    pub symlinks: Vec<(String, String)>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloCellFilteringResult {
    pub log_main: String,
    pub empty_drops_requested: bool,
    pub empty_drops: Option<SoloEmptyDropsCrResult>,
    pub output_results: Option<SoloOutputResults>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloEmptyDropsCrResult {
    pub log_main: String,
    pub early_return: Option<String>,
    pub feat_det_n: u32,
    pub min_umi: u32,
    pub candidate_first: u32,
    pub candidate_last: u32,
    pub extra_cells: u32,
    pub p_values: Vec<(u32, f64, f64)>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloFeatureProcessRecordsResult {
    pub created_directories: Vec<String>,
    pub files: std::collections::BTreeMap<String, String>,
    pub symlinks: Vec<(String, String)>,
    pub log_main: String,
    pub stats_output: Option<SoloStatsOutput>,
    pub cell_filtering: Option<SoloCellFilteringResult>,
    pub returned_after_quant_transcript: bool,
    pub count_cb_gene_umi_called: bool,
    pub quant_transcript_called: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloProcessAndOutputResult {
    pub files: std::collections::BTreeMap<String, String>,
    pub symlinks: Vec<(String, String)>,
    pub created_directories: Vec<String>,
    pub log_stdout: String,
    pub log_main: String,
    pub feature_results: Vec<SoloFeatureProcessRecordsResult>,
    pub returned_after_barcode_output: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloFeatureQuantTranscriptResult {
    pub files: std::collections::BTreeMap<String, String>,
    pub log_main: String,
    pub stdout: String,
    pub cluster_expression: std::collections::BTreeMap<u32, Vec<f64>>,
    pub returned_no_cluster_file: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloConstructorCellFilteringResult {
    pub solo: Solo,
    pub cell_filtering: Option<SoloCellFilteringResult>,
    pub log_stdout: String,
    pub log_main: String,
    pub exited: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloStatsOutput {
    pub files: std::collections::BTreeMap<String, String>,
}
