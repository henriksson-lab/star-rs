#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SpliceGraph` at STAR/source/SpliceGraph.h:12."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpliceGraph {
    pub super_trome: SuperTranscriptome,
    pub ra: Option<ReadAlign>,
    pub super_tr_seed_count: Vec<u16>,
    pub scoring_matrix: Vec<Vec<i32>>,
    pub score_two_columns: [Vec<i32>; 2],
    pub direction_matrix: Vec<u8>,
    pub sj_dindex: Vec<u32>,
    pub gap_penalty: i8,
    pub match_score: i8,
    pub mismatch_penalty: i8,
    pub align_info: SpliceGraphAlignInfo,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpliceGraphAlignInfo {
    pub n_map: u32,
    pub n_mm: u32,
    pub n_i: u32,
    pub n_d: u32,
    pub n_sj: u32,
    pub a_start: [u32; 2],
    pub a_end: [u32; 2],
}

#[doc = "Original `SpliceGraph::SpliceGraph` at STAR/source/SpliceGraph.cpp:8. Args: superTrome: SuperTranscriptome, P: Parameters, RA: ReadAlign"]
pub fn splicegraph_l8_splicegraph_splicegraph(
    super_trome: crate::super_transcriptome::SuperTranscriptome,
    ra: Option<crate::read_align::ReadAlign>,
) -> crate::splice_graph::SpliceGraph {
    let scoring_rows = super_trome.sj_donor_nmax as usize + 2;
    crate::splice_graph::SpliceGraph {
        super_tr_seed_count: vec![0; 2 * super_trome.n as usize],
        scoring_matrix: vec![vec![0; SPLICEGRAPH_MAX_SEQ_LENGTH]; scoring_rows],
        score_two_columns: [
            vec![0; SPLICEGRAPH_MAX_SEQ_LENGTH],
            vec![0; SPLICEGRAPH_MAX_SEQ_LENGTH],
        ],
        direction_matrix: Vec::new(),
        sj_dindex: vec![0; super_trome.sj_donor_nmax as usize],
        gap_penalty: -1,
        match_score: 1,
        mismatch_penalty: -1,
        super_trome,
        ra,
        ..Default::default()
    }
}

#[doc = "Original `SpliceGraph::~SpliceGraph` at STAR/source/SpliceGraph.cpp:28. Args: "]
pub fn splicegraph_l28_splicegraph_splicegraph(
    splice_graph: &mut crate::splice_graph::SpliceGraph,
) {
    splice_graph.scoring_matrix.clear();
    splice_graph.score_two_columns[0].clear();
    splice_graph.score_two_columns[1].clear();
    splice_graph.direction_matrix.clear();
    splice_graph.sj_dindex.clear();
    splice_graph.super_tr_seed_count.clear();
}
