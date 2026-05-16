#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `Quantifications` at STAR/source/Quantifications.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Quantifications {
    pub gene_counts: QuantificationGeneCounts,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct QuantificationGeneCounts {
    pub n_ge: u32,
    pub n_type: i32,
    pub c_multi: u64,
    pub c_ambig: Vec<u64>,
    pub c_none: Vec<u64>,
    pub g_count: Vec<Vec<u64>>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct QuantTranscriptomeBamRequest {
    pub transcript: Transcript,
    pub n_align_t: u64,
    pub i_align_t: u64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct QuantTranscriptomeResult {
    pub n_align_t: u64,
    pub align_t: Vec<Transcript>,
    pub bam_requests: Vec<QuantTranscriptomeBamRequest>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericBamOutputRequest {
    pub al1: Transcript,
    pub al2: Transcript,
    pub i_tr: u64,
    pub chim_n: u64,
    pub is_best_chim_align: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericDetectionOldOutputResult {
    pub chim_n: u64,
    pub chim_sam: String,
    pub chim_junction: String,
    pub bam_requests: Vec<ChimericBamOutputRequest>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignGenomeTransformResult {
    pub al_best: Transcript,
    pub al_mult: Vec<Transcript>,
    pub al_n: u64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AlignBamRequest {
    pub transcript: Transcript,
    pub n_tr_out: u64,
    pub i_tr_out: u64,
    pub tr_chr_start: u64,
    pub mate_chr: u64,
    pub mate_start: u64,
    pub mate_strand: i8,
    pub align_type: i32,
    pub mate_map: Option<[bool; 2]>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AlignBamResult {
    pub n_lines: u32,
    pub records: [Vec<u8>; 2],
    pub record_sizes: [u32; 2],
    pub sam_flags: [u16; 2],
    pub mapq: [i32; 2],
    pub n_cigar: [u32; 2],
    pub signal_records: [Option<SignalFromBamRecord>; 2],
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericBamOutputResult {
    pub chim_represent: i32,
    pub chim_type: i32,
    pub representative_request_index: i32,
    pub supplementary_request_index: i32,
    pub bam_requests: Vec<AlignBamRequest>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericDetectionRequest {
    pub detector: String,
    pub n_w: u64,
    pub read_length: Vec<u64>,
    pub max_non_chim_align_score: i32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericDetectionResult {
    pub request: Option<ChimericDetectionRequest>,
    pub old_output_requested: bool,
    pub old_output: Option<ChimericDetectionOldOutputResult>,
    pub mult_output: Option<ChimericDetectionMultResult>,
    pub mult_map_select_requested: bool,
    pub mapped_filter_requested: bool,
    pub pe_tr_chim: Vec<Transcript>,
    pub chim_record: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericDetectionMultResult {
    pub chim_record: bool,
    pub chim_n: u64,
    pub chim_score_best: i32,
    pub min_score_to_consider: i32,
    pub chim_junction: String,
    pub bam_outputs: Vec<ChimericBamOutputResult>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeOverlapMergeMapResult {
    pub map_one_read_requested: bool,
    pub chimeric_detection: Option<ChimericDetectionResult>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WriteSamResult {
    pub out_bam_bytes: u64,
    pub unmap_type: i32,
    pub mate_mapped: [bool; 2],
    pub sam: String,
    pub bam_requests: Vec<AlignBamRequest>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct OutputAlignmentsResult {
    pub out_bam_bytes: u64,
    pub unmap_type: i32,
    pub out_filter_by_sjout_pass: bool,
    pub mate_mapped: [bool; 2],
    pub sam: String,
    pub write_sam: WriteSamResult,
    pub quant_transcriptome: Option<QuantTranscriptomeResult>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignOneReadResult {
    pub status: i32,
    pub map_one_read_requested: bool,
    pub splice_graph_log: String,
    pub pe_overlap: PeOverlapMergeMapResult,
    pub chimeric_detection: Option<ChimericDetectionResult>,
    pub wasp: WaspMapResult,
    pub output_alignments: Option<OutputAlignmentsResult>,
    pub tr_mult: Vec<Transcript>,
    pub aligns_gen_out: ReadAlignGenomeTransformResult,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WaspMapRequest {
    pub alleles: Vec<u8>,
    pub read1: [Vec<u8>; 3],
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WaspMapOutcome {
    pub unmap_type: i32,
    pub n_tr: u64,
    pub align: Transcript,
    pub transformed: Option<ReadAlignGenomeTransformResult>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WaspMapResult {
    pub wasp_type: i32,
    pub requests: Vec<WaspMapRequest>,
}

#[doc = "Original `Quantifications::Quantifications` at STAR/source/Quantifications.cpp:3. Args: nGeIn: uint32"]
pub fn quantifications_l3_quantifications_quantifications(
    n_ge_in: u32,
) -> crate::quantifications::Quantifications {
    let n_type = 3usize;
    crate::quantifications::Quantifications {
        gene_counts: crate::quantifications::QuantificationGeneCounts {
            n_ge: n_ge_in,
            n_type: n_type as i32,
            c_multi: 0,
            c_ambig: vec![0; n_type],
            c_none: vec![0; n_type],
            g_count: vec![vec![0; n_ge_in as usize]; n_type],
        },
    }
}

#[doc = "Original `Quantifications::addQuants` at STAR/source/Quantifications.cpp:25. Args: quantsIn: Quantifications"]
pub fn quantifications_l25_quantifications_addquants(
    quantifications: &mut crate::quantifications::Quantifications,
    quants_in: &crate::quantifications::Quantifications,
) {
    quantifications.gene_counts.c_multi += quants_in.gene_counts.c_multi;
    for itype in 0..quantifications.gene_counts.n_type as usize {
        quantifications.gene_counts.c_ambig[itype] += quants_in.gene_counts.c_ambig[itype];
        quantifications.gene_counts.c_none[itype] += quants_in.gene_counts.c_none[itype];
        for ii in 0..quantifications.gene_counts.n_ge as usize {
            quantifications.gene_counts.g_count[itype][ii] +=
                quants_in.gene_counts.g_count[itype][ii];
        }
    }
}
