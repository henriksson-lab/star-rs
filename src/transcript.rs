#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `Transcript` at STAR/source/Transcript.h:10."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Transcript {
    pub exons: Vec<[u32; crate::include_define::EX_SIZE]>,
    pub cigar: Vec<[u32; 2]>,
    pub canon_sj: Vec<i32>,
    pub sj_annot: Vec<u8>,
    pub shift_sj: Vec<[u32; 2]>,
    pub sj_str: Vec<u8>,
    pub sj_yes: bool,
    pub intron_motifs: [u32; 3],
    pub sj_motif_strand: u32,
    pub n_exons: u32,
    pub l_read: u32,
    pub read_length: Vec<u32>,
    pub read_nmates: u32,
    pub read_length_original: Vec<u32>,
    pub read_length_pair_original: u32,
    pub c_start: u32,
    pub r_start: u32,
    pub ro_start: u32,
    pub ro_str: u32,
    pub r_length: u32,
    pub mapped_length: u32,
    pub g_start: u32,
    pub g_length: u32,
    pub chr: u32,
    pub str_: u32,
    pub i_frag: i32,
    pub primary_flag: bool,
    pub max_score: i32,
    pub n_match: u32,
    pub n_mm: u32,
    pub n_gap: u32,
    pub l_gap: u32,
    pub l_del: u32,
    pub l_ins: u32,
    pub n_del: u32,
    pub n_ins: u32,
    pub n_unique: u32,
    pub n_anchor: u32,
    pub extend_l: u32,
    pub haplo_type: u32,
    pub var_ind: Vec<i32>,
    pub var_gen_coord: Vec<u32>,
    pub var_read_coord: Vec<u32>,
    pub var_allele: Vec<u8>,
    pub i_read: u64,
    pub read_name: String,
}

#[doc = "Original `Transcript::Transcript` at STAR/source/Transcript.cpp:3. Args: "]
pub fn transcript_l3_transcript_transcript() -> crate::transcript::Transcript {
    let mut tr = crate::transcript::Transcript::default();
    transcript_l8_transcript_reset(&mut tr);
    tr
}

#[doc = "Original `Transcript::reset` at STAR/source/Transcript.cpp:8. Args: "]
pub fn transcript_l8_transcript_reset(tr: &mut crate::transcript::Transcript) {
    tr.extend_l = 0;
    tr.primary_flag = false;
    tr.r_start = 0;
    tr.ro_start = 0;
    tr.r_length = 0;
    tr.g_start = 0;
    tr.g_length = 0;
    tr.max_score = 0;
    tr.n_match = 0;
    tr.n_mm = 0;
    tr.n_gap = 0;
    tr.l_gap = 0;
    tr.l_del = 0;
    tr.l_ins = 0;
    tr.n_del = 0;
    tr.n_ins = 0;
    tr.n_unique = 0;
    tr.n_anchor = 0;
}

#[doc = "Original `Transcript::add` at STAR/source/Transcript.cpp:28. Args: trIn: Transcript"]
pub fn transcript_l28_transcript_add(
    tr: &mut crate::transcript::Transcript,
    tr_in: &crate::transcript::Transcript,
) {
    tr.max_score += tr_in.max_score;
    tr.n_match += tr_in.n_match;
    tr.n_mm += tr_in.n_mm;
    tr.n_gap += tr_in.n_gap;
    tr.l_gap += tr_in.l_gap;
    tr.l_del += tr_in.l_del;
    tr.n_del += tr_in.n_del;
    tr.l_ins += tr_in.l_ins;
    tr.n_ins += tr_in.n_ins;
    tr.n_unique += tr_in.n_unique;
}

#[doc = "Original `Transcript::extractSpliceJunctions` at STAR/source/Transcript.cpp:38. Args: sjOut: vector<array<uint64,2>>, annotYes: bool"]
pub fn transcript_l38_transcript_extractsplicejunctions(
    tr: &crate::transcript::Transcript,
    sj_out: &mut Vec<[u64; 2]>,
    annot_yes: &mut bool,
) {
    *annot_yes = true;
    for iex in 0..tr.n_exons.saturating_sub(1) as usize {
        if tr.canon_sj[iex] >= 0 {
            let sj0 = tr.exons[iex][EX_G] as u64 + tr.exons[iex][EX_L] as u64;
            let sj1 = tr.exons[iex + 1][EX_G] as u64 - sj0;
            sj_out.push([sj0, sj1]);
            if tr.sj_annot[iex] == 0 {
                *annot_yes = false;
            }
        }
    }
}

#[doc = "Original `Transcript::chrStartLengthExtended` at STAR/source/Transcript.cpp:53. Args: "]
pub fn transcript_l53_transcript_chrstartlengthextended(
    tr: &crate::transcript::Transcript,
) -> u64 {
    let start1 = tr.c_start as u64 - tr.exons[0][EX_R] as u64;
    let length1 = tr.exons[tr.n_exons as usize - 1][EX_G] as u64 + tr.l_read as u64
        - tr.exons[tr.n_exons as usize - 1][EX_R] as u64
        - tr.exons[0][EX_G] as u64
        + tr.exons[0][EX_R] as u64;
    (start1 << 32) | length1
}
