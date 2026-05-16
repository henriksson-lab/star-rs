#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `Transcript` at STAR/source/Transcript.h:10."]
#[derive(Clone, Debug, PartialEq)]
pub struct Transcript {
    /// Per-exon coordinates. STAR's C++ uses `uint exons[MAX_N_EXONS][EX_SIZE]`
    /// — a fixed-size stack array. We match that layout so `clone()` is a
    /// bounded memcpy instead of `MAX_N_EXONS` Vec allocations. The logical
    /// length is `n_exons`; entries beyond that are unused.
    pub exons: [[u64; crate::include_define::EX_SIZE]; crate::include_define::MAX_N_EXONS],
    pub cigar: Vec<[u32; 2]>,
    pub canon_sj: [i32; crate::include_define::MAX_N_EXONS],
    pub sj_annot: [u8; crate::include_define::MAX_N_EXONS],
    pub shift_sj: [[u64; 2]; crate::include_define::MAX_N_EXONS],
    pub sj_str: [u8; crate::include_define::MAX_N_EXONS],
    pub sj_yes: bool,
    pub intron_motifs: [u64; 3],
    pub sj_motif_strand: u8,
    pub n_exons: u64,
    pub l_read: u64,
    pub read_length: [u64; crate::include_define::MAX_N_MATES],
    pub read_nmates: u64,
    pub read_length_original: [u64; crate::include_define::MAX_N_MATES],
    pub read_length_pair_original: u64,
    pub c_start: u64,
    pub r_start: u64,
    pub ro_start: u64,
    pub ro_str: u64,
    pub r_length: u64,
    pub mapped_length: u64,
    pub g_start: u64,
    pub g_length: u64,
    pub chr: u64,
    pub str_: u64,
    pub i_frag: i32,
    pub primary_flag: bool,
    pub max_score: i32,
    pub n_match: u64,
    pub n_mm: u64,
    pub n_gap: u64,
    pub l_gap: u64,
    pub l_del: u64,
    pub l_ins: u64,
    pub n_del: u64,
    pub n_ins: u64,
    pub n_unique: u64,
    pub n_anchor: u64,
    pub extend_l: u64,
    pub haplo_type: u32,
    pub var_ind: Vec<i32>,
    pub var_gen_coord: Vec<i32>,
    pub var_read_coord: Vec<i32>,
    pub var_allele: Vec<u8>,
    pub i_read: u64,
    pub read_name: String,
}

impl Default for Transcript {
    fn default() -> Self {
        Self {
            exons: [[0; crate::include_define::EX_SIZE]; crate::include_define::MAX_N_EXONS],
            cigar: Vec::new(),
            canon_sj: [0; crate::include_define::MAX_N_EXONS],
            sj_annot: [0; crate::include_define::MAX_N_EXONS],
            shift_sj: [[0; 2]; crate::include_define::MAX_N_EXONS],
            sj_str: [0; crate::include_define::MAX_N_EXONS],
            sj_yes: false,
            intron_motifs: [0; 3],
            sj_motif_strand: 0,
            n_exons: 0,
            l_read: 0,
            read_length: [0; crate::include_define::MAX_N_MATES],
            read_nmates: 0,
            read_length_original: [0; crate::include_define::MAX_N_MATES],
            read_length_pair_original: 0,
            c_start: 0,
            r_start: 0,
            ro_start: 0,
            ro_str: 0,
            r_length: 0,
            mapped_length: 0,
            g_start: 0,
            g_length: 0,
            chr: 0,
            str_: 0,
            i_frag: 0,
            primary_flag: false,
            max_score: 0,
            n_match: 0,
            n_mm: 0,
            n_gap: 0,
            l_gap: 0,
            l_del: 0,
            l_ins: 0,
            n_del: 0,
            n_ins: 0,
            n_unique: 0,
            n_anchor: 0,
            extend_l: 0,
            haplo_type: 0,
            var_ind: Vec::new(),
            var_gen_coord: Vec::new(),
            var_read_coord: Vec::new(),
            var_allele: Vec::new(),
            i_read: 0,
            read_name: String::new(),
        }
    }
}

impl Transcript {
    /// Field-wise copy. STAR's C++ does a struct value-copy here (fixed-size
    /// arrays inline). With our matching fixed-size arrays this is a memcpy
    /// for those fields; the few remaining Vec/String fields use `clone_from`
    /// to reuse existing capacity.
    pub fn copy_from(&mut self, other: &Self) {
        self.exons = other.exons;
        self.cigar.clone_from(&other.cigar);
        self.canon_sj = other.canon_sj;
        self.sj_annot = other.sj_annot;
        self.shift_sj = other.shift_sj;
        self.sj_str = other.sj_str;
        self.sj_yes = other.sj_yes;
        self.intron_motifs = other.intron_motifs;
        self.sj_motif_strand = other.sj_motif_strand;
        self.n_exons = other.n_exons;
        self.l_read = other.l_read;
        self.read_length = other.read_length;
        self.read_nmates = other.read_nmates;
        self.read_length_original = other.read_length_original;
        self.read_length_pair_original = other.read_length_pair_original;
        self.c_start = other.c_start;
        self.r_start = other.r_start;
        self.ro_start = other.ro_start;
        self.ro_str = other.ro_str;
        self.r_length = other.r_length;
        self.mapped_length = other.mapped_length;
        self.g_start = other.g_start;
        self.g_length = other.g_length;
        self.chr = other.chr;
        self.str_ = other.str_;
        self.i_frag = other.i_frag;
        self.primary_flag = other.primary_flag;
        self.max_score = other.max_score;
        self.n_match = other.n_match;
        self.n_mm = other.n_mm;
        self.n_gap = other.n_gap;
        self.l_gap = other.l_gap;
        self.l_del = other.l_del;
        self.l_ins = other.l_ins;
        self.n_del = other.n_del;
        self.n_ins = other.n_ins;
        self.n_unique = other.n_unique;
        self.n_anchor = other.n_anchor;
        self.extend_l = other.extend_l;
        self.haplo_type = other.haplo_type;
        self.var_ind.clone_from(&other.var_ind);
        self.var_gen_coord.clone_from(&other.var_gen_coord);
        self.var_read_coord.clone_from(&other.var_read_coord);
        self.var_allele.clone_from(&other.var_allele);
        self.i_read = other.i_read;
        self.read_name.clone_from(&other.read_name);
    }
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
pub fn transcript_l53_transcript_chrstartlengthextended(tr: &crate::transcript::Transcript) -> u64 {
    let start1 = tr.c_start as u64 - tr.exons[0][EX_R] as u64;
    let length1 = tr.exons[tr.n_exons as usize - 1][EX_G] as u64 + tr.l_read as u64
        - tr.exons[tr.n_exons as usize - 1][EX_R] as u64
        - tr.exons[0][EX_G] as u64
        + tr.exons[0][EX_R] as u64;
    (start1 << 32) | length1
}
