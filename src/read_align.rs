#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ReadAlign` at STAR/source/ReadAlign.h:22."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlign {
    pub chim_record: bool,
    pub wasp_type: i32,
    pub unmap_type: i32,
    pub i_read: u64,
    pub stats_ra: Stats,
    pub tr_best: Transcript,
    pub l_read: u64,
    pub out_filter_mismatch_nmax_total: u64,
    pub read_nmates: u64,
    pub map_marker: u64,
    pub n_a: u64,
    pub n_p: u64,
    pub n_split: u64,
    pub n_w: u64,
    pub n_wall: u64,
    pub n_tr: u64,
    pub revert_strand: bool,
    pub tr_mult: Vec<Transcript>,
    pub n_win_tr: Vec<u64>,
    pub tr_all: Vec<Vec<Transcript>>,
    pub n_um: [u64; 2],
    pub stored_lmin: u64,
    pub uniq_lmax: u64,
    pub uniq_lmax_ind: u64,
    pub mult_lmax: u64,
    pub mult_lmax_n: u64,
    pub mult_nmin_l: u64,
    pub mult_nmin: u64,
    pub mult_nmax: u64,
    pub mult_nmax_l: u64,
    pub chim_n: u64,
    pub chim_str: u64,
    pub chim_repeat0: u64,
    pub chim_repeat1: u64,
    pub chim_j0: u64,
    pub chim_j1: u64,
    pub chim_motif: i32,
    pub tr_chim: Vec<Transcript>,
    pub max_score_mate: Vec<i32>,
    pub read_length: Vec<u64>,
    pub read_length_original: Vec<u64>,
    pub read_length_pair_original: u64,
    pub read_name: String,
    pub i_read_all: u64,
    pub read_filter: i32,
    pub read_files_index: u32,
    pub read_file_type: i32,
    pub out_bam_bytes: u64,
    pub rng_mult_order_seed: u64,
    pub rng_uniform_real_0_to_1: [f64; 2],
    pub align_tr_all: Vec<Transcript>,
    pub split_r: [Vec<u64>; 3],
    pub tr_array: Vec<Transcript>,
    pub tr_array_pointer_n: usize,
    pub tr_init: Box<Transcript>,
    pub aligns_gen_out_al_mult: Vec<Transcript>,
    pub read0: Vec<Vec<u8>>,
    pub qual0: Vec<Vec<u8>>,
    pub read_name_mates: Vec<Vec<u8>>,
    pub read_name_extra: Vec<String>,
    pub out_bam_one_align_nbytes: Vec<u64>,
    pub out_bam_one_align: Vec<Vec<u8>>,
    pub chunk_out_chim_junction_opened: bool,
    pub solo_read: SoloRead,
    pub spl_graph_present: bool,
    pub qual_hist: Vec<Vec<u64>>,
    pub pe_ov: ReadAlignPeOverlap,
    pub read1: [Vec<u8>; 3],
    pub clip_mates: Vec<Vec<ClipMate>>,
    pub mates_cigar: Vec<String>,
    pub win_bin: [Vec<u32>; 2],
    pub win_bin_touched: Vec<[u32; 2]>,
    pub wc: Vec<[u64; crate::include_define::WC_SIZE]>,
    pub wa: Vec<Vec<[u64; crate::include_define::WA_SIZE]>>,
    pub n_wa: Vec<u64>,
    pub n_wap: Vec<u64>,
    pub wal_rec: Vec<u64>,
    pub w_last_anchor: Vec<u64>,
    pub wa_incl: Vec<bool>,
    pub score_seed_best: Vec<i32>,
    pub score_seed_best_mm: Vec<u64>,
    pub score_seed_best_ind: Vec<u64>,
    pub seed_chain: Vec<u64>,
    pub pc: Vec<[u64; crate::include_define::PC_SIZE]>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignPeOverlap {
    pub yes: bool,
    pub n_ov: u64,
    pub ov_s: u64,
    pub mate_start: [u64; 2],
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StoredAlign {
    pub i_dir: u64,
    pub shift: u64,
    pub n_rep: u64,
    pub l: u64,
    pub ind_start_end: [u64; 2],
    pub i_frag: u64,
}

#[doc = "Original `ReadAlign::ReadAlign` at STAR/source/ReadAlign.cpp:6. Args: Pin: Parameters, genomeIn: Genome, TrIn: Transcriptome, iChunk: int"]
pub fn readalign_l6_readalign_readalign(
    p: &crate::parameters_chimeric::Parameters,
    genome: &crate::genome::Genome,
    _tr_in: Option<&crate::transcriptome::Transcriptome>,
    i_chunk: i32,
) -> crate::read_align::ReadAlign {
    let read_nmates = p.read_nmates as u64;
    let align_windows = p.align_windows_per_read_nmax as usize;
    let seed_per_window = p.seed_per_window_nmax as usize;
    let align_transcripts = p.align_transcripts_per_read_nmax as usize;
    let read_nends = p.read_nends as usize;
    let mut clip_mates = Vec::new();
    if !p.p_clip.adapter_type.is_empty()
        && p.p_clip.in_[0].n.len() >= read_nends
        && p.p_clip.in_[1].n.len() >= read_nends
    {
        parametersclip_initialize_l84_parametersclip_initializeclipmates(
            &p.p_clip,
            &mut clip_mates,
        );
    }

    let mut read_align = crate::read_align::ReadAlign {
        read_nmates,
        rng_mult_order_seed: (p.run_rng_seed as u64).wrapping_mul(i_chunk as u64 + 1),
        rng_uniform_real_0_to_1: [0.0, 1.0],
        align_tr_all: if p.quant_tr_sam_yes {
            vec![crate::transcript::Transcript::default(); align_transcripts]
        } else {
            Vec::new()
        },
        win_bin: [
            vec![UINT_WIN_BIN_MAX; p.win_bin_n as usize],
            vec![UINT_WIN_BIN_MAX; p.win_bin_n as usize],
        ],
        win_bin_touched: Vec::new(),
        split_r: [
            vec![0; p.max_nsplit as usize],
            vec![0; p.max_nsplit as usize],
            vec![0; p.max_nsplit as usize],
        ],
        pc: vec![[0; PC_SIZE]; p.seed_per_read_nmax as usize],
        wc: vec![[0; WC_SIZE]; align_windows],
        n_wa: vec![0; align_windows],
        n_wap: vec![0; align_windows],
        wal_rec: vec![0; align_windows],
        w_last_anchor: vec![0; align_windows],
        wa: vec![vec![[0; WA_SIZE]; seed_per_window]; align_windows],
        wa_incl: vec![false; seed_per_window],
        score_seed_best: vec![0; seed_per_window],
        score_seed_best_mm: vec![0; seed_per_window],
        score_seed_best_ind: vec![0; seed_per_window],
        seed_chain: vec![0; seed_per_window],
        tr_all: vec![Vec::new(); align_windows + 1],
        n_win_tr: vec![0; align_windows],
        tr_array: vec![crate::transcript::Transcript::default(); align_transcripts],
        tr_array_pointer_n: align_transcripts,
        tr_init: Box::new(crate::transcript::Transcript::default()),
        aligns_gen_out_al_mult: if genome.genome_out.conv_yes {
            vec![crate::transcript::Transcript::default(); p.out_sam_mult_nmax as usize]
        } else {
            Vec::new()
        },
        read0: vec![vec![0; DEF_READ_SEQ_LENGTH_MAX + 1]; read_nends],
        qual0: vec![vec![0; DEF_READ_SEQ_LENGTH_MAX + 1]; read_nends],
        read_name_mates: vec![vec![0; DEF_READ_NAME_LENGTH_MAX]; read_nends],
        read_name_extra: vec![String::new(); read_nends],
        read_name: String::new(),
        read1: [
            vec![0; DEF_READ_SEQ_LENGTH_MAX + 1],
            vec![0; DEF_READ_SEQ_LENGTH_MAX + 1],
            vec![0; DEF_READ_SEQ_LENGTH_MAX + 1],
        ],
        qual_hist: vec![vec![0; 256]; read_nends],
        out_bam_one_align_nbytes: vec![0; read_nmates as usize + 2],
        out_bam_one_align: vec![vec![0; 100_000]; read_nmates as usize + 2],
        chunk_out_chim_junction_opened: true,
        solo_read: soloread_l3_soloread_soloread(p, i_chunk),
        clip_mates,
        max_score_mate: vec![0; read_nmates as usize],
        read_length: vec![0; read_nends],
        read_length_original: vec![0; read_nends],
        spl_graph_present: p.p_ge.g_type_string == "SuperTranscriptome",
        ..Default::default()
    };
    readalign_l113_readalign_resetn(&mut read_align);
    read_align
}

#[doc = "Original `ReadAlign::resetN` at STAR/source/ReadAlign.cpp:113. Args: "]
pub fn readalign_l113_readalign_resetn(read_align: &mut crate::read_align::ReadAlign) {
    read_align.map_marker = 0;
    read_align.n_a = 0;
    read_align.n_p = 0;
    read_align.n_w = 0;
    read_align.n_tr = 0;
    read_align.n_um[0] = 0;
    read_align.n_um[1] = 0;
    read_align.stored_lmin = 0;
    read_align.uniq_lmax = 0;
    read_align.uniq_lmax_ind = 0;
    read_align.mult_lmax = 0;
    read_align.mult_lmax_n = 0;
    read_align.mult_nmin_l = 0;
    read_align.mult_nmin = 0;
    read_align.mult_nmax = 0;
    read_align.mult_nmax_l = 0;
    read_align.chim_n = 0;

    for ii in 0..read_align.read_nmates as usize {
        read_align.max_score_mate[ii] = 0;
    }
}
