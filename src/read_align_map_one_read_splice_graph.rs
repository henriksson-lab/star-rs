#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::mapOneReadSpliceGraph` at STAR/source/ReadAlign_mapOneReadSpliceGraph.cpp:6. Args: "]
pub fn readalign_maponereadsplicegraph_l6_readalign_maponereadsplicegraph(
    read_align: &mut crate::read_align::ReadAlign,
    splice_graph: &mut crate::splice_graph::SpliceGraph,
    map_gen: &crate::genome::Genome,
    out_filter_match_nmin: u32,
    seed_multimap_nmax: u32,
) -> String {
    if read_align.l_read < out_filter_match_nmin as u64 {
        read_align.map_marker = MARKER_READ_TOO_SHORT as u64;
        read_align.tr_best.r_length = 0;
        read_align.n_w = 0;
        return String::new();
    }

    readalign_l113_readalign_resetn(read_align);

    let mut tr_init = transcript_l3_transcript_transcript();
    tr_init.chr = 0;
    tr_init.str_ = 0;
    tr_init.ro_str = 0;
    tr_init.c_start = 0;
    tr_init.g_length = 0;
    tr_init.i_read = read_align.i_read_all;
    tr_init.l_read = read_align.l_read;
    tr_init.read_name = read_align.read_name.clone();
    tr_init.n_exons = 0;
    for (i, &v) in read_align
        .read_length_original
        .iter()
        .take(crate::include_define::MAX_N_MATES)
        .enumerate()
    {
        tr_init.read_length_original[i] = v;
    }
    tr_init.read_length_pair_original = read_align.read_length_pair_original;
    for (i, &v) in read_align
        .read_length
        .iter()
        .take(crate::include_define::MAX_N_MATES)
        .enumerate()
    {
        tr_init.read_length[i] = v;
    }
    tr_init.read_nmates = read_align.read_nmates;
    read_align.tr_best = tr_init;
    splice_graph.ra = Some(read_align.clone());

    let l_read = read_align.l_read as usize;
    let out = splicegraph_findsupertr_l5_splicegraph_findsupertr(
        splice_graph,
        &read_align.read1[0][..l_read],
        &read_align.read1[2][..l_read],
        read_align.l_read,
        &read_align.read_name,
        map_gen,
        seed_multimap_nmax,
    );

    if let Some(ra_out) = splice_graph.ra.as_ref() {
        read_align.tr_best = ra_out.tr_best.clone();
        read_align.n_w = ra_out.n_w;
        read_align.n_win_tr = ra_out.n_win_tr.clone();
        read_align.tr_all = ra_out.tr_all.clone();
    }

    out
}
