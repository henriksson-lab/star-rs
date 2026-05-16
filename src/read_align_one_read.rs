#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::oneRead` at STAR/source/ReadAlign_oneRead.cpp:8. Args: "]
pub fn readalign_oneread_l8_readalign_oneread(
    read_align: &mut crate::read_align::ReadAlign,
    read_in_streams: &mut [&mut dyn std::io::BufRead],
    p: &mut crate::parameters_chimeric::Parameters,
    map_gen: &crate::genome::Genome,
    transcriptome: &mut crate::transcriptome::Transcriptome,
    mut splice_graph: Option<&mut crate::splice_graph::SpliceGraph>,
    mapped_standard_ra: Option<&crate::read_align::ReadAlign>,
    mapped_pe_merge_ra: Option<&crate::read_align::ReadAlign>,
    pe_merge_ra: &mut crate::read_align::ReadAlign,
    wasp_ra: &mut crate::read_align::ReadAlign,
    gen_out: Option<&crate::genome::Genome>,
    wasp_remap_outcomes: &[crate::quantifications::WaspMapOutcome],
    chunk_out_sj: &mut crate::out_sj::OutSJ,
    chunk_out_sj1: &mut crate::out_sj::OutSJ,
    chunk_out_filter_by_sjout_files: &mut [String],
    chunk_out_unmapped_reads_stream: &mut [String],
    out_sam_stream: &mut String,
    primary_pick_fraction: f64,
    chim_detector_result: Option<bool>,
) -> Result<crate::quantifications::ReadAlignOneReadResult, String> {
    let mut result = crate::quantifications::ReadAlignOneReadResult::default();
    let read_nends = p.read_nends as usize;
    if read_in_streams.len() < read_nends {
        return Err("EXITING because of FATAL ERROR: not enough read input streams\n".to_string());
    }
    if read_align.read_length.len() < read_nends.max(2) {
        read_align.read_length.resize(read_nends.max(2), 0);
    }
    if read_align.read_length_original.len() < read_nends.max(2) {
        read_align.read_length_original.resize(read_nends.max(2), 0);
    }
    if read_align.read0.len() < read_nends {
        read_align.read0.resize(read_nends, Vec::new());
    }
    if read_align.qual0.len() < read_nends {
        read_align.qual0.resize(read_nends, Vec::new());
    }
    if read_align.read_name_mates.len() < read_nends {
        read_align
            .read_name_mates
            .resize(read_nends, vec![0; DEF_READ_NAME_LENGTH_MAX]);
    }
    if read_align.read_name_extra.len() < read_nends {
        read_align.read_name_extra.resize(read_nends, String::new());
    }
    if read_align.clip_mates.len() < read_nends {
        read_align.clip_mates.resize(
            read_nends,
            vec![crate::clip_mate::ClipMate::default(); 2],
        );
    }
    if read_align.qual_hist.len() < p.read_nmates as usize {
        read_align
            .qual_hist
            .resize(p.read_nmates as usize, vec![0; 256]);
    }

    let mut read_status = vec![0i32; read_nends];
    let mut read0 = vec![String::new(); read_nends];
    let mut qual0 = vec![String::new(); read_nends];
    let mut read_name_mates = vec![String::new(); read_nends];
    let mut read_filter = read_align.read_filter as u8;

    for im in 0..read_nends {
        if read_align.clip_mates[im].len() < 2 {
            read_align.clip_mates[im].resize(2, crate::clip_mate::ClipMate::default());
        }
        read_status[im] = readload_l4_readload(
            read_in_streams[im],
            p,
            &mut read_align.read_length[im],
            &mut read_align.read_length_original[im],
            &mut read_name_mates[im],
            &mut read0[im],
            &mut read_align.read1[im],
            &mut qual0[im],
            &mut read_align.clip_mates[im],
            &mut read_align.i_read_all,
            &mut read_align.read_files_index,
            &mut read_filter,
            &mut read_align.read_name_extra[im],
        )?;
        if read_status[im] != read_status[0] {
            return Err("EXITING because of FATAL ERROR: read files are not consistent, reached the end of the one before the other one\nSOLUTION: Check you your input files: they may be corrupted\n".to_string());
        }
    }
    read_align.read_filter = read_filter as i32;

    if read_status[0] == -1 {
        result.status = -1;
        return Ok(result);
    }

    read_align.read_name = read_name_mates.first().cloned().unwrap_or_default();
    for im in 0..read_nends {
        read_align.read0[im] = read0[im].as_bytes().to_vec();
        read_align.qual0[im] = qual0[im].as_bytes().to_vec();
        read_align.read_name_mates[im] = read_name_mates[im].as_bytes().to_vec();
    }

    if p.out_filter_by_sjout_stage != 2 {
        for im in 0..p.read_nmates as usize {
            let left = read_align.clip_mates[im][0].clipped_n as usize;
            let right = read_align.clip_mates[im][1].clipped_n as usize;
            let end = (read_align.read_length_original[im] as usize).saturating_sub(right);
            for ix in left..end {
                if let Some(q) = qual0[im].as_bytes().get(ix) {
                    read_align.qual_hist[im][*q as usize] += 1;
                }
            }
        }
    }

    if p.read_nmates == 2 {
        read_align.l_read = read_align.read_length[0] + read_align.read_length[1] + 1;
        read_align.read_length_pair_original =
            read_align.read_length_original[0] + read_align.read_length_original[1] + 1;
        if read_align.l_read as usize > DEF_READ_SEQ_LENGTH_MAX {
            return Err(format!(
                "EXITING because of FATAL ERROR in reads input: Lread of the pair = {}   while DEF_readSeqLengthMax={}\nRead Name={}\nSOLUTION: increase DEF_readSeqLengthMax in IncludeDefine.h and re-compile STAR\n",
                read_align.l_read,
                DEF_READ_SEQ_LENGTH_MAX,
                read_name_mates.first().map(String::as_str).unwrap_or("")
            ));
        }

        let l0 = read_align.read_length[0] as usize;
        let l1 = read_align.read_length[1] as usize;
        if read_align.read1[0].len() < read_align.l_read as usize {
            read_align.read1[0].resize(read_align.l_read as usize, 0);
        }
        read_align.read1[0][l0] = MARK_FRAG_SPACER_BASE;
        let (mate0, mate1_and_rest) = read_align.read1.split_at_mut(1);
        sequencefuns_l4_complementseqnumbers(
            &mate1_and_rest[0],
            &mut mate0[0][l0 + 1..l0 + 1 + l1],
            read_align.read_length[1],
        );
        for ii in 0..read_align.read_length[1] / 2 {
            let left = ii as usize + l0 + 1;
            let right = read_align.l_read as usize - ii as usize - 1;
            read_align.read1[0].swap(right, left);
        }
    } else {
        read_align.l_read = read_align.read_length[0];
        read_align.read_length_pair_original = read_align.read_length_original[0];
        read_align.read_length[1] = 0;
    }

    read_align.read_file_type = read_status[0];
    if read_align.read1[1].len() < read_align.l_read as usize {
        read_align.read1[1].resize(read_align.l_read as usize, 0);
    }
    let (read_forward, read_reverse_and_rest) = read_align.read1.split_at_mut(1);
    sequencefuns_l4_complementseqnumbers(
        &read_forward[0],
        &mut read_reverse_and_rest[0],
        read_align.l_read,
    );
    if read_align.read1[2].len() < read_align.l_read as usize {
        read_align.read1[2].resize(read_align.l_read as usize, 0);
    }
    for ii in 0..read_align.l_read as usize {
        read_align.read1[2][read_align.l_read as usize - ii - 1] = read_align.read1[1][ii];
    }

    read_align.stats_ra.read_n += 1;
    read_align.stats_ra.read_bases += (read_align.read_length[0] + read_align.read_length[1]) as u32;
    let read_bases = read_align.read_length[0] + read_align.read_length[1];
    read_align.out_filter_mismatch_nmax_total = std::cmp::min(
        p.out_filter_mismatch_nmax as u64,
        (p.out_filter_mismatch_nover_read_lmax * read_bases as f64) as u64,
    );

    if p.p_ge.g_type_string == "SpliceGraph" || p.p_ge.g_type_string == "SuperTranscriptome" {
        let mut default_splice_graph;
        let splice_graph = if let Some(splice_graph) = splice_graph.as_deref_mut() {
            splice_graph
        } else {
            default_splice_graph = crate::splice_graph::SpliceGraph::default();
            &mut default_splice_graph
        };
        result.splice_graph_log =
            readalign_maponereadsplicegraph_l6_readalign_maponereadsplicegraph(
                read_align,
                splice_graph,
                map_gen,
                p.out_filter_match_nmin,
                p.seed_multimap_nmax,
            );
    } else if let Some(mapped) = mapped_standard_ra {
        let read_state = (
            read_align.read_length.clone(),
            read_align.read_length_original.clone(),
            read_align.read_length_pair_original,
            read_align.l_read,
            read_align.read_name.clone(),
            read_align.i_read_all,
            read_align.read_filter,
            read_align.read_files_index,
            read_align.read_file_type,
            read_align.read0.clone(),
            read_align.qual0.clone(),
            read_align.read_name_mates.clone(),
            read_align.read_name_extra.clone(),
            read_align.read1.clone(),
            read_align.clip_mates.clone(),
            read_align.qual_hist.clone(),
            read_align.stats_ra.clone(),
            read_align.out_filter_mismatch_nmax_total,
        );
        *read_align = mapped.clone();
        read_align.read_length = read_state.0;
        read_align.read_length_original = read_state.1;
        read_align.read_length_pair_original = read_state.2;
        read_align.l_read = read_state.3;
        read_align.read_name = read_state.4;
        read_align.i_read_all = read_state.5;
        read_align.read_filter = read_state.6;
        read_align.read_files_index = read_state.7;
        read_align.read_file_type = read_state.8;
        read_align.read0 = read_state.9;
        read_align.qual0 = read_state.10;
        read_align.read_name_mates = read_state.11;
        read_align.read_name_extra = read_state.12;
        read_align.read1 = read_state.13;
        read_align.clip_mates = read_state.14;
        read_align.qual_hist = read_state.15;
        read_align.stats_ra = read_state.16;
        read_align.out_filter_mismatch_nmax_total = read_state.17;
    } else {
        result.map_one_read_requested = true;
        readalign_maponeread_l6_readalign_maponeread(read_align, map_gen, p)?;
    }

    result.pe_overlap = readalign_peoverlapmergemap_l4_readalign_peoverlapmergemap(
        read_align,
        pe_merge_ra,
        mapped_pe_merge_ra,
        p,
        map_gen,
        &read_align.tr_init.clone(),
        &map_gen.g,
        p.p_ge.sjdb_score,
        p.score_ins_base,
        p.score_ins_open,
        p.score_del_base,
        p.score_del_open,
        p.score_gap_noncan,
        p.score_gap,
        p.score_gap_gcag,
        p.score_gap_atac,
        p.score_genomic_length_log2scale,
        chim_detector_result,
    )?;

    let rng_snapshot = read_align.rng_uniform_real_0_to_1;
    result.tr_mult = readalign_multmapselect_inner(
        &mut read_align.n_tr,
        &mut read_align.tr_best,
        read_align.l_read,
        read_align.n_w,
        &read_align.n_win_tr,
        map_gen,
        &read_align.tr_all,
        p.out_filter_multimap_score_range,
        p.out_filter_multimap_nmax,
        p.out_multimapper_order_random,
        p.out_sam_mult_nmax_is_limited,
        &p.out_sam_primary_flag,
        &rng_snapshot,
    )?;
    if p.read_nmates == 2 {
        expand_collapsed_same_locus_paired_multimaps(
            read_align,
            map_gen,
            &mut result.tr_mult,
            read_align.tr_best.max_score,
            p.out_multimapper_order_random,
            p.out_sam_mult_nmax_is_limited,
            &p.out_sam_primary_flag,
            &rng_snapshot,
        );
    }
    readalign_mappedfilter_l3_readalign_mappedfilter(
        read_align,
        p.out_filter_score_min,
        p.out_filter_score_min_over_lread,
        p.out_filter_match_nmin,
        p.out_filter_match_nmin_over_lread,
        p.out_filter_mismatch_nover_lmax,
        p.out_filter_multimap_nmax,
    );

    result.aligns_gen_out = readalign_transformgenome_l5_readalign_transformgenome(
        read_align,
        map_gen,
        &mut result.tr_mult,
        p.out_filter_multimap_nmax,
        read_align.tr_best.max_score,
        p.out_multimapper_order_random,
        p.out_sam_mult_nmax_is_limited,
        &p.out_sam_primary_flag,
        &rng_snapshot,
    );

    if !read_align.pe_ov.yes {
        result.chimeric_detection =
            Some(readalign_chimericdetection_l16_readalign_chimericdetection(
                read_align,
                p,
                map_gen,
                chim_detector_result,
            )?);
    }

    if p.p_ch.out_bam && read_align.chim_record {
        result.status = 0;
        return Ok(result);
    }

    result.wasp = readalign_waspmap_l3_readalign_waspmap(
        read_align,
        wasp_ra,
        p,
        map_gen,
        gen_out,
        Some(&result.aligns_gen_out),
        wasp_remap_outcomes,
    );

    let mut solo_read = std::mem::take(&mut read_align.solo_read);
    let output = readalign_outputalignments_l5_readalign_outputalignments(
        read_align,
        p,
        map_gen,
        transcriptome,
        &mut solo_read,
        &result.tr_mult,
        &result.aligns_gen_out,
        chunk_out_sj,
        chunk_out_sj1,
        chunk_out_filter_by_sjout_files,
        chunk_out_unmapped_reads_stream,
        &read_name_mates,
        &read_align.read_name_extra.clone(),
        &read0,
        &qual0,
        read_align.read_file_type,
        out_sam_stream,
        primary_pick_fraction,
    )?;
    read_align.solo_read = solo_read;
    read_align.out_bam_bytes += output.out_bam_bytes;
    result.output_alignments = Some(output);
    result.status = 0;
    Ok(result)
}
