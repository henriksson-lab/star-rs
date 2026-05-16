#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::peOverlapMergeMap` at STAR/source/ReadAlign_peOverlapMergeMap.cpp:4. Args: "]
pub fn readalign_peoverlapmergemap_l4_readalign_peoverlapmergemap(
    read_align: &mut crate::read_align::ReadAlign,
    pe_merge_ra: &mut crate::read_align::ReadAlign,
    mapped_pe_merge_ra: Option<&crate::read_align::ReadAlign>,
    p: &crate::parameters_chimeric::Parameters,
    map_gen: &crate::genome::Genome,
    tr_init: &crate::transcript::Transcript,
    genome: &[u8],
    sjdb_score: i32,
    score_ins_base: i32,
    score_ins_open: i32,
    score_del_base: i32,
    score_del_open: i32,
    score_gap_noncan: i32,
    score_gap: i32,
    score_gap_gcag: i32,
    score_gap_atac: i32,
    score_genomic_length_log2scale: f64,
    chim_detector_result: Option<bool>,
) -> Result<crate::quantifications::PeOverlapMergeMapResult, String> {
    let mut result = crate::quantifications::PeOverlapMergeMapResult::default();
    if p.pe_overlap_nbases_min == 0 || p.read_nmates != 2 {
        read_align.pe_ov.yes = false;
        return Ok(result);
    }

    readalign_waspmap_l115_readalign_copyread(pe_merge_ra, read_align);
    readalign_peoverlapmergemap_l79_readalign_pemergemates(
        pe_merge_ra,
        p.pe_overlap_mmp,
        p.pe_overlap_nbases_min,
    );
    read_align.pe_ov = pe_merge_ra.pe_ov.clone();
    read_align.pe_ov.yes = false;

    if read_align.pe_ov.n_ov == 0 {
        return Ok(result);
    }

    result.map_one_read_requested = true;
    if let Some(mapped) = mapped_pe_merge_ra {
        *pe_merge_ra = mapped.clone();
    } else {
        readalign_maponeread_l6_readalign_maponeread(pe_merge_ra, map_gen, p)?;
    }

    if pe_merge_ra.n_w == 0 {
        return Ok(result);
    }

    let pe_score = read_align.tr_best.max_score;
    pe_merge_ra.pe_ov = read_align.pe_ov.clone();
    readalign_peoverlapmergemap_l266_readalign_peoverlapsetope(
        read_align,
        pe_merge_ra,
        tr_init,
        genome,
        sjdb_score,
        score_ins_base,
        score_ins_open,
        score_del_base,
        score_del_open,
        score_gap_noncan,
        score_gap,
        score_gap_gcag,
        score_gap_atac,
        score_genomic_length_log2scale,
    );

    let chim_result = readalign_chimericdetectionpemerged_l5_readalign_chimericdetectionpemerged(
        read_align,
        pe_merge_ra,
        p,
        map_gen,
        chim_detector_result,
    )?;
    result.chimeric_detection = Some(chim_result);

    if pe_score <= read_align.tr_best.max_score || read_align.chim_record {
        read_align.pe_ov.yes = true;
    }

    Ok(result)
}

#[doc = "Original `ReadAlign::peMergeMates` at STAR/source/ReadAlign_peOverlapMergeMap.cpp:79. Args: "]
pub fn readalign_peoverlapmergemap_l79_readalign_pemergemates(
    ra: &mut crate::read_align::ReadAlign,
    pe_overlap_mmp: f64,
    pe_overlap_nbases_min: u32,
) {
    let s1 = sequencefuns_l317_localsearchnismm(
        &ra.read1[0],
        ra.read_length[0],
        &ra.read1[0][ra.read_length[0] as usize + 1..],
        ra.read_length[1],
        pe_overlap_mmp,
    );
    let s0 = sequencefuns_l317_localsearchnismm(
        &ra.read1[0][ra.read_length[0] as usize + 1..],
        ra.read_length[1],
        &ra.read1[0],
        ra.read_length[0],
        pe_overlap_mmp,
    );

    let o1 = ra.read_length[1].min(ra.read_length[0] - s1);
    let o0 = ra.read_length[0].min(ra.read_length[1] - s0);

    ra.pe_ov.n_ov = o0.max(o1);
    if ra.pe_ov.n_ov < pe_overlap_nbases_min {
        ra.pe_ov.n_ov = 0;
        return;
    }

    if o1 >= o0 {
        ra.pe_ov.mate_start[0] = 0;
        ra.pe_ov.mate_start[1] = s1;
        if o1 < ra.read_length[1] {
            let src_start = ra.read_length[0] as usize + 1 + o1 as usize;
            let src_end = ra.read_length[0] as usize + 1 + ra.read_length[1] as usize;
            let dst_start = ra.read_length[0] as usize;
            ra.read1[0].copy_within(src_start..src_end, dst_start);
        }
    } else {
        ra.pe_ov.mate_start[1] = 0;
        ra.pe_ov.mate_start[0] = s0;
        let l_read = ra.l_read as usize;
        let read0_len = ra.read_length[0] as usize;
        let read1_len = ra.read_length[1] as usize;
        if ra.read1[0].len() < l_read + read0_len {
            ra.read1[0].resize(l_read + read0_len, 0);
        }
        ra.read1[0].copy_within(0..read0_len, l_read);
        ra.read1[0].copy_within(read0_len + 1..read0_len + 1 + read1_len, 0);
        if o0 < ra.read_length[0] {
            let src_start = l_read + o0 as usize;
            let src_end = l_read + read0_len;
            let dst_start = read1_len;
            ra.read1[0].copy_within(src_start..src_end, dst_start);
        }
    }

    ra.l_read = ra.l_read - ra.pe_ov.n_ov - 1;
    if ra.read_length.len() < 2 {
        ra.read_length.resize(2, 0);
    }
    if ra.read_length_original.len() < 2 {
        ra.read_length_original.resize(2, 0);
    }
    ra.read_length[0] = ra.l_read;
    ra.read_length[1] = 0;
    ra.read_length_original[0] = ra.l_read;
    ra.read_length_original[1] = 0;
    ra.read_nmates = 1;

    for mate in 1..=2 {
        if ra.read1[mate].len() < ra.l_read as usize {
            ra.read1[mate].resize(ra.l_read as usize, 0);
        }
    }
    let read0_prefix = ra.read1[0][..ra.l_read as usize].to_vec();
    sequencefuns_l4_complementseqnumbers(&read0_prefix, &mut ra.read1[1], ra.l_read);
    for ii in 0..ra.l_read as usize {
        ra.read1[2][ra.l_read as usize - ii - 1] = ra.read1[1][ii];
    }
}

#[doc = "Original `Transcript::peOverlapSEtoPE` at STAR/source/ReadAlign_peOverlapMergeMap.cpp:136. Args: mateStart: uint, t: Transcript"]
pub fn readalign_peoverlapmergemap_l136_transcript_peoverlapsetope(
    transcript: &mut crate::transcript::Transcript,
    mate_start: &[u32; 2],
    t: &crate::transcript::Transcript,
) {
    let m_len = [
        transcript.read_length[t.str_ as usize],
        transcript.read_length[1 - t.str_ as usize],
    ];
    let m_sta2 = [0_u32, m_len[0] + 1];
    let mut m_sta = [mate_start[0], mate_start[1]];
    if t.str_ == 1 {
        for ii in 0..2 {
            m_sta[ii] = t.l_read - transcript.read_length[ii] - m_sta[ii];
        }
        m_sta.swap(0, 1);
    }
    let m_end = [m_sta[0] + m_len[0], m_sta[1] + m_len[1]];

    transcript.n_exons = 0;
    for imate in 0..2 {
        for iex in 0..t.n_exons as usize {
            if t.exons[iex][EX_R] >= m_end[imate]
                || t.exons[iex][EX_R] + t.exons[iex][EX_L] <= m_sta[imate]
            {
                continue;
            }

            let out_exon = transcript.n_exons as usize;
            if transcript.exons.len() <= out_exon {
                transcript.exons.resize(out_exon + 1, [0; EX_SIZE]);
            }
            if transcript.canon_sj.len() <= out_exon {
                transcript.canon_sj.resize(out_exon + 1, 0);
            }
            if transcript.sj_annot.len() <= out_exon {
                transcript.sj_annot.resize(out_exon + 1, 0);
            }
            if transcript.sj_str.len() <= out_exon {
                transcript.sj_str.resize(out_exon + 1, 0);
            }
            if transcript.shift_sj.len() <= out_exon {
                transcript.shift_sj.resize(out_exon + 1, [0; 2]);
            }

            transcript.exons[out_exon][EX_IFRAG] = if imate == 0 { t.str_ } else { 1 - t.str_ };
            transcript.exons[out_exon][EX_SJA] = t.exons[iex][EX_SJA];
            if iex < t.n_exons as usize - 1 {
                transcript.canon_sj[out_exon] = t.canon_sj[iex];
                transcript.sj_annot[out_exon] = t.sj_annot[iex];
                transcript.sj_str[out_exon] = t.sj_str[iex];
                transcript.shift_sj[out_exon] = t.shift_sj[iex];
            }

            if t.exons[iex][EX_R] >= m_sta[imate] {
                transcript.exons[out_exon][EX_G] = t.exons[iex][EX_G];
                transcript.exons[out_exon][EX_L] = t.exons[iex][EX_L];
                transcript.exons[out_exon][EX_R] =
                    t.exons[iex][EX_R] - m_sta[imate] + m_sta2[imate];
            } else {
                transcript.exons[out_exon][EX_R] = m_sta2[imate];
                let delta = m_sta[imate] - t.exons[iex][EX_R];
                transcript.exons[out_exon][EX_L] = t.exons[iex][EX_L] - delta;
                transcript.exons[out_exon][EX_G] = t.exons[iex][EX_G] + delta;
            }

            if t.exons[iex][EX_R] + t.exons[iex][EX_L] > m_end[imate] {
                transcript.exons[out_exon][EX_L] -=
                    t.exons[iex][EX_R] + t.exons[iex][EX_L] - m_end[imate];
            }

            transcript.n_exons += 1;
            if transcript.n_exons > MAX_N_EXONS {
                transcript.max_score = 0;
                transcript.n_exons = 0;
                return;
            }
        }

        let last = transcript.n_exons as usize - 1;
        transcript.canon_sj[last] = -3;
        transcript.sj_annot[last] = 0;
        transcript.sj_str[last] = 0;
        transcript.shift_sj[last] = [0, 0];
    }

    transcript.intron_motifs = t.intron_motifs;
    transcript.sj_motif_strand = t.sj_motif_strand;
    transcript.chr = t.chr;
    transcript.str_ = t.str_;
    transcript.ro_str = t.ro_str;
    transcript.g_start = t.g_start;
    transcript.g_length = t.g_length;
    transcript.c_start = t.c_start;

    transcript.r_length = 0;
    for iex in 0..transcript.n_exons as usize {
        transcript.r_length += transcript.exons[iex][EX_L];
    }
    transcript.mapped_length = transcript.r_length;
    transcript.r_start = transcript.exons[0][EX_R];
    transcript.ro_start = if transcript.ro_str == 0 {
        transcript.r_start
    } else {
        transcript.l_read - transcript.r_start - transcript.r_length
    };

    transcript.n_gap = t.n_gap;
    transcript.l_gap = t.l_gap;
    transcript.n_del = t.n_del;
    transcript.n_ins = t.n_ins;
    transcript.l_del = t.n_del;
    transcript.l_ins = t.l_ins;
    transcript.n_unique = t.n_unique;
    transcript.n_anchor = t.n_anchor;
}

#[doc = "Original `ReadAlign::peOverlapSEtoPE` at STAR/source/ReadAlign_peOverlapMergeMap.cpp:266. Args: seRA: ReadAlign"]
pub fn readalign_peoverlapmergemap_l266_readalign_peoverlapsetope(
    ra: &mut crate::read_align::ReadAlign,
    se_ra: &crate::read_align::ReadAlign,
    tr_init: &crate::transcript::Transcript,
    g: &[u8],
    sjdb_score: i32,
    score_ins_base: i32,
    score_ins_open: i32,
    score_del_base: i32,
    score_del_open: i32,
    score_gap_noncan: i32,
    score_gap: i32,
    score_gap_gcag: i32,
    score_gap_atac: i32,
    score_genomic_length_log2scale: f64,
) {
    let mut best_score = -10 * ra.l_read as i32;
    ra.tr_best = tr_init.clone();
    ra.tr_all.clear();
    ra.n_win_tr.clear();

    for i_w in 0..se_ra.n_w as usize {
        let mut out_window = Vec::<crate::transcript::Transcript>::new();
        let n_win_tr = se_ra
            .n_win_tr
            .get(i_w)
            .copied()
            .unwrap_or_else(|| se_ra.tr_all.get(i_w).map_or(0, |w| w.len() as u64))
            as usize;

        if let Some(se_window) = se_ra.tr_all.get(i_w) {
            for se_tr in se_window.iter().take(n_win_tr) {
                let mut tr = tr_init.clone();
                readalign_peoverlapmergemap_l136_transcript_peoverlapsetope(
                    &mut tr,
                    &ra.pe_ov.mate_start,
                    se_tr,
                );
                if tr.n_exons == 0 {
                    continue;
                }

                transcript_alignscore_l4_transcript_alignscore(
                    &mut tr,
                    &ra.read1[0],
                    &ra.read1[2],
                    g,
                    sjdb_score,
                    score_ins_base,
                    score_ins_open,
                    score_del_base,
                    score_del_open,
                    score_gap_noncan,
                    score_gap,
                    score_gap_gcag,
                    score_gap_atac,
                    score_genomic_length_log2scale,
                );

                out_window.push(tr);
                let last = out_window.len() - 1;
                if out_window[last].max_score > out_window[0].max_score {
                    out_window.swap(last, 0);
                }
            }
        }

        if !out_window.is_empty() {
            ra.n_win_tr.push(out_window.len() as u64);
            if out_window[0].max_score > best_score {
                ra.tr_best = out_window[0].clone();
                best_score = ra.tr_best.max_score;
            }
            ra.tr_all.push(out_window);
        }
    }

    ra.n_w = ra.tr_all.len() as u32;
}

#[doc = "Original `ReadAlign::peOverlapChimericSEtoPE` at STAR/source/ReadAlign_peOverlapMergeMap.cpp:308. Args: seTrIn1: Transcript, seTrIn2: Transcript, peTrOut1: Transcript, peTrOut2: Transcript"]
pub fn readalign_peoverlapmergemap_l308_readalign_peoverlapchimericsetope(
    tr_init: &crate::transcript::Transcript,
    mate_start: &[u32; 2],
    read_length_original: &[u32],
    se_tr_in1: &crate::transcript::Transcript,
    se_tr_in2: &crate::transcript::Transcript,
    pe_tr_out1: &mut crate::transcript::Transcript,
    pe_tr_out2: &mut crate::transcript::Transcript,
) {
    let mut temp_tr_chim = [tr_init.clone(), tr_init.clone()];
    readalign_peoverlapmergemap_l136_transcript_peoverlapsetope(
        &mut temp_tr_chim[0],
        mate_start,
        se_tr_in1,
    );
    readalign_peoverlapmergemap_l136_transcript_peoverlapsetope(
        &mut temp_tr_chim[1],
        mate_start,
        se_tr_in2,
    );

    let mut seg_len = [[0u32; 2]; 2];
    let mut seg_ex = [0u32; 2];
    let mut i1 = 0usize;
    let mut i2 = 0usize;
    let mut pos_of_junction_in_read = 0u32;

    for ii in 0..2 {
        for iex in 0..temp_tr_chim[ii].n_exons as usize {
            if temp_tr_chim[ii].exons[iex][EX_IFRAG] == temp_tr_chim[ii].exons[0][EX_IFRAG] {
                seg_len[ii][0] += temp_tr_chim[ii].exons[iex][EX_L];
                seg_ex[ii] = iex as u32;
            } else {
                seg_len[ii][1] += temp_tr_chim[ii].exons[iex][EX_L];
            }
        }

        let read_len0 = read_length_original[temp_tr_chim[ii].exons[0][EX_IFRAG] as usize];
        let read_len1 = read_length_original[1 - temp_tr_chim[ii].exons[0][EX_IFRAG] as usize];
        for jj in 0..2 {
            let cur_pos = if temp_tr_chim[ii].exons[jj][EX_R] > read_len0 {
                read_len0 + read_len1 + 1 - temp_tr_chim[ii].exons[jj][EX_R]
            } else {
                temp_tr_chim[ii].exons[jj][EX_R]
            };
            if seg_len[ii][jj] < seg_len[i1][i2]
                || (seg_len[ii][jj] == seg_len[i1][i2] && cur_pos > pos_of_junction_in_read)
            {
                pos_of_junction_in_read = cur_pos;
                i1 = ii;
                i2 = jj;
            }
        }
    }

    if i2 == 1 {
        temp_tr_chim[i1].n_exons = seg_ex[i1] + 1;
        temp_tr_chim[i1]
            .exons
            .truncate(temp_tr_chim[i1].n_exons as usize);
    } else {
        let start = seg_ex[i1] as usize + 1;
        let old_n = temp_tr_chim[i1].n_exons as usize;
        let new_n = old_n.saturating_sub(start);
        for iex in 0..new_n {
            let iex1 = iex + start;
            temp_tr_chim[i1].exons[iex] = temp_tr_chim[i1].exons[iex1];
            if temp_tr_chim[i1].canon_sj.len() > iex1 {
                temp_tr_chim[i1].canon_sj[iex] = temp_tr_chim[i1].canon_sj[iex1];
            }
            if temp_tr_chim[i1].sj_annot.len() > iex1 {
                temp_tr_chim[i1].sj_annot[iex] = temp_tr_chim[i1].sj_annot[iex1];
            }
            if temp_tr_chim[i1].sj_str.len() > iex1 {
                temp_tr_chim[i1].sj_str[iex] = temp_tr_chim[i1].sj_str[iex1];
            }
            if temp_tr_chim[i1].shift_sj.len() > iex1 {
                temp_tr_chim[i1].shift_sj[iex] = temp_tr_chim[i1].shift_sj[iex1];
            }
        }
        temp_tr_chim[i1].n_exons = new_n as u32;
        temp_tr_chim[i1].exons.truncate(new_n);
    }

    *pe_tr_out1 = temp_tr_chim[0].clone();
    *pe_tr_out2 = temp_tr_chim[1].clone();
}
