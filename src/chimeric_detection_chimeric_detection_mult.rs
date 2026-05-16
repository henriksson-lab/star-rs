#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `chimericAlignScore` at STAR/source/ChimericDetection_chimericDetectionMult.cpp:6. Args: seg1: ChimericSegment, seg2: ChimericSegment"]
pub fn chimericdetection_chimericdetectionmult_l6_chimericalignscore(
    seg1: &crate::chimeric_segment::ChimericSegment,
    seg2: &crate::chimeric_segment::ChimericSegment,
) -> i32 {
    let mut chim_score = 0;
    let chim_overlap = if seg2.ro_s > seg1.ro_s {
        if seg2.ro_s > seg1.ro_e {
            0
        } else {
            seg1.ro_e - seg2.ro_s + 1
        }
    } else if seg2.ro_e < seg1.ro_s {
        0
    } else {
        seg2.ro_e - seg1.ro_s + 1
    };
    let diff_mates = (seg1.ro_e < seg1.align.read_length[0]
        && seg2.ro_s >= seg1.align.read_length[0])
        || (seg2.ro_e < seg1.align.read_length[0] && seg1.ro_s >= seg1.align.read_length[0]);

    if seg1.ro_e > seg1.segment_min + seg1.ro_s + chim_overlap
        && seg2.ro_e > seg1.segment_min + seg2.ro_s + chim_overlap
        && (diff_mates
            || ((seg1.ro_e + seg1.segment_read_gap_max + 1) >= seg2.ro_s
                && (seg2.ro_e + seg1.segment_read_gap_max + 1) >= seg1.ro_s))
    {
        chim_score = seg1.align.max_score + seg2.align.max_score - chim_overlap as i32;
    }

    chim_score
}

#[doc = "Original `ChimericDetection::chimericDetectionMult` at STAR/source/ChimericDetection_chimericDetectionMult.cpp:23. Args: nW: uint, readLength: uint, maxNonChimAlignScore: int, PEunmergedRA: ReadAlign"]
pub fn chimericdetection_chimericdetectionmult_l23_chimericdetection_chimericdetectionmult(
    chim_det: &mut crate::chimeric_detection::ChimericDetection,
    n_w: u32,
    read_length: &[u32],
    max_non_chim_align_score: i32,
    pe_unmerged_ra: Option<&crate::read_align::ReadAlign>,
    pe_tr_init: Option<&crate::transcript::Transcript>,
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
) -> crate::quantifications::ChimericDetectionMultResult {
    chim_det.chim_aligns.clear();

    let mut result = crate::quantifications::ChimericDetectionMultResult::default();
    let p_ch = &chim_det.p.p_ch;
    let mut chim_score_best = 0_i32;
    let mut best_chim_align = 0usize;
    let max_possible_align_score = read_length.first().copied().unwrap_or(0) as i32
        + read_length.get(1).copied().unwrap_or(0) as i32;
    let mut min_score_to_consider = p_ch.score_min;
    if max_non_chim_align_score >= min_score_to_consider {
        min_score_to_consider = max_non_chim_align_score + 1;
    }
    if max_possible_align_score - p_ch.score_drop_max > min_score_to_consider {
        min_score_to_consider = max_possible_align_score - p_ch.score_drop_max;
    }

    for i_w1 in 0..n_w as usize {
        let n_win1 = chim_det
            .n_win_tr
            .get(i_w1)
            .copied()
            .unwrap_or_else(|| chim_det.tr_all.get(i_w1).map_or(0, |w| w.len() as u32))
            as usize;
        for i_a1 in 0..n_win1 {
            let Some(align1) = chim_det.tr_all.get(i_w1).and_then(|w| w.get(i_a1)) else {
                continue;
            };
            let seg1 = chimericsegment_l3_chimericsegment_chimericsegment(p_ch, align1);
            if !chimericsegment_l19_chimericsegment_segmentcheck(&seg1) {
                continue;
            }

            for i_w2 in i_w1..n_w as usize {
                let n_win2 = chim_det
                    .n_win_tr
                    .get(i_w2)
                    .copied()
                    .unwrap_or_else(|| chim_det.tr_all.get(i_w2).map_or(0, |w| w.len() as u32))
                    as usize;
                let i_a2_start = if i_w1 != i_w2 { 0 } else { i_a1 + 1 };
                for i_a2 in i_a2_start..n_win2 {
                    let Some(align2) = chim_det.tr_all.get(i_w2).and_then(|w| w.get(i_a2)) else {
                        continue;
                    };
                    let seg2 = chimericsegment_l3_chimericsegment_chimericsegment(p_ch, align2);
                    if !chimericsegment_l19_chimericsegment_segmentcheck(&seg2) {
                        continue;
                    }
                    if seg1.str_ != 0 && seg2.str_ != 0 && seg2.str_ != seg1.str_ {
                        continue;
                    }

                    let chim_score =
                        chimericdetection_chimericdetectionmult_l6_chimericalignscore(&seg1, &seg2);
                    if chim_score >= min_score_to_consider {
                        let mut ch_al = chimericalign_l3_chimericalign_chimericalign(
                            seg1.clone(),
                            seg2,
                            chim_score,
                            p_ch.junction_overhang_min,
                        );
                        if !chimericalign_l17_chimericalign_chimericcheck(&ch_al) {
                            continue;
                        }
                        chimericalign_chimericstitching_l3_chimericalign_chimericstitching(
                            &mut ch_al,
                            &chim_det.out_gen.g,
                            [&chim_det.read1[0], &chim_det.read1[1]],
                            p_ch,
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
                        if ch_al.chim_score >= min_score_to_consider {
                            chim_det.chim_aligns.push(ch_al);
                            let last = chim_det.chim_aligns.len() - 1;
                            if chim_det.chim_aligns[last].chim_score > chim_score_best {
                                chim_score_best = chim_det.chim_aligns[last].chim_score;
                                best_chim_align = last;
                                let score_range_floor =
                                    chim_score_best - p_ch.multimap_score_range as i32;
                                if score_range_floor > min_score_to_consider {
                                    min_score_to_consider = score_range_floor;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    result.chim_score_best = chim_score_best;
    result.min_score_to_consider = min_score_to_consider;
    if chim_score_best == 0 {
        return result;
    }

    let mut chim_n = 0_u32;
    for ch_al in &chim_det.chim_aligns {
        if ch_al.chim_score >= min_score_to_consider {
            chim_n += 1;
        }
    }
    result.chim_n = chim_n;
    if chim_n > p_ch.multimap_nmax {
        return result;
    }

    let mut i_tr = 0_u32;
    for i in 0..chim_det.chim_aligns.len() {
        if chim_det.chim_aligns[i].chim_score >= min_score_to_consider {
            if p_ch.out_junctions {
                if let Some(ra) = chim_det.ra.as_ref().or(pe_unmerged_ra) {
                    result.chim_junction.push_str(
                        &chimericalign_chimericjunctionoutput_l4_chimericalign_chimericjunctionoutput(
                            &chim_det.chim_aligns[i],
                            &chim_det.out_gen,
                            &chim_det.p,
                            &ra.read_name,
                            ra.read_files_index,
                            chim_n,
                            max_non_chim_align_score,
                            pe_unmerged_ra.is_some(),
                            chim_score_best,
                            max_possible_align_score,
                            None,
                        ),
                    );
                }
            }

            if p_ch.out_bam {
                let mut al1 = chim_det.chim_aligns[i].al1.clone();
                let mut al2 = chim_det.chim_aligns[i].al2.clone();
                let ra_for_bam = if let Some(pe_ra) = pe_unmerged_ra {
                    let tr_init = pe_tr_init.unwrap_or(&pe_ra.tr_best);
                    readalign_peoverlapmergemap_l308_readalign_peoverlapchimericsetope(
                        tr_init,
                        &pe_ra.pe_ov.mate_start,
                        &pe_ra.read_length_original,
                        &chim_det.chim_aligns[i].al1,
                        &chim_det.chim_aligns[i].al2,
                        &mut al1,
                        &mut al2,
                    );
                    pe_ra
                } else if let Some(ra) = chim_det.ra.as_ref() {
                    ra
                } else {
                    i_tr += 1;
                    continue;
                };
                result.bam_outputs.push(
                    chimericalign_chimericbamoutput_l7_chimericalign_chimericbamoutput(
                        &al1,
                        &al2,
                        ra_for_bam,
                        &chim_det.out_gen,
                        i_tr,
                        chim_n,
                        i == best_chim_align,
                        &chim_det.p,
                    ),
                );
            }
            i_tr += 1;
        }
    }

    result.chim_record = chim_n > 0;
    result
}
