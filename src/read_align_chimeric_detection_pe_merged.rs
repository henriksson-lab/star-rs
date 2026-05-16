#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::chimericDetectionPEmerged` at STAR/source/ReadAlign_chimericDetectionPEmerged.cpp:5. Args: seRA: ReadAlign"]
pub fn readalign_chimericdetectionpemerged_l5_readalign_chimericdetectionpemerged(
    read_align: &mut crate::read_align::ReadAlign,
    se_ra: &mut crate::read_align::ReadAlign,
    p: &crate::parameters_chimeric::Parameters,
    map_gen: &crate::genome::Genome,
    detector_result: Option<bool>,
) -> Result<crate::quantifications::ChimericDetectionResult, String> {
    read_align.chim_record = false;
    if p.p_ch.segment_min == 0 {
        return Ok(crate::quantifications::ChimericDetectionResult::default());
    }

    let mut result = crate::quantifications::ChimericDetectionResult::default();
    if p.p_ch.multimap_nmax == 0 {
        result.mult_map_select_requested = true;
        result.mapped_filter_requested = true;
        if detector_result.is_none() {
            let rng_snapshot = se_ra.rng_uniform_real_0_to_1;
            let _tr_mult = readalign_multmapselect_inner(
                &mut se_ra.n_tr,
                &mut se_ra.tr_best,
                se_ra.l_read,
                se_ra.n_w,
                &se_ra.n_win_tr,
                map_gen,
                &se_ra.tr_all,
                p.out_filter_multimap_score_range,
                p.out_filter_multimap_nmax,
                p.out_multimapper_order_random,
                p.out_sam_mult_nmax_is_limited,
                &p.out_sam_primary_flag,
                &rng_snapshot,
            )?;
            readalign_mappedfilter_l3_readalign_mappedfilter(
                se_ra,
                p.out_filter_score_min,
                p.out_filter_score_min_over_lread,
                p.out_filter_match_nmin,
                p.out_filter_match_nmin_over_lread,
                p.out_filter_mismatch_nover_lmax,
                p.out_filter_multimap_nmax,
            );
        }
        result.request = Some(crate::quantifications::ChimericDetectionRequest {
            detector: "chimericDetectionOld".to_string(),
            n_w: se_ra.n_w as u32,
            read_length: se_ra.read_length.iter().map(|&v| v as u32).collect(),
            max_non_chim_align_score: se_ra.tr_best.max_score,
        });
        read_align.chim_record = if let Some(chim_record) = detector_result {
            chim_record
        } else {
            readalign_chimericdetectionold_l7_readalign_chimericdetectionold(se_ra, p, map_gen)
        };
        if !read_align.chim_record {
            result.chim_record = false;
            return Ok(result);
        }

        if se_ra.tr_chim.len() >= 2 {
            let mut pe_tr_out1 = crate::transcript::Transcript::default();
            let mut pe_tr_out2 = crate::transcript::Transcript::default();
            readalign_peoverlapmergemap_l308_readalign_peoverlapchimericsetope(
                &read_align.tr_best,
                &read_align.pe_ov.mate_start,
                &read_align.read_length_original,
                &se_ra.tr_chim[0],
                &se_ra.tr_chim[1],
                &mut pe_tr_out1,
                &mut pe_tr_out2,
            );
            read_align.tr_chim = vec![pe_tr_out1, pe_tr_out2];
            result.pe_tr_chim = read_align.tr_chim.clone();
        }
        result.old_output_requested = true;
        if detector_result.is_none()
            && read_align.tr_chim.len() >= 2
            && read_align.tr_chim[0].n_exons > 0
            && read_align.tr_chim[1].n_exons > 0
        {
            let mut tr_chim = [read_align.tr_chim[0].clone(), read_align.tr_chim[1].clone()];
            let read0_strings: Vec<String> = read_align
                .read0
                .iter()
                .map(|read| String::from_utf8_lossy(read).into_owned())
                .collect();
            let qual0_strings: Vec<String> = read_align
                .qual0
                .iter()
                .map(|qual| String::from_utf8_lossy(qual).into_owned())
                .collect();
            let old_output =
                readalign_chimericdetectionoldoutput_l5_readalign_chimericdetectionoldoutput(
                    read_align.chim_record,
                    &mut tr_chim,
                    read_align,
                    p,
                    map_gen,
                    &read0_strings,
                    read_align.read_file_type,
                    &qual0_strings,
                    &read_align.read_name_extra,
                    p.read_files_in.len(),
                    read_align.chim_j0,
                    read_align.chim_j1,
                    read_align.chim_motif,
                    read_align.chim_repeat0,
                    read_align.chim_repeat1,
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
                )?;
            read_align.tr_chim[0] = tr_chim[0].clone();
            read_align.tr_chim[1] = tr_chim[1].clone();
            result.old_output = Some(old_output);
        }
    } else {
        let read_len_sum: u64 = read_align.read_length.iter().take(2).copied().sum();
        if read_align.tr_best.max_score
            <= read_len_sum as i32 - p.p_ch.nonchim_score_drop_min as i32
        {
            result.request = Some(crate::quantifications::ChimericDetectionRequest {
                detector: "chimericDetectionMult".to_string(),
                n_w: se_ra.n_w as u32,
                read_length: se_ra.read_length.iter().map(|&v| v as u32).collect(),
                max_non_chim_align_score: se_ra.tr_best.max_score,
            });
            if let Some(chim_record) = detector_result {
                read_align.chim_record = chim_record;
            } else {
                let n_win_tr: Vec<u32> = se_ra
                    .n_win_tr
                    .iter()
                    .map(|&n| n.min(u32::MAX as u64) as u32)
                    .collect();
                let mut chim_det = chimericdetection_l3_chimericdetection_chimericdetection(
                    p.clone(),
                    se_ra.tr_all.clone(),
                    n_win_tr,
                    [se_ra.read1[0].clone(), se_ra.read1[1].clone()],
                    map_gen.clone(),
                    p.p_ch.out_junctions,
                    se_ra.clone(),
                );
                chim_det.n_w = se_ra.n_w as u32;
                let pe_ra_snapshot = read_align.clone();
                let mult_output =
                    chimericdetection_chimericdetectionmult_l23_chimericdetection_chimericdetectionmult(
                        &mut chim_det,
                        se_ra.n_w,
                        &se_ra.read_length,
                        se_ra.tr_best.max_score,
                        Some(&pe_ra_snapshot),
                        None,
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
                    );
                read_align.chim_record = mult_output.chim_record;
                result.mult_output = Some(mult_output);
            }
        }
    }

    if read_align.chim_record {
        read_align.stats_ra.chimeric_all += 1;
    }
    result.chim_record = read_align.chim_record;
    Ok(result)
}
