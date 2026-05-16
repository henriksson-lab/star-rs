#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::chimericDetection` at STAR/source/ReadAlign_chimericDetection.cpp:16. Args: "]
pub fn readalign_chimericdetection_l16_readalign_chimericdetection(
    read_align: &mut crate::read_align::ReadAlign,
    p: &crate::parameters_chimeric::Parameters,
    map_gen: &crate::genome::Genome,
    detector_result: Option<bool>,
) -> Result<crate::quantifications::ChimericDetectionResult, String> {
    read_align.chim_record = false;

    if p.p_ch.segment_min == 0 {
        return Ok(crate::quantifications::ChimericDetectionResult {
            chim_record: read_align.chim_record,
            ..Default::default()
        });
    }
    if p.out_filter_by_sjout_stage > 1 {
        return Ok(crate::quantifications::ChimericDetectionResult {
            chim_record: read_align.chim_record,
            ..Default::default()
        });
    }

    let mut result = crate::quantifications::ChimericDetectionResult::default();
    if p.p_ch.multimap_nmax == 0 {
        result.request = Some(crate::quantifications::ChimericDetectionRequest {
            detector: "chimericDetectionOld".to_string(),
            n_w: read_align.n_w,
            read_length: read_align.read_length.clone(),
            max_non_chim_align_score: read_align.tr_best.max_score,
        });
        read_align.chim_record = if let Some(chim_record) = detector_result {
            chim_record
        } else {
            readalign_chimericdetectionold_l7_readalign_chimericdetectionold(read_align, p, map_gen)
        };
        result.old_output_requested = true;
        if read_align.tr_chim.len() < 2 {
            read_align
                .tr_chim
                .resize(2, crate::transcript::Transcript::default());
        }
        let mut tr_chim = [read_align.tr_chim[0].clone(), read_align.tr_chim[1].clone()];
        let old_output = if detector_result.is_none()
            && (!read_align.chim_record || (tr_chim[0].n_exons > 0 && tr_chim[1].n_exons > 0))
        {
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
            )?
        } else {
            crate::quantifications::ChimericDetectionOldOutputResult::default()
        };
        read_align.tr_chim[0] = tr_chim[0].clone();
        read_align.tr_chim[1] = tr_chim[1].clone();
        result.old_output = Some(old_output);
    } else {
        let read_len_sum: u32 = read_align.read_length.iter().take(2).copied().sum();
        if read_align.tr_best.max_score
            <= read_len_sum as i32 - p.p_ch.nonchim_score_drop_min as i32
        {
            result.request = Some(crate::quantifications::ChimericDetectionRequest {
                detector: "chimericDetectionMult".to_string(),
                n_w: read_align.n_w,
                read_length: read_align.read_length.clone(),
                max_non_chim_align_score: read_align.tr_best.max_score,
            });
            if let Some(chim_record) = detector_result {
                read_align.chim_record = chim_record;
            } else {
                let n_win_tr: Vec<u32> = read_align
                    .n_win_tr
                    .iter()
                    .map(|&n| n.min(u32::MAX as u64) as u32)
                    .collect();
                let mut chim_det = chimericdetection_l3_chimericdetection_chimericdetection(
                    p.clone(),
                    read_align.tr_all.clone(),
                    n_win_tr,
                    [read_align.read1[0].clone(), read_align.read1[1].clone()],
                    map_gen.clone(),
                    p.p_ch.out_junctions,
                    read_align.clone(),
                );
                chim_det.n_w = read_align.n_w;
                let mult_output =
                    chimericdetection_chimericdetectionmult_l23_chimericdetection_chimericdetectionmult(
                        &mut chim_det,
                        read_align.n_w,
                        &read_align.read_length,
                        read_align.tr_best.max_score,
                        None,
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
