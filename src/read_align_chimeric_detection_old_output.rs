#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::chimericDetectionOldOutput` at STAR/source/ReadAlign_chimericDetectionOldOutput.cpp:5. Args: "]
pub fn readalign_chimericdetectionoldoutput_l5_readalign_chimericdetectionoldoutput(
    chim_record: bool,
    tr_chim: &mut [crate::transcript::Transcript; 2],
    read_align: &crate::read_align::ReadAlign,
    p: &crate::parameters_chimeric::Parameters,
    map_gen: &crate::genome::Genome,
    read0: &[String],
    read_file_type: i32,
    qual0: &[String],
    read_name_extra: &[String],
    read_files_in_n: usize,
    chim_j0: u32,
    chim_j1: u32,
    chim_motif: i32,
    chim_repeat0: u32,
    chim_repeat1: u32,
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
) -> Result<crate::quantifications::ChimericDetectionOldOutputResult, String> {
    let mut result = crate::quantifications::ChimericDetectionOldOutputResult::default();
    if !chim_record {
        return Ok(result);
    }

    transcript_alignscore_l4_transcript_alignscore(
        &mut tr_chim[0],
        &read_align.read1[0],
        &read_align.read1[2],
        &map_gen.g,
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
    transcript_alignscore_l4_transcript_alignscore(
        &mut tr_chim[1],
        &read_align.read1[0],
        &read_align.read1[2],
        &map_gen.g,
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

    if p.p_ch.out_bam {
        result
            .bam_requests
            .push(crate::quantifications::ChimericBamOutputRequest {
                al1: tr_chim[0].clone(),
                al2: tr_chim[1].clone(),
                i_tr: 0,
                chim_n: 1,
                is_best_chim_align: true,
            });
    }

    if p.p_ch.out_sam_old {
        result.chim_n = 2;
        let tr0_paired = tr_chim[0].exons[0][EX_IFRAG]
            != tr_chim[0].exons[tr_chim[0].n_exons as usize - 1][EX_IFRAG];
        let tr1_paired = tr_chim[1].exons[0][EX_IFRAG]
            != tr_chim[1].exons[tr_chim[1].n_exons as usize - 1][EX_IFRAG];
        if tr0_paired {
            tr_chim[0].primary_flag = true;
            tr_chim[1].primary_flag = false;
        } else if tr1_paired {
            tr_chim[1].primary_flag = true;
            tr_chim[0].primary_flag = false;
        } else if tr_chim[0].exons[0][EX_IFRAG] != tr_chim[1].exons[0][EX_IFRAG] {
            tr_chim[0].primary_flag = true;
            tr_chim[1].primary_flag = true;
        } else {
            let chim_represent = if tr_chim[0].max_score > tr_chim[1].max_score {
                0
            } else {
                1
            };
            tr_chim[chim_represent].primary_flag = true;
            tr_chim[1 - chim_represent].primary_flag = false;
        }

        for i_tr in 0..result.chim_n as usize {
            let (mate_chr, mate_start, mate_strand) = if p.read_nmates == 2 {
                let other = 1 - i_tr;
                let mut iex = 0usize;
                if tr_chim[other].exons[0][EX_IFRAG]
                    != tr_chim[other].exons[tr_chim[other].n_exons as usize - 1][EX_IFRAG]
                {
                    while iex < tr_chim[other].n_exons as usize {
                        if tr_chim[other].exons[iex][EX_IFRAG] != tr_chim[i_tr].exons[0][EX_IFRAG] {
                            break;
                        }
                        iex += 1;
                    }
                }

                (
                    tr_chim[other].chr,
                    tr_chim[other].exons[iex][EX_G],
                    if tr_chim[other].str_ != tr_chim[other].exons[iex][EX_IFRAG] {
                        1
                    } else {
                        0
                    },
                )
            } else {
                (u32::MAX, u32::MAX, -1)
            };

            readalign_outputtranscriptsam_l5_readalign_outputtranscriptsam(
                &tr_chim[i_tr],
                result.chim_n,
                i_tr as u32,
                mate_chr,
                mate_start,
                mate_strand,
                -1,
                None,
                &mut result.chim_sam,
                read_align.read_filter as u8,
                &read_align.read_name,
                read0,
                read_file_type,
                qual0,
                p,
                read_align.read_files_index,
                read_name_extra,
                read_align.l_read,
                &read_align.read_length,
                &read_align.read_length_original,
                &read_align.clip_mates,
                &read_align.read1,
                map_gen,
            )?;
        }
    }

    if p.p_ch.out_junctions {
        result.chim_junction.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            map_gen.chr_name[tr_chim[0].chr as usize],
            chim_j0 as u64 - map_gen.chr_start[tr_chim[0].chr as usize] + 1,
            if tr_chim[0].str_ == 0 { "+" } else { "-" },
            map_gen.chr_name[tr_chim[1].chr as usize],
            chim_j1 as u64 - map_gen.chr_start[tr_chim[1].chr as usize] + 1,
            if tr_chim[1].str_ == 0 { "+" } else { "-" },
            chim_motif,
            chim_repeat0,
            chim_repeat1,
            read_align.read_name.get(1..).unwrap_or(""),
            tr_chim[0].exons[0][EX_G] as u64 - map_gen.chr_start[tr_chim[0].chr as usize] + 1,
            readalign_outputtranscriptcigarp_l4_readalign_outputtranscriptcigarp(
                read_align,
                &tr_chim[0],
                read_files_in_n,
            ),
            tr_chim[1].exons[0][EX_G] as u64 - map_gen.chr_start[tr_chim[1].chr as usize] + 1,
            readalign_outputtranscriptcigarp_l4_readalign_outputtranscriptcigarp(
                read_align,
                &tr_chim[1],
                read_files_in_n,
            )
        ));
        if p.out_sam_attr_present.rg {
            result.chim_junction.push_str(&format!(
                "\t{}",
                p.out_sam_attr_rg[read_align.read_files_index as usize]
            ));
        }
        result.chim_junction.push('\n');
    }

    Ok(result)
}
