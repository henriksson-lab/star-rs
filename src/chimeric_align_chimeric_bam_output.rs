#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ChimericAlign::chimericBAMoutput` at STAR/source/ChimericAlign_chimericBAMoutput.cpp:7. Args: al1: Transcript, al2: Transcript, RA: ReadAlign, iTr: uint, chimN: uint, isBestChimAlign: bool, P: Parameters"]
pub fn chimericalign_chimericbamoutput_l7_chimericalign_chimericbamoutput(
    al1: &crate::transcript::Transcript,
    al2: &crate::transcript::Transcript,
    read_align: &crate::read_align::ReadAlign,
    map_gen: &crate::genome::Genome,
    i_tr: u32,
    chim_n: u32,
    is_best_chim_align: bool,
    p: &crate::parameters_chimeric::Parameters,
) -> crate::quantifications::ChimericBamOutputResult {
    let mut tr_chim = [al1.clone(), al2.clone()];
    let mut result = crate::quantifications::ChimericBamOutputResult {
        chim_represent: -999,
        ..Default::default()
    };

    if tr_chim[0].exons[0][EX_IFRAG] != tr_chim[0].exons[tr_chim[0].n_exons as usize - 1][EX_IFRAG]
    {
        result.chim_represent = 0;
        result.chim_type = 1;
    } else if tr_chim[1].exons[0][EX_IFRAG]
        != tr_chim[1].exons[tr_chim[1].n_exons as usize - 1][EX_IFRAG]
    {
        result.chim_represent = 1;
        result.chim_type = 1;
    } else if tr_chim[0].exons[0][EX_IFRAG] != tr_chim[1].exons[0][EX_IFRAG] {
        result.chim_represent = -1;
        result.chim_type = 2;
    } else {
        result.chim_represent = if tr_chim[0].max_score > tr_chim[1].max_score {
            0
        } else {
            1
        };
        result.chim_type = 3;
    }

    result.representative_request_index = -1;
    result.supplementary_request_index = -1;

    for itr in 0..2usize {
        tr_chim[itr].primary_flag = is_best_chim_align;
        let mate_chr: u32;
        let mate_start_abs: u32;
        let mate_strand: i8;
        let align_type: i32;

        if result.chim_type == 2 {
            mate_chr = tr_chim[1 - itr].chr;
            mate_start_abs = tr_chim[1 - itr].exons[0][EX_G];
            mate_strand = i8::from(tr_chim[1 - itr].str_ != tr_chim[1 - itr].exons[0][EX_IFRAG]);
            align_type = -10;
        } else {
            mate_chr = u32::MAX;
            mate_start_abs = u32::MAX;
            mate_strand = 0;
            if result.chim_represent == itr as i32 {
                align_type = -10;
                result.representative_request_index = result.bam_requests.len() as i32;
            } else {
                align_type = if p.p_ch.out_bam_hard_clip {
                    if (itr as u32) == tr_chim[itr].str_ {
                        -12
                    } else {
                        -11
                    }
                } else {
                    -13
                };
                result.supplementary_request_index = result.bam_requests.len() as i32;
                if result.chim_type == 1 {
                    let chim_represent = result.chim_represent as usize;
                    let mut iex = 0usize;
                    while iex < tr_chim[chim_represent].n_exons.saturating_sub(1) as usize {
                        if tr_chim[chim_represent].exons[iex][EX_IFRAG]
                            != tr_chim[itr].exons[0][EX_IFRAG]
                        {
                            break;
                        }
                        iex += 1;
                    }
                    let mate_chr1 = tr_chim[chim_represent].chr;
                    let mate_start1 = tr_chim[chim_represent].exons[iex][EX_G];
                    let mate_strand1 = i8::from(
                        tr_chim[chim_represent].str_
                            != tr_chim[chim_represent].exons[iex][EX_IFRAG],
                    );
                    let tr_chr_start = map_gen.chr_start[tr_chim[itr].chr as usize];
                    let mate_start = mate_start1.wrapping_sub(
                        map_gen.chr_start[mate_chr1.min(map_gen.n_chr_real - 1) as usize] as u32,
                    );
                    result
                        .bam_requests
                        .push(crate::quantifications::AlignBamRequest {
                            transcript: tr_chim[itr].clone(),
                            n_tr_out: chim_n,
                            i_tr_out: i_tr,
                            tr_chr_start,
                            mate_chr: mate_chr1,
                            mate_start,
                            mate_strand: mate_strand1,
                            align_type,
                            mate_map: None,
                        });
                    continue;
                }
            }
        }

        let tr_chr_start = map_gen.chr_start[tr_chim[itr].chr as usize];
        let mate_chr_index = if mate_chr < map_gen.n_chr_real {
            mate_chr
        } else {
            0
        };
        let mate_start =
            mate_start_abs.wrapping_sub(map_gen.chr_start[mate_chr_index as usize] as u32);
        result
            .bam_requests
            .push(crate::quantifications::AlignBamRequest {
                transcript: tr_chim[itr].clone(),
                n_tr_out: chim_n,
                i_tr_out: i_tr,
                tr_chr_start,
                mate_chr,
                mate_start,
                mate_strand,
                align_type,
                mate_map: None,
            });
    }

    let _ = read_align;
    result
}
