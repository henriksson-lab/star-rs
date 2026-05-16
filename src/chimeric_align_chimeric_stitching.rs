#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ChimericAlign::chimericStitching` at STAR/source/ChimericAlign_chimericStitching.cpp:3. Args: genSeq: char, Read1: char"]
pub fn chimericalign_chimericstitching_l3_chimericalign_chimericstitching(
    chim: &mut crate::chimeric_align::ChimericAlign,
    gen_seq: &[u8],
    read1: [&[u8]; 2],
    p_ch: &crate::parameters_chimeric::ParametersChimeric,
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
    if chim.stitching_done {
        return;
    }

    chim.stitching_done = true;
    chim.chim_str = std::cmp::max(chim.seg1.str_ as i32, chim.seg2.str_ as i32);
    chim.chim_repeat1 = 0;
    chim.chim_repeat2 = 0;
    chim.chim_j1 = 0;
    chim.chim_j2 = 0;
    chim.chim_motif = 0;

    let ex1 = chim.ex1 as usize;
    let ex2 = chim.ex2 as usize;

    if chim.al1.exons[ex1][EX_IFRAG] < chim.al2.exons[ex2][EX_IFRAG] {
        chim.chim_motif = -1;
        chim.chim_j1 = if chim.al1.str_ == 1 {
            chim.al1.exons[ex1][EX_G] - 1
        } else {
            chim.al1.exons[ex1][EX_G] + chim.al1.exons[ex1][EX_L]
        };
        chim.chim_j2 = if chim.al2.str_ == 0 {
            chim.al2.exons[ex2][EX_G] - 1
        } else {
            chim.al2.exons[ex2][EX_G] + chim.al2.exons[ex2][EX_L]
        };
    } else {
        let ro_start0 = if chim.al1.str_ == 0 {
            chim.al1.exons[ex1][EX_R]
        } else {
            chim.al1.l_read - chim.al1.exons[ex1][EX_R] - chim.al1.exons[ex1][EX_L]
        };
        let ro_start1 = if chim.al2.str_ == 0 {
            chim.al2.exons[ex2][EX_R]
        } else {
            chim.al1.l_read - chim.al2.exons[ex2][EX_R] - chim.al2.exons[ex2][EX_L]
        };

        let mut j_rbest = 0_u64;
        let mut j_score = 0_i32;
        let mut j_score_best = -999999_i32;
        let j_rmax0 = ro_start1 + chim.al2.exons[ex2][EX_L];
        let j_rmax = if j_rmax0 > ro_start0 {
            j_rmax0 - ro_start0 - 1
        } else {
            0
        };

        let mut j_r = 0_u64;
        while j_r < j_rmax {
            if chim.al1.read_length.first().copied().unwrap_or(0) == j_r {
                j_r += 1;
            }

            let b_r = read1[0][(ro_start0 + j_r) as usize];
            let mut b0;
            let mut b1;
            if chim.al1.str_ == 0 {
                b0 = gen_seq[(chim.al1.exons[ex1][EX_G] + j_r) as usize];
            } else {
                b0 = gen_seq
                    [(chim.al1.exons[ex1][EX_G] + chim.al1.exons[ex1][EX_L] - 1 - j_r) as usize];
                if b0 < 4 {
                    b0 = 3 - b0;
                }
            }

            if chim.al2.str_ == 0 {
                b1 = gen_seq[(chim.al2.exons[ex2][EX_G] - ro_start1 + ro_start0 + j_r) as usize];
            } else {
                b1 = gen_seq[(chim.al2.exons[ex2][EX_G] + chim.al2.exons[ex2][EX_L] - 1 + ro_start1
                    - ro_start0
                    - j_r) as usize];
                if b1 < 4 {
                    b1 = 3 - b1;
                }
            }

            if (p_ch.filter_genomic_n && (b0 > 3 || b1 > 3)) || b_r > 3 {
                chim.chim_score = 0;
                return;
            }

            let mut b01;
            let mut b02;
            let mut b11;
            let mut b12;
            if chim.al1.str_ == 0 {
                b01 = gen_seq[(chim.al1.exons[ex1][EX_G] + j_r + 1) as usize];
                b02 = gen_seq[(chim.al1.exons[ex1][EX_G] + j_r + 2) as usize];
            } else {
                b01 = gen_seq[(chim.al1.exons[ex1][EX_G] + chim.al1.exons[ex1][EX_L] - 1 - j_r - 1)
                    as usize];
                if b01 < 4 {
                    b01 = 3 - b01;
                }
                b02 = gen_seq[(chim.al1.exons[ex1][EX_G] + chim.al1.exons[ex1][EX_L] - 1 - j_r - 2)
                    as usize];
                if b02 < 4 {
                    b02 = 3 - b02;
                }
            }
            if chim.al2.str_ == 0 {
                b11 =
                    gen_seq[(chim.al2.exons[ex2][EX_G] - ro_start1 + ro_start0 + j_r - 1) as usize];
                b12 = gen_seq[(chim.al2.exons[ex2][EX_G] - ro_start1 + ro_start0 + j_r) as usize];
            } else {
                b11 = gen_seq[(chim.al2.exons[ex2][EX_G] + chim.al2.exons[ex2][EX_L] - 1
                    + ro_start1
                    - ro_start0
                    - j_r
                    + 1) as usize];
                if b11 < 4 {
                    b11 = 3 - b11;
                }
                b12 = gen_seq[(chim.al2.exons[ex2][EX_G] + chim.al2.exons[ex2][EX_L] - 1
                    + ro_start1
                    - ro_start0
                    - j_r) as usize];
                if b12 < 4 {
                    b12 = 3 - b12;
                }
            }

            let mut j_motif = 0_i32;
            if b01 == 2 && b02 == 3 && b11 == 0 && b12 == 2 {
                if chim.chim_str != 2 {
                    j_motif = 1;
                }
            } else if b01 == 1 && b02 == 3 && b11 == 0 && b12 == 1 && chim.chim_str != 1 {
                j_motif = 2;
            }

            if b_r == b0 && b_r != b1 {
                j_score += 1;
            } else if b_r != b0 && b_r == b1 {
                j_score -= 1;
            }

            let j_score_j = if j_motif == 0 {
                j_score + p_ch.score_junction_non_gtag
            } else {
                j_score
            };
            if j_score_j > j_score_best || (j_score_j == j_score_best && j_motif > 0) {
                chim.chim_motif = j_motif;
                j_rbest = j_r;
                j_score_best = j_score_j;
            }

            j_r += 1;
        }

        if chim.al1.str_ == 1 {
            chim.al1.exons[ex1][EX_R] += chim.al1.exons[ex1][EX_L] - j_rbest - 1;
            chim.al1.exons[ex1][EX_G] += chim.al1.exons[ex1][EX_L] - j_rbest - 1;
            chim.al1.exons[ex1][EX_L] = j_rbest + 1;
            chim.chim_j1 = chim.al1.exons[ex1][EX_G] - 1;
        } else {
            chim.al1.exons[ex1][EX_L] = j_rbest + 1;
            chim.chim_j1 = chim.al1.exons[ex1][EX_G] + chim.al1.exons[ex1][EX_L];
        }

        if chim.al2.str_ == 0 {
            chim.al2.exons[ex2][EX_R] += ro_start0 + j_rbest + 1 - ro_start1;
            chim.al2.exons[ex2][EX_G] += ro_start0 + j_rbest + 1 - ro_start1;
            chim.al2.exons[ex2][EX_L] =
                ro_start1 + chim.al2.exons[ex2][EX_L] - ro_start0 - j_rbest - 1;
            chim.chim_j2 = chim.al2.exons[ex2][EX_G] - 1;
        } else {
            chim.al2.exons[ex2][EX_L] =
                ro_start1 + chim.al2.exons[ex2][EX_L] - ro_start0 - j_rbest - 1;
            chim.chim_j2 = chim.al2.exons[ex2][EX_G] + chim.al2.exons[ex2][EX_L];
        }

        let mut repeat = 0_u64;
        while repeat < 100 {
            let mut b0 = if chim.al1.str_ == 0 {
                gen_seq[(chim.chim_j1 as u64 + repeat) as usize]
            } else {
                gen_seq[(chim.chim_j1 as u64 - repeat) as usize]
            };
            if chim.al1.str_ != 0 && b0 < 4 {
                b0 = 3 - b0;
            }
            let mut b1 = if chim.al2.str_ == 0 {
                gen_seq[(chim.chim_j2 as u64 + 1 + repeat) as usize]
            } else {
                gen_seq[(chim.chim_j2 as u64 - 1 - repeat) as usize]
            };
            if chim.al2.str_ != 0 && b1 < 4 {
                b1 = 3 - b1;
            }
            if b0 != b1 {
                break;
            }
            repeat += 1;
        }
        chim.chim_repeat2 = repeat;

        repeat = 0;
        while repeat < 100 {
            let mut b0 = if chim.al1.str_ == 0 {
                gen_seq[(chim.chim_j1 as u64 - 1 - repeat) as usize]
            } else {
                gen_seq[(chim.chim_j1 as u64 + 1 + repeat) as usize]
            };
            if chim.al1.str_ != 0 && b0 < 4 {
                b0 = 3 - b0;
            }
            let mut b1 = if chim.al2.str_ == 0 {
                gen_seq[(chim.chim_j2 as u64 - repeat) as usize]
            } else {
                gen_seq[(chim.chim_j2 as u64 + repeat) as usize]
            };
            if chim.al2.str_ != 0 && b1 < 4 {
                b1 = 3 - b1;
            }
            if b0 != b1 {
                break;
            }
            repeat += 1;
        }
        chim.chim_repeat1 = repeat;
    }

    if chim.chim_motif >= 0
        && (chim.al1.exons[ex1][EX_L] < p_ch.junction_overhang_min
            || chim.al2.exons[ex2][EX_L] < p_ch.junction_overhang_min)
    {
        chim.chim_score = 0;
        return;
    }

    let score1 = transcript_alignscore_l4_transcript_alignscore(
        &mut chim.al1,
        read1[0],
        read1[1],
        gen_seq,
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
    let score2 = transcript_alignscore_l4_transcript_alignscore(
        &mut chim.al2,
        read1[0],
        read1[1],
        gen_seq,
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
    chim.chim_score = score1
        + score2
        + if chim.chim_motif == 0 {
            p_ch.score_junction_non_gtag
        } else {
            0
        };
}
