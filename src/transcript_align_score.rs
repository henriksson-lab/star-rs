#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Transcript::alignScore` at STAR/source/Transcript_alignScore.cpp:4. Args: Read1: char, G: char, P: Parameters"]
pub fn transcript_alignscore_l4_transcript_alignscore(
    transcript: &mut crate::transcript::Transcript,
    read1_0: &[u8],
    read1_2: &[u8],
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
) -> i32 {
    transcript.max_score = 0;
    transcript.n_mm = 0;
    transcript.n_match = 0;

    if transcript.n_exons == 0 {
        return transcript.max_score;
    }

    let r = if transcript.ro_str == 0 {
        read1_0
    } else {
        read1_2
    };
    for iex in 0..transcript.n_exons as usize {
        for ii in 0..transcript.exons[iex][EX_L] as usize {
            let r1 = r[ii + transcript.exons[iex][EX_R] as usize];
            let g1 = g[ii + transcript.exons[iex][EX_G] as usize];
            if r1 > 3 || g1 > 3 {
            } else if r1 == g1 {
                transcript.max_score += 1;
                transcript.n_match += 1;
            } else {
                transcript.n_mm += 1;
                transcript.max_score -= 1;
            }
        }
    }

    for iex in 0..transcript.n_exons as usize - 1 {
        if transcript.sj_annot[iex] == 1 {
            transcript.max_score += sjdb_score;
        } else {
            match transcript.canon_sj[iex] {
                -3 => {}
                -2 => {
                    transcript.max_score += (transcript.exons[iex + 1][EX_R]
                        - transcript.exons[iex][EX_R]
                        - transcript.exons[iex][EX_L])
                        as i32
                        * score_ins_base
                        + score_ins_open;
                }
                -1 => {
                    transcript.max_score += (transcript.exons[iex + 1][EX_G]
                        - transcript.exons[iex][EX_G]
                        - transcript.exons[iex][EX_L])
                        as i32
                        * score_del_base
                        + score_del_open;
                }
                0 => transcript.max_score += score_gap_noncan + score_gap,
                1 | 2 => transcript.max_score += score_gap,
                3 | 4 => transcript.max_score += score_gap_gcag + score_gap,
                5 | 6 => transcript.max_score += score_gap_atac + score_gap,
                _ => {}
            }
        }
    }

    if score_genomic_length_log2scale != 0.0 {
        let last = transcript.n_exons as usize - 1;
        let genomic_length = std::cmp::max(
            1,
            transcript.exons[last][EX_G] + transcript.exons[last][EX_L] - transcript.exons[0][EX_G],
        );
        transcript.max_score +=
            ((genomic_length as f64).log2() * score_genomic_length_log2scale - 0.5).ceil() as i32;
    }

    transcript.max_score
}
