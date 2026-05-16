#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Transcript::generateCigarP` at STAR/source/Transcript_generateCigarP.cpp:4. Args: "]
pub fn transcript_generatecigarp_l4_transcript_generatecigarp(
    tr: &crate::transcript::Transcript,
) -> String {
    let mut cigar = String::new();
    let mut left_mate = 0usize;
    if tr.read_nmates > 1 {
        left_mate = tr.str_ as usize;
    }

    let read_len_left = tr.read_length_original[left_mate];
    let mut trim_l = tr.exons[0][EX_R].wrapping_sub(if tr.exons[0][EX_R] < read_len_left {
        0
    } else {
        read_len_left + 1
    });
    if trim_l > 0 {
        cigar.push_str(&format!("{}S", trim_l));
    }

    for ii in 0..tr.n_exons as usize {
        if ii > 0 {
            let prev_end_g = tr.exons[ii - 1][EX_G] + tr.exons[ii - 1][EX_L];
            let gap_g = tr.exons[ii][EX_G].wrapping_sub(prev_end_g);

            if tr.exons[ii][EX_G] >= prev_end_g {
                if tr.canon_sj[ii - 1] == -3 {
                    let s1 = read_len_left - (tr.exons[ii - 1][EX_R] + tr.exons[ii - 1][EX_L]);
                    let s2 = tr.exons[ii][EX_R] - (read_len_left + 1);
                    if s1 > 0 {
                        cigar.push_str(&format!("{}S", s1));
                    }
                    cigar.push_str(&format!("{}p", gap_g));
                    if s2 > 0 {
                        cigar.push_str(&format!("{}S", s2));
                    }
                } else {
                    let gap_r =
                        tr.exons[ii][EX_R] - tr.exons[ii - 1][EX_R] - tr.exons[ii - 1][EX_L];
                    if gap_r > 0 {
                        cigar.push_str(&format!("{}I", gap_r));
                    }
                    if tr.canon_sj[ii - 1] >= 0 || tr.sj_annot[ii - 1] == 1 {
                        cigar.push_str(&format!("{}N", gap_g));
                    } else if gap_g > 0 {
                        cigar.push_str(&format!("{}D", gap_g));
                    }
                }
            } else {
                cigar.push_str(&format!(
                    "-{}p",
                    (tr.exons[ii - 1][EX_G] + tr.exons[ii - 1][EX_L]) - tr.exons[ii][EX_G]
                ));
            }
        }
        cigar.push_str(&format!("{}M", tr.exons[ii][EX_L]));
    }

    let last = tr.n_exons as usize - 1;
    trim_l = (if tr.exons[last][EX_R] < read_len_left {
        read_len_left
    } else {
        tr.read_length_pair_original
    }) - tr.exons[last][EX_R]
        - tr.exons[last][EX_L];
    if trim_l > 0 {
        cigar.push_str(&format!("{}S", trim_l));
    }
    cigar
}
