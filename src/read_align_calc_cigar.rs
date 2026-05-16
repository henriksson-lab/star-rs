#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::calcCIGAR` at STAR/source/ReadAlign_calcCIGAR.cpp:3. Args: trOut: Transcript, nMates: uint, iExMate: uint, leftMate: uint"]
pub fn readalign_calccigar_l3_readalign_calccigar(
    read_align: &mut crate::read_align::ReadAlign,
    tr_out: &crate::transcript::Transcript,
    n_mates: u64,
    i_ex_mate: u64,
    left_mate: u64,
) {
    read_align.mates_cigar.clear();
    for imate in 0..n_mates {
        let i_ex1 = if imate == 0 { 0 } else { i_ex_mate + 1 };
        let i_ex2 = if imate == 0 {
            i_ex_mate
        } else {
            tr_out.n_exons - 1
        };
        let mate = tr_out.exons[i_ex1 as usize][EX_IFRAG] as usize;
        let str_ = tr_out.str_ as usize;

        let mut cigar = String::new();

        let trim_l = if str_ == 0 && mate == 0 {
            read_align.clip_mates[mate][0].clipped_n
        } else if str_ == 0 && mate == 1 {
            read_align.clip_mates[mate][1].clipped_n
        } else if str_ == 1 && mate == 0 {
            read_align.clip_mates[mate][1].clipped_n
        } else {
            read_align.clip_mates[mate][0].clipped_n
        };

        let left_mate_usize = left_mate as usize;
        let trim_l1 = trim_l as u64 + tr_out.exons[i_ex1 as usize][EX_R]
            - if tr_out.exons[i_ex1 as usize][EX_R] < read_align.read_length[left_mate_usize] {
                0
            } else {
                read_align.read_length[left_mate_usize] + 1
            };
        if trim_l1 > 0 {
            cigar.push_str(&format!("{}S", trim_l1));
        }

        for ii in i_ex1..=i_ex2 {
            let ii_usize = ii as usize;
            if ii > i_ex1 {
                let prev = ii_usize - 1;
                let gap_g = tr_out.exons[ii_usize][EX_G]
                    - (tr_out.exons[prev][EX_G] + tr_out.exons[prev][EX_L]);
                let gap_r = tr_out.exons[ii_usize][EX_R]
                    - tr_out.exons[prev][EX_R]
                    - tr_out.exons[prev][EX_L];
                if gap_r > 0 {
                    cigar.push_str(&format!("{}I", gap_r));
                }
                if tr_out.canon_sj[prev] >= 0 || tr_out.sj_annot[prev] == 1 {
                    cigar.push_str(&format!("{}N", gap_g));
                } else if gap_g > 0 {
                    cigar.push_str(&format!("{}D", gap_g));
                }
            }
            cigar.push_str(&format!("{}M", tr_out.exons[ii_usize][EX_L]));
        }

        let i_ex1_usize = i_ex1 as usize;
        let i_ex2_usize = i_ex2 as usize;
        let trim_r1 = (if tr_out.exons[i_ex1_usize][EX_R] < read_align.read_length[left_mate_usize]
        {
            read_align.read_length_original[left_mate_usize]
        } else {
            read_align.read_length[left_mate_usize] + 1 + read_align.read_length_original[mate]
        }) - tr_out.exons[i_ex2_usize][EX_R]
            - tr_out.exons[i_ex2_usize][EX_L]
            - trim_l as u64;
        if trim_r1 > 0 {
            cigar.push_str(&format!("{}S", trim_r1));
        }
        read_align.mates_cigar.push(cigar);
    }
}
