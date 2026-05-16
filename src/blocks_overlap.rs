#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `blocksOverlap` at STAR/source/blocksOverlap.cpp:3. Args: t1: Transcript, t2: Transcript"]
pub fn blocksoverlap_l3_blocksoverlap(
    t1: &crate::transcript::Transcript,
    t2: &crate::transcript::Transcript,
) -> u64 {
    let mut i1 = 0usize;
    let mut i2 = 0usize;
    let mut n_overlap = 0;
    while i1 < t1.n_exons as usize && i2 < t2.n_exons as usize {
        let rs1 = t1.exons[i1][EX_R];
        let rs2 = t2.exons[i2][EX_R];
        let re1 = t1.exons[i1][EX_R] + t1.exons[i1][EX_L];
        let re2 = t2.exons[i2][EX_R] + t2.exons[i2][EX_L];
        let gs1 = t1.exons[i1][EX_G];
        let gs2 = t2.exons[i2][EX_G];

        if rs1 >= re2 {
            i2 += 1;
        } else if rs2 >= re1 {
            i1 += 1;
        } else if gs1.wrapping_sub(rs1) != gs2.wrapping_sub(rs2) {
            if re1 >= re2 {
                i2 += 1;
            }
            if re2 >= re1 {
                i1 += 1;
            }
        } else {
            n_overlap += re1.min(re2) - rs1.max(rs2);
            if re1 >= re2 {
                i2 += 1;
            }
            if re2 >= re1 {
                i1 += 1;
            }
        }
    }
    n_overlap
}
