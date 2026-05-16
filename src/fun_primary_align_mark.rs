#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `funPrimaryAlignMark` at STAR/source/funPrimaryAlignMark.cpp:3. Args: trMult: Transcript, nTr: uint64, P: Parameters, maxScore: int, rngUniformReal0to1: std::uniform_real_distribution<double>, rngMultOrder: std::mt19937"]
pub fn funprimaryalignmark_l3_funprimaryalignmark(
    tr_mult: &mut [crate::transcript::Transcript],
    n_tr: u64,
    max_score: i32,
    out_multimapper_order_random: bool,
    out_sam_mult_nmax_is_limited: bool,
    out_sam_primary_flag: &str,
    rng_uniform_real_0_to_1: &[f64],
) {
    if n_tr == 1 {
        tr_mult[0].primary_flag = true;
    } else {
        let mut nbest = 0usize;
        if out_multimapper_order_random || out_sam_mult_nmax_is_limited {
            for itr in 0..n_tr as usize {
                if tr_mult[itr].max_score == max_score {
                    tr_mult.swap(itr, nbest);
                    nbest += 1;
                }
            }
        }

        if out_multimapper_order_random {
            let mut irng = 0usize;
            for itr in (1..nbest).rev() {
                let rand1 = (rng_uniform_real_0_to_1[irng] * itr as f64 + 0.5) as usize;
                irng += 1;
                tr_mult.swap(itr, rand1);
            }
            for itr in (1..(n_tr as usize - nbest)).rev() {
                let rand1 = (rng_uniform_real_0_to_1[irng] * itr as f64 + 0.5) as usize;
                irng += 1;
                tr_mult.swap(nbest + itr, nbest + rand1);
            }
        }

        if out_sam_primary_flag == "AllBestScore" {
            for itr in 0..n_tr as usize {
                if tr_mult[itr].max_score == max_score {
                    tr_mult[itr].primary_flag = true;
                }
            }
        } else if out_multimapper_order_random || out_sam_mult_nmax_is_limited {
            tr_mult[0].primary_flag = true;
        } else {
            tr_mult[0].primary_flag = true;
        }
    }
}
