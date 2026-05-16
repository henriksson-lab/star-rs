#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ClipMate::initialize` at STAR/source/ClipMate_initialize.cpp:5. Args: Nin: uint32, adSeqIn: string, NafterAdin: uint32, adMMpIn: double"]
pub fn clipmate_initialize_l5_clipmate_initialize(
    clip_mate: &mut crate::clip_mate::ClipMate,
    n_in: u32,
    ad_seq_in: &str,
    n_after_ad_in: u32,
    ad_mmp_in: f64,
) {
    clip_mate.n = n_in;

    clip_mate.ad_seq = ad_seq_in.to_string();
    if clip_mate.ad_seq == "-" {
        clip_mate.ad_seq.clear();
        clip_mate.ad_seq_num.clear();
    } else if clip_mate.ad_seq == "polyA" {
        clip_mate.ad_seq_num.clear();
        clip_mate.ad_seq_num.resize(DEF_READ_SEQ_LENGTH_MAX, 0);
    } else {
        clip_mate.ad_seq_num.resize(clip_mate.ad_seq.len(), 0);
        sequencefuns_l131_convertnucleotidestonumbers(
            clip_mate.ad_seq.as_bytes(),
            &mut clip_mate.ad_seq_num,
            clip_mate.ad_seq.len() as u64,
        );
    }

    if clip_mate.n == 0 && clip_mate.ad_seq.is_empty() {
        clip_mate.type_ = -1;
    }

    if clip_mate.type_ == 10 {
        clip_mate.cr4 = Some(clipcr4_l3_clipcr4_clipcr4());
    }

    clip_mate.n_after_ad = n_after_ad_in;
    clip_mate.ad_mmp = ad_mmp_in;
}
