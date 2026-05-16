#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ClipMate::clip` at STAR/source/ClipMate_clip.cpp:5. Args: Lread: uint, seqNum: char"]
pub fn clipmate_clip_l5_clipmate_clip(
    clip_mate: &mut crate::clip_mate::ClipMate,
    l_read: &mut u64,
    seq_num: &mut [u8],
) -> u32 {
    clip_mate.clipped_n = 0;

    if clip_mate.type_ < 0 {
        return 0;
    }

    let l_read_old = *l_read;

    if clip_mate.n > 0 {
        if *l_read > clip_mate.n as u64 {
            *l_read -= clip_mate.n as u64;
            clip_mate.clipped_n += clip_mate.n;
            if clip_mate.type_ == 0 {
                seq_num.copy_within(clip_mate.n as usize..l_read_old as usize, 0);
            }
        } else {
            *l_read = 0;
            clip_mate.clipped_n = l_read_old as u32;
        }
    }

    if !clip_mate.ad_seq.is_empty() {
        match clip_mate.type_ {
            1 => {
                clip_mate.clipped_ad_n = (*l_read as u32)
                    - sequencefuns_l293_localsearch(
                        seq_num,
                        *l_read as u32,
                        &clip_mate.ad_seq_num,
                        clip_mate.ad_seq_num.len() as u32,
                        clip_mate.ad_mmp,
                    );
            }
            10 => {
                clip_mate.clipped_ad_n = std::cmp::min(clip_mate.clipped_info as u32, *l_read as u32);
                seq_num.copy_within(clip_mate.clipped_ad_n as usize..*l_read as usize, 0);
            }
            11 => {
                clip_mate.clipped_ad_n = clipcr4_l82_clipcr4_polytail3p(seq_num, *l_read as u32);
            }
            _ => {}
        }

        *l_read -= clip_mate.clipped_ad_n as u64;
        clip_mate.clipped_n += clip_mate.clipped_ad_n;
    }

    if clip_mate.n_after_ad > 0 {
        if *l_read > clip_mate.n_after_ad as u64 {
            *l_read -= clip_mate.n_after_ad as u64;
            clip_mate.clipped_n += clip_mate.n_after_ad;
            if clip_mate.type_ == 0 {
                seq_num.copy_within(
                    clip_mate.n_after_ad as usize..(*l_read + clip_mate.n_after_ad as u64) as usize,
                    0,
                );
            }
        } else {
            *l_read = 0;
            clip_mate.clipped_n = l_read_old as u32;
        }
    }

    clip_mate.clipped_n
}
