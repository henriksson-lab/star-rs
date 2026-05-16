#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ClipMate::clipChunk` at STAR/source/ClipMate_clipChunk.cpp:7. Args: chArr: char, chSize: uint64"]
pub fn clipmate_clipchunk_l7_clipmate_clipchunk(
    clip_mate: &mut crate::clip_mate::ClipMate,
    ch_arr: &mut [u8],
    ch_size: u64,
) -> Result<(), String> {
    if clip_mate.type_ != 10 {
        return Ok(());
    }

    let Some(cr4) = clip_mate.cr4.as_mut() else {
        return Err("ClipMate::clipChunk requires CR4 state for type=10".to_string());
    };

    let mut ch_a1 = 0usize;
    let ch_size = ch_size as usize;
    let mut ch_good = true;
    while ch_good && ch_a1 < ch_size {
        let mut db_n1 = cr4.db_n as i32;
        for idb in 0..cr4.db_n as usize {
            let name_end_rel = clipmate_clipchunk_l55_findchar(&ch_arr[ch_a1..], b'\n');
            ch_a1 += name_end_rel + 1;

            let seq_end_rel = clipmate_clipchunk_l55_findchar(&ch_arr[ch_a1..], b'\n');
            let ch_a2 = ch_a1 + seq_end_rel;
            let r_l = (ch_a2 - ch_a1) as u32;

            clipcr4_l43_clipcr4_opalfilloneseq(cr4, idb as u32, &ch_arr[ch_a1..ch_a2], r_l);
            cr4.store_clip[idb] = (ch_a2 + 1) as u32;

            ch_a1 = ch_a2 + 3 + r_l as usize + 1;
            if ch_a1 > ch_size {
                ch_good = false;
                db_n1 = idb as i32 + 1;
                break;
            }
            if ch_a1 == ch_size {
                ch_good = false;
                db_n1 = idb as i32 + 1;
                break;
            }
        }

        clipcr4_l72_clipcr4_opalalign(
            cr4,
            &clip_mate.ad_seq_num,
            clip_mate.ad_seq_num.len() as u32,
            db_n1,
        );

        for idb in 0..db_n1 as usize {
            let l = cr4.opal_res[idb].end_location_target + 1;
            let s = cr4.opal_res[idb].score;
            let l0 = s < 20 || (s == 20 && l > 26) || (s == 21 && l > 30);
            let offset = cr4.store_clip[idb] as usize;
            if offset >= ch_arr.len() {
                return Err("ClipMate::clipChunk storeClip offset is out of bounds".to_string());
            }
            ch_arr[offset] = if l0 { 0 } else { l as u8 };
        }
    }

    Ok(())
}

#[doc = "Original `findChar` at STAR/source/ClipMate_clipChunk.cpp:55. Args: arr: char, c: char"]
pub fn clipmate_clipchunk_l55_findchar(arr: &[u8], c: u8) -> usize {
    let mut index = 0;
    while arr[index] != c {
        index += 1;
    }
    index
}
