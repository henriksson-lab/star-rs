#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::addBAMtags` at STAR/source/SoloFeature_addBAMtags.cpp:5. Args: bam0: char, size0: uint32, bam1: char"]
pub fn solofeature_addbamtags_l5_solofeature_addbamtags(
    bam0: &[u8],
    size0: u32,
    read_info: &[crate::solo_feature::SoloFeatureReadInfo],
    cb_wlstr: &[String],
    umi_l: u32,
) -> Vec<u8> {
    let size0_usize = size0 as usize;
    let mut iread_bytes = [0u8; 8];
    iread_bytes.copy_from_slice(&bam0[size0_usize..size0_usize + 8]);
    let iread = (u64::from_ne_bytes(iread_bytes) >> 32) as usize;

    let mut cb = "-".to_string();
    let mut umi = "-".to_string();
    if read_info[iread].cb + 1 != 0 {
        cb = cb_wlstr[read_info[iread].cb as usize].clone();
    }
    if read_info[iread].umi.wrapping_add(1) != 0 {
        umi = sequencefuns_l267_convertnuclint64tostring(read_info[iread].umi, umi_l);
    }

    let mut bam1 = vec![0u8; size0_usize + cb.len() + umi.len() + 16];
    bam1[..size0_usize].copy_from_slice(&bam0[..size0_usize]);
    let mut size1 = size0_usize;
    size1 += bamfunctions_l124_bamattrarraywrite(cb.as_str(), b"CB", &mut bam1[size1..]) as usize;
    size1 += bamfunctions_l124_bamattrarraywrite(umi.as_str(), b"UB", &mut bam1[size1..]) as usize;
    let block_len = (size1 as u32 - std::mem::size_of::<u32>() as u32).to_ne_bytes();
    bam1[..4].copy_from_slice(&block_len);
    bam1.truncate(size1);
    bam1
}
