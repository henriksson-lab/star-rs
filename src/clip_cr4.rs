#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ClipCR4` at STAR/source/ClipCR4.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ClipCR4 {
    pub db_n: u32,
    pub score_matrix: Vec<i32>,
    pub read_len: u32,
    pub alphabet_length: i32,
    pub gap_open: i32,
    pub gap_ext: i32,
    pub db_seq_arr: Vec<u8>,
    pub db_seqs_len: Vec<i32>,
    pub store_clip: Vec<u32>,
    pub opal_res: Vec<OpalSearchResult>,
}

#[doc = "Original `ClipCR4::ClipCR4` at STAR/source/ClipCR4.cpp:3. Args: "]
pub fn clipcr4_l3_clipcr4_clipcr4() -> crate::clip_cr4::ClipCR4 {
    let db_n = 64u32;
    let read_len = 91u32;
    crate::clip_cr4::ClipCR4 {
        db_n,
        score_matrix: vec![
            1, -2, -2, -2, -2, -2, 1, -2, -2, -2, -2, -2, 1, -2, -2, -2, -2, -2, 1, -2, -2, -2, -2,
            -2, 0,
        ],
        read_len,
        alphabet_length: 5,
        gap_open: 2,
        gap_ext: 2,
        db_seq_arr: vec![0; db_n as usize * read_len as usize],
        db_seqs_len: vec![read_len as i32; db_n as usize],
        store_clip: vec![0; db_n as usize],
        opal_res: vec![crate::opal::opal::OpalSearchResult::default(); db_n as usize],
    }
}

#[doc = "Original `ClipCR4::opalFillOneSeq` at STAR/source/ClipCR4.cpp:43. Args: idb: uint32, seq: char, seqL: uint32"]
pub fn clipcr4_l43_clipcr4_opalfilloneseq(
    clip: &mut crate::clip_cr4::ClipCR4,
    idb: u32,
    seq: &[u8],
    seq_l: u32,
) {
    let min_len = std::cmp::min(seq_l, clip.read_len) as usize;
    let offset = idb as usize * clip.read_len as usize;
    let copy_len = std::cmp::min(min_len, seq.len());
    clip.db_seq_arr[offset..offset + copy_len].copy_from_slice(&seq[..copy_len]);
    if copy_len < min_len {
        clip.db_seq_arr[offset + copy_len..offset + min_len].fill(4);
    }
    for ib in 0..min_len {
        clip.db_seq_arr[offset + ib] = match clip.db_seq_arr[offset + ib] {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        };
    }
    if seq_l < clip.read_len {
        let fill_start = offset + seq_l as usize;
        let fill_end = offset + clip.read_len as usize;
        clip.db_seq_arr[fill_start..fill_end].fill(4);
    }
}

#[doc = "Original `ClipCR4::opalAlign` at STAR/source/ClipCR4.cpp:72. Args: query: uint8_t, queryLen: uint32, dbN1: int"]
pub fn clipcr4_l72_clipcr4_opalalign(
    clip: &mut crate::clip_cr4::ClipCR4,
    query: &[u8],
    query_len: u32,
    db_n1: i32,
) -> i32 {
    for idb in 0..db_n1 as usize {
        opal_l1553_opalinitsearchresult(&mut clip.opal_res[idb]);
    }

    let read_len = clip.read_len as usize;
    let db_storage: Vec<Vec<u8>> = (0..db_n1 as usize)
        .map(|idb| {
            let start = idb * read_len;
            clip.db_seq_arr[start..start + read_len].to_vec()
        })
        .collect();
    let db_slices: Vec<&[u8]> = db_storage.iter().map(|seq| seq.as_slice()).collect();

    opal_l1437_opalsearchdatabase(
        query,
        query_len as i32,
        &db_slices,
        db_n1,
        &clip.db_seqs_len,
        clip.gap_open,
        clip.gap_ext,
        &clip.score_matrix,
        clip.alphabet_length,
        &mut clip.opal_res[..db_n1 as usize],
        OPAL_SEARCH_SCORE_END,
        OPAL_MODE_OV,
        OPAL_OVERFLOW_BUCKETS,
    )
}

#[doc = "Original `ClipCR4::polyTail3p` at STAR/source/ClipCR4.cpp:82. Args: seq: char, seqLen: uint32"]
pub fn clipcr4_l82_clipcr4_polytail3p(seq: &[u8], seq_len: u32) -> u32 {
    if seq_len < 20 {
        return 0;
    }

    let mut ib1 = seq_len - 1;
    let mut score = 0i32;
    let mut score1 = 0i32;
    for ib in 1..=seq_len {
        if seq[(seq_len - ib) as usize] == 0 {
            score += 1;
            if score * 10 >= ib as i32 * 7 {
                ib1 = ib;
                score1 = score;
            }
        } else {
            score -= 2;
            if ib as i32 - score > 27 {
                break;
            }
        }
    }
    if score1 < 20 {
        ib1 = 0;
    }
    ib1
}
