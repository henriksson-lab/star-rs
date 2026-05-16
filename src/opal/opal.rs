#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `Cell` at STAR/source/opal/opal.cpp:1206."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Cell {
    pub h: i64,
    pub e: i64,
    pub f: i64,
}

#[doc = "Original struct `CellEH` at STAR/source/opal/opal.cpp:156."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CellEH {
    pub h: (),
    pub e: (),
}

#[doc = "Original struct `OpalSearchResult` at STAR/source/opal/opal.h:47."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OpalSearchResult {
    pub score_set: i32,
    pub score: i32,
    pub end_location_target: i32,
    pub end_location_query: i32,
    pub start_location_target: i32,
    pub start_location_query: i32,
    pub alignment: Option<Vec<u8>>,
    pub alignment_length: i32,
}

#[doc = "Original class `Simd` at STAR/source/opal/opal.cpp:551."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Simd {}

#[doc = "Original struct `Simd<char>` at STAR/source/opal/opal.cpp:554."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Simd_char_ {
    pub numseqs: i64,
    pub satarthm: bool,
}

#[doc = "Original struct `Simd<int>` at STAR/source/opal/opal.cpp:580."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Simd_int_ {
    pub numseqs: i64,
    pub satarthm: bool,
}

#[doc = "Original struct `Simd<short>` at STAR/source/opal/opal.cpp:567."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Simd_short_ {
    pub numseqs: i64,
    pub satarthm: bool,
}

#[doc = "Original struct `SimdSW` at STAR/source/opal/opal.cpp:95."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SimdSW {}

#[doc = "Original struct `SimdSW<char>` at STAR/source/opal/opal.cpp:98."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SimdSW_char_ {
    pub numseqs: i64,
    pub satarthm: bool,
    pub negrange: bool,
}

#[doc = "Original struct `SimdSW<int>` at STAR/source/opal/opal.cpp:126."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SimdSW_int_ {
    pub numseqs: i64,
    pub satarthm: bool,
    pub negrange: bool,
}

#[doc = "Original struct `SimdSW<short>` at STAR/source/opal/opal.cpp:112."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SimdSW_short_ {
    pub numseqs: i64,
    pub satarthm: bool,
    pub negrange: bool,
}

#[doc = "Original `simdIsAllZeroes` at STAR/source/opal/opal.cpp:87. Args: a: __mxxxi"]
pub fn opal_l87_simdisallzeroes(a: &[u8]) -> i32 {
    if a.iter().all(|value| *value == 0) {
        1
    } else {
        0
    }
}

#[doc = "Original `add` at STAR/source/opal/opal.cpp:103. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l103_add(a: [i8; 32], b: [i8; 32]) -> [i8; 32] {
    let mut out = [0i8; 32];
    for ii in 0..32 {
        out[ii] = a[ii].saturating_add(b[ii]);
    }
    out
}

#[doc = "Original `sub` at STAR/source/opal/opal.cpp:104. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l104_sub(a: [i8; 32], b: [i8; 32]) -> [i8; 32] {
    let mut out = [0i8; 32];
    for ii in 0..32 {
        out[ii] = a[ii].saturating_sub(b[ii]);
    }
    out
}

#[doc = "Original `min` at STAR/source/opal/opal.cpp:105. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l105_min(a: [i8; 32], b: [i8; 32]) -> [i8; 32] {
    let mut out = [0i8; 32];
    for ii in 0..32 {
        out[ii] = if (a[ii] as u8) < (b[ii] as u8) {
            a[ii]
        } else {
            b[ii]
        };
    }
    out
}

#[doc = "Original `max` at STAR/source/opal/opal.cpp:106. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l106_max(a: [i8; 32], b: [i8; 32]) -> [i8; 32] {
    let mut out = [0i8; 32];
    for ii in 0..32 {
        out[ii] = if (a[ii] as u8) > (b[ii] as u8) {
            a[ii]
        } else {
            b[ii]
        };
    }
    out
}

#[doc = "Original `cmpgt` at STAR/source/opal/opal.cpp:107. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l107_cmpgt(a: [i8; 32], b: [i8; 32]) -> [i8; 32] {
    let mut out = [0i8; 32];
    for ii in 0..32 {
        out[ii] = if a[ii] > b[ii] { -1 } else { 0 };
    }
    out
}

#[doc = "Original `set1` at STAR/source/opal/opal.cpp:108. Args: a: int"]
pub fn opal_l108_set1(a: i32) -> [i8; 32] {
    [a as i8; 32]
}

#[doc = "Original `add` at STAR/source/opal/opal.cpp:117. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l117_add(a: [i16; 16], b: [i16; 16]) -> [i16; 16] {
    let mut out = [0i16; 16];
    for ii in 0..16 {
        out[ii] = a[ii].saturating_add(b[ii]);
    }
    out
}

#[doc = "Original `sub` at STAR/source/opal/opal.cpp:118. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l118_sub(a: [i16; 16], b: [i16; 16]) -> [i16; 16] {
    let mut out = [0i16; 16];
    for ii in 0..16 {
        out[ii] = a[ii].saturating_sub(b[ii]);
    }
    out
}

#[doc = "Original `min` at STAR/source/opal/opal.cpp:119. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l119_min(a: [i16; 16], b: [i16; 16]) -> [i16; 16] {
    let mut out = [0i16; 16];
    for ii in 0..16 {
        out[ii] = a[ii].min(b[ii]);
    }
    out
}

#[doc = "Original `max` at STAR/source/opal/opal.cpp:120. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l120_max(a: [i16; 16], b: [i16; 16]) -> [i16; 16] {
    let mut out = [0i16; 16];
    for ii in 0..16 {
        out[ii] = a[ii].max(b[ii]);
    }
    out
}

#[doc = "Original `cmpgt` at STAR/source/opal/opal.cpp:121. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l121_cmpgt(a: [i16; 16], b: [i16; 16]) -> [i16; 16] {
    let mut out = [0i16; 16];
    for ii in 0..16 {
        out[ii] = if a[ii] > b[ii] { -1 } else { 0 };
    }
    out
}

#[doc = "Original `set1` at STAR/source/opal/opal.cpp:122. Args: a: int"]
pub fn opal_l122_set1(a: i32) -> [i16; 16] {
    [a as i16; 16]
}

#[doc = "Original `add` at STAR/source/opal/opal.cpp:131. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l131_add(a: [i32; 8], b: [i32; 8]) -> [i32; 8] {
    let mut out = [0i32; 8];
    for ii in 0..8 {
        out[ii] = a[ii].wrapping_add(b[ii]);
    }
    out
}

#[doc = "Original `sub` at STAR/source/opal/opal.cpp:132. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l132_sub(a: [i32; 8], b: [i32; 8]) -> [i32; 8] {
    let mut out = [0i32; 8];
    for ii in 0..8 {
        out[ii] = a[ii].wrapping_sub(b[ii]);
    }
    out
}

#[doc = "Original `min` at STAR/source/opal/opal.cpp:133. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l133_min(a: [i32; 8], b: [i32; 8]) -> [i32; 8] {
    let mut out = [0i32; 8];
    for ii in 0..8 {
        out[ii] = a[ii].min(b[ii]);
    }
    out
}

#[doc = "Original `max` at STAR/source/opal/opal.cpp:134. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l134_max(a: [i32; 8], b: [i32; 8]) -> [i32; 8] {
    let mut out = [0i32; 8];
    for ii in 0..8 {
        out[ii] = a[ii].max(b[ii]);
    }
    out
}

#[doc = "Original `cmpgt` at STAR/source/opal/opal.cpp:135. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l135_cmpgt(a: [i32; 8], b: [i32; 8]) -> [i32; 8] {
    let mut out = [0i32; 8];
    for ii in 0..8 {
        out[ii] = if a[ii] > b[ii] { -1 } else { 0 };
    }
    out
}

#[doc = "Original `set1` at STAR/source/opal/opal.cpp:136. Args: a: int"]
pub fn opal_l136_set1(a: i32) -> [i32; 8] {
    [a; 8]
}

#[doc = "Original `print_mmxxxi` at STAR/source/opal/opal.cpp:148. Args: mm: __mxxxi"]
pub fn opal_l148_print_mmxxxi<T: std::fmt::Display>(mm: &[T]) -> String {
    let mut out = String::new();
    for value in mm {
        out.push_str(&format!("{} ", value));
    }
    out
}

#[doc = "Original `searchDatabaseSW_` at STAR/source/opal/opal.cpp:167. Args: query: unsigned char, queryLength: int, db: unsigned char, dbLength: int, dbSeqLengths: int, gapOpen: int, gapExt: int, scoreMatrix: int, alphabetLength: int, results: OpalSearchResult, searchType: int, calculated: bool, overflowMethod: int"]
pub fn opal_l167_searchdatabasesw(
    query: &[u8],
    query_length: i32,
    db: &[&[u8]],
    db_length: i32,
    db_seq_lengths: &[i32],
    gap_open: i32,
    gap_ext: i32,
    score_matrix: &[i32],
    alphabet_length: i32,
    results: &mut [crate::opal::opal::OpalSearchResult],
    search_type: i32,
    calculated: &[bool],
    _overflow_method: i32,
) -> i32 {
    for i in 0..db_length as usize {
        if calculated[i] {
            continue;
        }
        let target = db[i];
        let target_length = db_seq_lengths[i] as usize;
        let mut prev_h = vec![0i32; query_length as usize];
        let mut prev_e = vec![0i32; query_length as usize];
        let mut max_h = 0i32;
        let mut best_score = i32::MIN;
        let mut best_row = -1i32;
        let mut best_column = -1i32;

        for c in 0..target_length {
            let mut u_f = 0i32;
            let mut u_h = 0i32;
            let mut ul_h = 0i32;
            let mut rows_with_improvement = Vec::new();
            let mut curr_h = vec![0i32; query_length as usize];
            let mut curr_e = vec![0i32; query_length as usize];

            for r in 0..query_length as usize {
                let e = (prev_h[r] - gap_open).max(prev_e[r] - gap_ext);
                let f = (u_h - gap_open).max(u_f - gap_ext);
                let score =
                    score_matrix[query[r] as usize * alphabet_length as usize + target[c] as usize];
                let h = 0.max(e.max(f.max(ul_h + score)));
                if search_type != OPAL_SEARCH_SCORE && h > max_h {
                    rows_with_improvement.push(r as i32);
                }
                max_h = max_h.max(h);
                u_f = f;
                u_h = h;
                ul_h = prev_h[r];
                curr_e[r] = e;
                curr_h[r] = h;
            }

            if search_type != OPAL_SEARCH_SCORE {
                for r in rows_with_improvement {
                    let h = curr_h[r as usize];
                    if h > best_score {
                        best_score = h;
                        best_row = r;
                        best_column = c as i32;
                    }
                }
            }

            prev_h = curr_h;
            prev_e = curr_e;
        }

        opal_l1565_opalsearchresultsetscore(&mut results[i], max_h);
        if search_type != OPAL_SEARCH_SCORE {
            results[i].end_location_query = best_row;
            results[i].end_location_target = best_column;
        } else {
            results[i].end_location_query = -1;
            results[i].end_location_target = -1;
        }
    }

    0
}

#[doc = "Original `loadNextSequence` at STAR/source/opal/opal.cpp:474. Args: nextDbSeqIdx: int, dbLength: int, currDbSeqIdx: int, currDbSeqPos: unsigned char, currDbSeqLength: int, db: unsigned char, dbSeqLengths: int, calculated: bool, numEndedDbSeqs: int"]
pub fn opal_l474_loadnextsequence<'a>(
    next_db_seq_idx: &mut i32,
    db_length: i32,
    curr_db_seq_idx: &mut i32,
    curr_db_seq_pos: &mut Option<&'a [u8]>,
    curr_db_seq_length: &mut i32,
    db: &'a [&'a [u8]],
    db_seq_lengths: &[i32],
    calculated: &[bool],
    num_ended_db_seqs: &mut i32,
) -> bool {
    while *next_db_seq_idx < db_length && calculated[*next_db_seq_idx as usize] {
        *next_db_seq_idx += 1;
        *num_ended_db_seqs += 1;
    }
    if *next_db_seq_idx < db_length {
        *curr_db_seq_idx = *next_db_seq_idx;
        *curr_db_seq_pos = Some(db[*next_db_seq_idx as usize]);
        *curr_db_seq_length = db_seq_lengths[*next_db_seq_idx as usize];
        *next_db_seq_idx += 1;
        true
    } else {
        *curr_db_seq_idx = -1;
        *curr_db_seq_length = -1;
        *curr_db_seq_pos = None;
        false
    }
}

#[doc = "Original `searchDatabaseSW` at STAR/source/opal/opal.cpp:498. Args: query: unsigned char, queryLength: int, db: unsigned char, dbLength: int, dbSeqLengths: int, gapOpen: int, gapExt: int, scoreMatrix: int, alphabetLength: int, results: OpalSearchResult, searchType: int, skip: bool, overflowMethod: int"]
pub fn opal_l498_searchdatabasesw(
    query: &[u8],
    query_length: i32,
    db: &[&[u8]],
    db_length: i32,
    db_seq_lengths: &[i32],
    gap_open: i32,
    gap_ext: i32,
    score_matrix: &[i32],
    alphabet_length: i32,
    results: &mut [crate::opal::opal::OpalSearchResult],
    search_type: i32,
    skip: Option<&[bool]>,
    overflow_method: i32,
) -> i32 {
    let chunk_size = if overflow_method == OPAL_OVERFLOW_BUCKETS {
        1024
    } else {
        db_length
    };
    let mut start_idx = 0;
    while start_idx < db_length {
        let db_length_ = if start_idx + chunk_size >= db_length {
            db_length - start_idx
        } else {
            chunk_size
        };
        let mut calculated = vec![false; db_length_ as usize];
        for i in 0..db_length_ as usize {
            calculated[i] = skip.map(|s| s[start_idx as usize + i]).unwrap_or(false);
        }
        let start = start_idx as usize;
        let end = start + db_length_ as usize;
        let result_code = opal_l167_searchdatabasesw(
            query,
            query_length,
            &db[start..end],
            db_length_,
            &db_seq_lengths[start..end],
            gap_open,
            gap_ext,
            score_matrix,
            alphabet_length,
            &mut results[start..end],
            search_type,
            &calculated,
            overflow_method,
        );
        if result_code != 0 {
            return result_code;
        }
        start_idx += chunk_size;
    }

    0
}

#[doc = "Original `add` at STAR/source/opal/opal.cpp:558. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l558_add(a: [i8; 32], b: [i8; 32]) -> [i8; 32] {
    let mut out = [0i8; 32];
    for ii in 0..32 {
        out[ii] = a[ii].saturating_add(b[ii]);
    }
    out
}

#[doc = "Original `sub` at STAR/source/opal/opal.cpp:559. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l559_sub(a: [i8; 32], b: [i8; 32]) -> [i8; 32] {
    let mut out = [0i8; 32];
    for ii in 0..32 {
        out[ii] = a[ii].saturating_sub(b[ii]);
    }
    out
}

#[doc = "Original `min` at STAR/source/opal/opal.cpp:560. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l560_min(a: [i8; 32], b: [i8; 32]) -> [i8; 32] {
    let mut out = [0i8; 32];
    for ii in 0..32 {
        out[ii] = a[ii].min(b[ii]);
    }
    out
}

#[doc = "Original `max` at STAR/source/opal/opal.cpp:561. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l561_max(a: [i8; 32], b: [i8; 32]) -> [i8; 32] {
    let mut out = [0i8; 32];
    for ii in 0..32 {
        out[ii] = a[ii].max(b[ii]);
    }
    out
}

#[doc = "Original `cmpgt` at STAR/source/opal/opal.cpp:562. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l562_cmpgt(a: [i8; 32], b: [i8; 32]) -> [i8; 32] {
    let mut out = [0i8; 32];
    for ii in 0..32 {
        out[ii] = if a[ii] > b[ii] { -1 } else { 0 };
    }
    out
}

#[doc = "Original `set1` at STAR/source/opal/opal.cpp:563. Args: a: int"]
pub fn opal_l563_set1(a: i32) -> [i8; 32] {
    [a as i8; 32]
}

#[doc = "Original `add` at STAR/source/opal/opal.cpp:571. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l571_add(a: [i16; 16], b: [i16; 16]) -> [i16; 16] {
    let mut out = [0i16; 16];
    for ii in 0..16 {
        out[ii] = a[ii].saturating_add(b[ii]);
    }
    out
}

#[doc = "Original `sub` at STAR/source/opal/opal.cpp:572. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l572_sub(a: [i16; 16], b: [i16; 16]) -> [i16; 16] {
    let mut out = [0i16; 16];
    for ii in 0..16 {
        out[ii] = a[ii].saturating_sub(b[ii]);
    }
    out
}

#[doc = "Original `min` at STAR/source/opal/opal.cpp:573. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l573_min(a: [i16; 16], b: [i16; 16]) -> [i16; 16] {
    let mut out = [0i16; 16];
    for ii in 0..16 {
        out[ii] = a[ii].min(b[ii]);
    }
    out
}

#[doc = "Original `max` at STAR/source/opal/opal.cpp:574. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l574_max(a: [i16; 16], b: [i16; 16]) -> [i16; 16] {
    let mut out = [0i16; 16];
    for ii in 0..16 {
        out[ii] = a[ii].max(b[ii]);
    }
    out
}

#[doc = "Original `cmpgt` at STAR/source/opal/opal.cpp:575. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l575_cmpgt(a: [i16; 16], b: [i16; 16]) -> [i16; 16] {
    let mut out = [0i16; 16];
    for ii in 0..16 {
        out[ii] = if a[ii] > b[ii] { -1 } else { 0 };
    }
    out
}

#[doc = "Original `set1` at STAR/source/opal/opal.cpp:576. Args: a: int"]
pub fn opal_l576_set1(a: i32) -> [i16; 16] {
    [a as i16; 16]
}

#[doc = "Original `add` at STAR/source/opal/opal.cpp:584. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l584_add(a: [i32; 8], b: [i32; 8]) -> [i32; 8] {
    let mut out = [0i32; 8];
    for ii in 0..8 {
        out[ii] = a[ii].wrapping_add(b[ii]);
    }
    out
}

#[doc = "Original `sub` at STAR/source/opal/opal.cpp:585. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l585_sub(a: [i32; 8], b: [i32; 8]) -> [i32; 8] {
    let mut out = [0i32; 8];
    for ii in 0..8 {
        out[ii] = a[ii].wrapping_sub(b[ii]);
    }
    out
}

#[doc = "Original `min` at STAR/source/opal/opal.cpp:586. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l586_min(a: [i32; 8], b: [i32; 8]) -> [i32; 8] {
    let mut out = [0i32; 8];
    for ii in 0..8 {
        out[ii] = a[ii].min(b[ii]);
    }
    out
}

#[doc = "Original `max` at STAR/source/opal/opal.cpp:587. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l587_max(a: [i32; 8], b: [i32; 8]) -> [i32; 8] {
    let mut out = [0i32; 8];
    for ii in 0..8 {
        out[ii] = a[ii].max(b[ii]);
    }
    out
}

#[doc = "Original `cmpgt` at STAR/source/opal/opal.cpp:588. Args: a: __mxxxi, b: __mxxxi"]
pub fn opal_l588_cmpgt(a: [i32; 8], b: [i32; 8]) -> [i32; 8] {
    let mut out = [0i32; 8];
    for ii in 0..8 {
        out[ii] = if a[ii] > b[ii] { -1 } else { 0 };
    }
    out
}

#[doc = "Original `set1` at STAR/source/opal/opal.cpp:589. Args: a: int"]
pub fn opal_l589_set1(a: i32) -> [i32; 8] {
    [a; 8]
}

#[doc = "Original `searchDatabase_` at STAR/source/opal/opal.cpp:597. Args: query: unsigned char, queryLength: int, db: unsigned char, dbLength: int, dbSeqLengths: int, gapOpen: int, gapExt: int, scoreMatrix: int, alphabetLength: int, results: OpalSearchResult, searchType: int, calculated: bool, overflowMethod: int"]
pub fn opal_l597_searchdatabase(
    query: &[u8],
    query_length: i32,
    db: &[&[u8]],
    db_length: i32,
    db_seq_lengths: &[i32],
    gap_open: i32,
    gap_ext: i32,
    score_matrix: &[i32],
    alphabet_length: i32,
    results: &mut [crate::opal::opal::OpalSearchResult],
    search_type: i32,
    calculated: &[bool],
    _overflow_method: i32,
    mode: i32,
) -> i32 {
    if mode != OPAL_MODE_NW && mode != OPAL_MODE_HW && mode != OPAL_MODE_OV {
        return OPAL_ERR_INVALID_MODE;
    }

    for i in 0..db_length as usize {
        if calculated[i] {
            continue;
        }

        let target = db[i];
        let target_length = db_seq_lengths[i];
        let lower_score_bound = i32::MIN + gap_open.max(gap_ext);
        let mut initial_column = vec![0i32; query_length as usize];
        let initial_e = vec![lower_score_bound; query_length as usize];
        for r in 0..query_length as usize {
            initial_column[r] = -gap_open - r as i32 * gap_ext;
        }
        let mut prev_h = initial_column;
        let mut prev_e = initial_e;
        let mut h = i32::MIN;
        let mut max_score = i32::MIN;
        let mut max_score_column = -1i32;
        let mut last_column = vec![lower_score_bound; query_length as usize];

        for c in 0..target_length as usize {
            let mut curr_h = vec![lower_score_bound; query_length as usize];
            let mut curr_e = vec![lower_score_bound; query_length as usize];
            let mut u_f = lower_score_bound;
            let mut u_h = -gap_open - c as i32 * gap_ext;
            let mut ul_h = if c == 0 { 0 } else { u_h + gap_ext };

            for r in 0..query_length as usize {
                let e = (prev_h[r] - gap_open).max(prev_e[r] - gap_ext);
                let f = (u_h - gap_open).max(u_f - gap_ext);
                let score =
                    score_matrix[query[r] as usize * alphabet_length as usize + target[c] as usize];
                h = e.max(f.max(ul_h + score));
                u_f = f;
                u_h = h;
                ul_h = prev_h[r];
                curr_h[r] = h;
                curr_e[r] = e;
            }

            if mode == OPAL_MODE_HW || mode == OPAL_MODE_OV {
                if h > max_score {
                    max_score = h;
                    max_score_column = c as i32;
                }
            }
            last_column = curr_h.clone();
            prev_h = curr_h;
            prev_e = curr_e;
        }

        if mode == OPAL_MODE_NW {
            opal_l1565_opalsearchresultsetscore(&mut results[i], h);
            if search_type != OPAL_SEARCH_SCORE {
                results[i].end_location_query = query_length - 1;
                results[i].end_location_target = target_length - 1;
            } else {
                results[i].end_location_query = -1;
                results[i].end_location_target = -1;
            }
        } else if mode == OPAL_MODE_HW {
            opal_l1565_opalsearchresultsetscore(&mut results[i], max_score);
            if search_type != OPAL_SEARCH_SCORE {
                results[i].end_location_query = query_length - 1;
                results[i].end_location_target = max_score_column;
            } else {
                results[i].end_location_query = -1;
                results[i].end_location_target = -1;
            }
        } else {
            for r in 0..query_length as usize {
                if last_column[r] > max_score {
                    max_score = last_column[r];
                }
            }
            opal_l1565_opalsearchresultsetscore(&mut results[i], max_score);
            if search_type != OPAL_SEARCH_SCORE {
                let mut row = 0i32;
                while row < query_length && last_column[row as usize] != max_score {
                    row += 1;
                }
                results[i].end_location_query = row;
                results[i].end_location_target = target_length - 1;
            } else {
                results[i].end_location_query = -1;
                results[i].end_location_target = -1;
            }
        }
    }

    0
}

#[doc = "Original `searchDatabase` at STAR/source/opal/opal.cpp:986. Args: query: unsigned char, queryLength: int, db: unsigned char, dbLength: int, dbSeqLengths: int, gapOpen: int, gapExt: int, scoreMatrix: int, alphabetLength: int, results: OpalSearchResult, searchType: int, skip: bool, overflowMethod: int"]
pub fn opal_l986_searchdatabase(
    query: &[u8],
    query_length: i32,
    db: &[&[u8]],
    db_length: i32,
    db_seq_lengths: &[i32],
    gap_open: i32,
    gap_ext: i32,
    score_matrix: &[i32],
    alphabet_length: i32,
    results: &mut [crate::opal::opal::OpalSearchResult],
    search_type: i32,
    skip: Option<&[bool]>,
    overflow_method: i32,
    mode: i32,
) -> i32 {
    if mode == OPAL_MODE_SW {
        return opal_l498_searchdatabasesw(
            query,
            query_length,
            db,
            db_length,
            db_seq_lengths,
            gap_open,
            gap_ext,
            score_matrix,
            alphabet_length,
            results,
            search_type,
            skip,
            overflow_method,
        );
    }

    if mode != OPAL_MODE_NW && mode != OPAL_MODE_HW && mode != OPAL_MODE_OV {
        return OPAL_ERR_INVALID_MODE;
    }

    let chunk_size = if overflow_method == OPAL_OVERFLOW_BUCKETS {
        1024
    } else {
        db_length
    };
    let mut start_idx = 0;
    while start_idx < db_length {
        let db_length_ = if start_idx + chunk_size >= db_length {
            db_length - start_idx
        } else {
            chunk_size
        };
        let start = start_idx as usize;
        let end = start + db_length_ as usize;
        let mut calculated = vec![false; db_length_ as usize];
        for i in 0..db_length_ as usize {
            calculated[i] = skip.map(|s| s[start + i]).unwrap_or(false);
        }

        let result_code = opal_l597_searchdatabase(
            query,
            query_length,
            &db[start..end],
            db_length_,
            &db_seq_lengths[start..end],
            gap_open,
            gap_ext,
            score_matrix,
            alphabet_length,
            &mut results[start..end],
            search_type,
            &calculated,
            overflow_method,
            mode,
        );
        if result_code != 0 {
            return result_code;
        }
        start_idx += chunk_size;
    }

    0
}

#[doc = "Original `arrayMax` at STAR/source/opal/opal.cpp:1031. Args: array: T, length: int"]
pub fn opal_l1031_arraymax<T: Copy + Default + PartialOrd>(array: &[T], length: i32) -> T {
    let length = if length <= 0 {
        0
    } else {
        std::cmp::min(length as usize, array.len())
    };
    if length == 0 {
        return T::default();
    }
    let mut max_element_idx = 0usize;
    for i in 1..length {
        if array[i] > array[max_element_idx] {
            max_element_idx = i;
        }
    }
    array[max_element_idx]
}

#[doc = "Original `gapPenalty` at STAR/source/opal/opal.cpp:1048. Args: length: int, gapOpen: int, gapExt: int"]
pub fn opal_l1048_gappenalty(length: i32, gap_open: i32, gap_ext: i32) -> i32 {
    if length > 0 {
        gap_open + gap_ext * (length - 1)
    } else {
        0
    }
}

#[doc = "Original `calculateBottomBandBorderOV` at STAR/source/opal/opal.cpp:1059. Args: k: int, Q: int, T: int, Go: int, Ge: int, M: int"]
pub fn opal_l1059_calculatebottombandborderov(
    k: i32,
    q: i32,
    t: i32,
    go: i32,
    ge: i32,
    m: i32,
) -> i32 {
    let mut border = 0;
    border = border.max((q - t).min(-1 * (k + go - ge - m * t) / ge));
    let border_candidate = -1 * (k - m * q + go - ge) / (ge + m);
    if border_candidate > q - t {
        border = border.max(border_candidate);
    }
    border.min(q - 1)
}

#[doc = "Original `calculateTopBandBorderHW` at STAR/source/opal/opal.cpp:1074. Args: k: int, Q: int, T: int, Go: int, Ge: int, M: int"]
pub fn opal_l1074_calculatetopbandborderhw(
    k: i32,
    q: i32,
    t: i32,
    go: i32,
    ge: i32,
    m: i32,
) -> i32 {
    let mut border = 0;
    border = border.max((t - q).min(-1 * (k - m * q + go) / ge + 1));
    let border_candidate = -1 * (k - t * m + 2 * go + ge * (q - t - 2)) / (2 * ge + m);
    if border_candidate > t - q {
        border = border.max(border_candidate);
    }
    border.min(t - 1)
}

#[doc = "Original `calculateBottomBandBorderHW` at STAR/source/opal/opal.cpp:1089. Args: k: int, Q: int, T: int, Go: int, Ge: int, M: int"]
pub fn opal_l1089_calculatebottombandborderhw(
    k: i32,
    q: i32,
    t: i32,
    go: i32,
    ge: i32,
    m: i32,
) -> i32 {
    let mut border = 0;
    let border_candidate = -1 * (k + go - ge - q * m) / (ge + m);
    if border_candidate >= q - t {
        border = border.max(border_candidate);
    }
    if -2 * go - ge * (q - t - 2) + m * t >= k {
        border = border.max(q - t - 1);
    }
    border.min(q - 1)
}

#[doc = "Original `calculateBottomBandBorderNW` at STAR/source/opal/opal.cpp:1106. Args: k: int, Q: int, T: int, Go: int, Ge: int, M: int"]
pub fn opal_l1106_calculatebottombandbordernw(
    k: i32,
    q: i32,
    t: i32,
    go: i32,
    ge: i32,
    m: i32,
) -> i32 {
    let mut border = 0;
    let border_candidate = -1 * (k + 2 * go - m * q + ge * (t - q - 2)) / (2 * ge + m);
    if border_candidate > q - t {
        border = border.max(border_candidate);
    }
    if q - t <= -1 * (k + go - m * t - ge) / ge {
        border = border.max(q - t);
    }
    if -2 * go - ge * (q - t - 2) + m * t >= k {
        border = border.max(q - t - 1);
    }
    border.min(q - 1)
}

#[doc = "Original `calculateBandBorders` at STAR/source/opal/opal.cpp:1153. Args: k: int, mode: int, Q: int, T: int, Go: int, Ge: int, M: int"]
pub fn opal_l1153_calculatebandborders(
    k: i32,
    mode: i32,
    q: i32,
    t: i32,
    go: i32,
    ge: i32,
    m: i32,
) -> (i32, i32) {
    if mode == OPAL_MODE_OV || mode == OPAL_MODE_SW {
        if m * q.min(t) >= k {
            (
                opal_l1059_calculatebottombandborderov(k, q, t, go, ge, m),
                opal_l1059_calculatebottombandborderov(k, t, q, go, ge, m),
            )
        } else {
            (-1, -1)
        }
    } else if mode == OPAL_MODE_HW {
        if m * q.min(t) - opal_l1048_gappenalty(q - q.min(t), go, ge) >= k {
            (
                opal_l1089_calculatebottombandborderhw(k, q, t, go, ge, m),
                opal_l1074_calculatetopbandborderhw(k, q, t, go, ge, m),
            )
        } else {
            (-1, -1)
        }
    } else if mode == OPAL_MODE_NW {
        if m * q.min(t) - opal_l1048_gappenalty((q - t).abs(), go, ge) >= k {
            (
                opal_l1106_calculatebottombandbordernw(k, q, t, go, ge, m),
                opal_l1106_calculatebottombandbordernw(k, t, q, go, ge, m),
            )
        } else {
            (-1, -1)
        }
    } else {
        unreachable!()
    }
}

#[doc = "Original `createReverseCopy` at STAR/source/opal/opal.cpp:1188. Args: seq: unsigned char, length: int"]
pub fn opal_l1188_createreversecopy(seq: &[u8], length: i32) -> Vec<u8> {
    let length = if length <= 0 {
        0
    } else {
        std::cmp::min(length as usize, seq.len())
    };
    let mut r_seq = vec![0; length];
    for i in 0..length {
        r_seq[i] = seq[length - i - 1];
    }
    r_seq
}

#[doc = "Original `revertArray` at STAR/source/opal/opal.cpp:1197. Args: array: T, length: int"]
pub fn opal_l1197_revertarray<T: Copy>(array: &mut [T], length: i32) {
    let length = if length <= 0 {
        0
    } else {
        std::cmp::min(length as usize, array.len())
    };
    for i in 0..(length / 2) {
        let tmp = array[i];
        array[i] = array[length - 1 - i];
        array[length - 1 - i] = tmp;
    }
}

#[doc = "Original `findAlignment` at STAR/source/opal/opal.cpp:1238. Args: query: unsigned char, queryLength: int, target: unsigned char, targetLength: int, gapOpen: int, gapExt: int, scoreMatrix: int, alphabetLength: int, scoreLimit: int, result: OpalSearchResult, mode: int"]
pub fn opal_l1238_findalignment(
    query: &[u8],
    query_length: i32,
    target: &[u8],
    target_length: i32,
    gap_open: i32,
    gap_ext: i32,
    score_matrix: &[i32],
    alphabet_length: i32,
    score_limit: i32,
    result: &mut crate::opal::opal::OpalSearchResult,
    mode: i32,
) {
    #[derive(Clone, Copy)]
    struct Cell {
        h: i32,
        e: i32,
        f: i32,
    }
    #[derive(Clone, Copy)]
    enum Field {
        H,
        E,
        F,
    }

    let band_borders = opal_l1153_calculatebandborders(
        score_limit,
        mode,
        query_length,
        target_length,
        gap_open,
        gap_ext,
        opal_l1031_arraymax(score_matrix, alphabet_length * alphabet_length),
    );
    assert!(band_borders.0 >= 0 && band_borders.0 < query_length);
    assert!(band_borders.1 >= 0 && band_borders.1 < target_length);

    let lower_score_bound = i32::MIN + gap_open.max(gap_ext);
    let mut matrix = vec![
        vec![
            Cell {
                h: lower_score_bound,
                e: lower_score_bound,
                f: lower_score_bound,
            };
            query_length as usize
        ];
        target_length as usize
    ];
    let mut initial_column = vec![
        Cell {
            h: 0,
            e: lower_score_bound,
            f: lower_score_bound,
        };
        query_length as usize
    ];
    for r in 0..query_length as usize {
        initial_column[r].h = -gap_open - r as i32 * gap_ext;
    }

    let mut prev_column = initial_column.clone();
    let mut max_score = i32::MIN;
    let mut h = i32::MIN;
    let mut c = 0;
    while c < target_length && max_score < score_limit {
        let r_band_start = 0.max(c - band_borders.1);
        let r_band_end = (query_length - 1).min(c + band_borders.0);

        let (mut u_f, mut u_h, mut ul_h) = if r_band_start == 0 {
            let u_h = -gap_open - c * gap_ext;
            (
                lower_score_bound,
                u_h,
                if c == 0 { 0 } else { u_h + gap_ext },
            )
        } else {
            (
                lower_score_bound,
                lower_score_bound,
                prev_column[r_band_start as usize - 1].h,
            )
        };

        for r in r_band_start..=r_band_end {
            let r_usize = r as usize;
            let e = (prev_column[r_usize].h - gap_open).max(prev_column[r_usize].e - gap_ext);
            let f = (u_h - gap_open).max(u_f - gap_ext);
            let score = score_matrix
                [query[r_usize] as usize * alphabet_length as usize + target[c as usize] as usize];
            h = e.max(f.max(ul_h + score));

            if mode == OPAL_MODE_SW || (mode == OPAL_MODE_OV && c == target_length - 1) {
                max_score = max_score.max(h);
            }

            u_f = f;
            u_h = h;
            ul_h = prev_column[r_usize].h;

            matrix[c as usize][r_usize].h = h;
            matrix[c as usize][r_usize].e = e;
            matrix[c as usize][r_usize].f = f;
        }

        for r in 0..r_band_start as usize {
            matrix[c as usize][r] = Cell {
                h: lower_score_bound,
                e: lower_score_bound,
                f: lower_score_bound,
            };
        }
        for r in r_band_end as usize + 1..query_length as usize {
            matrix[c as usize][r] = Cell {
                h: lower_score_bound,
                e: lower_score_bound,
                f: lower_score_bound,
            };
        }

        if mode == OPAL_MODE_HW || mode == OPAL_MODE_OV {
            max_score = max_score.max(h);
        }
        prev_column = matrix[c as usize].clone();
        c += 1;
    }
    let last_column_idx = c - 1;

    result.start_location_target = 0;
    result.start_location_query = 0;
    result.score_set = 1;
    match mode {
        OPAL_MODE_NW => {
            opal_l1565_opalsearchresultsetscore(result, h);
            result.end_location_target = target_length - 1;
            result.end_location_query = query_length - 1;
        }
        OPAL_MODE_HW => {
            opal_l1565_opalsearchresultsetscore(result, max_score);
            result.end_location_target = last_column_idx;
            result.end_location_query = query_length - 1;
        }
        OPAL_MODE_SW | OPAL_MODE_OV => {
            opal_l1565_opalsearchresultsetscore(result, max_score);
            result.end_location_target = last_column_idx;
            let mut r = 0;
            while r < query_length && matrix[last_column_idx as usize][r as usize].h != max_score {
                r += 1;
            }
            assert!(r < query_length);
            result.end_location_query = r;
        }
        _ => unreachable!(),
    }

    let mut alignment = Vec::with_capacity(
        (result.end_location_query + result.end_location_target).max(0) as usize,
    );
    let mut r_idx = result.end_location_query;
    let mut c_idx = result.end_location_target;
    let mut field = Field::H;
    while r_idx >= 0 && c_idx >= 0 {
        let cell = matrix[c_idx as usize][r_idx as usize];
        match field {
            Field::H => {
                if cell.h == cell.e {
                    field = Field::E;
                } else if cell.h == cell.f {
                    field = Field::F;
                } else {
                    alignment.push(if query[r_idx as usize] == target[c_idx as usize] {
                        OPAL_ALIGN_MATCH
                    } else {
                        OPAL_ALIGN_MISMATCH
                    });
                    c_idx -= 1;
                    r_idx -= 1;
                }
            }
            Field::E => {
                field = if c_idx > 0
                    && cell.e == matrix[c_idx as usize - 1][r_idx as usize].h - gap_open
                {
                    Field::H
                } else {
                    Field::E
                };
                alignment.push(OPAL_ALIGN_INS);
                c_idx -= 1;
            }
            Field::F => {
                field = if r_idx > 0
                    && cell.f == matrix[c_idx as usize][r_idx as usize - 1].h - gap_open
                {
                    Field::H
                } else {
                    Field::F
                };
                alignment.push(OPAL_ALIGN_DEL);
                r_idx -= 1;
            }
        }
    }
    while r_idx >= 0 {
        alignment.push(OPAL_ALIGN_DEL);
        r_idx -= 1;
    }
    while c_idx >= 0 {
        alignment.push(OPAL_ALIGN_INS);
        c_idx -= 1;
    }
    assert_eq!(r_idx, -1);
    assert_eq!(c_idx, -1);
    let alignment_length = alignment.len() as i32;
    opal_l1197_revertarray(&mut alignment, alignment_length);
    result.alignment_length = alignment.len() as i32;
    result.alignment = Some(alignment);
}

#[doc = "Original `opalSearchDatabase` at STAR/source/opal/opal.cpp:1437. Args: query: unsigned char, queryLength: int, db: unsigned char, dbLength: int, dbSeqLengths: int, gapOpen: int, gapExt: int, scoreMatrix: int, alphabetLength: int, results: OpalSearchResult, searchType: int, mode: int, overflowMethod: int"]
pub fn opal_l1437_opalsearchdatabase(
    query: &[u8],
    query_length: i32,
    db: &[&[u8]],
    db_length: i32,
    db_seq_lengths: &[i32],
    gap_open: i32,
    gap_ext: i32,
    score_matrix: &[i32],
    alphabet_length: i32,
    results: &mut [crate::opal::opal::OpalSearchResult],
    search_type: i32,
    mode: i32,
    overflow_method: i32,
) -> i32 {
    let mut skip = vec![false; db_length as usize];
    for i in 0..db_length as usize {
        skip[i] = opal_l1561_opalsearchresultisempty(&results[i]) == 0
            && (search_type == OPAL_SEARCH_SCORE
                || (results[i].end_location_query >= 0 && results[i].end_location_target >= 0));
    }

    let status = if mode == OPAL_MODE_SW {
        opal_l498_searchdatabasesw(
            query,
            query_length,
            db,
            db_length,
            db_seq_lengths,
            gap_open,
            gap_ext,
            score_matrix,
            alphabet_length,
            results,
            search_type,
            Some(&skip),
            overflow_method,
        )
    } else if mode == OPAL_MODE_NW || mode == OPAL_MODE_HW || mode == OPAL_MODE_OV {
        for i in 0..db_length as usize {
            if skip[i] {
                continue;
            }
            let target = db[i];
            let target_length = db_seq_lengths[i];
            let lower_score_bound = i32::MIN + gap_open.max(gap_ext);
            let mut initial_column = vec![0i32; query_length as usize];
            let initial_e = vec![lower_score_bound; query_length as usize];
            for r in 0..query_length as usize {
                initial_column[r] = -gap_open - r as i32 * gap_ext;
            }
            let mut prev_h = initial_column;
            let mut prev_e = initial_e;
            let mut h = i32::MIN;
            let mut max_score = i32::MIN;
            let mut max_score_column = -1i32;
            let mut last_column = vec![lower_score_bound; query_length as usize];

            for c in 0..target_length as usize {
                let mut curr_h = vec![lower_score_bound; query_length as usize];
                let mut curr_e = vec![lower_score_bound; query_length as usize];
                let mut u_f = lower_score_bound;
                let mut u_h = -gap_open - c as i32 * gap_ext;
                let mut ul_h = if c == 0 { 0 } else { u_h + gap_ext };

                for r in 0..query_length as usize {
                    let e = (prev_h[r] - gap_open).max(prev_e[r] - gap_ext);
                    let f = (u_h - gap_open).max(u_f - gap_ext);
                    let score = score_matrix
                        [query[r] as usize * alphabet_length as usize + target[c] as usize];
                    h = e.max(f.max(ul_h + score));
                    u_f = f;
                    u_h = h;
                    ul_h = prev_h[r];
                    curr_h[r] = h;
                    curr_e[r] = e;
                }

                if mode == OPAL_MODE_HW || mode == OPAL_MODE_OV {
                    if h > max_score {
                        max_score = h;
                        max_score_column = c as i32;
                    }
                }
                last_column = curr_h.clone();
                prev_h = curr_h;
                prev_e = curr_e;
            }

            if mode == OPAL_MODE_NW {
                opal_l1565_opalsearchresultsetscore(&mut results[i], h);
                if search_type != OPAL_SEARCH_SCORE {
                    results[i].end_location_query = query_length - 1;
                    results[i].end_location_target = target_length - 1;
                }
            } else if mode == OPAL_MODE_HW {
                opal_l1565_opalsearchresultsetscore(&mut results[i], max_score);
                if search_type != OPAL_SEARCH_SCORE {
                    results[i].end_location_query = query_length - 1;
                    results[i].end_location_target = max_score_column;
                }
            } else {
                for r in 0..query_length as usize {
                    if last_column[r] > max_score {
                        max_score = last_column[r];
                    }
                }
                opal_l1565_opalsearchresultsetscore(&mut results[i], max_score);
                if search_type != OPAL_SEARCH_SCORE {
                    let mut row = 0i32;
                    while row < query_length && last_column[row as usize] != max_score {
                        row += 1;
                    }
                    results[i].end_location_query = row;
                    results[i].end_location_target = target_length - 1;
                }
            }
            if search_type != OPAL_SEARCH_SCORE {
                if mode == OPAL_MODE_NW {
                    results[i].end_location_query = query_length - 1;
                    results[i].end_location_target = target_length - 1;
                }
            } else {
                results[i].end_location_query = -1;
                results[i].end_location_target = -1;
            }
        }
        0
    } else {
        OPAL_ERR_INVALID_MODE
    };
    if status != 0 {
        return status;
    }

    if search_type == OPAL_SEARCH_ALIGNMENT {
        let r_query = opal_l1188_createreversecopy(query, query_length);
        for i in 0..db_length as usize {
            if mode == OPAL_MODE_SW && results[i].score == 0 {
                results[i].alignment = None;
                results[i].alignment_length = 0;
                results[i].start_location_query = -1;
                results[i].start_location_target = -1;
                results[i].end_location_query = -1;
                results[i].end_location_target = -1;
            } else {
                let align_query_length = results[i].end_location_query + 1;
                let align_query = &r_query
                    [query_length as usize - align_query_length as usize..query_length as usize];
                let align_target_length = results[i].end_location_target + 1;
                let align_target = opal_l1188_createreversecopy(db[i], align_target_length);
                let mut result = crate::opal::opal::OpalSearchResult::default();
                opal_l1238_findalignment(
                    align_query,
                    align_query_length,
                    &align_target,
                    align_target_length,
                    gap_open,
                    gap_ext,
                    score_matrix,
                    alphabet_length,
                    results[i].score,
                    &mut result,
                    mode,
                );
                assert_eq!(results[i].score, result.score);
                results[i].start_location_query =
                    align_query_length - result.end_location_query - 1;
                results[i].start_location_target =
                    align_target_length - result.end_location_target - 1;
                results[i].alignment_length = result.alignment_length;
                results[i].alignment = result.alignment;
                if let Some(alignment) = results[i].alignment.as_mut() {
                    opal_l1197_revertarray(alignment, results[i].alignment_length);
                }
            }
        }
    } else {
        for i in 0..db_length as usize {
            results[i].alignment = None;
            results[i].alignment_length = -1;
            results[i].start_location_query = -1;
            results[i].start_location_target = -1;
        }
    }

    0
}

#[doc = "Original `opalSearchDatabaseCharSW` at STAR/source/opal/opal.cpp:1526. Args: query: unsigned char, queryLength: int, db: unsigned char, dbLength: int, dbSeqLengths: int, gapOpen: int, gapExt: int, scoreMatrix: int, alphabetLength: int, results: OpalSearchResult"]
pub fn opal_l1526_opalsearchdatabasecharsw(
    query: &[u8],
    query_length: i32,
    db: &[&[u8]],
    db_length: i32,
    db_seq_lengths: &[i32],
    gap_open: i32,
    gap_ext: i32,
    score_matrix: &[i32],
    alphabet_length: i32,
    results: &mut [crate::opal::opal::OpalSearchResult],
) -> i32 {
    let calculated = vec![false; db_length as usize];
    let result_code = opal_l498_searchdatabasesw(
        query,
        query_length,
        db,
        db_length,
        db_seq_lengths,
        gap_open,
        gap_ext,
        score_matrix,
        alphabet_length,
        results,
        OPAL_SEARCH_SCORE,
        Some(&calculated),
        OPAL_OVERFLOW_SIMPLE,
    );
    for i in 0..db_length as usize {
        if !calculated[i] && opal_l1561_opalsearchresultisempty(&results[i]) != 0 {
            results[i].score = -1;
            results[i].score_set = 0;
        }
    }
    result_code
}

#[doc = "Original `opalInitSearchResult` at STAR/source/opal/opal.cpp:1553. Args: result: OpalSearchResult"]
pub fn opal_l1553_opalinitsearchresult(result: &mut crate::opal::opal::OpalSearchResult) {
    result.score_set = 0;
    result.start_location_target = -1;
    result.start_location_query = -1;
    result.end_location_target = -1;
    result.end_location_query = -1;
    result.alignment = None;
    result.alignment_length = 0;
}

#[doc = "Original `opalSearchResultIsEmpty` at STAR/source/opal/opal.cpp:1561. Args: result: OpalSearchResult"]
pub fn opal_l1561_opalsearchresultisempty(
    result: &crate::opal::opal::OpalSearchResult,
) -> i32 {
    if result.score_set == 0 { 1 } else { 0 }
}

#[doc = "Original `opalSearchResultSetScore` at STAR/source/opal/opal.cpp:1565. Args: result: OpalSearchResult, score: int"]
pub fn opal_l1565_opalsearchresultsetscore(
    result: &mut crate::opal::opal::OpalSearchResult,
    score: i32,
) {
    result.score_set = 1;
    result.score = score;
}
