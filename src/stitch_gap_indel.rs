#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `stitchGapIndel` at STAR/source/stitchGapIndel.cpp:4. Args: rAend: uint, gAend: uint, rBstart: uint, gBstart: uint, L: uint, gapStart: uint, gapEnd: uint, R: char, G: char, P: Parameters, iRbest: uint, nMM: uint"]
pub fn stitchgapindel_l4_stitchgapindel(
    r_aend: u32,
    g_aend: u32,
    r_bstart: u32,
    g_bstart: u32,
    l: u32,
    gap_start: u32,
    gap_end: u32,
    r: &[u8],
    g: &[u8],
    score_del_base: i32,
    score_del_open: i32,
    i_rbest: &mut u32,
    n_mm: &mut u32,
) -> i32 {
    let gap_length = gap_end - gap_start + 1;
    let in_del =
        (g_bstart - g_aend - 1) as i64 - gap_length as i64 - (r_bstart - r_aend - 1) as i64;

    if in_del == 0 {
        return -1;
    }

    if in_del > 0 {
        let mut score2 = 0;
        let mut score2best = -1;
        *i_rbest = 0;
        for i_r in 1..r_bstart - r_aend {
            let mut i_g1 = g_aend + i_r;
            let mut i_g2 = i_g1 + in_del as u32;
            if i_g1 >= gap_start {
                i_g1 += gap_length;
            }
            if i_g2 >= gap_start {
                i_g2 += gap_length;
            }

            if r[(r_aend + i_r) as usize] == g[i_g1 as usize]
                && r[(r_aend + i_r) as usize] != g[i_g2 as usize]
            {
                score2 += 1;
            } else if r[(r_aend + i_r) as usize] != g[i_g1 as usize]
                && r[(r_aend + i_r) as usize] == g[i_g2 as usize]
            {
                score2 -= 1;
            }

            if score2 > score2best {
                score2best = score2;
                *i_rbest = i_r;
            }
        }

        *n_mm = 0;
        score2 = l as i32 - in_del as i32 * score_del_base - score_del_open;
        for i_r in 1..r_bstart - r_aend {
            let mut i_g = g_aend + i_r;
            if i_r > *i_rbest {
                i_g += in_del as u32;
            }
            if i_g >= gap_start {
                i_g += gap_length;
            }

            if r[(r_aend + i_r) as usize] == g[i_g as usize] {
                score2 += 1;
            } else if r[(r_aend + i_r) as usize] != g[i_g as usize]
                && r[(r_aend + i_r) as usize] < 4
                && g[i_g as usize] < 4
            {
                score2 -= 1;
                *n_mm += 1;
            }
        }
        score2
    } else {
        -1
    }
}
