#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `extendAlign` at STAR/source/extendAlign.cpp:6. Args: R: char, G: char, rStart: uint, gStart: uint, dR: int, dG: int, L: uint, Lprev: uint, nMMprev: uint, nMMmax: uint, pMMmax: double, extendToEnd: bool, trA: Transcript"]
pub fn extendalign_l6_extendalign(
    r: &[u8],
    g: &[u8],
    r_start: u64,
    g_start: u64,
    d_r: i32,
    d_g: i32,
    l: u64,
    l_prev: u64,
    n_mm_prev: u64,
    n_mm_max: u64,
    p_mm_max: f64,
    extend_to_end: bool,
    tr_a: &mut crate::transcript::Transcript,
) -> bool {
    let mut score = 0i32;
    let mut n_match = 0u64;
    let mut n_mm = 0u64;
    tr_a.max_score = 0;

    if (d_r != 1 && d_r != -1) || (d_g != 1 && d_g != -1) {
        return extendalign_indexed_fallback(
            r,
            g,
            r_start,
            g_start,
            d_r,
            d_g,
            l,
            l_prev,
            n_mm_prev,
            n_mm_max,
            p_mm_max,
            extend_to_end,
            tr_a,
        );
    }

    let r_start_usize = r_start as usize;
    let g_start_usize = g_start as usize;

    let r_safe = if r_start_usize >= r.len() {
        0
    } else if d_r > 0 {
        (r.len() - r_start_usize) as u64
    } else {
        r_start + 1
    };
    let g_safe = if g_start_usize >= g.len() {
        0
    } else if d_g > 0 {
        (g.len() - g_start_usize) as u64
    } else {
        g_start + 1
    };

    if extend_to_end {
        let mut i_ext = 0u64;
        let l_safe = l.min(r_safe).min(g_safe);
        let mut r_pos = r_start_usize;
        let mut g_pos = g_start_usize;
        while i_ext < l {
            if i_ext >= g_safe || unsafe { *g.get_unchecked(g_pos) } == GENOME_SPACING_CHAR {
                tr_a.extend_l = 0;
                tr_a.max_score = -999_999_999;
                tr_a.n_match = 0;
                tr_a.n_mm = n_mm_max + 1;
                return true;
            }
            if i_ext >= l_safe {
                break;
            }
            let rb = unsafe { *r.get_unchecked(r_pos) };
            if rb == MARK_FRAG_SPACER_BASE {
                break;
            }
            let gb = unsafe { *g.get_unchecked(g_pos) };
            if rb > 3 || gb > 3 {
                i_ext += 1;
                if d_r > 0 {
                    r_pos += 1;
                } else {
                    r_pos = r_pos.wrapping_sub(1);
                }
                if d_g > 0 {
                    g_pos += 1;
                } else {
                    g_pos = g_pos.wrapping_sub(1);
                }
                continue;
            }

            if gb == rb {
                n_match += 1;
                score += 1;
            } else {
                n_mm += 1;
                score -= 1;
            }
            i_ext += 1;
            if d_r > 0 {
                r_pos += 1;
            } else {
                r_pos = r_pos.wrapping_sub(1);
            }
            if d_g > 0 {
                g_pos += 1;
            } else {
                g_pos = g_pos.wrapping_sub(1);
            }
        }

        if i_ext > 0 {
            tr_a.extend_l = i_ext;
            tr_a.max_score = score;
            tr_a.n_match = n_match;
            tr_a.n_mm = n_mm;
            return true;
        }
        return false;
    }

    let l_safe = l.min(r_safe).min(g_safe);
    let mut r_pos = r_start_usize;
    let mut g_pos = g_start_usize;
    for i in 0..l_safe {
        let rb = unsafe { *r.get_unchecked(r_pos) };
        let gb = unsafe { *g.get_unchecked(g_pos) };
        if gb == GENOME_SPACING_CHAR || rb == MARK_FRAG_SPACER_BASE {
            break;
        }
        if rb > 3 || gb > 3 {
            if d_r > 0 {
                r_pos += 1;
            } else {
                r_pos = r_pos.wrapping_sub(1);
            }
            if d_g > 0 {
                g_pos += 1;
            } else {
                g_pos = g_pos.wrapping_sub(1);
            }
            continue;
        }

        if gb == rb {
            n_match += 1;
            score += 1;
            if score > tr_a.max_score
                && (n_mm + n_mm_prev) as f64
                    <= (p_mm_max * (l_prev + i + 1) as f64).min(n_mm_max as f64)
            {
                tr_a.extend_l = i + 1;
                tr_a.max_score = score;
                tr_a.n_match = n_match;
                tr_a.n_mm = n_mm;
            }
        } else {
            if (n_mm + n_mm_prev) as f64 >= (p_mm_max * (l_prev + l) as f64).min(n_mm_max as f64) {
                break;
            }
            n_mm += 1;
            score -= 1;
        }
        if d_r > 0 {
            r_pos += 1;
        } else {
            r_pos = r_pos.wrapping_sub(1);
        }
        if d_g > 0 {
            g_pos += 1;
        } else {
            g_pos = g_pos.wrapping_sub(1);
        }
    }

    tr_a.extend_l > 0
}

fn extendalign_indexed_fallback(
    r: &[u8],
    g: &[u8],
    r_start: u64,
    g_start: u64,
    d_r: i32,
    d_g: i32,
    l: u64,
    l_prev: u64,
    n_mm_prev: u64,
    n_mm_max: u64,
    p_mm_max: f64,
    extend_to_end: bool,
    tr_a: &mut crate::transcript::Transcript,
) -> bool {
    let mut score = 0i32;
    let mut n_match = 0u64;
    let mut n_mm = 0u64;
    tr_a.max_score = 0;

    if extend_to_end {
        let mut i_ext = 0u64;
        while i_ext < l {
            let i_s = d_r as i64 * i_ext as i64;
            let i_g = d_g as i64 * i_ext as i64;
            let r_pos = r_start as i64 + i_s;
            let g_pos = g_start as i64 + i_g;

            if g_pos < 0 || g_pos as usize >= g.len() || g[g_pos as usize] == GENOME_SPACING_CHAR {
                tr_a.extend_l = 0;
                tr_a.max_score = -999_999_999;
                tr_a.n_match = 0;
                tr_a.n_mm = n_mm_max + 1;
                return true;
            }
            if r_pos < 0 || r_pos as usize >= r.len() {
                break;
            }
            if r[r_pos as usize] == MARK_FRAG_SPACER_BASE {
                break;
            }
            if r[r_pos as usize] > 3 || g[g_pos as usize] > 3 {
                i_ext += 1;
                continue;
            }

            if g[g_pos as usize] == r[r_pos as usize] {
                n_match += 1;
                score += 1;
            } else {
                n_mm += 1;
                score -= 1;
            }
            i_ext += 1;
        }

        if i_ext > 0 {
            tr_a.extend_l = i_ext;
            tr_a.max_score = score;
            tr_a.n_match = n_match;
            tr_a.n_mm = n_mm;
            return true;
        }
        return false;
    }

    for i in 0..l {
        let i_s = d_r as i64 * i as i64;
        let i_g = d_g as i64 * i as i64;
        let r_pos = r_start as i64 + i_s;
        let g_pos = g_start as i64 + i_g;

        if g_pos < 0
            || g_pos as usize >= g.len()
            || r_pos < 0
            || r_pos as usize >= r.len()
            || g[g_pos as usize] == GENOME_SPACING_CHAR
            || r[r_pos as usize] == MARK_FRAG_SPACER_BASE
        {
            break;
        }
        if r[r_pos as usize] > 3 || g[g_pos as usize] > 3 {
            continue;
        }

        if g[g_pos as usize] == r[r_pos as usize] {
            n_match += 1;
            score += 1;
            if score > tr_a.max_score
                && (n_mm + n_mm_prev) as f64
                    <= (p_mm_max * (l_prev + i + 1) as f64).min(n_mm_max as f64)
            {
                tr_a.extend_l = i + 1;
                tr_a.max_score = score;
                tr_a.n_match = n_match;
                tr_a.n_mm = n_mm;
            }
        } else {
            if (n_mm + n_mm_prev) as f64 >= (p_mm_max * (l_prev + l) as f64).min(n_mm_max as f64) {
                break;
            }
            n_mm += 1;
            score -= 1;
        }
    }

    tr_a.extend_l > 0
}
