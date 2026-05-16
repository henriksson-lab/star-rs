#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::assignAlignToWindow` at STAR/source/ReadAlign_assignAlignToWindow.cpp:6. Args: a1: uint, aLength: uint, aStr: uint, aNrep: uint, aFrag: uint, aRstart: uint, aAnchor: bool, sjA: uint"]
pub fn readalign_assignaligntowindow_l6_readalign_assignaligntowindow(
    read_align: &mut crate::read_align::ReadAlign,
    a1: u64,
    a_length: u64,
    a_str: u64,
    a_nrep: u64,
    a_frag: u64,
    a_rstart: u64,
    a_anchor: bool,
    sj_a: u64,
    win_bin_nbits: u32,
    seed_per_window_nmax: u32,
) -> Result<(), String> {
    let i_w = read_align.win_bin[a_str as usize][(a1 >> win_bin_nbits) as usize];

    if i_w >= UINT_WIN_BIN_MAX || (!a_anchor && a_length < read_align.wal_rec[i_w as usize]) {
        return Ok(());
    }
    let i_w_usize = i_w as usize;

    let mut i_a = 0usize;
    while i_a < read_align.n_wa[i_w_usize] as usize {
        let wa = read_align.wa[i_w_usize][i_a];
        if a_frag == wa[WA_I_FRAG]
            && wa[WA_SJ_A] == sj_a
            && a1 + wa[WA_R_START] == wa[WA_G_START] + a_rstart
            && ((a_rstart >= wa[WA_R_START] && a_rstart < wa[WA_R_START] + wa[WA_LENGTH])
                || (a_rstart + a_length >= wa[WA_R_START]
                    && a_rstart + a_length < wa[WA_R_START] + wa[WA_LENGTH]))
        {
            break;
        }
        i_a += 1;
    }
    if i_a < read_align.n_wa[i_w_usize] as usize {
        if a_length > read_align.wa[i_w_usize][i_a][WA_LENGTH] {
            let mut i_a0 = 0usize;
            while i_a0 < read_align.n_wa[i_w_usize] as usize {
                if i_a0 != i_a && a_rstart < read_align.wa[i_w_usize][i_a0][WA_R_START] {
                    break;
                }
                i_a0 += 1;
            }
            if i_a0 > i_a {
                i_a0 -= 1;
            }

            if i_a0 < i_a {
                for i_a1 in (i_a0 + 1..=i_a).rev() {
                    read_align.wa[i_w_usize][i_a1] = read_align.wa[i_w_usize][i_a1 - 1];
                }
            } else if i_a0 > i_a {
                for i_a1 in i_a..i_a0 {
                    read_align.wa[i_w_usize][i_a1] = read_align.wa[i_w_usize][i_a1 + 1];
                }
            }

            read_align.wa[i_w_usize][i_a0][WA_R_START] = a_rstart;
            read_align.wa[i_w_usize][i_a0][WA_LENGTH] = a_length;
            read_align.wa[i_w_usize][i_a0][WA_G_START] = a1;
            read_align.wa[i_w_usize][i_a0][WA_N_REP] = a_nrep;
            read_align.wa[i_w_usize][i_a0][WA_ANCHOR] = u64::from(a_anchor);
            read_align.wa[i_w_usize][i_a0][WA_I_FRAG] = a_frag;
            read_align.wa[i_w_usize][i_a0][WA_SJ_A] = sj_a;
        }
        return Ok(());
    }

    if read_align.n_wa[i_w_usize] == seed_per_window_nmax as u64 {
        read_align.wal_rec[i_w_usize] = read_align.l_read + 1;
        for i_a in 0..read_align.n_wa[i_w_usize] as usize {
            if read_align.wa[i_w_usize][i_a][WA_ANCHOR] != 1 {
                read_align.wal_rec[i_w_usize] = std::cmp::min(
                    read_align.wal_rec[i_w_usize],
                    read_align.wa[i_w_usize][i_a][WA_LENGTH],
                );
            }
        }

        if read_align.wal_rec[i_w_usize] == read_align.l_read + 1 {
            read_align.map_marker = MARKER_TOO_MANY_ANCHORS_PER_WINDOW as u64;
            read_align.n_w = 0;
            return Ok(());
        }

        if !a_anchor && a_length < read_align.wal_rec[i_w_usize] {
            return Ok(());
        }

        let mut i_a1 = 0usize;
        for i_a in 0..read_align.n_wa[i_w_usize] as usize {
            if read_align.wa[i_w_usize][i_a][WA_ANCHOR] == 1
                || read_align.wa[i_w_usize][i_a][WA_LENGTH] > read_align.wal_rec[i_w_usize]
            {
                read_align.wa[i_w_usize][i_a1] = read_align.wa[i_w_usize][i_a];
                i_a1 += 1;
            }
        }
        read_align.n_wa[i_w_usize] = i_a1 as u64;

        if !a_anchor && a_length <= read_align.wal_rec[i_w_usize] {
            read_align.n_wap[i_w_usize] = 0;
        }
    }

    if a_anchor || a_length > read_align.wal_rec[i_w_usize] {
        if read_align.n_wa[i_w_usize] >= seed_per_window_nmax as u64 {
            return Err("BUG: iA>=P.seedPerWindowNmax in stitchPieces, exiting".to_string());
        }

        let mut i_a = 0usize;
        while i_a < read_align.n_wa[i_w_usize] as usize {
            if a_rstart < read_align.wa[i_w_usize][i_a][WA_R_START] {
                break;
            }
            i_a += 1;
        }
        for i_a1 in (i_a + 1..=read_align.n_wa[i_w_usize] as usize).rev() {
            read_align.wa[i_w_usize][i_a1] = read_align.wa[i_w_usize][i_a1 - 1];
        }

        read_align.wa[i_w_usize][i_a][WA_R_START] = a_rstart;
        read_align.wa[i_w_usize][i_a][WA_LENGTH] = a_length;
        read_align.wa[i_w_usize][i_a][WA_G_START] = a1;
        read_align.wa[i_w_usize][i_a][WA_N_REP] = a_nrep;
        read_align.wa[i_w_usize][i_a][WA_ANCHOR] = u64::from(a_anchor);
        read_align.wa[i_w_usize][i_a][WA_I_FRAG] = a_frag;
        read_align.wa[i_w_usize][i_a][WA_SJ_A] = sj_a;

        read_align.n_wa[i_w_usize] += 1;
        read_align.n_wap[i_w_usize] += 1;
        if a_anchor && read_align.w_last_anchor[i_w_usize] < i_a as u64 {
            read_align.w_last_anchor[i_w_usize] = i_a as u64;
        }
    }
    Ok(())
}
