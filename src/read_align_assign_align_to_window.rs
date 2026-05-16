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
) {
    let i_w = unsafe {
        *read_align
            .win_bin
            .get_unchecked(a_str as usize)
            .get_unchecked((a1 >> win_bin_nbits) as usize)
    };

    if i_w >= UINT_WIN_BIN_MAX || (!a_anchor && a_length < read_align.wal_rec[i_w as usize]) {
        return;
    }
    let i_w_usize = i_w as usize;
    let wa = &mut read_align.wa[i_w_usize];
    let mut n_wa = read_align.n_wa[i_w_usize] as usize;
    let mut n_wap = read_align.n_wap[i_w_usize];
    let mut wal_rec = read_align.wal_rec[i_w_usize];
    let new_wa = [
        a_length,
        a_rstart,
        a1,
        a_nrep,
        u64::from(a_anchor),
        a_frag,
        sj_a,
    ];

    let mut i_a = 0usize;
    while i_a < n_wa {
        let wai = unsafe { *wa.get_unchecked(i_a) };
        if a_frag == wai[WA_I_FRAG]
            && wai[WA_SJ_A] == sj_a
            && a1 + wai[WA_R_START] == wai[WA_G_START] + a_rstart
            && ((a_rstart >= wai[WA_R_START] && a_rstart < wai[WA_R_START] + wai[WA_LENGTH])
                || (a_rstart + a_length >= wai[WA_R_START]
                    && a_rstart + a_length < wai[WA_R_START] + wai[WA_LENGTH]))
        {
            break;
        }
        i_a += 1;
    }
    if i_a < n_wa {
        if a_length > wa[i_a][WA_LENGTH] {
            let mut i_a0 = 0usize;
            while i_a0 < n_wa {
                if i_a0 != i_a && a_rstart < unsafe { wa.get_unchecked(i_a0) }[WA_R_START] {
                    break;
                }
                i_a0 += 1;
            }
            if i_a0 > i_a {
                i_a0 -= 1;
            }

            if i_a0 < i_a {
                for i_a1 in (i_a0 + 1..=i_a).rev() {
                    unsafe {
                        *wa.get_unchecked_mut(i_a1) = *wa.get_unchecked(i_a1 - 1);
                    }
                }
            } else if i_a0 > i_a {
                for i_a1 in i_a..i_a0 {
                    unsafe {
                        *wa.get_unchecked_mut(i_a1) = *wa.get_unchecked(i_a1 + 1);
                    }
                }
            }

            unsafe {
                *wa.get_unchecked_mut(i_a0) = new_wa;
            }
        }
        return;
    }

    if n_wa == seed_per_window_nmax as usize {
        wal_rec = read_align.l_read + 1;
        for i_a in 0..n_wa {
            let row = unsafe { wa.get_unchecked(i_a) };
            if row[WA_ANCHOR] != 1 {
                wal_rec = std::cmp::min(wal_rec, row[WA_LENGTH]);
            }
        }

        if wal_rec == read_align.l_read + 1 {
            read_align.wal_rec[i_w_usize] = wal_rec;
            read_align.map_marker = MARKER_TOO_MANY_ANCHORS_PER_WINDOW as u64;
            read_align.n_w = 0;
            return;
        }

        if !a_anchor && a_length < wal_rec {
            read_align.wal_rec[i_w_usize] = wal_rec;
            return;
        }

        let mut i_a1 = 0usize;
        for i_a in 0..n_wa {
            let row = unsafe { *wa.get_unchecked(i_a) };
            if row[WA_ANCHOR] == 1 || row[WA_LENGTH] > wal_rec {
                unsafe {
                    *wa.get_unchecked_mut(i_a1) = row;
                }
                i_a1 += 1;
            }
        }
        n_wa = i_a1;

        if !a_anchor && a_length <= wal_rec {
            n_wap = 0;
        }
    }

    if a_anchor || a_length > wal_rec {
        if n_wa >= seed_per_window_nmax as usize {
            panic!("BUG: iA>=P.seedPerWindowNmax in stitchPieces, exiting");
        }

        let mut i_a = 0usize;
        while i_a < n_wa {
            if a_rstart < unsafe { wa.get_unchecked(i_a) }[WA_R_START] {
                break;
            }
            i_a += 1;
        }
        for i_a1 in (i_a + 1..=n_wa).rev() {
            unsafe {
                *wa.get_unchecked_mut(i_a1) = *wa.get_unchecked(i_a1 - 1);
            }
        }

        unsafe {
            *wa.get_unchecked_mut(i_a) = new_wa;
        }

        n_wa += 1;
        if n_wap == 0 {
            read_align.n_wap_touched.push(i_w_usize);
        }
        n_wap += 1;
        if a_anchor && read_align.w_last_anchor[i_w_usize] < i_a as u64 {
            read_align.w_last_anchor[i_w_usize] = i_a as u64;
        }
    }
    read_align.n_wa[i_w_usize] = n_wa as u64;
    read_align.n_wap[i_w_usize] = n_wap;
    read_align.wal_rec[i_w_usize] = wal_rec;
}
