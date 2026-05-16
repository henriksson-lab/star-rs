#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::createExtendWindowsWithAlign` at STAR/source/ReadAlign_createExtendWindowsWithAlign.cpp:7. Args: a1: uint, aStr: uint"]
pub fn readalign_createextendwindowswithalign_l7_readalign_createextendwindowswithalign(
    read_align: &mut crate::read_align::ReadAlign,
    map_gen: &crate::genome::Genome,
    p: &crate::parameters_chimeric::Parameters,
    a1: u32,
    a_str: u32,
) -> i32 {
    let a_bin = a1 >> p.win_bin_nbits;
    let mut i_bin_left = a_bin;
    let mut i_bin_right = a_bin;
    let strand = a_str as usize;
    let mut i_win = u32::MAX;
    let mut i_win_right = u32::MAX;

    if read_align.win_bin[strand][a_bin as usize] >= UINT_WIN_BIN_MAX {
        let mut flag_merge_left = false;
        if a_bin > 0 {
            let mut i_bin;
            let left_limit = if a_bin > p.win_anchor_dist_nbins {
                a_bin - p.win_anchor_dist_nbins
            } else {
                0
            };
            i_bin = a_bin - 1;
            while i_bin >= left_limit {
                if read_align.win_bin[strand][i_bin as usize] < UINT_WIN_BIN_MAX {
                    flag_merge_left = true;
                    break;
                }
                if i_bin == left_limit {
                    break;
                }
                i_bin -= 1;
            }
            flag_merge_left = flag_merge_left
                && map_gen.chr_bin[(i_bin >> p.win_bin_chr_nbits) as usize]
                    == map_gen.chr_bin[(a_bin >> p.win_bin_chr_nbits) as usize];
            if flag_merge_left {
                i_win = read_align.win_bin[strand][i_bin as usize];
                i_bin_left = read_align.wc[i_win as usize][WC_G_START];
                for ii in i_bin + 1..=a_bin {
                    readalign_set_win_bin(read_align, strand, ii, i_win);
                }
            }
        }

        let mut flag_merge_right = false;
        if a_bin + 1 < p.win_bin_n {
            let mut i_bin;
            let right_limit = std::cmp::min(a_bin + p.win_anchor_dist_nbins + 1, p.win_bin_n);
            i_bin = a_bin + 1;
            while i_bin < right_limit {
                if read_align.win_bin[strand][i_bin as usize] < UINT_WIN_BIN_MAX {
                    flag_merge_right = true;
                    break;
                }
                i_bin += 1;
            }

            flag_merge_right = flag_merge_right
                && map_gen.chr_bin[(i_bin >> p.win_bin_chr_nbits) as usize]
                    == map_gen.chr_bin[(a_bin >> p.win_bin_chr_nbits) as usize];
            if flag_merge_right {
                while i_bin + 1 < p.win_bin_n
                    && read_align.win_bin[strand][i_bin as usize]
                        == read_align.win_bin[strand][(i_bin + 1) as usize]
                {
                    i_bin += 1;
                }
                i_bin_right = i_bin;
                i_win_right = read_align.win_bin[strand][i_bin as usize];
                if !flag_merge_left {
                    i_win = read_align.win_bin[strand][i_bin as usize];
                }
                for ii in a_bin..=i_bin {
                    readalign_set_win_bin(read_align, strand, ii, i_win);
                }
            }
        }

        if !flag_merge_left && !flag_merge_right {
            i_win = read_align.n_w;
            readalign_set_win_bin(read_align, strand, a_bin, i_win);
            let i_win_usize = i_win as usize;
            if read_align.wc.len() <= i_win_usize {
                read_align.wc.resize(i_win_usize + 1, [0; WC_SIZE]);
            }
            read_align.wc[i_win_usize][WC_CHR] =
                map_gen.chr_bin[(a_bin >> p.win_bin_chr_nbits) as usize];
            read_align.wc[i_win_usize][WC_STR] = a_str;
            read_align.wc[i_win_usize][WC_G_END] = a_bin;
            read_align.wc[i_win_usize][WC_G_START] = a_bin;
            read_align.n_w += 1;
            if read_align.n_w >= p.align_windows_per_read_nmax {
                read_align.n_w = p.align_windows_per_read_nmax - 1;
                return EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS;
            }
        } else {
            let i_win_usize = i_win as usize;
            read_align.wc[i_win_usize][WC_G_START] = i_bin_left;
            read_align.wc[i_win_usize][WC_G_END] = i_bin_right;
            if flag_merge_left && flag_merge_right {
                let i_win_right_usize = i_win_right as usize;
                read_align.wc[i_win_right_usize][WC_G_START] = 1;
                read_align.wc[i_win_right_usize][WC_G_END] = 0;
            }
        }
    }
    0
}
