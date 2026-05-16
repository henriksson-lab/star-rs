#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::stitchPieces` at STAR/source/ReadAlign_stitchPieces.cpp:12. Args: R: char, Lread: uint"]
pub fn readalign_reset_touched_win_bins(
    read_align: &mut crate::read_align::ReadAlign,
    win_bin_n: usize,
) {
    if read_align.win_bin[0].len() < win_bin_n {
        read_align.win_bin[0].resize(win_bin_n, UINT_WIN_BIN_MAX);
    }
    if read_align.win_bin[1].len() < win_bin_n {
        read_align.win_bin[1].resize(win_bin_n, UINT_WIN_BIN_MAX);
    }

    for i in 0..read_align.win_bin_touched.len() {
        let [strand, bin] = read_align.win_bin_touched[i];
        let strand = strand as usize;
        let bin = bin as usize;
        if strand < 2 && bin < read_align.win_bin[strand].len() {
            read_align.win_bin[strand][bin] = UINT_WIN_BIN_MAX;
        }
    }
    read_align.win_bin_touched.clear();
}

pub fn readalign_set_win_bin(
    read_align: &mut crate::read_align::ReadAlign,
    strand: usize,
    bin: u32,
    value: u32,
) {
    let bin_usize = bin as usize;
    if read_align.win_bin[strand][bin_usize] >= UINT_WIN_BIN_MAX {
        read_align.win_bin_touched.push([strand as u32, bin]);
    }
    read_align.win_bin[strand][bin_usize] = value;
}

pub fn readalign_stitchpieces_l12_readalign_stitchpieces(
    read_align: &mut crate::read_align::ReadAlign,
    r: [&[u8]; 3],
    l_read: u64,
    map_gen: &crate::genome::Genome,
    p: &crate::parameters_chimeric::Parameters,
) -> Result<(), String> {
    let win_bin_n = p.win_bin_n as usize;
    readalign_reset_touched_win_bins(read_align, win_bin_n);

    let align_windows = p.align_windows_per_read_nmax as usize;
    let seed_per_window = p.seed_per_window_nmax as usize;
    if read_align.wc.len() < align_windows {
        read_align.wc.resize(align_windows, [0; WC_SIZE]);
    }
    if read_align.n_wa.len() < align_windows {
        read_align.n_wa.resize(align_windows, 0);
    }
    if read_align.n_wap.len() < align_windows {
        read_align.n_wap.resize(align_windows, 0);
    }
    if read_align.wal_rec.len() < align_windows {
        read_align.wal_rec.resize(align_windows, 0);
    }
    if read_align.w_last_anchor.len() < align_windows {
        read_align.w_last_anchor.resize(align_windows, 0);
    }
    if read_align.wa.len() < align_windows {
        read_align.wa.resize(align_windows, Vec::new());
    }
    for i_w in 0..align_windows {
        if read_align.wa[i_w].len() < seed_per_window {
            read_align.wa[i_w].resize(seed_per_window, [0; WA_SIZE]);
        }
    }
    if read_align.wa_incl.len() < seed_per_window {
        read_align.wa_incl.resize(seed_per_window, false);
    }
    if read_align.tr_all.len() < align_windows {
        read_align.tr_all.resize(align_windows, Vec::new());
    }
    if read_align.n_win_tr.len() < align_windows {
        read_align.n_win_tr.resize(align_windows, 0);
    }

    read_align.n_w = 0;
    for i_p in 0..read_align.n_p as usize {
        if read_align.pc[i_p][PC_NREP] <= p.win_anchor_multimap_nmax as u64 {
            let a_dir = read_align.pc[i_p][PC_DIR] as u32;
            let a_length = read_align.pc[i_p][PC_LENGTH] as u32;

            for i_sa in read_align.pc[i_p][PC_SASTART]..=read_align.pc[i_p][PC_SAEND] {
                let mut a1 = genome_sa_index_value(map_gen, i_sa);
                let mut a_str = (a1 >> map_gen.gstrand_bit) as u32;
                a1 &= map_gen.gstrand_mask as u64;

                if a_dir == 1 && a_str == 0 {
                    a_str = 1;
                } else if a_dir == 0 && a_str == 1 {
                    let Some(converted) = map_gen.n_genome.checked_sub(a_length as u64 + a1) else {
                        continue;
                    };
                    a1 = converted;
                } else if a_dir == 1 && a_str == 1 {
                    a_str = 0;
                    let Some(converted) = map_gen.n_genome.checked_sub(a_length as u64 + a1) else {
                        continue;
                    };
                    a1 = converted;
                }

                if read_align.revert_strand {
                    a_str = 1 - a_str;
                }

                if a1 >= map_gen.sj_gstart {
                    let mut a1_d = 0;
                    let mut a_length_d = 0;
                    let mut a1_a = 0;
                    let mut a_length_a = 0;
                    let mut sj1 = 0;
                    if sjalignsplit_l3_sjalignsplit(
                        a1,
                        a_length as u64,
                        map_gen,
                        &mut a1_d,
                        &mut a_length_d,
                        &mut a1_a,
                        &mut a_length_a,
                        &mut sj1,
                    ) {
                        let add_status =
                            readalign_createextendwindowswithalign_l7_readalign_createextendwindowswithalign(
                                read_align,
                                map_gen,
                                p,
                                a1_d as u32,
                                a_str as u32,
                            );
                        if add_status == EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS {
                            break;
                        }
                        let add_status =
                            readalign_createextendwindowswithalign_l7_readalign_createextendwindowswithalign(
                                read_align,
                                map_gen,
                                p,
                                a1_a as u32,
                                a_str as u32,
                            );
                        if add_status == EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS {
                            break;
                        }
                    }
                } else {
                    let add_status =
                        readalign_createextendwindowswithalign_l7_readalign_createextendwindowswithalign(
                            read_align,
                            map_gen,
                            p,
                            a1 as u32,
                            a_str as u32,
                        );
                    if add_status == EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS {
                        break;
                    }
                }
            }
        }
    }

    for i_win in 0..read_align.n_w as usize {
        if read_align.wc[i_win][WC_G_START] <= read_align.wc[i_win][WC_G_END] {
            let mut wb = read_align.wc[i_win][WC_G_START];
            let mut ii = 0;
            while ii < p.win_flank_nbins
                && wb > 0
                && map_gen.chr_bin[((wb - 1) >> p.win_bin_chr_nbits) as usize] as u64
                    == read_align.wc[i_win][WC_CHR]
            {
                wb -= 1;
                readalign_set_win_bin(
                    read_align,
                    read_align.wc[i_win][WC_STR] as usize,
                    wb as u32,
                    i_win as u32,
                );
                ii += 1;
            }
            read_align.wc[i_win][WC_G_START] = wb;

            wb = read_align.wc[i_win][WC_G_END];
            ii = 0;
            while ii < p.win_flank_nbins
                && wb + 1 < p.win_bin_n as u64
                && map_gen.chr_bin[((wb + 1) >> p.win_bin_chr_nbits) as usize] as u64
                    == read_align.wc[i_win][WC_CHR]
            {
                wb += 1;
                readalign_set_win_bin(
                    read_align,
                    read_align.wc[i_win][WC_STR] as usize,
                    wb as u32,
                    i_win as u32,
                );
                ii += 1;
            }
            read_align.wc[i_win][WC_G_END] = wb;
        }
        read_align.n_wa[i_win] = 0;
        read_align.wal_rec[i_win] = 0;
        read_align.w_last_anchor[i_win] = u64::MAX;
    }
    read_align.n_wall = read_align.n_w;

    for i_p in 0..read_align.n_p as usize {
        let a_nrep = read_align.pc[i_p][PC_NREP];
        let a_frag = read_align.pc[i_p][PC_IFRAG];
        let a_length = read_align.pc[i_p][PC_LENGTH];
        let a_dir = read_align.pc[i_p][PC_DIR];
        let a_anchor = read_align.pc[i_p][PC_NREP] <= p.win_anchor_multimap_nmax as u64;

        for ii in 0..read_align.n_w as usize {
            read_align.n_wap[ii] = 0;
        }

        for i_sa in read_align.pc[i_p][PC_SASTART]..=read_align.pc[i_p][PC_SAEND] {
            let mut a1 = genome_sa_index_value(map_gen, i_sa);
            let mut a_str = a1 >> map_gen.gstrand_bit;
            a1 &= map_gen.gstrand_mask as u64;
            let mut a_rstart = read_align.pc[i_p][PC_R_START];

            if a_dir == 1 && a_str == 0 {
                a_str = 1;
                let Some(converted_rstart) = l_read.checked_sub(a_length + a_rstart) else {
                    continue;
                };
                a_rstart = converted_rstart;
            } else if a_dir == 0 && a_str == 1 {
                let Some(converted_rstart) = l_read.checked_sub(a_length + a_rstart) else {
                    continue;
                };
                let Some(converted) = map_gen.n_genome.checked_sub(a_length + a1) else {
                    continue;
                };
                a_rstart = converted_rstart;
                a1 = converted;
            } else if a_dir == 1 && a_str == 1 {
                a_str = 0;
                let Some(converted) = map_gen.n_genome.checked_sub(a_length + a1) else {
                    continue;
                };
                a1 = converted;
            }

            if read_align.revert_strand {
                a_str = 1 - a_str;
            }

            if a1 >= map_gen.sj_gstart {
                let mut a1_d = 0;
                let mut a_length_d = 0;
                let mut a1_a = 0;
                let mut a_length_a = 0;
                let mut isj1 = 0;
                if sjalignsplit_l3_sjalignsplit(
                    a1,
                    a_length,
                    map_gen,
                    &mut a1_d,
                    &mut a_length_d,
                    &mut a1_a,
                    &mut a_length_a,
                    &mut isj1,
                ) {
                    readalign_assignaligntowindow_l6_readalign_assignaligntowindow(
                        read_align,
                        a1_d,
                        a_length_d,
                        a_str,
                        a_nrep,
                        a_frag,
                        a_rstart,
                        a_anchor,
                        isj1,
                        p.win_bin_nbits,
                        p.seed_per_window_nmax,
                    )?;
                    readalign_assignaligntowindow_l6_readalign_assignaligntowindow(
                        read_align,
                        a1_a,
                        a_length_a,
                        a_str,
                        a_nrep,
                        a_frag,
                        a_rstart + a_length_d,
                        a_anchor,
                        isj1,
                        p.win_bin_nbits,
                        p.seed_per_window_nmax,
                    )?;
                } else {
                    continue;
                }
            } else {
                readalign_assignaligntowindow_l6_readalign_assignaligntowindow(
                    read_align,
                    a1,
                    a_length,
                    a_str,
                    a_nrep,
                    a_frag,
                    a_rstart,
                    a_anchor,
                    u64::MAX,
                    p.win_bin_nbits,
                    p.seed_per_window_nmax,
                )?;
            }
        }
    }
    read_align.tr_best = (*read_align.tr_init).clone();
    let mut i_w1 = 0usize;
    let mut tr_n_total = 0u64;

    for i_w in 0..read_align.n_w as usize {
        if read_align.n_wa[i_w] == 0 {
            continue;
        }

        if read_align.w_last_anchor[i_w] < read_align.n_wa[i_w] {
            let last_anchor = read_align.w_last_anchor[i_w] as usize;
            read_align.wa[i_w][last_anchor][WA_ANCHOR] = 2;
        }

        for ii in 0..read_align.n_wa[i_w] as usize {
            read_align.wa_incl[ii] = false;
        }

        let mut tr_a = (*read_align.tr_init).clone();
        tr_a.chr = read_align.wc[i_w][WC_CHR];
        tr_a.str_ = read_align.wc[i_w][WC_STR];
        tr_a.ro_str = if read_align.revert_strand {
            1 - tr_a.str_
        } else {
            tr_a.str_
        };
        tr_a.max_score = 0;

        if tr_n_total + p.align_transcripts_per_window_nmax as u64
            >= p.align_transcripts_per_read_nmax as u64
        {
            break;
        }

        let mut w_tr = std::mem::take(&mut read_align.tr_all[i_w1]);
        w_tr.clear();
        w_tr.push(crate::transcript::Transcript::default());
        let mut n_win_tr = 0u64;
        let read_index = if tr_a.ro_str == 0 { 0 } else { 2 };
        let read_bases: &[u8] = if read_index < r.len() {
            &r[read_index]
        } else if tr_a.ro_str as usize % 2 < r.len() {
            &r[tr_a.ro_str as usize % 2]
        } else {
            &[]
        };
        let wa_row = &read_align.wa[i_w];
        let wa_incl = &mut read_align.wa_incl;
        let out_filter_mismatch_nmax_total = read_align.out_filter_mismatch_nmax_total;
        let read_length = &read_align.read_length;
        let max_score_mate = &mut read_align.max_score_mate;
        stitchwindowaligns_l8_stitchwindowaligns(
            0,
            read_align.n_wa[i_w],
            0,
            wa_incl,
            0,
            0,
            tr_a,
            l_read,
            wa_row,
            read_bases,
            map_gen,
            p,
            &mut w_tr,
            &mut n_win_tr,
            out_filter_mismatch_nmax_total,
            read_length,
            max_score_mate,
        )?;

        if n_win_tr == 0 {
            continue;
        }

        if w_tr[0].max_score > read_align.tr_best.max_score
            || (w_tr[0].max_score == read_align.tr_best.max_score
                && w_tr[0].g_length < read_align.tr_best.g_length)
        {
            read_align.tr_best = w_tr[0].clone();
        }

        if read_align.tr_all.len() <= i_w1 {
            read_align.tr_all.resize(i_w1 + 1, Vec::new());
        }
        read_align.tr_all[i_w1] = w_tr;
        if read_align.n_win_tr.len() <= i_w1 {
            read_align.n_win_tr.resize(i_w1 + 1, 0);
        }
        read_align.n_win_tr[i_w1] = n_win_tr;
        tr_n_total += n_win_tr;
        i_w1 += 1;
    }

    read_align.n_w = i_w1 as u64;
    read_align.n_tr = tr_n_total;

    if read_align.tr_best.max_score == 0 {
        read_align.map_marker = MARKER_NO_GOOD_WINDOW as u64;
        read_align.n_w = 0;
    }

    Ok(())
}
