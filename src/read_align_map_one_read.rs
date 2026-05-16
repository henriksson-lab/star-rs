#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::mapOneRead` at STAR/source/ReadAlign_mapOneRead.cpp:6. Args: "]
pub fn readalign_maponeread_l6_readalign_maponeread(
    read_align: &mut crate::read_align::ReadAlign,
    map_gen: &crate::genome::Genome,
    p: &crate::parameters_chimeric::Parameters,
) -> Result<i32, String> {
    read_align.revert_strand = false;

    if read_align.l_read > 0 {
        if read_align.split_r[0].len() < p.max_nsplit as usize {
            read_align.split_r[0].resize(p.max_nsplit as usize, 0);
        }
        if read_align.split_r[1].len() < p.max_nsplit as usize {
            read_align.split_r[1].resize(p.max_nsplit as usize, 0);
        }
        if read_align.split_r[2].len() < p.max_nsplit as usize {
            read_align.split_r[2].resize(p.max_nsplit as usize, 0);
        }
        read_align.n_split = sequencefuns_l411_qualitysplit(
            &read_align.read1[0],
            read_align.l_read,
            p.max_nsplit as u64,
            p.seed_split_min as u64,
            &mut read_align.split_r,
        );
    } else {
        read_align.n_split = 0;
    }
    readalign_l113_readalign_resetn(read_align);
    transcript_l8_transcript_reset(&mut read_align.tr_init);
    read_align.tr_init.chr = 0;
    read_align.tr_init.str_ = 0;
    read_align.tr_init.ro_str = 0;
    read_align.tr_init.c_start = 0;
    read_align.tr_init.g_length = 0;
    read_align.tr_init.i_read = read_align.i_read;
    read_align.tr_init.l_read = read_align.l_read;
    read_align.tr_init.n_exons = 0;
    // Copy per-mate read lengths into the fixed-size array (matches C++ STAR layout).
    for (i, &v) in read_align
        .read_length_original
        .iter()
        .take(crate::include_define::MAX_N_MATES)
        .enumerate()
    {
        read_align.tr_init.read_length_original[i] = v;
    }
    read_align.tr_init.read_length_pair_original = read_align.read_length_pair_original;
    for (i, &v) in read_align
        .read_length
        .iter()
        .take(crate::include_define::MAX_N_MATES)
        .enumerate()
    {
        read_align.tr_init.read_length[i] = v;
    }
    read_align.tr_init.read_nmates = read_align.read_nmates;
    read_align
        .tr_init
        .read_name
        .clone_from(&read_align.read_name);
    // tr_best mirrors tr_init's initial state (C++: `trBest = trInit` pointer-aliasing).
    // Reuse tr_best's existing Vec capacity instead of allocating fresh.
    {
        let crate::read_align::ReadAlign {
            tr_init, tr_best, ..
        } = read_align;
        tr_best.copy_from(tr_init);
    }

    let seed_search_start_lmax = std::cmp::min(
        p.seed_search_start_lmax as u64,
        (p.seed_search_start_lmax_over_lread * read_align.l_read.saturating_sub(1) as f64) as u64,
    );

    let mut stored_aligns = std::mem::take(&mut read_align.stored_aligns_buf);
    stored_aligns.clear();
    for ip in 0..read_align.n_split as usize {
        let n_start =
            if p.seed_search_start_lmax > 0 && seed_search_start_lmax < read_align.split_r[1][ip] {
                read_align.split_r[1][ip] / seed_search_start_lmax + 1
            } else {
                1
            };
        let l_start = read_align.split_r[1][ip] / n_start;
        let mut flag_dir_map = true;

        for i_dir in 0..2_u64 {
            for istart in 0..n_start {
                if flag_dir_map || istart > 0 {
                    let mut l_mapped = 0_u64;
                    while istart * l_start + l_mapped + (p.seed_map_min as u64)
                        < read_align.split_r[1][ip]
                    {
                        let shift = if i_dir == 0 {
                            read_align.split_r[0][ip] + istart * l_start + l_mapped
                        } else {
                            read_align.split_r[0][ip] + read_align.split_r[1][ip]
                                - istart * l_start
                                - 1
                                - l_mapped
                        };
                        let seed_length = read_align.split_r[1][ip] - l_mapped - istart * l_start;
                        let mut l = 0;
                        readalign_maxmappablelength2strands_l5_readalign_maxmappablelength2strands(
                            map_gen,
                            &map_gen.p_ge,
                            [&read_align.read1[0], &read_align.read1[1]],
                            shift,
                            seed_length,
                            i_dir,
                            0,
                            map_gen.n_sa.saturating_sub(1),
                            &mut l,
                            read_align.split_r[2][ip],
                            &mut stored_aligns,
                        )?;
                        for stored in stored_aligns.drain(..) {
                            readalign_storealigns_l10_readalign_storealigns(
                                read_align,
                                p,
                                stored.i_dir,
                                stored.shift,
                                stored.n_rep,
                                stored.l,
                                stored.ind_start_end,
                                stored.i_frag,
                            )?;
                        }
                        if i_dir == 0
                            && istart == 0
                            && l_mapped == 0
                            && shift + l == read_align.split_r[1][ip]
                        {
                            flag_dir_map = false;
                        }
                        l_mapped += l;
                        if l == 0 {
                            break;
                        }
                    }
                }

                if p.seed_search_lmax > 0 {
                    let shift = if i_dir == 0 {
                        read_align.split_r[0][ip] + istart * l_start
                    } else {
                        read_align.split_r[0][ip] + read_align.split_r[1][ip] - istart * l_start - 1
                    };
                    let seed_length = std::cmp::min(
                        p.seed_search_lmax as u64,
                        if i_dir == 0 {
                            read_align.split_r[0][ip] + read_align.split_r[1][ip] - shift
                        } else {
                            shift + 1
                        },
                    );
                    let mut l = 0;
                    readalign_maxmappablelength2strands_l5_readalign_maxmappablelength2strands(
                        map_gen,
                        &map_gen.p_ge,
                        [&read_align.read1[0], &read_align.read1[1]],
                        shift,
                        seed_length,
                        i_dir,
                        0,
                        map_gen.n_sa.saturating_sub(1),
                        &mut l,
                        read_align.split_r[2][ip],
                        &mut stored_aligns,
                    )?;
                    for stored in stored_aligns.drain(..) {
                        readalign_storealigns_l10_readalign_storealigns(
                            read_align,
                            p,
                            stored.i_dir,
                            stored.shift,
                            stored.n_rep,
                            stored.l,
                            stored.ind_start_end,
                            stored.i_frag,
                        )?;
                    }
                }
            }
        }
    }

    if read_align.l_read < p.out_filter_match_nmin as u64 {
        read_align.map_marker = MARKER_READ_TOO_SHORT as u64;
        read_align.tr_best.r_length = 0;
        read_align.n_w = 0;
    } else if read_align.n_split == 0 {
        read_align.map_marker = MARKER_NO_GOOD_PIECES as u64;
        read_align.tr_best.r_length = read_align.split_r[1][0];
        read_align.n_w = 0;
    } else if read_align.n_split > 0 && read_align.n_a == 0 {
        read_align.map_marker = MARKER_ALL_PIECES_EXCEED_SEED_MULTIMAP_NMAX as u64;
        read_align.tr_best.r_length = read_align.mult_nmin_l;
        read_align.n_w = 0;
    } else if read_align.n_split > 0 && read_align.n_a > 0 {
        let read1_ptrs = [
            read_align.read1[0].as_ptr(),
            read_align.read1[1].as_ptr(),
            read_align.read1[2].as_ptr(),
        ];
        let read1_lens = [
            read_align.read1[0].len(),
            read_align.read1[1].len(),
            read_align.read1[2].len(),
        ];
        let r = unsafe {
            [
                std::slice::from_raw_parts(read1_ptrs[0], read1_lens[0]),
                std::slice::from_raw_parts(read1_ptrs[1], read1_lens[1]),
                std::slice::from_raw_parts(read1_ptrs[2], read1_lens[2]),
            ]
        };
        readalign_stitchpieces_l12_readalign_stitchpieces(
            read_align,
            r,
            read_align.l_read,
            map_gen,
            p,
        )?;
    }

    read_align.stored_aligns_buf = stored_aligns;
    Ok(0)
}
