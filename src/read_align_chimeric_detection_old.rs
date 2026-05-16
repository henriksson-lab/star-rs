#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::chimericDetectionOld` at STAR/source/ReadAlign_chimericDetectionOld.cpp:7. Args: "]
pub fn readalign_chimericdetectionold_l7_readalign_chimericdetectionold(
    read_align: &mut crate::read_align::ReadAlign,
    p: &crate::parameters_chimeric::Parameters,
    map_gen: &crate::genome::Genome,
) -> bool {
    // Defer the tr_best clone until after the cheap early-exit checks below.
    // Most reads fail the segment_min check at lines 14-29 and never reach the
    // point where we need an owned copy.
    if read_align.n_tr > p.p_ch.main_segment_mult_nmax as u64 && read_align.n_tr != 2 {
        return false;
    }

    {
        let tr_best = &read_align.tr_best;
        if !(p.p_ch.segment_min > 0
            && tr_best.r_length >= p.p_ch.segment_min as u64
            && (tr_best.exons[tr_best.n_exons as usize - 1][EX_R]
                + tr_best.exons[tr_best.n_exons as usize - 1][EX_L]
                + p.p_ch.segment_min as u64
                <= read_align.l_read
                || tr_best.exons[0][EX_R] >= p.p_ch.segment_min as u64)
            && tr_best.intron_motifs[0] == 0
            && (tr_best.intron_motifs[1] == 0 || tr_best.intron_motifs[2] == 0))
        {
            return false;
        }
    }

    let tr_best = read_align.tr_best.clone();
    let mut chim_score_best = 0_i32;
    let mut chim_score_next = 0_i32;
    if read_align.tr_chim.len() < 2 {
        read_align
            .tr_chim
            .resize(2, crate::transcript::Transcript::default());
    }
    read_align.tr_chim[0].clone_from(&tr_best);
    let mut tr_chim1 = crate::transcript::Transcript::default();
    let mut tr_chim1_set = false;

    let mut ro_start1 = if tr_best.str_ == 0 {
        tr_best.exons[0][EX_R]
    } else {
        read_align.l_read
            - tr_best.exons[tr_best.n_exons as usize - 1][EX_R]
            - tr_best.exons[tr_best.n_exons as usize - 1][EX_L]
    };
    let mut ro_end1 = if tr_best.str_ == 0 {
        tr_best.exons[tr_best.n_exons as usize - 1][EX_R]
            + tr_best.exons[tr_best.n_exons as usize - 1][EX_L]
            - 1
    } else {
        read_align.l_read - tr_best.exons[0][EX_R] - 1
    };
    let read_len0 = read_align.read_length.first().copied().unwrap_or(0);
    let read_len1 = read_align.read_length.get(1).copied().unwrap_or(0);
    if ro_start1 > read_len0 {
        ro_start1 -= 1;
    }
    if ro_end1 > read_len0 {
        ro_end1 -= 1;
    }

    let mut chim_str_best = 0_u64;
    read_align.chim_str = if tr_best.intron_motifs[1] == 0 && tr_best.intron_motifs[2] == 0 {
        0
    } else if (tr_best.str_ == 0) == (tr_best.intron_motifs[1] > 0) {
        1
    } else {
        2
    };

    for i_w in 0..read_align.n_w as usize {
        let n_win_tr = read_align
            .n_win_tr
            .get(i_w)
            .copied()
            .unwrap_or_else(|| read_align.tr_all.get(i_w).map_or(0, |w| w.len() as u64))
            as usize;
        for i_wt in 0..n_win_tr {
            let Some(tr) = read_align
                .tr_all
                .get(i_w)
                .and_then(|w| w.get(i_wt))
                .cloned()
            else {
                continue;
            };
            let best_is_window_first = read_align
                .tr_all
                .get(i_w)
                .and_then(|w| w.first())
                .is_some_and(|first| *first == tr_best);
            if !best_is_window_first && i_wt > 0 {
                break;
            }
            if best_is_window_first && i_wt == 0 {
                continue;
            }
            if tr.intron_motifs[0] > 0 {
                continue;
            }
            let chim_str1 = if tr.intron_motifs[1] == 0 && tr.intron_motifs[2] == 0 {
                0
            } else if (tr.str_ == 0) == (tr.intron_motifs[1] > 0) {
                1
            } else {
                2
            };
            if read_align.chim_str != 0 && chim_str1 != 0 && read_align.chim_str != chim_str1 {
                continue;
            }

            let mut ro_start2 = if tr.str_ == 0 {
                tr.exons[0][EX_R]
            } else {
                read_align.l_read
                    - tr.exons[tr.n_exons as usize - 1][EX_R]
                    - tr.exons[tr.n_exons as usize - 1][EX_L]
            };
            let mut ro_end2 = if tr.str_ == 0 {
                tr.exons[tr.n_exons as usize - 1][EX_R] + tr.exons[tr.n_exons as usize - 1][EX_L]
                    - 1
            } else {
                read_align.l_read - tr.exons[0][EX_R] - 1
            };
            if ro_start2 > read_len0 {
                ro_start2 -= 1;
            }
            if ro_end2 > read_len0 {
                ro_end2 -= 1;
            }

            let chim_overlap = if ro_start2 > ro_start1 {
                if ro_start2 > ro_end1 {
                    0
                } else {
                    ro_end1 - ro_start2 + 1
                }
            } else if ro_end2 < ro_start1 {
                0
            } else {
                ro_end2 - ro_start1 + 1
            };
            let diff_mates = (ro_end1 < read_len0 && ro_start2 >= read_len0)
                || (ro_end2 < read_len0 && ro_start1 >= read_len0);

            if ro_end1 > p.p_ch.segment_min as u64 + ro_start1 + chim_overlap
                && ro_end2 > p.p_ch.segment_min as u64 + ro_start2 + chim_overlap
                && (diff_mates
                    || ((ro_end1 + p.p_ch.segment_read_gap_max as u64 + 1) >= ro_start2
                        && (ro_end2 + p.p_ch.segment_read_gap_max as u64 + 1) >= ro_start1))
            {
                let chim_score = tr_best.max_score + tr.max_score - chim_overlap as i32;
                let overlap1 = if i_wt > 0 && chim_score_best > 0 {
                    blocksoverlap_l3_blocksoverlap(&read_align.tr_chim[1], &tr)
                } else {
                    0
                };

                if chim_score > chim_score_best {
                    read_align.tr_chim[1] = tr.clone();
                    tr_chim1 = tr.clone();
                    tr_chim1_set = true;
                    if overlap1 == 0 {
                        chim_score_next = chim_score_best;
                    }
                    chim_score_best = chim_score;
                    read_align.tr_chim[1].ro_start = if read_align.tr_chim[1].ro_str == 0 {
                        read_align.tr_chim[1].r_start
                    } else {
                        read_align.l_read
                            - read_align.tr_chim[1].r_start
                            - read_align.tr_chim[1].r_length
                    };
                    read_align.tr_chim[1].c_start = read_align.tr_chim[1].g_start
                        - map_gen.chr_start[read_align.tr_chim[1].chr as usize];
                    chim_str_best = chim_str1;
                } else if chim_score > chim_score_next && overlap1 == 0 {
                    chim_score_next = chim_score;
                }
            }
        }
    }

    if !(chim_score_best >= p.p_ch.score_min
        && chim_score_best + p.p_ch.score_drop_max >= (read_len0 + read_len1) as i32)
    {
        return false;
    }

    if read_align.n_tr > p.p_ch.main_segment_mult_nmax as u64 {
        if !tr_chim1_set
            || (read_align.tr_mult.first() != Some(&tr_chim1)
                && read_align.tr_mult.get(1) != Some(&tr_chim1))
        {
            return false;
        }
    }

    if read_align.chim_str == 0 {
        read_align.chim_str = chim_str_best;
    }
    read_align.chim_n = 0;
    if chim_score_next + p.p_ch.score_separation >= chim_score_best {
        return false;
    }
    if read_align.tr_chim[0].ro_start > read_align.tr_chim[1].ro_start {
        read_align.tr_chim.swap(0, 1);
    }

    let e0 = if read_align.tr_chim[0].str_ == 1 {
        0
    } else {
        read_align.tr_chim[0].n_exons - 1
    } as usize;
    let e1 = if read_align.tr_chim[1].str_ == 0 {
        0
    } else {
        read_align.tr_chim[1].n_exons - 1
    } as usize;

    read_align.chim_repeat0 = 0;
    read_align.chim_repeat1 = 0;
    read_align.chim_j0 = 0;
    read_align.chim_j1 = 0;
    read_align.chim_motif = 0;
    read_align.chim_n = 2;

    if read_align.tr_chim[0].exons[e0][EX_IFRAG] > read_align.tr_chim[1].exons[e1][EX_IFRAG] {
        return false;
    } else if read_align.tr_chim[0].exons[e0][EX_IFRAG] < read_align.tr_chim[1].exons[e1][EX_IFRAG]
    {
        read_align.chim_n = 2;
        read_align.chim_motif = -1;
        read_align.chim_j0 = if read_align.tr_chim[0].str_ == 1 {
            read_align.tr_chim[0].exons[e0][EX_G] - 1
        } else {
            read_align.tr_chim[0].exons[e0][EX_G] + read_align.tr_chim[0].exons[e0][EX_L]
        };
        read_align.chim_j1 = if read_align.tr_chim[1].str_ == 0 {
            read_align.tr_chim[1].exons[e1][EX_G] - 1
        } else {
            read_align.tr_chim[1].exons[e1][EX_G] + read_align.tr_chim[1].exons[e1][EX_L]
        };
    } else {
        if !(read_align.tr_chim[0].exons[e0][EX_L] >= p.p_ch.junction_overhang_min as u64
            && read_align.tr_chim[1].exons[e1][EX_L] >= p.p_ch.junction_overhang_min as u64)
        {
            return false;
        }
        let ro_start0 = if read_align.tr_chim[0].str_ == 0 {
            read_align.tr_chim[0].exons[e0][EX_R]
        } else {
            read_align.l_read
                - read_align.tr_chim[0].exons[e0][EX_R]
                - read_align.tr_chim[0].exons[e0][EX_L]
        };
        let ro_start1 = if read_align.tr_chim[1].str_ == 0 {
            read_align.tr_chim[1].exons[e1][EX_R]
        } else {
            read_align.l_read
                - read_align.tr_chim[1].exons[e1][EX_R]
                - read_align.tr_chim[1].exons[e1][EX_L]
        };

        let mut j_rbest = 0_u64;
        let mut j_score = 0_i32;
        let mut j_score_best = -999999_i32;
        let j_rmax = if ro_start1 + read_align.tr_chim[1].exons[e1][EX_L] > ro_start0 {
            ro_start1 + read_align.tr_chim[1].exons[e1][EX_L] - ro_start0 - 1
        } else {
            0
        };
        let mut j_r = 0_u64;
        while j_r < j_rmax {
            if j_r == read_len0 {
                j_r += 1;
            }
            let b_r = read_align.read1[0][(ro_start0 + j_r) as usize];
            let mut b0;
            let mut b1;
            if read_align.tr_chim[0].str_ == 0 {
                b0 = map_gen.g[(read_align.tr_chim[0].exons[e0][EX_G] + j_r) as usize];
            } else {
                b0 = map_gen.g[(read_align.tr_chim[0].exons[e0][EX_G]
                    + read_align.tr_chim[0].exons[e0][EX_L]
                    - 1
                    - j_r) as usize];
                if b0 < 4 {
                    b0 = 3 - b0;
                }
            }
            if read_align.tr_chim[1].str_ == 0 {
                b1 = map_gen.g[(read_align.tr_chim[1].exons[e1][EX_G] - ro_start1 + ro_start0 + j_r)
                    as usize];
            } else {
                b1 = map_gen.g[(read_align.tr_chim[1].exons[e1][EX_G]
                    + read_align.tr_chim[1].exons[e1][EX_L]
                    - 1
                    + ro_start1
                    - ro_start0
                    - j_r) as usize];
                if b1 < 4 {
                    b1 = 3 - b1;
                }
            }
            if (p.p_ch.filter_genomic_n && (b0 > 3 || b1 > 3)) || b_r > 3 {
                read_align.chim_n = 0;
                break;
            }

            let mut b01;
            let mut b02;
            let mut b11;
            let mut b12;
            if read_align.tr_chim[0].str_ == 0 {
                b01 = map_gen.g[(read_align.tr_chim[0].exons[e0][EX_G] + j_r + 1) as usize];
                b02 = map_gen.g[(read_align.tr_chim[0].exons[e0][EX_G] + j_r + 2) as usize];
            } else {
                b01 = map_gen.g[(read_align.tr_chim[0].exons[e0][EX_G]
                    + read_align.tr_chim[0].exons[e0][EX_L]
                    - 1
                    - j_r
                    - 1) as usize];
                if b01 < 4 {
                    b01 = 3 - b01;
                }
                b02 = map_gen.g[(read_align.tr_chim[0].exons[e0][EX_G]
                    + read_align.tr_chim[0].exons[e0][EX_L]
                    - 1
                    - j_r
                    - 2) as usize];
                if b02 < 4 {
                    b02 = 3 - b02;
                }
            }
            if read_align.tr_chim[1].str_ == 0 {
                b11 =
                    map_gen.g[(read_align.tr_chim[1].exons[e1][EX_G] - ro_start1 + ro_start0 + j_r
                        - 1) as usize];
                b12 = map_gen.g[(read_align.tr_chim[1].exons[e1][EX_G] - ro_start1
                    + ro_start0
                    + j_r) as usize];
            } else {
                b11 = map_gen.g[(read_align.tr_chim[1].exons[e1][EX_G]
                    + read_align.tr_chim[1].exons[e1][EX_L]
                    - 1
                    + ro_start1
                    - ro_start0
                    - j_r
                    + 1) as usize];
                if b11 < 4 {
                    b11 = 3 - b11;
                }
                b12 = map_gen.g[(read_align.tr_chim[1].exons[e1][EX_G]
                    + read_align.tr_chim[1].exons[e1][EX_L]
                    - 1
                    + ro_start1
                    - ro_start0
                    - j_r) as usize];
                if b12 < 4 {
                    b12 = 3 - b12;
                }
            }
            let mut j_motif = 0_i32;
            if b01 == 2 && b02 == 3 && b11 == 0 && b12 == 2 {
                if read_align.chim_str != 2 {
                    j_motif = 1;
                }
            } else if b01 == 1 && b02 == 3 && b11 == 0 && b12 == 1 {
                if read_align.chim_str != 1 {
                    j_motif = 2;
                }
            }
            if b_r == b0 && b_r != b1 {
                j_score += 1;
            } else if b_r != b0 && b_r == b1 {
                j_score -= 1;
            }
            let j_score_j = if j_motif == 0 {
                j_score + p.p_ch.score_junction_non_gtag
            } else {
                j_score
            };
            if j_score_j > j_score_best || (j_score_j == j_score_best && j_motif > 0) {
                read_align.chim_motif = j_motif;
                j_rbest = j_r;
                j_score_best = j_score_j;
            }
            j_r += 1;
        }
        if read_align.chim_n == 0 {
            return false;
        }
        if read_align.chim_motif == 0 {
            chim_score_best += 1 + p.p_ch.score_junction_non_gtag;
            if !(chim_score_best >= p.p_ch.score_min
                && chim_score_best + p.p_ch.score_drop_max >= (read_len0 + read_len1) as i32)
            {
                return false;
            }
        }

        if read_align.tr_chim[0].str_ == 1 {
            read_align.tr_chim[0].exons[e0][EX_R] +=
                read_align.tr_chim[0].exons[e0][EX_L] - j_rbest - 1;
            read_align.tr_chim[0].exons[e0][EX_G] +=
                read_align.tr_chim[0].exons[e0][EX_L] - j_rbest - 1;
            read_align.tr_chim[0].exons[e0][EX_L] = j_rbest + 1;
            read_align.chim_j0 = read_align.tr_chim[0].exons[e0][EX_G] - 1;
        } else {
            read_align.tr_chim[0].exons[e0][EX_L] = j_rbest + 1;
            read_align.chim_j0 =
                read_align.tr_chim[0].exons[e0][EX_G] + read_align.tr_chim[0].exons[e0][EX_L];
        }

        if read_align.tr_chim[1].str_ == 0 {
            read_align.tr_chim[1].exons[e1][EX_R] += ro_start0 + j_rbest + 1 - ro_start1;
            read_align.tr_chim[1].exons[e1][EX_G] += ro_start0 + j_rbest + 1 - ro_start1;
            read_align.tr_chim[1].exons[e1][EX_L] =
                ro_start1 + read_align.tr_chim[1].exons[e1][EX_L] - ro_start0 - j_rbest - 1;
            read_align.chim_j1 = read_align.tr_chim[1].exons[e1][EX_G] - 1;
        } else {
            read_align.tr_chim[1].exons[e1][EX_L] =
                ro_start1 + read_align.tr_chim[1].exons[e1][EX_L] - ro_start0 - j_rbest - 1;
            read_align.chim_j1 =
                read_align.tr_chim[1].exons[e1][EX_G] + read_align.tr_chim[1].exons[e1][EX_L];
        }

        let mut repeat = 0_u64;
        while repeat < 100 {
            let mut b0 = if read_align.tr_chim[0].str_ == 0 {
                map_gen.g[(read_align.chim_j0 + repeat) as usize]
            } else {
                map_gen.g[(read_align.chim_j0 - repeat) as usize]
            };
            if read_align.tr_chim[0].str_ != 0 && b0 < 4 {
                b0 = 3 - b0;
            }
            let mut b1 = if read_align.tr_chim[1].str_ == 0 {
                map_gen.g[(read_align.chim_j1 + 1 + repeat) as usize]
            } else {
                map_gen.g[(read_align.chim_j1 - 1 - repeat) as usize]
            };
            if read_align.tr_chim[1].str_ != 0 && b1 < 4 {
                b1 = 3 - b1;
            }
            if b0 != b1 {
                break;
            }
            repeat += 1;
        }
        read_align.chim_repeat1 = repeat;

        repeat = 0;
        while repeat < 100 {
            let mut b0 = if read_align.tr_chim[0].str_ == 0 {
                map_gen.g[(read_align.chim_j0 - 1 - repeat) as usize]
            } else {
                map_gen.g[(read_align.chim_j0 + 1 + repeat) as usize]
            };
            if read_align.tr_chim[0].str_ != 0 && b0 < 4 {
                b0 = 3 - b0;
            }
            let mut b1 = if read_align.tr_chim[1].str_ == 0 {
                map_gen.g[(read_align.chim_j1 - repeat) as usize]
            } else {
                map_gen.g[(read_align.chim_j1 + repeat) as usize]
            };
            if read_align.tr_chim[1].str_ != 0 && b1 < 4 {
                b1 = 3 - b1;
            }
            if b0 != b1 {
                break;
            }
            repeat += 1;
        }
        read_align.chim_repeat0 = repeat;
    }

    let distance = if read_align.tr_chim[0].str_ == 0 {
        read_align.chim_j1.wrapping_sub(read_align.chim_j0) as u64 + 1
    } else {
        read_align.chim_j0.wrapping_sub(read_align.chim_j1) as u64 + 1
    };
    if read_align.tr_chim[0].str_ != read_align.tr_chim[1].str_
        || read_align.tr_chim[0].chr != read_align.tr_chim[1].chr
        || distance
            > if read_align.chim_motif >= 0 {
                p.align_intron_max as u64
            } else {
                p.align_mates_gap_max as u64
            }
    {
        if read_align.chim_motif >= 0
            && (read_align.tr_chim[0].exons[e0][EX_L]
                < p.p_ch.junction_overhang_min as u64 + read_align.chim_repeat0
                || read_align.tr_chim[1].exons[e1][EX_L]
                    < p.p_ch.junction_overhang_min as u64 + read_align.chim_repeat1)
        {
            return false;
        }
        return true;
    }

    false
}
