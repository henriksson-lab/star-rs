#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `stitchAlignToTranscript` at STAR/source/stitchAlignToTranscript.cpp:9. Args: rAend: uint, gAend: uint, rBstart: uint, gBstart: uint, L: uint, iFragB: uint, sjAB: uint, P: Parameters, R: char, mapGen: Genome, trA: Transcript, outFilterMismatchNmaxTotal: uint"]
pub fn stitchaligntotranscript_l9_stitchaligntotranscript(
    r_aend: u64,
    g_aend: u64,
    mut r_bstart: u64,
    mut g_bstart: u64,
    mut l: u64,
    i_frag_b: u64,
    sj_ab: u64,
    p: &crate::parameters_chimeric::Parameters,
    r: &[u8],
    map_gen: &crate::genome::Genome,
    tr_a: &mut crate::transcript::Transcript,
    out_filter_mismatch_nmax_total: u64,
) -> i32 {
    if tr_a.n_exons >= MAX_N_EXONS {
        return -1_000_010;
    }

    if tr_a.exons.len() <= tr_a.n_exons as usize {
        tr_a.exons.resize(tr_a.n_exons as usize + 1, [0; EX_SIZE]);
    }
    if tr_a.canon_sj.len() <= tr_a.n_exons as usize {
        tr_a.canon_sj.resize(tr_a.n_exons as usize + 1, 0);
    }
    if tr_a.sj_annot.len() <= tr_a.n_exons as usize {
        tr_a.sj_annot.resize(tr_a.n_exons as usize + 1, 0);
    }
    if tr_a.shift_sj.len() <= tr_a.n_exons as usize {
        tr_a.shift_sj.resize(tr_a.n_exons as usize + 1, [0; 2]);
    }
    if tr_a.sj_str.len() <= tr_a.n_exons as usize {
        tr_a.sj_str.resize(tr_a.n_exons as usize + 1, 0);
    }

    let g = &map_gen.g;
    let mut score = 0_i32;
    let last = tr_a.n_exons as usize - 1;

    if sj_ab != u64::MAX
        && tr_a.exons[last][EX_SJA] == sj_ab
        && tr_a.exons[last][EX_IFRAG] == i_frag_b
        && r_bstart == r_aend + 1
        && g_aend + 1 < g_bstart
    {
        let sj = sj_ab as usize;
        if map_gen.sjdb_motif[sj] == 0
            && (l <= map_gen.sjdb_shift_right[sj] as u64
                || tr_a.exons[last][EX_L] <= map_gen.sjdb_shift_left[sj] as u64)
        {
            return -1_000_006;
        }
        let out = tr_a.n_exons as usize;
        tr_a.exons[out][EX_L] = l;
        tr_a.exons[out][EX_R] = r_bstart;
        tr_a.exons[out][EX_G] = g_bstart;
        tr_a.canon_sj[last] = map_gen.sjdb_motif[sj] as i32;
        tr_a.shift_sj[last][0] = map_gen.sjdb_shift_left[sj] as u64;
        tr_a.shift_sj[last][1] = map_gen.sjdb_shift_right[sj] as u64;
        tr_a.sj_annot[last] = 1;
        tr_a.sj_str[last] = map_gen.sjdb_strand[sj];
        tr_a.n_exons += 1;
        tr_a.n_match += l;
        for _ii in r_bstart..r_bstart + l {
            score += 1;
        }
        score += p.p_ge.sjdb_score;
    } else {
        tr_a.sj_annot[last] = 0;
        tr_a.sj_str[last] = 0;

        if tr_a.exons[last][EX_IFRAG] == i_frag_b {
            let g_bend = g_bstart.wrapping_add(l).wrapping_sub(1);
            let r_bend = r_bstart.wrapping_add(l).wrapping_sub(1);

            if r_bend <= r_aend {
                return -1_000_001;
            }
            if g_bend <= g_aend && tr_a.exons[last][EX_IFRAG] == i_frag_b {
                return -1_000_002;
            }

            if r_bstart <= r_aend {
                g_bstart += r_aend - r_bstart + 1;
                r_bstart = r_aend + 1;
                l = r_bend.wrapping_sub(r_bstart).wrapping_add(1);
            }

            for _ii in r_bstart..=r_bend {
                score += 1;
            }

            let g_gap = g_bstart as i32 - g_aend as i32 - 1;
            let r_gap = r_bstart as i32 - r_aend as i32 - 1;
            let mut n_match = l;
            let mut n_mm = 0_u64;
            let mut del = 0_u64;
            let mut ins = 0_u64;
            let mut n_ins = 0_u64;
            let mut n_del = 0_u64;
            let mut j_r = 0_i32;
            let mut j_can = 999_i32;
            let g_bstart1 = g_bstart as i32 - r_gap - 1;

            if g_gap == 0 && r_gap == 0 {
            } else if g_gap > 0 && r_gap > 0 && r_gap == g_gap {
                for ii in 1..=r_gap {
                    let gp = (g_aend as i32 + ii) as usize;
                    let rp = (r_aend as i32 + ii) as usize;
                    if g[gp] < 4 && r[rp] < 4 {
                        if r[rp] == g[gp] {
                            score += 1;
                            n_match += 1;
                        } else {
                            score -= 1;
                            n_mm += 1;
                        }
                    }
                }
            } else if g_gap > r_gap {
                n_del = 1;
                del = (g_gap - r_gap) as u64;
                if del > p.align_intron_max as u64 && p.align_intron_max > 0 {
                    return -1_000_003;
                }

                let mut score1 = 0_i32;
                let mut j_r1 = 1_i32;
                loop {
                    j_r1 -= 1;
                    let rp_i = r_aend as i32 + j_r1;
                    let gb_i = g_bstart1 + j_r1;
                    let ga_i = g_aend as i32 + j_r1;
                    if rp_i < 0
                        || gb_i < 0
                        || ga_i < 0
                        || rp_i as usize >= r.len()
                        || gb_i as usize >= g.len()
                        || ga_i as usize >= g.len()
                    {
                        break;
                    }
                    let rp = rp_i as usize;
                    let gb = gb_i as usize;
                    let ga = ga_i as usize;
                    if r[rp] != g[gb] && g[gb] < 4 && r[rp] == g[ga] {
                        score1 -= 1;
                    }
                    if !(score1 + p.score_stitch_sj_shift >= 0
                        && tr_a.exons[last][EX_L] as i32 + j_r1 > 1)
                    {
                        break;
                    }
                }

                let mut max_score2 = -999_999_i32;
                score1 = 0;
                let mut j_pen = 0_i32;
                while j_r1 < r_bend as i32 - r_aend as i32 {
                    let rp_i = r_aend as i32 + j_r1;
                    let ga_i = g_aend as i32 + j_r1;
                    let gb_i = g_bstart1 + j_r1;
                    if rp_i < 0
                        || ga_i < 0
                        || gb_i < 0
                        || rp_i as usize >= r.len()
                        || ga_i as usize >= g.len()
                        || gb_i as usize >= g.len()
                    {
                        break;
                    }
                    let rp = rp_i as usize;
                    let ga = ga_i as usize;
                    let gb = gb_i as usize;
                    if r[rp] == g[ga] && r[rp] != g[gb] {
                        score1 += 1;
                    }
                    if r[rp] != g[ga] && r[rp] == g[gb] {
                        score1 -= 1;
                    }

                    let mut j_can1 = -1_i32;
                    let mut j_pen1 = 0_i32;
                    let mut score2 = score1;
                    if del >= p.align_intron_min as u64 {
                        let d1_i = g_aend as i32 + j_r1 + 1;
                        let d2_i = g_aend as i32 + j_r1 + 2;
                        let a1_i = g_bstart1 + j_r1 - 1;
                        let a2_i = g_bstart1 + j_r1;
                        if d1_i >= 0
                            && d2_i >= 0
                            && a1_i >= 0
                            && a2_i >= 0
                            && (d1_i as usize) < g.len()
                            && (d2_i as usize) < g.len()
                            && (a1_i as usize) < g.len()
                            && (a2_i as usize) < g.len()
                        {
                            let d1 = g[d1_i as usize];
                            let d2 = g[d2_i as usize];
                            let a1 = g[a1_i as usize];
                            let a2 = g[a2_i as usize];
                            if d1 == 2 && d2 == 3 && a1 == 0 && a2 == 2 {
                                j_can1 = 1;
                            } else if d1 == 1 && d2 == 3 && a1 == 0 && a2 == 1 {
                                j_can1 = 2;
                            } else if d1 == 2 && d2 == 1 && a1 == 0 && a2 == 2 {
                                j_can1 = 3;
                                j_pen1 = p.score_gap_gcag;
                            } else if d1 == 1 && d2 == 3 && a1 == 2 && a2 == 1 {
                                j_can1 = 4;
                                j_pen1 = p.score_gap_gcag;
                            } else if d1 == 0 && d2 == 3 && a1 == 0 && a2 == 1 {
                                j_can1 = 5;
                                j_pen1 = p.score_gap_atac;
                            } else if d1 == 2 && d2 == 3 && a1 == 0 && a2 == 3 {
                                j_can1 = 6;
                                j_pen1 = p.score_gap_atac;
                            } else {
                                j_can1 = 0;
                                j_pen1 = p.score_gap_noncan;
                            }
                        } else {
                            j_can1 = 0;
                            j_pen1 = p.score_gap_noncan;
                        }
                        score2 += j_pen1;
                    }

                    if max_score2 < score2 {
                        max_score2 = score2;
                        j_r = j_r1;
                        j_can = j_can1;
                        j_pen = j_pen1;
                    }
                    j_r1 += 1;
                }

                let mut jj_l = 0_u64;
                let mut jj_r = 0_u64;
                while g_aend as i32 + j_r >= jj_l as i32
                    && g[(g_aend as i32 - jj_l as i32 + j_r) as usize]
                        == g[(g_bstart1 - jj_l as i32 + j_r) as usize]
                    && g[(g_aend as i32 - jj_l as i32 + j_r) as usize] < 4
                    && jj_l <= MAX_SJ_REPEAT_SEARCH as u64
                {
                    jj_l += 1;
                }
                while {
                    let ga = g_aend as i64 + jj_r as i64 + j_r as i64 + 1;
                    let gb = g_bstart1 as i64 + jj_r as i64 + j_r as i64 + 1;
                    ga >= 0
                        && gb >= 0
                        && (ga as u64) < map_gen.n_genome
                        && (gb as usize) < g.len()
                        && g[ga as usize] == g[gb as usize]
                        && g[ga as usize] < 4
                        && jj_r <= MAX_SJ_REPEAT_SEARCH as u64
                } {
                    jj_r += 1;
                }

                if j_can <= 0 {
                    j_r -= jj_l as i32;
                    if tr_a.exons[last][EX_L] as i32 + j_r < 1 {
                        return -1_000_005;
                    }
                    jj_r += jj_l;
                    jj_l = 0;
                }

                for ii in std::cmp::min(1, j_r + 1)..=std::cmp::max(r_gap, j_r) {
                    let g1 = if ii <= j_r {
                        g_aend as i32 + ii
                    } else {
                        g_bstart1 + ii
                    } as usize;
                    let rp = (r_aend as i32 + ii) as usize;
                    if g[g1] < 4 && r[rp] < 4 {
                        if r[rp] == g[g1] {
                            if ii >= 1 && ii <= r_gap {
                                score += 1;
                                n_match = n_match.wrapping_add(1);
                            }
                        } else {
                            score -= 1;
                            n_mm = n_mm.wrapping_add(1);
                            if ii < 1 || ii > r_gap {
                                score -= 1;
                                n_match = n_match.wrapping_sub(1);
                            }
                        }
                    }
                }

                if map_gen.sjdb_n > 0 {
                    let j_s = (g_aend as i64 + j_r as i64 + 1) as u64;
                    let j_e = (g_bstart1 + j_r) as i64 as u64;
                    let sjdb_ind = binarysearch2_l3_binarysearch2(
                        j_s,
                        j_e,
                        &map_gen.sjdb_start,
                        &map_gen.sjdb_end,
                        map_gen.sjdb_n as i32,
                    );
                    if sjdb_ind < 0 {
                        if del >= p.align_intron_min as u64 {
                            score += p.score_gap + j_pen;
                        } else {
                            score += del as i32 * p.score_del_base + p.score_del_open;
                            j_can = -1;
                            tr_a.sj_annot[last] = 0;
                        }
                    } else {
                        let sj = sjdb_ind as usize;
                        j_can = map_gen.sjdb_motif[sj] as i32;
                        if map_gen.sjdb_motif[sj] == 0 {
                            if l <= map_gen.sjdb_shift_left[sj] as u64
                                || tr_a.exons[last][EX_L] <= map_gen.sjdb_shift_left[sj] as u64
                            {
                                return -1_000_006;
                            }
                            j_r += map_gen.sjdb_shift_left[sj] as i32;
                            if r_aend as i32 + j_r >= r_bend as i32 {
                                return -1_000_006;
                            }
                            jj_l = map_gen.sjdb_shift_left[sj] as u64;
                            jj_r = map_gen.sjdb_shift_right[sj] as u64;
                        }
                        tr_a.sj_annot[last] = 1;
                        tr_a.sj_str[last] = map_gen.sjdb_strand[sj];
                        score += p.p_ge.sjdb_score;
                    }
                } else if del >= p.align_intron_min as u64 {
                    score += p.score_gap + j_pen;
                } else {
                    score += del as i32 * p.score_del_base + p.score_del_open;
                    j_can = -1;
                    tr_a.sj_annot[last] = 0;
                }

                tr_a.shift_sj[last][0] = jj_l;
                tr_a.shift_sj[last][1] = jj_r;
                tr_a.canon_sj[last] = j_can;
                if tr_a.sj_annot[last] == 0 {
                    tr_a.sj_str[last] = if j_can > 0 { (2 - j_can % 2) as u8 } else { 0 };
                }
            } else if r_gap > g_gap {
                ins = (r_gap - g_gap) as u64;
                n_ins = 1;
                if g_gap == 0 {
                    j_r = 0;
                } else if g_gap < 0 {
                    j_r = 0;
                    for _ii in 0..-g_gap {
                        score -= 1;
                    }
                } else {
                    let mut score1 = 0_i32;
                    let mut max_score1 = 0_i32;
                    for j_r1 in 1..=g_gap {
                        let gp = (g_aend as i32 + j_r1) as usize;
                        if g[gp] < 4 {
                            let r1 = (r_aend as i32 + j_r1) as usize;
                            let r2 = (r_aend as i32 + ins as i32 + j_r1) as usize;
                            score1 += if r[r1] == g[gp] { 1 } else { -1 };
                            score1 += if r[r2] == g[gp] { -1 } else { 1 };
                        }
                        if score1 > max_score1
                            || (score1 == max_score1 && p.align_insertion_flush.flush_right)
                        {
                            max_score1 = score1;
                            j_r = j_r1;
                        }
                    }
                    for ii in 1..=g_gap {
                        let r1 =
                            (r_aend as i32 + ii + if ii <= j_r { 0 } else { ins as i32 }) as usize;
                        let gp = (g_aend as i32 + ii) as usize;
                        if g[gp] < 4 && r[r1] < 4 {
                            if r[r1] == g[gp] {
                                score += 1;
                                n_match += 1;
                            } else {
                                score -= 1;
                                n_mm += 1;
                            }
                        }
                    }
                }

                if p.align_insertion_flush.flush_right {
                    while j_r < r_bend as i32 - r_aend as i32 - ins as i32 {
                        let rp = (r_aend as i32 + j_r + 1) as usize;
                        let gp = (g_aend as i32 + j_r + 1) as usize;
                        if r[rp] != g[gp] || g[gp] == 4 {
                            break;
                        }
                        j_r += 1;
                    }
                    if j_r == r_bend as i32 - r_aend as i32 - ins as i32 {
                        return -1_000_009;
                    }
                }
                score += ins as i32 * p.score_ins_base + p.score_ins_open;
                j_can = -2;
            }

            if tr_a.n_mm + n_mm <= out_filter_mismatch_nmax_total
                && (j_can < 0
                    || (j_can < 7
                        && n_mm
                            <= p.align_sj_stitch_mismatch_nmax
                                .get(((j_can + 1) / 2) as usize)
                                .copied()
                                .unwrap_or(0) as u64))
            {
                tr_a.n_mm = tr_a.n_mm.wrapping_add(n_mm);
                tr_a.n_match = tr_a.n_match.wrapping_add(n_match);
                if del >= p.align_intron_min as u64 {
                    tr_a.n_gap += n_del;
                    tr_a.l_gap += del;
                } else {
                    tr_a.n_del += n_del;
                    tr_a.l_del += del;
                }

                if del == 0 && ins == 0 {
                    tr_a.exons[last][EX_L] += r_bend - r_aend;
                } else if del > 0 {
                    tr_a.exons[last][EX_L] = (tr_a.exons[last][EX_L] as i32 + j_r) as u64;
                    let out = tr_a.n_exons as usize;
                    tr_a.exons[out][EX_L] = (r_bend as i32 - r_aend as i32 - j_r) as u64;
                    tr_a.exons[out][EX_R] = (r_aend as i32 + j_r + 1) as u64;
                    tr_a.exons[out][EX_G] = (g_bstart1 + j_r + 1) as u64;
                    tr_a.n_exons += 1;
                } else if ins > 0 {
                    tr_a.n_ins += n_ins;
                    tr_a.l_ins += ins;
                    tr_a.exons[last][EX_L] = (tr_a.exons[last][EX_L] as i32 + j_r) as u64;
                    let out = tr_a.n_exons as usize;
                    tr_a.exons[out][EX_L] =
                        (r_bend as i32 - r_aend as i32 - j_r - ins as i32) as u64;
                    tr_a.exons[out][EX_R] = (r_aend as i32 + j_r + ins as i32 + 1) as u64;
                    tr_a.exons[out][EX_G] = (g_aend as i32 + 1 + j_r) as u64;
                    tr_a.canon_sj[last] = -2;
                    tr_a.sj_annot[last] = 0;
                    tr_a.n_exons += 1;
                }
            } else {
                return -1_000_007;
            }
        } else if g_bstart + tr_a.exons[0][EX_R] + p.align_ends_protrude.n_bases_max as u64
            >= tr_a.exons[0][EX_G]
            || tr_a.exons[0][EX_G] < tr_a.exons[0][EX_R]
        {
            if p.align_mates_gap_max > 0
                && g_bstart
                    > tr_a.exons[last][EX_G] + tr_a.exons[last][EX_L] + p.align_mates_gap_max as u64
            {
                return -1_000_004;
            }

            for _ii in r_bstart..r_bstart + l {
                score += 1;
            }

            let mut tr_extend = crate::transcript::Transcript::default();
            transcript_l8_transcript_reset(&mut tr_extend);
            if extendalign_l6_extendalign(
                r,
                g,
                r_aend + 1,
                g_aend + 1,
                1,
                1,
                DEF_READ_SEQ_LENGTH_MAX as u64,
                tr_a.n_match,
                tr_a.n_mm,
                out_filter_mismatch_nmax_total,
                p.out_filter_mismatch_nover_lmax,
                p.align_ends_type.ext[tr_a.exons[last][EX_IFRAG] as usize][1],
                &mut tr_extend,
            ) {
                transcript_l28_transcript_add(tr_a, &tr_extend);
                score += tr_extend.max_score;
                tr_a.exons[last][EX_L] += tr_extend.extend_l;
            }

            let out = tr_a.n_exons as usize;
            tr_a.exons[out][EX_R] = r_bstart;
            tr_a.exons[out][EX_G] = g_bstart;
            tr_a.exons[out][EX_L] = l;
            tr_a.n_match += l;

            transcript_l8_transcript_reset(&mut tr_extend);
            let extlen = if p.align_ends_type.ext[i_frag_b as usize][1] {
                DEF_READ_SEQ_LENGTH_MAX as u64
            } else {
                g_bstart
                    .wrapping_sub(tr_a.exons[0][EX_G])
                    .wrapping_add(tr_a.exons[0][EX_R])
            };
            if extendalign_l6_extendalign(
                r,
                g,
                r_bstart.wrapping_sub(1),
                g_bstart.wrapping_sub(1),
                -1,
                -1,
                extlen,
                tr_a.n_match,
                tr_a.n_mm,
                out_filter_mismatch_nmax_total,
                p.out_filter_mismatch_nover_lmax,
                p.align_ends_type.ext[i_frag_b as usize][1],
                &mut tr_extend,
            ) {
                transcript_l28_transcript_add(tr_a, &tr_extend);
                score += tr_extend.max_score;
                tr_a.exons[out][EX_R] = tr_a.exons[out][EX_R].wrapping_sub(tr_extend.extend_l);
                tr_a.exons[out][EX_G] = tr_a.exons[out][EX_G].wrapping_sub(tr_extend.extend_l);
                tr_a.exons[out][EX_L] += tr_extend.extend_l;
            }

            tr_a.canon_sj[last] = -3;
            tr_a.sj_annot[last] = 0;
            tr_a.n_exons += 1;
        } else {
            return -1_000_008;
        }
    }

    let last_new = tr_a.n_exons as usize - 1;
    tr_a.exons[last_new][EX_IFRAG] = i_frag_b;
    tr_a.exons[last_new][EX_SJA] = sj_ab;
    score
}
