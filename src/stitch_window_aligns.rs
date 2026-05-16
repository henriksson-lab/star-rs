#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `stitchWindowAligns` at STAR/source/stitchWindowAligns.cpp:8. Args: iA: uint, nA: uint, Score: int, WAincl: bool, tR2: uint, tG2: uint, trA: Transcript, Lread: uint, WA: uiWA, R: char, mapGen: Genome, P: Parameters, wTr: Transcript, nWinTr: uint, RA: ReadAlign"]
pub fn stitchwindowaligns_l8_stitchwindowaligns(
    i_a: u32,
    n_a: u32,
    score: i32,
    wa_incl: &mut [bool],
    t_r2: u32,
    t_g2: u32,
    tr_a: crate::transcript::Transcript,
    l_read: u32,
    wa: &[[u32; WA_SIZE]],
    r: &[u8],
    map_gen: &crate::genome::Genome,
    p: &crate::parameters_chimeric::Parameters,
    w_tr: &mut Vec<crate::transcript::Transcript>,
    n_win_tr: &mut u32,
    out_filter_mismatch_nmax_total: u32,
    read_length: &[u32],
    max_score_mate: &mut Vec<i32>,
) -> Result<(), String> {
    if i_a >= n_a && t_r2 == 0 {
        return Ok(());
    }

    if i_a >= n_a {
        let mut tr_a = tr_a;
        let mut score = score;
        let mut t_r2 = t_r2;
        let mut t_g2 = t_g2;
        let mut tr_step = crate::transcript::Transcript::default();
        let v_order = if tr_a.ro_str == 0 { [0, 1] } else { [1, 0] };

        for ord in v_order {
            match ord {
                0 => {
                    if tr_a.r_start > 0
                        && tr_a.g_start > 0
                        && t_r2 >= tr_a.r_start
                        && tr_a.n_exons > 0
                    {
                        transcript_l8_transcript_reset(&mut tr_step);
                        let imate = tr_a.exons[0][EX_IFRAG] as usize;
                        let ext = p
                            .align_ends_type
                            .ext
                            .get(imate)
                            .and_then(|row| row.get((tr_a.str_ != imate as u32) as usize))
                            .copied()
                            .unwrap_or(false);
                        if extendalign_l6_extendalign(
                            r,
                            &map_gen.g,
                            tr_a.r_start - 1,
                            tr_a.g_start - 1,
                            -1,
                            -1,
                            tr_a.r_start,
                            t_r2 - tr_a.r_start + 1,
                            tr_a.n_mm,
                            out_filter_mismatch_nmax_total,
                            p.out_filter_mismatch_nover_lmax,
                            ext,
                            &mut tr_step,
                        ) {
                            transcript_l28_transcript_add(&mut tr_a, &tr_step);
                            score += tr_step.max_score;
                            tr_a.exons[0][EX_R] = tr_a.r_start - tr_step.extend_l;
                            tr_a.r_start = tr_a.exons[0][EX_R];
                            tr_a.exons[0][EX_G] = tr_a.g_start - tr_step.extend_l;
                            tr_a.g_start = tr_a.exons[0][EX_G];
                            tr_a.exons[0][EX_L] += tr_step.extend_l;
                        }
                    }
                }
                1 => {
                    if t_r2 < l_read && t_r2 >= tr_a.r_start && tr_a.n_exons > 0 {
                        transcript_l8_transcript_reset(&mut tr_step);
                        let last = tr_a.n_exons as usize - 1;
                        let imate = tr_a.exons[last][EX_IFRAG] as usize;
                        let ext = p
                            .align_ends_type
                            .ext
                            .get(imate)
                            .and_then(|row| row.get((imate as u32 == tr_a.str_) as usize))
                            .copied()
                            .unwrap_or(false);
                        if extendalign_l6_extendalign(
                            r,
                            &map_gen.g,
                            t_r2 + 1,
                            t_g2 + 1,
                            1,
                            1,
                            l_read - t_r2 - 1,
                            t_r2 - tr_a.r_start + 1,
                            tr_a.n_mm,
                            out_filter_mismatch_nmax_total,
                            p.out_filter_mismatch_nover_lmax,
                            ext,
                            &mut tr_step,
                        ) {
                            transcript_l28_transcript_add(&mut tr_a, &tr_step);
                            score += tr_step.max_score;
                            t_r2 += tr_step.extend_l;
                            t_g2 += tr_step.extend_l;
                            tr_a.exons[last][EX_L] += tr_step.extend_l;
                        }
                    }
                }
                _ => {}
            }
        }

        if tr_a.n_exons == 0 {
            return Ok(());
        }
        if !map_gen.chr_bin.is_empty() {
            let bin = (tr_a.exons[0][EX_G] >> map_gen.p_ge.g_chr_bin_nbits) as usize;
            if let Some(chr) = map_gen.chr_bin.get(bin) {
                tr_a.chr = *chr;
            }
        }
        let chr = tr_a.chr as usize;
        if !p.align_soft_clip_at_reference_ends_yes && chr < map_gen.chr_start.len() {
            let chr_end =
                map_gen.chr_start[chr] + map_gen.chr_length.get(chr).copied().unwrap_or(0);
            let last = tr_a.n_exons as usize - 1;
            if (tr_a.exons[last][EX_G] as u64 + l_read as u64 - tr_a.exons[last][EX_R] as u64)
                > chr_end
                || (tr_a.exons[0][EX_G] as u64)
                    < map_gen.chr_start[chr] + tr_a.exons[0][EX_R] as u64
            {
                return Ok(());
            }
        }

        tr_a.r_length = 0;
        for exon in tr_a.exons.iter().take(tr_a.n_exons as usize) {
            tr_a.r_length += exon[EX_L];
        }
        tr_a.g_length = t_g2.wrapping_add(1).wrapping_sub(tr_a.g_start);

        for isj in 0..tr_a.n_exons.saturating_sub(1) as usize {
            if tr_a.canon_sj.get(isj).copied().unwrap_or(-1) >= 0 {
                if tr_a.sj_annot.get(isj).copied().unwrap_or(0) == 1 {
                    let left_bad = tr_a.exons[isj][EX_L] < p.align_sjdb_overhang_min
                        && (isj == 0
                            || tr_a.canon_sj.get(isj - 1).copied().unwrap_or(-1) == -3
                            || (tr_a.sj_annot.get(isj - 1).copied().unwrap_or(0) == 0
                                && tr_a.canon_sj.get(isj - 1).copied().unwrap_or(-1) >= 0));
                    let right_bad = tr_a.exons[isj + 1][EX_L] < p.align_sjdb_overhang_min
                        && (isj == tr_a.n_exons as usize - 2
                            || tr_a.canon_sj.get(isj + 1).copied().unwrap_or(-1) == -3
                            || (tr_a.sj_annot.get(isj + 1).copied().unwrap_or(0) == 0
                                && tr_a.canon_sj.get(isj + 1).copied().unwrap_or(-1) >= 0));
                    if left_bad || right_bad {
                        return Ok(());
                    }
                } else if tr_a.exons[isj][EX_L]
                    < p.align_sj_overhang_min + tr_a.shift_sj.get(isj).copied().unwrap_or([0; 2])[0]
                    || tr_a.exons[isj + 1][EX_L]
                        < p.align_sj_overhang_min
                            + tr_a.shift_sj.get(isj).copied().unwrap_or([0; 2])[1]
                {
                    return Ok(());
                }
            }
        }
        if tr_a.n_exons > 1
            && tr_a
                .sj_annot
                .get(tr_a.n_exons as usize - 2)
                .copied()
                .unwrap_or(0)
                == 1
            && tr_a.exons[tr_a.n_exons as usize - 1][EX_L] < p.align_sjdb_overhang_min
        {
            return Ok(());
        }

        let mut sj_n = 0_u32;
        tr_a.intron_motifs = [0; 3];
        tr_a.sj_yes = false;
        for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
            if tr_a.canon_sj.get(iex).copied().unwrap_or(-1) >= 0 {
                sj_n += 1;
                let sj_str = tr_a.sj_str.get(iex).copied().unwrap_or(0) as usize;
                if sj_str < tr_a.intron_motifs.len() {
                    tr_a.intron_motifs[sj_str] += 1;
                }
                tr_a.sj_yes = true;
            }
        }
        tr_a.sj_motif_strand = if tr_a.intron_motifs[1] > 0 && tr_a.intron_motifs[2] == 0 {
            1
        } else if tr_a.intron_motifs[1] == 0 && tr_a.intron_motifs[2] > 0 {
            2
        } else {
            0
        };
        if tr_a.intron_motifs[1] > 0
            && tr_a.intron_motifs[2] > 0
            && p.out_filter_intron_strands == "RemoveInconsistentStrands"
        {
            return Ok(());
        }
        if sj_n > 0 && tr_a.sj_motif_strand == 0 && p.out_sam_strand_field_type == 1 {
            return Ok(());
        }
        if p.out_filter_intron_motifs == "RemoveNoncanonical" {
            for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
                if tr_a.canon_sj.get(iex).copied().unwrap_or(-1) == 0 {
                    return Ok(());
                }
            }
        } else if p.out_filter_intron_motifs == "RemoveNoncanonicalUnannotated" {
            for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
                if tr_a.canon_sj.get(iex).copied().unwrap_or(-1) == 0
                    && tr_a.sj_annot.get(iex).copied().unwrap_or(0) == 0
                {
                    return Ok(());
                }
            }
        } else if p.out_filter_intron_motifs != "None" && !p.out_filter_intron_motifs.is_empty() {
            return Err(format!(
                "EXITING because of FATAL INPUT error: unrecognized value of --outFilterIntronMotifs={}\nSOLUTION: re-run STAR with --outFilterIntronMotifs = None -OR- RemoveNoncanonical -OR- RemoveNoncanonicalUnannotated\n",
                p.out_filter_intron_motifs
            ));
        }

        let mut nsj = 0_u32;
        let mut exl = 0_u32;
        for iex in 0..tr_a.n_exons as usize {
            exl += tr_a.exons[iex][EX_L];
            if iex == tr_a.n_exons as usize - 1
                || tr_a.canon_sj.get(iex).copied().unwrap_or(-1) == -3
            {
                let mate = tr_a.exons[iex][EX_IFRAG] as usize;
                let mate_len = read_length.get(mate).copied().unwrap_or(l_read);
                if nsj > 0
                    && (exl < p.align_spliced_mate_map_lmin
                        || exl
                            < (p.align_spliced_mate_map_lmin_over_lmate * mate_len as f64) as u32)
                {
                    return Ok(());
                }
                exl = 0;
                nsj = 0;
            } else if tr_a.canon_sj.get(iex).copied().unwrap_or(-1) >= 0 {
                nsj += 1;
            }
        }

        if p.out_filter_by_sjout_stage == 2 {
            for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
                if tr_a.canon_sj.get(iex).copied().unwrap_or(-1) >= 0
                    && tr_a.sj_annot.get(iex).copied().unwrap_or(0) == 0
                {
                    let js = tr_a.exons[iex][EX_G] + tr_a.exons[iex][EX_L];
                    let je = tr_a.exons[iex + 1][EX_G] - 1;
                    if binarysearch2_l3_binarysearch2(
                        js,
                        je,
                        &p.sj_novel_start,
                        &p.sj_novel_end,
                        p.sj_novel_n as i32,
                    ) < 0
                    {
                        return Ok(());
                    }
                }
            }
        }

        if tr_a.exons[0][EX_IFRAG] != tr_a.exons[tr_a.n_exons as usize - 1][EX_IFRAG] {
            let last = tr_a.n_exons as usize - 1;
            if tr_a.exons[last][EX_G] + tr_a.exons[last][EX_L] <= tr_a.exons[0][EX_G] {
                return Ok(());
            }
            let mut iex_m2 = tr_a.n_exons as usize;
            for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
                if tr_a.canon_sj.get(iex).copied().unwrap_or(-1) == -3 {
                    iex_m2 = iex + 1;
                    break;
                }
            }
            if iex_m2 > 0
                && iex_m2 < tr_a.n_exons as usize
                && tr_a.exons[iex_m2 - 1][EX_G] + tr_a.exons[iex_m2 - 1][EX_L]
                    > tr_a.exons[iex_m2][EX_G]
            {
                if tr_a.exons[0][EX_G]
                    > tr_a.exons[iex_m2][EX_G]
                        + tr_a.exons[0][EX_R]
                        + p.align_ends_protrude.n_bases_max as u32
                    || tr_a.exons[iex_m2 - 1][EX_G] + tr_a.exons[iex_m2 - 1][EX_L]
                        > tr_a.exons[last][EX_G] + l_read - tr_a.exons[last][EX_R]
                            + p.align_ends_protrude.n_bases_max as u32
                {
                    return Ok(());
                }
                let mut iex1 = 1usize;
                let mut iex2 = iex_m2 + 1;
                while iex1 < iex_m2 {
                    if tr_a.exons[iex1][EX_G]
                        >= tr_a.exons[iex2 - 1][EX_G] + tr_a.exons[iex2 - 1][EX_L]
                    {
                        break;
                    }
                    iex1 += 1;
                }
                while iex1 < iex_m2 && iex2 < tr_a.n_exons as usize {
                    if tr_a.canon_sj.get(iex1 - 1).copied().unwrap_or(-1) < 0 {
                        iex1 += 1;
                        continue;
                    }
                    if tr_a.canon_sj.get(iex2 - 1).copied().unwrap_or(-1) < 0 {
                        iex2 += 1;
                        continue;
                    }
                    if tr_a.exons[iex1][EX_G] != tr_a.exons[iex2][EX_G]
                        || tr_a.exons[iex1 - 1][EX_G] + tr_a.exons[iex1 - 1][EX_L]
                            != tr_a.exons[iex2 - 1][EX_G] + tr_a.exons[iex2 - 1][EX_L]
                    {
                        return Ok(());
                    }
                    iex1 += 1;
                    iex2 += 1;
                }
            }
        }

        if p.score_genomic_length_log2scale != 0.0 {
            let last = tr_a.n_exons as usize - 1;
            let glen = tr_a.exons[last][EX_G] + tr_a.exons[last][EX_L] - tr_a.exons[0][EX_G];
            score += ((glen as f64).log2() * p.score_genomic_length_log2scale - 0.5).ceil() as i32;
            score = score.max(0);
        }
        tr_a.ro_start = if tr_a.ro_str == 0 {
            tr_a.r_start
        } else {
            l_read - tr_a.r_start - tr_a.r_length
        };
        tr_a.max_score = score;
        if tr_a.exons[0][EX_IFRAG] == tr_a.exons[tr_a.n_exons as usize - 1][EX_IFRAG] {
            tr_a.i_frag = tr_a.exons[0][EX_IFRAG] as i32;
            let mate = tr_a.i_frag as usize;
            if max_score_mate.len() <= mate {
                max_score_mate.resize(mate + 1, i32::MIN);
            }
            max_score_mate[mate] = max_score_mate[mate].max(score);
        } else {
            tr_a.i_frag = -1;
        }
        score += transcript_variationadjust_l4_transcript_variationadjust(&mut tr_a, map_gen, r);
        tr_a.max_score = score;

        let best_score = w_tr.first().map(|tr| tr.max_score).unwrap_or(i32::MIN);
        let mate_ok = tr_a.i_frag >= 0
            && max_score_mate
                .get(tr_a.i_frag as usize)
                .map(|mate_score| {
                    score.saturating_add(p.out_filter_multimap_score_range) >= *mate_score
                })
                .unwrap_or(false);
        if score.saturating_add(p.out_filter_multimap_score_range) >= best_score
            || mate_ok
            || p.p_ch.segment_min > 0
        {
            tr_a.mapped_length = tr_a
                .exons
                .iter()
                .take(tr_a.n_exons as usize)
                .map(|ex| ex[EX_L])
                .sum();
            let mut i_tr = 0usize;
            while i_tr < *n_win_tr as usize {
                let n_overlap = blocksoverlap_l3_blocksoverlap(&tr_a, &w_tr[i_tr]);
                let u_new = tr_a.mapped_length - n_overlap;
                let u_old = w_tr[i_tr].mapped_length - n_overlap;
                if u_new == 0 && score < w_tr[i_tr].max_score {
                    break;
                } else if u_old == 0 {
                    let n_win = *n_win_tr as usize;
                    w_tr[i_tr..n_win].rotate_left(1);
                    *n_win_tr -= 1;
                } else if u_old > 0 && (u_new > 0 || score >= w_tr[i_tr].max_score) {
                    i_tr += 1;
                }
            }
            if i_tr == *n_win_tr as usize {
                let mut insert_at = 0usize;
                while insert_at < *n_win_tr as usize {
                    if score > w_tr[insert_at].max_score
                        || (score == w_tr[insert_at].max_score
                            && tr_a.g_length < w_tr[insert_at].g_length)
                    {
                        break;
                    }
                    insert_at += 1;
                }
                if w_tr.len() <= *n_win_tr as usize {
                    w_tr.push(crate::transcript::Transcript::default());
                }
                let n_win = *n_win_tr as usize;
                w_tr[insert_at..=n_win].rotate_right(1);
                w_tr[insert_at] = tr_a;
                let max_win = if p.align_transcripts_per_window_nmax > 0 {
                    p.align_transcripts_per_window_nmax
                } else {
                    p.align_transcripts_per_read_nmax
                };
                if max_win == 0 || *n_win_tr < max_win {
                    *n_win_tr += 1;
                }
            }
        }
        return Ok(());
    }

    let skip_allowed = wa[i_a as usize][WA_ANCHOR] != 2 || tr_a.n_anchor > 0;
    let mut d_score = 0_i32;
    let (mut tr_ai, tr_a_skip) = if skip_allowed {
        (tr_a.clone(), Some(tr_a))
    } else {
        (tr_a, None)
    };
    if tr_ai.n_exons > 0 {
        d_score = stitchaligntotranscript_l9_stitchaligntotranscript(
            t_r2,
            t_g2,
            wa[i_a as usize][WA_R_START],
            wa[i_a as usize][WA_G_START],
            wa[i_a as usize][WA_LENGTH],
            wa[i_a as usize][WA_I_FRAG],
            wa[i_a as usize][WA_SJ_A],
            p,
            r,
            map_gen,
            &mut tr_ai,
            out_filter_mismatch_nmax_total,
        );
    } else {
        tr_ai.exons.resize(1, [0; EX_SIZE]);
        tr_ai.canon_sj.resize(1, -1);
        tr_ai.sj_annot.resize(1, 0);
        tr_ai.shift_sj.resize(1, [0; 2]);
        tr_ai.sj_str.resize(1, 0);
        tr_ai.exons[0][EX_R] = wa[i_a as usize][WA_R_START];
        tr_ai.r_start = tr_ai.exons[0][EX_R];
        tr_ai.exons[0][EX_G] = wa[i_a as usize][WA_G_START];
        tr_ai.g_start = tr_ai.exons[0][EX_G];
        tr_ai.exons[0][EX_L] = wa[i_a as usize][WA_LENGTH];
        tr_ai.exons[0][EX_IFRAG] = wa[i_a as usize][WA_I_FRAG];
        tr_ai.exons[0][EX_SJA] = wa[i_a as usize][WA_SJ_A];
        tr_ai.n_exons = 1;
        d_score += wa[i_a as usize][WA_LENGTH] as i32;
        tr_ai.n_match = wa[i_a as usize][WA_LENGTH];
        for incl in wa_incl.iter_mut().take(n_a as usize) {
            *incl = false;
        }
    }

    if d_score > -1_000_000 {
        wa_incl[i_a as usize] = true;
        if wa[i_a as usize][WA_N_REP] == 1 {
            tr_ai.n_unique += 1;
        }
        if wa[i_a as usize][WA_ANCHOR] > 0 {
            tr_ai.n_anchor += 1;
        }
        stitchwindowaligns_l8_stitchwindowaligns(
            i_a + 1,
            n_a,
            score + d_score,
            wa_incl,
            wa[i_a as usize][WA_R_START]
                .wrapping_add(wa[i_a as usize][WA_LENGTH])
                .wrapping_sub(1),
            wa[i_a as usize][WA_G_START]
                .wrapping_add(wa[i_a as usize][WA_LENGTH])
                .wrapping_sub(1),
            tr_ai,
            l_read,
            wa,
            r,
            map_gen,
            p,
            w_tr,
            n_win_tr,
            out_filter_mismatch_nmax_total,
            read_length,
            max_score_mate,
        )?;
    }

    if let Some(tr_a) = tr_a_skip {
        wa_incl[i_a as usize] = false;
        stitchwindowaligns_l8_stitchwindowaligns(
            i_a + 1,
            n_a,
            score,
            wa_incl,
            t_r2,
            t_g2,
            tr_a,
            l_read,
            wa,
            r,
            map_gen,
            p,
            w_tr,
            n_win_tr,
            out_filter_mismatch_nmax_total,
            read_length,
            max_score_mate,
        )?;
    }
    Ok(())
}
