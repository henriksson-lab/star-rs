#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::stitchWindowSeeds` at STAR/source/ReadAlign_stitchWindowSeeds.cpp:12. Args: iW: uint, iWrec: uint, WAexcl: bool, R: char"]
pub fn readalign_stitchwindowseeds_l12_readalign_stitchwindowseeds(
    read_align: &mut crate::read_align::ReadAlign,
    i_w: u32,
    i_wrec: u32,
    wa_excl: Option<&[bool]>,
    r: &[u8],
    map_gen: &crate::genome::Genome,
    p: &crate::parameters_chimeric::Parameters,
) -> Result<bool, String> {
    let i_w_usize = i_w as usize;
    let n_wa = read_align.n_wa[i_w_usize] as usize;
    if read_align.score_seed_best.len() < n_wa {
        read_align.score_seed_best.resize(n_wa, 0);
        read_align.score_seed_best_mm.resize(n_wa, 0);
        read_align.score_seed_best_ind.resize(n_wa, 0);
        read_align.seed_chain.resize(n_wa, 0);
    }
    if read_align.wa_incl.len() < n_wa {
        read_align.wa_incl.resize(n_wa, false);
    }
    for v in read_align.wa_incl.iter_mut().take(n_wa) {
        *v = false;
    }

    let mut tr_a1;
    for i_s1 in 0..n_wa {
        read_align.score_seed_best[i_s1] = 0;
        read_align.score_seed_best_mm[i_s1] = 0;
        read_align.score_seed_best_ind[i_s1] = u64::MAX;
        if wa_excl
            .map(|x| x.get(i_s1).copied().unwrap_or(false))
            .unwrap_or(false)
        {
            continue;
        }
        for i_s2 in 0..=i_s1 {
            if wa_excl
                .map(|x| x.get(i_s1).copied().unwrap_or(false))
                .unwrap_or(false)
            {
                continue;
            }
            tr_a1 = (*read_align.tr_init).clone();
            let score2;
            if i_s2 < i_s1 {
                tr_a1.n_exons = 1;
                tr_a1.n_mm = read_align.score_seed_best_mm[i_s2];
                if tr_a1.exons.is_empty() {
                    tr_a1.exons.resize(1, [0; EX_SIZE]);
                }
                if tr_a1.canon_sj.is_empty() {
                    tr_a1.canon_sj.resize(1, 0);
                }
                if tr_a1.sj_annot.is_empty() {
                    tr_a1.sj_annot.resize(1, 0);
                }
                if tr_a1.shift_sj.is_empty() {
                    tr_a1.shift_sj.resize(1, [0; 2]);
                }
                if tr_a1.sj_str.is_empty() {
                    tr_a1.sj_str.resize(1, 0);
                }
                let wa2 = read_align.wa[i_w_usize][i_s2];
                let wa1 = read_align.wa[i_w_usize][i_s1];
                tr_a1.exons[0][EX_R] = wa2[WA_R_START];
                tr_a1.exons[0][EX_G] = wa2[WA_G_START];
                tr_a1.exons[0][EX_L] = wa2[WA_LENGTH];
                tr_a1.exons[0][EX_IFRAG] = wa2[WA_I_FRAG];
                tr_a1.exons[0][EX_SJA] = wa2[WA_SJ_A];
                score2 = stitchaligntotranscript_l9_stitchaligntotranscript(
                    wa2[WA_R_START] + wa2[WA_LENGTH] - 1,
                    wa2[WA_G_START] + wa2[WA_LENGTH] - 1,
                    wa1[WA_R_START],
                    wa1[WA_G_START],
                    wa1[WA_LENGTH],
                    wa1[WA_I_FRAG],
                    wa1[WA_SJ_A],
                    p,
                    r,
                    map_gen,
                    &mut tr_a1,
                    read_align.out_filter_mismatch_nmax_total,
                );

                if p.out_filter_by_sjout_stage == 2 && tr_a1.n_exons > 1 {
                    let iex = 0usize;
                    if tr_a1.canon_sj[iex] >= 0 && tr_a1.sj_annot[iex] == 0 {
                        let j_s = (tr_a1.exons[iex][EX_G] + tr_a1.exons[iex][EX_L]) as u32;
                        let j_e = (tr_a1.exons[iex + 1][EX_G] - 1) as u32;
                        if binarysearch2_l3_binarysearch2(
                            j_s,
                            j_e,
                            &p.sj_novel_start,
                            &p.sj_novel_end,
                            p.sj_novel_n as i32,
                        ) < 0
                        {
                            return Ok(false);
                        }
                    }
                }

                let exon_long_enough = tr_a1.exons[0][EX_L]
                    >= if tr_a1.sj_annot[0] == 0 {
                        p.align_sj_overhang_min as u64
                    } else {
                        p.align_sjdb_overhang_min as u64
                    };
                if exon_long_enough
                    && score2 > 0
                    && score2 + read_align.score_seed_best[i_s2] > read_align.score_seed_best[i_s1]
                {
                    read_align.score_seed_best[i_s1] = score2 + read_align.score_seed_best[i_s2];
                    read_align.score_seed_best_mm[i_s1] = tr_a1.n_mm;
                    read_align.score_seed_best_ind[i_s1] = i_s2 as u64;
                }
            } else {
                let wa1 = read_align.wa[i_w_usize][i_s1];
                let mut score2_local = wa1[WA_LENGTH] as i32;
                if wa1[WA_R_START] > 0
                    && extendalign_l6_extendalign(
                        r,
                        &map_gen.g,
                        wa1[WA_R_START] - 1,
                        wa1[WA_G_START] - 1,
                        -1,
                        -1,
                        wa1[WA_R_START],
                        100_000,
                        0,
                        read_align.out_filter_mismatch_nmax_total,
                        p.out_filter_mismatch_nover_lmax,
                        p.align_ends_type.ext[wa1[WA_I_FRAG] as usize]
                            [read_align.tr_best.str_ as usize],
                        &mut tr_a1,
                    )
                {
                    score2_local += tr_a1.max_score;
                }
                let exon_long_enough =
                    wa1[WA_LENGTH] + tr_a1.extend_l >= p.align_sj_overhang_min as u64;
                if exon_long_enough && score2_local > read_align.score_seed_best[i_s1] {
                    read_align.score_seed_best[i_s1] = score2_local;
                    read_align.score_seed_best_ind[i_s1] = i_s1 as u64;
                }
            }
        }
    }

    let mut score_best = 0_i32;
    let mut score_best_ind = 0_u64;
    for i_s1 in 0..n_wa {
        tr_a1 = (*read_align.tr_init).clone();
        let wa1 = read_align.wa[i_w_usize][i_s1];
        let t_r2 = wa1[WA_R_START] + wa1[WA_LENGTH];
        let t_g2 = wa1[WA_G_START] + wa1[WA_LENGTH];
        if t_r2 < read_align.l_read - 1
            && extendalign_l6_extendalign(
                r,
                &map_gen.g,
                t_r2,
                t_g2,
                1,
                1,
                read_align.l_read - t_r2,
                100_000,
                read_align.score_seed_best_mm[i_s1],
                read_align.out_filter_mismatch_nmax_total,
                p.out_filter_mismatch_nover_lmax,
                p.align_ends_type.ext[wa1[WA_I_FRAG] as usize]
                    [1 - read_align.tr_best.str_ as usize],
                &mut tr_a1,
            )
        {
            read_align.score_seed_best[i_s1] += tr_a1.max_score;
        }

        let exon_long_enough = wa1[WA_LENGTH] + tr_a1.extend_l >= p.align_sj_overhang_min as u64;
        if exon_long_enough && read_align.score_seed_best[i_s1] > score_best {
            score_best = read_align.score_seed_best[i_s1];
            score_best_ind = i_s1 as u64;
        }
    }

    let mut seed_n = 0usize;
    loop {
        read_align.seed_chain[seed_n] = score_best_ind;
        seed_n += 1;
        read_align.wa_incl[score_best_ind as usize] = true;
        let prev = read_align.score_seed_best_ind[score_best_ind as usize];
        if score_best_ind > prev {
            score_best_ind = prev;
        } else {
            break;
        }
    }

    let mut tr_a = (*read_align.tr_init).clone();
    let first_seed = read_align.seed_chain[seed_n - 1] as usize;
    let first_wa = read_align.wa[i_w_usize][first_seed];
    let mut score = first_wa[WA_LENGTH] as i32;
    tr_a.max_score = score;
    tr_a.n_match = first_wa[WA_LENGTH];
    tr_a.n_mm = 0;
    if tr_a.exons.is_empty() {
        tr_a.exons.resize(1, [0; EX_SIZE]);
    }
    if tr_a.canon_sj.is_empty() {
        tr_a.canon_sj.resize(1, 0);
    }
    if tr_a.sj_annot.is_empty() {
        tr_a.sj_annot.resize(1, 0);
    }
    if tr_a.shift_sj.is_empty() {
        tr_a.shift_sj.resize(1, [0; 2]);
    }
    if tr_a.sj_str.is_empty() {
        tr_a.sj_str.resize(1, 0);
    }
    tr_a.exons[0][EX_R] = first_wa[WA_R_START];
    tr_a.r_start = first_wa[WA_R_START];
    tr_a.exons[0][EX_G] = first_wa[WA_G_START];
    tr_a.g_start = first_wa[WA_G_START];
    tr_a.exons[0][EX_L] = first_wa[WA_LENGTH];
    tr_a.exons[0][EX_IFRAG] = first_wa[WA_I_FRAG];
    tr_a.exons[0][EX_SJA] = first_wa[WA_SJ_A];
    tr_a.n_exons = 1;

    for i_sc in (1..seed_n).rev() {
        let i_s1 = read_align.seed_chain[i_sc] as usize;
        let i_s2 = read_align.seed_chain[i_sc - 1] as usize;
        let wa1 = read_align.wa[i_w_usize][i_s1];
        let wa2 = read_align.wa[i_w_usize][i_s2];
        let score_stitch = stitchaligntotranscript_l9_stitchaligntotranscript(
            wa1[WA_R_START] + wa1[WA_LENGTH] - 1,
            wa1[WA_G_START] + wa1[WA_LENGTH] - 1,
            wa2[WA_R_START],
            wa2[WA_G_START],
            wa2[WA_LENGTH],
            wa2[WA_I_FRAG],
            wa2[WA_SJ_A],
            p,
            r,
            map_gen,
            &mut tr_a,
            read_align.out_filter_mismatch_nmax_total,
        );
        score += score_stitch;
    }
    tr_a.max_score = score;

    tr_a1 = (*read_align.tr_init).clone();
    if tr_a.exons[0][EX_R] > 0
        && extendalign_l6_extendalign(
            r,
            &map_gen.g,
            tr_a.exons[0][EX_R] - 1,
            tr_a.exons[0][EX_G] - 1,
            -1,
            -1,
            tr_a.exons[0][EX_R],
            100_000,
            0,
            read_align.out_filter_mismatch_nmax_total,
            p.out_filter_mismatch_nover_lmax,
            p.align_ends_type.ext[tr_a.exons[0][EX_IFRAG] as usize][tr_a.str_ as usize],
            &mut tr_a1,
        )
    {
        transcript_l28_transcript_add(&mut tr_a, &tr_a1);
        tr_a.exons[0][EX_R] -= tr_a1.extend_l;
        tr_a.exons[0][EX_G] -= tr_a1.extend_l;
        tr_a.exons[0][EX_L] += tr_a1.extend_l;
        tr_a.r_start = tr_a.exons[0][EX_R];
        tr_a.g_start = tr_a.exons[0][EX_G];
    }

    let last_chain_seed = read_align.seed_chain[0] as usize;
    tr_a1 = (*read_align.tr_init).clone();
    let wa_last = read_align.wa[i_w_usize][last_chain_seed];
    let t_r2 = wa_last[WA_R_START] + wa_last[WA_LENGTH];
    let t_g2 = wa_last[WA_G_START] + wa_last[WA_LENGTH];
    if t_r2 < read_align.l_read
        && extendalign_l6_extendalign(
            r,
            &map_gen.g,
            t_r2,
            t_g2,
            1,
            1,
            read_align.l_read - t_r2,
            100_000,
            read_align.score_seed_best_mm[last_chain_seed],
            read_align.out_filter_mismatch_nmax_total,
            p.out_filter_mismatch_nover_lmax,
            p.align_ends_type.ext[tr_a.exons[tr_a.n_exons as usize - 1][EX_IFRAG] as usize]
                [1 - tr_a.str_ as usize],
            &mut tr_a1,
        )
    {
        transcript_l28_transcript_add(&mut tr_a, &tr_a1);
        let last = tr_a.n_exons as usize - 1;
        tr_a.exons[last][EX_L] += tr_a1.extend_l;
    }

    tr_a.r_length = 0;
    for isj in 0..tr_a.n_exons as usize {
        tr_a.r_length += tr_a.exons[isj][EX_L];
    }
    let last = tr_a.n_exons as usize - 1;
    tr_a.g_length = tr_a.exons[last][EX_G]
        .wrapping_add(1)
        .wrapping_sub(tr_a.g_start);
    tr_a.ro_start = if tr_a.ro_str == 0 {
        tr_a.r_start
    } else {
        read_align
            .l_read
            .wrapping_sub(tr_a.r_start)
            .wrapping_sub(tr_a.r_length)
    };

    if tr_a.exons[0][EX_IFRAG] == tr_a.exons[last][EX_IFRAG] {
        tr_a.i_frag = tr_a.exons[0][EX_IFRAG] as i32;
        let i_frag = tr_a.i_frag as usize;
        if read_align.max_score_mate.len() <= i_frag {
            read_align.max_score_mate.resize(i_frag + 1, 0);
        }
        read_align.max_score_mate[i_frag] = read_align.max_score_mate[i_frag].max(tr_a.max_score);
    } else {
        tr_a.i_frag = -1;
    }

    if p.score_genomic_length_log2scale != 0.0 {
        let genomic_length = tr_a.exons[last][EX_G] + tr_a.exons[last][EX_L] - tr_a.exons[0][EX_G];
        tr_a.max_score +=
            ((genomic_length as f64).log2() * p.score_genomic_length_log2scale - 0.5).ceil() as i32;
        tr_a.max_score = tr_a.max_score.max(0);
    }

    let mut sj_n = 0_u32;
    tr_a.intron_motifs = [0; 3];
    for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
        if tr_a.canon_sj[iex] >= 0 {
            sj_n += 1;
            tr_a.intron_motifs[tr_a.sj_str[iex] as usize] += 1;
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
        return Ok(false);
    }
    if sj_n > 0 && tr_a.sj_motif_strand == 0 && p.out_sam_strand_field_type == 1 {
        return Ok(false);
    }

    if p.out_filter_intron_motifs.is_empty() || p.out_filter_intron_motifs == "None" {
    } else if p.out_filter_intron_motifs == "RemoveNoncanonical" {
        for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
            if tr_a.canon_sj[iex] == 0 {
                return Ok(false);
            }
        }
    } else if p.out_filter_intron_motifs == "RemoveNoncanonicalUnannotated" {
        for iex in 0..tr_a.n_exons.saturating_sub(1) as usize {
            if tr_a.canon_sj[iex] == 0 && tr_a.sj_annot[iex] == 0 {
                return Ok(false);
            }
        }
    } else {
        return Err(format!(
            "EXITING because of FATAL INPUT error: unrecognized value of --outFilterIntronMotifs={}\nSOLUTION: re-run STAR with --outFilterIntronMotifs = None -OR- RemoveNoncanonical -OR- RemoveNoncanonicalUnannotated\n",
            p.out_filter_intron_motifs
        ));
    }

    let i_wrec_usize = i_wrec as usize;
    if read_align.tr_all.len() <= i_wrec_usize {
        read_align.tr_all.resize(i_wrec_usize + 1, Vec::new());
    }
    let record_index = if wa_excl.is_none() { 0 } else { 1 };
    if read_align.tr_all[i_wrec_usize].len() <= record_index {
        read_align.tr_all[i_wrec_usize].resize(
            record_index + 1,
            crate::transcript::Transcript::default(),
        );
    }
    read_align.tr_all[i_wrec_usize][record_index] = tr_a.clone();
    if read_align.n_win_tr.len() <= i_wrec_usize {
        read_align.n_win_tr.resize(i_wrec_usize + 1, 0);
    }
    read_align.n_win_tr[i_wrec_usize] = if wa_excl.is_none() { 1 } else { 2 };
    read_align.tr_best = tr_a;
    Ok(true)
}
