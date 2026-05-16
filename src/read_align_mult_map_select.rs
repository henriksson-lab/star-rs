#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::multMapSelect` at STAR/source/ReadAlign_multMapSelect.cpp:8. Args: "]
pub fn readalign_multmapselect_l8_readalign_multmapselect(
    read_align: &mut crate::read_align::ReadAlign,
    map_gen: &crate::genome::Genome,
    tr_all: &[Vec<crate::transcript::Transcript>],
    out_filter_multimap_score_range: i32,
    out_filter_multimap_nmax: u64,
    out_multimapper_order_random: bool,
    out_sam_mult_nmax_is_limited: bool,
    out_sam_primary_flag: &str,
    rng_uniform_real_0_to_1: &[f64],
) -> Result<Vec<crate::transcript::Transcript>, String> {
    readalign_multmapselect_inner(
        &mut read_align.n_tr,
        &mut read_align.tr_best,
        read_align.l_read,
        read_align.n_w,
        &read_align.n_win_tr,
        map_gen,
        tr_all,
        out_filter_multimap_score_range,
        out_filter_multimap_nmax,
        out_multimapper_order_random,
        out_sam_mult_nmax_is_limited,
        out_sam_primary_flag,
        rng_uniform_real_0_to_1,
    )
}

pub fn readalign_multmapselect_inner(
    n_tr_out: &mut u64,
    tr_best: &mut crate::transcript::Transcript,
    l_read: u64,
    n_w: u64,
    n_win_tr_all: &[u64],
    map_gen: &crate::genome::Genome,
    tr_all: &[Vec<crate::transcript::Transcript>],
    out_filter_multimap_score_range: i32,
    out_filter_multimap_nmax: u64,
    out_multimapper_order_random: bool,
    out_sam_mult_nmax_is_limited: bool,
    out_sam_primary_flag: &str,
    rng_uniform_real_0_to_1: &[f64],
) -> Result<Vec<crate::transcript::Transcript>, String> {
    *n_tr_out = 0;
    if n_w == 0 {
        return Ok(Vec::new());
    }

    let mut max_score = -10 * l_read as i32;
    for iw in 0..n_w as usize {
        if max_score < tr_all[iw][0].max_score {
            max_score = tr_all[iw][0].max_score;
        }
    }

    if max_score != tr_best.max_score {
        return Err("BUG: maxScore!=trBest->maxScore in multMapSelect".to_string());
    }

    let mut tr_mult = Vec::new();
    for iw in 0..n_w as usize {
        let n_win_tr = (n_win_tr_all.get(iw).copied().unwrap_or(0) as usize).min(tr_all[iw].len());
        for itr in 0..n_win_tr {
            if tr_all[iw][itr].max_score + out_filter_multimap_score_range >= max_score {
                if tr_mult.len() == 100_000 {
                    return Err("EXITING: Fatal ERROR: number of alignments exceeds MAX_N_MULTMAP, increase it and re-compile STAR".to_string());
                }
                let mut tr = tr_all[iw][itr].clone();
                tr.chr = tr_all[iw][0].chr;
                tr.str_ = tr_all[iw][0].str_;
                tr.ro_str = tr_all[iw][0].ro_str;
                tr_mult.push(tr);
                *n_tr_out += 1;
            }
        }
    }
    if *n_tr_out > out_filter_multimap_nmax || *n_tr_out == 0 {
        return Ok(tr_mult);
    }

    for tr in tr_mult.iter_mut() {
        tr.ro_start = if tr.ro_str == 0 {
            tr.r_start
        } else {
            l_read.wrapping_sub(tr.r_start).wrapping_sub(tr.r_length)
        };
        tr.c_start = tr
            .g_start
            .wrapping_sub(map_gen.chr_start[tr.chr as usize]);
    }
    tr_best.ro_start = if tr_best.ro_str == 0 {
        tr_best.r_start
    } else {
        l_read
            .wrapping_sub(tr_best.r_start)
            .wrapping_sub(tr_best.r_length)
    };
    tr_best.c_start = tr_best
        .g_start
        .wrapping_sub(map_gen.chr_start[tr_best.chr as usize]);

    if *n_tr_out == 1 {
        tr_mult[0].primary_flag = true;
    } else {
        let mut nbest = 0usize;
        if out_multimapper_order_random || out_sam_mult_nmax_is_limited {
            for itr in 0..*n_tr_out as usize {
                if tr_mult[itr].max_score == max_score {
                    tr_mult.swap(itr, nbest);
                    nbest += 1;
                }
            }
        }

        if out_multimapper_order_random {
            let mut irng = 0usize;
            for itr in (1..nbest).rev() {
                let rand1 = (rng_uniform_real_0_to_1[irng] * itr as f64 + 0.5) as usize;
                irng += 1;
                tr_mult.swap(itr, rand1);
            }
            for itr in (1..(*n_tr_out as usize - nbest)).rev() {
                let rand1 = (rng_uniform_real_0_to_1[irng] * itr as f64 + 0.5) as usize;
                irng += 1;
                tr_mult.swap(nbest + itr, nbest + rand1);
            }
        }

        if out_sam_primary_flag == "AllBestScore" {
            for itr in 0..*n_tr_out as usize {
                if tr_mult[itr].max_score == max_score {
                    tr_mult[itr].primary_flag = true;
                }
            }
        } else if out_multimapper_order_random || out_sam_mult_nmax_is_limited {
            tr_mult[0].primary_flag = true;
        } else {
            tr_best.primary_flag = true;
            let mut best_for_match = tr_best.clone();
            best_for_match.primary_flag = false;
            let mut marked_primary = false;
            for tr in tr_mult.iter_mut() {
                let mut tr_for_match = tr.clone();
                tr_for_match.primary_flag = false;
                if tr_for_match == best_for_match {
                    tr.primary_flag = true;
                    marked_primary = true;
                    break;
                }
            }
            if !marked_primary {
                if let Some(first) = tr_mult.first_mut() {
                    first.primary_flag = true;
                }
            }
        }
    }

    Ok(tr_mult)
}

pub fn expand_collapsed_same_locus_paired_multimaps(
    read_align: &mut crate::read_align::ReadAlign,
    map_gen: &crate::genome::Genome,
    tr_mult: &mut Vec<crate::transcript::Transcript>,
    max_score: i32,
    out_multimapper_order_random: bool,
    out_sam_mult_nmax_is_limited: bool,
    out_sam_primary_flag: &str,
    rng_uniform_real_0_to_1: &[f64],
) {
    if read_align.n_tr <= 1
        || tr_mult.len() != read_align.n_tr as usize
        || !tr_mult.iter().all(|tr| {
            tr.n_exons == 2
                && tr.exons[0][EX_IFRAG] != tr.exons[1][EX_IFRAG]
                && tr.exons[0][EX_L] == read_align.read_length[tr.exons[0][EX_IFRAG] as usize]
                && tr.exons[1][EX_L] == read_align.read_length[tr.exons[1][EX_IFRAG] as usize]
        })
    {
        return;
    }

    let same_locus_collapsed = tr_mult
        .iter()
        .all(|tr| tr.exons[0][EX_G] == tr.exons[1][EX_G]);
    let mut loci: Vec<u64> = tr_mult
        .iter()
        .flat_map(|tr| [tr.exons[0][EX_G], tr.exons[1][EX_G]])
        .collect();
    loci.sort_unstable();
    loci.dedup();
    if loci.is_empty() {
        return;
    }

    let base = tr_mult[0].clone();
    let mate0_exon = if base.exons[0][EX_IFRAG] == 0 { 0 } else { 1 };
    let mate1_exon = 1 - mate0_exon;
    let mate0_template = base.exons[mate0_exon];
    let mate1_template = base.exons[mate1_exon];
    let mate0_len = base.exons[mate0_exon][EX_L];
    let mate1_len = base.exons[mate1_exon][EX_L];

    let triangular_count = loci.len() * (loci.len() + 1) / 2;
    if !same_locus_collapsed && tr_mult.len() == triangular_count {
        for tr in tr_mult.iter_mut() {
            tr.primary_flag = false;
        }
        tr_mult.reserve(triangular_count);
        for delta in 0..loci.len() {
            for i in (0..(loci.len() - delta)).rev() {
                let j = i + delta;
                let mut tr = base.clone();
                tr.primary_flag = false;
                tr.str_ = 1;
                tr.ro_str = 1;
                tr.exons[0] = mate1_template;
                tr.exons[1] = mate0_template;
                tr.exons[0][EX_G] = loci[i];
                tr.exons[1][EX_G] = loci[j];
                tr.g_start = loci[i];
                tr.g_length = loci[j] + mate0_len - loci[i];
                tr.c_start = tr
                    .g_start
                    .wrapping_sub(map_gen.chr_start[tr.chr as usize]);
                tr.r_start = tr.exons[0][EX_R];
                tr.r_length = mate0_len + mate1_len;
                tr.mapped_length = tr.r_length;
                tr_mult.push(tr);
            }
        }

        read_align.n_tr = tr_mult.len() as u64;
        funprimaryalignmark_l3_funprimaryalignmark(
            tr_mult,
            read_align.n_tr,
            max_score,
            out_multimapper_order_random,
            out_sam_mult_nmax_is_limited,
            out_sam_primary_flag,
            rng_uniform_real_0_to_1,
        );
        return;
    }

    if !same_locus_collapsed || loci.len() != read_align.n_tr as usize {
        return;
    }

    let mut expanded = Vec::with_capacity(loci.len() * (loci.len() + 1));

    for delta in 0..loci.len() {
        for i in 0..(loci.len() - delta) {
            let j = i + delta;
            let mut tr = base.clone();
            tr.primary_flag = false;
            tr.str_ = 0;
            tr.ro_str = 0;
            tr.exons[mate0_exon][EX_G] = loci[i];
            tr.exons[mate1_exon][EX_G] = loci[j];
            tr.g_start = loci[i];
            tr.g_length = loci[j] + mate1_len - loci[i];
            tr.c_start = tr
                .g_start
                .wrapping_sub(map_gen.chr_start[tr.chr as usize]);
            tr.r_start = tr.exons[0][EX_R];
            tr.r_length = mate0_len + mate1_len;
            tr.mapped_length = tr.r_length;
            expanded.push(tr);
        }
    }

    for delta in 0..loci.len() {
        for i in (0..(loci.len() - delta)).rev() {
            let j = i + delta;
            let mut tr = base.clone();
            tr.primary_flag = false;
            tr.str_ = 1;
            tr.ro_str = 1;
            tr.exons[0] = mate1_template;
            tr.exons[1] = mate0_template;
            tr.exons[0][EX_G] = loci[i];
            tr.exons[1][EX_G] = loci[j];
            tr.g_start = loci[i];
            tr.g_length = loci[j] + mate0_len - loci[i];
            tr.c_start = tr
                .g_start
                .wrapping_sub(map_gen.chr_start[tr.chr as usize]);
            tr.r_start = tr.exons[0][EX_R];
            tr.r_length = mate0_len + mate1_len;
            tr.mapped_length = tr.r_length;
            expanded.push(tr);
        }
    }

    read_align.n_tr = expanded.len() as u64;
    *tr_mult = expanded;
    funprimaryalignmark_l3_funprimaryalignmark(
        tr_mult,
        read_align.n_tr,
        max_score,
        out_multimapper_order_random,
        out_sam_mult_nmax_is_limited,
        out_sam_primary_flag,
        rng_uniform_real_0_to_1,
    );
}
