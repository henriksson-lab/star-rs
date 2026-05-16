#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::transformGenome` at STAR/source/ReadAlign_transformGenome.cpp:5. Args: "]
pub fn readalign_transformgenome_l5_readalign_transformgenome(
    read_align: &crate::read_align::ReadAlign,
    map_gen: &crate::genome::Genome,
    tr_mult: &mut [crate::transcript::Transcript],
    out_filter_multimap_nmax: u64,
    max_score: i32,
    out_multimapper_order_random: bool,
    out_sam_mult_nmax_is_limited: bool,
    out_sam_primary_flag: &str,
    rng_uniform_real_0_to_1: &[f64],
) -> crate::quantifications::ReadAlignGenomeTransformResult {
    // STAR's C++ does `alignsGenOut.alBest = trBest` (pointer assignment).
    // The Rust port has to copy. Defer the copy past the early-exit check so
    // that reads which skip transformation (the common case when no genome
    // transformation is configured) don't pay for it.
    let mut aligns_gen_out = crate::quantifications::ReadAlignGenomeTransformResult::default();

    if !map_gen.genome_out.conv_yes
        || map_gen.p_ge.transform.type_ == 0
        || read_align.n_tr > out_filter_multimap_nmax
        || read_align.n_tr == 0
    {
        aligns_gen_out.al_best.copy_from(&read_align.tr_best);
        return aligns_gen_out;
    }

    aligns_gen_out.al_best.copy_from(&read_align.tr_best);

    for i_tr in 0..read_align.n_tr as usize {
        let is_best = tr_mult[i_tr] == read_align.tr_best;
        tr_mult[i_tr].haplo_type = if tr_mult[i_tr].chr < (map_gen.n_chr_real / 2) as u64 {
            1
        } else {
            2
        };

        let mut converted = tr_mult[i_tr].clone();
        if transcript_transformgenome_l4_transcript_transformgenome(
            &tr_mult[i_tr],
            map_gen,
            &mut converted,
        ) {
            if is_best {
                aligns_gen_out.al_best = converted.clone();
            }
            aligns_gen_out.al_mult.push(converted);
        }
    }

    if map_gen.p_ge.transform.type_ == 2 {
        let n_tr1 = aligns_gen_out.al_mult.len();
        let mut keep_tr = vec![true; n_tr1];

        for ia1 in 0..n_tr1 {
            if !keep_tr[ia1] {
                continue;
            }
            for ia2 in ia1 + 1..n_tr1 {
                if !keep_tr[ia1] {
                    continue;
                }

                let a1 = aligns_gen_out.al_mult[ia1].clone();
                let a2 = aligns_gen_out.al_mult[ia2].clone();
                let a1_end_ex = a1
                    .exons
                    .get(a1.n_exons as usize)
                    .unwrap_or_else(|| a1.exons.last().unwrap());
                let a2_end_ex = a2
                    .exons
                    .get(a2.n_exons as usize)
                    .unwrap_or_else(|| a2.exons.last().unwrap());

                if a1.chr == a2.chr
                    && a1.str_ == a2.str_
                    && a1.exons[0][EX_G].wrapping_sub(a1.exons[0][EX_R])
                        == a2.exons[0][EX_G].wrapping_sub(a2.exons[0][EX_R])
                    && a1_end_ex[EX_G] + a1_end_ex[EX_L] - a1_end_ex[EX_R]
                        == a2_end_ex[EX_G] + a2_end_ex[EX_L] - a2_end_ex[EX_R]
                {
                    aligns_gen_out.al_mult[ia1].haplo_type = 0;
                    aligns_gen_out.al_mult[ia2].haplo_type = 0;
                    if a1.max_score > a2.max_score {
                        keep_tr[ia2] = false;
                    } else {
                        keep_tr[ia1] = false;
                    }
                }
            }
        }

        let mut kept = Vec::new();
        for (ia1, keep) in keep_tr.into_iter().enumerate() {
            if keep {
                kept.push(aligns_gen_out.al_mult[ia1].clone());
            }
        }
        aligns_gen_out.al_mult = kept;
        aligns_gen_out.al_n = aligns_gen_out.al_mult.len() as u64;
    } else {
        aligns_gen_out.al_n = aligns_gen_out.al_mult.len() as u64;
    }

    funprimaryalignmark_l3_funprimaryalignmark(
        &mut aligns_gen_out.al_mult,
        aligns_gen_out.al_n as u64,
        max_score,
        out_multimapper_order_random,
        out_sam_mult_nmax_is_limited,
        out_sam_primary_flag,
        rng_uniform_real_0_to_1,
    );

    aligns_gen_out
}
