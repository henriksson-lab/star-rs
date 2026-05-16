#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `alignToTranscript` at STAR/source/Transcriptome_quantAlign.cpp:5. Args: aG: Transcript, trS1: uint, trStr1: uint8, exSE1: uint32, exLenCum1: uint32, exN1: uint16, aT: Transcript"]
pub fn transcriptome_quantalign_l5_aligntotranscript(
    a_g: &mut crate::transcript::Transcript,
    tr_s1: u32,
    tr_str1: u8,
    ex_se1: &[u32],
    ex_len_cum1: &[u32],
    ex_n1: u16,
    a_t: &mut crate::transcript::Transcript,
) -> i32 {
    let g1 = a_g.exons[0][EX_G] - tr_s1;
    let mut ex1 = servicefuns_l192_binarysearch1(g1, ex_se1, 2 * ex_n1 as u32);
    if ex1 >= 2 * ex_n1 as u32 {
        return 0;
    }

    if ex1 % 2 == 1 {
        if ex_se1[ex1 as usize] == g1 {
            ex1 -= 1;
        } else {
            return 0;
        }
    }
    ex1 /= 2;

    a_t.n_exons = 0;
    a_t.primary_flag = false;

    if a_g.canon_sj.len() < a_g.n_exons as usize {
        a_g.canon_sj.resize(a_g.n_exons as usize, 0);
    }
    a_g.canon_sj[a_g.n_exons as usize - 1] = -999;

    for iab in 0..a_g.n_exons as usize {
        if a_g.exons[iab][EX_G] + a_g.exons[iab][EX_L] > ex_se1[2 * ex1 as usize + 1] + tr_s1 + 1 {
            return 0;
        }

        if iab == 0 || a_g.canon_sj[iab - 1] < 0 {
            let at_exon = a_t.n_exons as usize;
            if a_t.exons.len() <= at_exon {
                a_t.exons.resize(at_exon + 1, [0; EX_SIZE]);
            }
            if a_t.canon_sj.len() <= at_exon {
                a_t.canon_sj.resize(at_exon + 1, 0);
            }
            a_t.exons[at_exon][EX_R] = a_g.exons[iab][EX_R];
            a_t.exons[at_exon][EX_G] =
                a_g.exons[iab][EX_G] - tr_s1 - ex_se1[2 * ex1 as usize] + ex_len_cum1[ex1 as usize];
            a_t.exons[at_exon][EX_L] = a_g.exons[iab][EX_L];
            a_t.exons[at_exon][EX_IFRAG] = a_g.exons[iab][EX_IFRAG];
            if a_t.n_exons > 0 {
                a_t.canon_sj[at_exon - 1] = a_g.canon_sj[iab - 1];
            }
            a_t.n_exons += 1;
        } else {
            let at_exon = a_t.n_exons as usize - 1;
            a_t.exons[at_exon][EX_L] += a_g.exons[iab][EX_L];
        }

        match a_g.canon_sj[iab] {
            -999 => {
                if tr_str1 == 2 {
                    let tr_length = ex_len_cum1[ex_n1 as usize - 1]
                        + ex_se1[2 * ex_n1 as usize - 1]
                        - ex_se1[2 * ex_n1 as usize - 2]
                        + 1;
                    for iex in 0..a_t.n_exons as usize {
                        a_t.exons[iex][EX_R] =
                            a_g.l_read - (a_t.exons[iex][EX_R] + a_t.exons[iex][EX_L]);
                        a_t.exons[iex][EX_G] =
                            tr_length - (a_t.exons[iex][EX_G] + a_t.exons[iex][EX_L]);
                    }
                    for iex in 0..a_t.n_exons as usize / 2 {
                        let iex2 = a_t.n_exons as usize - 1 - iex;
                        a_t.exons.swap(iex, iex2);
                    }
                    for iex in 0..(a_t.n_exons as usize - 1) / 2 {
                        let iex2 = a_t.n_exons as usize - 2 - iex;
                        a_t.canon_sj.swap(iex, iex2);
                    }
                }

                let n_at = a_t.n_exons as usize;
                a_t.sj_annot.resize(n_at, 0);
                a_t.shift_sj.resize(n_at, [0; 2]);
                a_t.sj_str.resize(n_at, 0);
                for iex in 0..n_at {
                    a_t.sj_annot[iex] = 0;
                    a_t.shift_sj[iex] = [0, 0];
                    a_t.sj_str[iex] = 0;
                }
                return 1;
            }
            -3 => {
                ex1 = servicefuns_l192_binarysearch1(
                    a_g.exons[iab + 1][EX_G] - tr_s1,
                    ex_se1,
                    2 * ex_n1 as u32,
                );
                if ex1 % 2 == 1 {
                    return 0;
                }
                ex1 /= 2;
            }
            -2 | -1 => {}
            _ => {
                if a_g.exons[iab][EX_G] + a_g.exons[iab][EX_L]
                    == ex_se1[2 * ex1 as usize + 1] + tr_s1 + 1
                    && a_g.exons[iab + 1][EX_G] == ex_se1[2 * (ex1 as usize + 1)] + tr_s1
                {
                    ex1 += 1;
                } else {
                    return 0;
                }
            }
        }
    }
    0
}

#[doc = "Original `Transcriptome::quantAlign` at STAR/source/Transcriptome_quantAlign.cpp:91. Args: aG: Transcript, aTall: Transcript"]
pub fn transcriptome_quantalign_l91_transcriptome_quantalign(
    transcriptome: &crate::transcriptome::Transcriptome,
    a_g: &mut crate::transcript::Transcript,
    a_tall: &mut Vec<crate::transcript::Transcript>,
) -> u32 {
    let mut n_atr = 0_u32;

    let mut tr1 = servicefuns_l239_binarysearch1a(
        a_g.exons[0][EX_G],
        &transcriptome.tr_s,
        transcriptome.n_tr as i32,
    );
    if tr1 == -1 {
        return 0;
    }

    let a_g_end = a_g.exons[a_g.n_exons as usize - 1][EX_G];

    tr1 += 1;
    loop {
        tr1 -= 1;
        let tr = tr1 as usize;
        if a_g_end <= transcriptome.tr_e[tr] {
            if a_tall.len() <= n_atr as usize {
                a_tall.resize(n_atr as usize + 1, Default::default());
            }
            let ex_i = transcriptome.tr_ex_i[tr] as usize;
            let ex_se_start = 2 * ex_i;
            let status = transcriptome_quantalign_l5_aligntotranscript(
                a_g,
                transcriptome.tr_s[tr],
                transcriptome.tr_str[tr],
                &transcriptome.ex_se[ex_se_start..],
                &transcriptome.ex_len_cum[ex_i..],
                transcriptome.tr_ex_n[tr],
                &mut a_tall[n_atr as usize],
            );
            if status == 1 {
                a_tall[n_atr as usize].chr = tr1 as u32;
                a_tall[n_atr as usize].str_ = if transcriptome.tr_str[tr] == 1 {
                    a_g.str_
                } else {
                    1 - a_g.str_
                };
                n_atr += 1;
            }
        }

        if !(transcriptome.tr_e_max[tr] >= a_g_end && tr1 > 0) {
            break;
        }
    }

    n_atr
}
