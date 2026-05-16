#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `alignToTranscript` at STAR/source/Transcriptome_classifyAlign.cpp:8. Args: aG: Transcript, trS1: uint64, exN1: uint16, exSE1: uint32, exLenCum1: uint32, distTrEnds: array<uint32,2>"]
pub fn transcriptome_classifyalign_l8_aligntotranscript(
    a_g: &crate::transcript::Transcript,
    tr_s1: u64,
    ex_n1: u16,
    ex_se1: &[u32],
    ex_len_cum1: &[u32],
    dist_tr_ends: &mut [u32; 2],
) -> i32 {
    let mut align_intronic = false;
    let mut align_exonic = false;
    let mut align_spans_exon_intr = false;
    let mut align_sj_concordant = true;

    let mut ex1 = 0_u32;
    let mut b_e = 0_u32;
    let mut e_e = 0_u32;
    let mut en_s = 0_u32;

    for iab in 0..a_g.n_exons {
        let b_e_prev = b_e;

        if u64::from(a_g.exons[iab as usize][EX_G]) < tr_s1 {
            return -1;
        }

        let b_s = (u64::from(a_g.exons[iab as usize][EX_G]) - tr_s1) as u32;
        b_e = b_s + a_g.exons[iab as usize][EX_L] - 1;

        if iab == 0 || a_g.canon_sj[iab as usize - 1] == -3 {
            if !servicefuns_l212_binarysearch_leleft(b_s, ex_se1, 2 * ex_n1 as u32, &mut ex1) {
                return -1;
            }
            ex1 /= 2;
        } else if a_g.canon_sj[iab as usize - 1] >= 0 {
            if b_e_prev == e_e && b_s == en_s {
                ex1 += 1;
            } else {
                align_sj_concordant = false;
                if !servicefuns_l212_binarysearch_leleft(b_s, ex_se1, 2 * ex_n1 as u32, &mut ex1) {
                    return -1;
                }
                ex1 /= 2;
            }
        }

        e_e = ex_se1[2 * ex1 as usize + 1];
        en_s = if ex1 + 1 < ex_n1 as u32 {
            ex_se1[2 * (ex1 as usize + 1)]
        } else {
            0
        };

        if b_s <= e_e {
            if b_e > e_e {
                align_spans_exon_intr = true;
            }
            align_exonic = true;

            if iab == 0 {
                dist_tr_ends[0] = ex_len_cum1[ex1 as usize] + b_s - ex_se1[2 * ex1 as usize];
            }
            let transcript_tail = if ex1 == ex_n1 as u32 - 1 {
                0
            } else {
                ex_se1[2 * ex_n1 as usize - 1] - ex_se1[2 * ex_n1 as usize - 2]
                    + 1
                    + ex_len_cum1[ex_n1 as usize - 1]
                    - ex_len_cum1[ex1 as usize + 1]
            };
            dist_tr_ends[1] = e_e.wrapping_sub(b_e).wrapping_add(transcript_tail);
        } else {
            if b_e >= en_s {
                align_spans_exon_intr = true;
            }
            align_intronic = true;
        }
    }

    if !align_sj_concordant {
        return -1;
    }

    if align_spans_exon_intr {
        ALIGN_VS_TRANSCRIPT_EXON_INTRON_SPAN
    } else if !align_intronic {
        ALIGN_VS_TRANSCRIPT_CONCORDANT
    } else if align_exonic {
        ALIGN_VS_TRANSCRIPT_EXON_INTRON
    } else {
        ALIGN_VS_TRANSCRIPT_INTRON
    }
}

#[doc = "Original `alignToTranscriptMinOverlap` at STAR/source/Transcriptome_classifyAlign.cpp:93. Args: aG: Transcript, trS1: uint, exSE1: uint32, exN1: uint16, minOverlapMinusOne: uint32"]
pub fn transcriptome_classifyalign_l93_aligntotranscriptminoverlap(
    a_g: &crate::transcript::Transcript,
    tr_s1: u32,
    ex_se1: &[u32],
    ex_n1: u16,
    min_overlap_minus_one: u32,
) -> i32 {
    let mut align_intronic = false;
    let mut align_exonic = false;
    let mut align_spans_exon_intr = false;
    let mut align_sj_concordant = true;

    let mut iab = 0_u32;
    while iab < a_g.n_exons {
        let b_s = a_g.exons[iab as usize][EX_G] - tr_s1;
        let ex1 = servicefuns_l192_binarysearch1(b_s, ex_se1, 2 * ex_n1 as u32) / 2;
        if ex1 == ex_n1 as u32 - 1 {
            align_exonic = true;
            break;
        }

        while iab < a_g.n_exons - 1
            && a_g.canon_sj[iab as usize] > -3
            && a_g.canon_sj[iab as usize] < 0
        {
            iab += 1;
        }

        let b_e = a_g.exons[iab as usize][EX_G] - tr_s1 + a_g.exons[iab as usize][EX_L] - 1;
        if b_e - b_s < min_overlap_minus_one {
            iab += 1;
            continue;
        }

        let e_e = ex_se1[2 * ex1 as usize + 1];
        let en_s = ex_se1[2 * ex1 as usize + 2];
        let en_e = ex_se1[2 * ex1 as usize + 3];

        if b_s + min_overlap_minus_one <= e_e {
            if b_e <= e_e + min_overlap_minus_one {
                align_exonic = true;
            } else {
                align_spans_exon_intr = true;
            }
        } else if b_s + min_overlap_minus_one < en_s {
            if b_e >= en_s + min_overlap_minus_one {
                align_spans_exon_intr = true;
            } else if b_e > e_e + min_overlap_minus_one {
                if en_s - e_e > 1_000_000 {
                    return -1;
                }
                align_intronic = true;
            }
        } else if b_e > en_e + min_overlap_minus_one {
            align_spans_exon_intr = true;
        } else if b_e >= en_s + min_overlap_minus_one {
            align_exonic = true;
        }

        if a_g.sj_yes && (align_intronic || align_spans_exon_intr) {
            align_sj_concordant = false;
            break;
        }

        iab += 1;
    }

    if !align_sj_concordant {
        return -1;
    }

    if align_spans_exon_intr {
        ALIGN_VS_TRANSCRIPT_EXON_INTRON_SPAN
    } else if !align_intronic {
        ALIGN_VS_TRANSCRIPT_CONCORDANT
    } else if align_exonic {
        ALIGN_VS_TRANSCRIPT_EXON_INTRON
    } else {
        ALIGN_VS_TRANSCRIPT_INTRON
    }
}

#[doc = "Original `Transcriptome::classifyAlign` at STAR/source/Transcriptome_classifyAlign.cpp:177. Args: alignG: Transcript, nAlignG: uint64, readAnnot: ReadAnnotations"]
pub fn transcriptome_classifyalign_l177_transcriptome_classifyalign(
    transcriptome: &crate::transcriptome::Transcriptome,
    p_solo: &crate::parameters_solo::ParametersSolo,
    align_g: &[crate::transcript::Transcript],
    read_annot: &mut crate::read_annotations::ReadAnnotations,
) {
    if read_annot.annot_features.len() <= SOLO_FEATURE_GENE as usize {
        read_annot
            .annot_features
            .resize(SOLO_FEATURE_GENE as usize + 1, Default::default());
    }

    let n_align_g = align_g.len();
    let ann_feat = &mut read_annot.annot_features[SOLO_FEATURE_GENE as usize];
    ann_feat.f_align.resize(n_align_g, Default::default());

    let mut re_ge = u32::MAX - 1;
    let mut re_ann = 0_u32;
    let velocyto_enabled = n_align_g == 1
        && (p_solo
            .feature_yes
            .get(SOLO_FEATURE_VELOCYTO as usize)
            .copied()
            .unwrap_or(false)
            || p_solo
                .feature_yes
                .get(SOLO_FEATURE_VELOCYTO_SIMPLE as usize)
                .copied()
                .unwrap_or(false));

    for (iag, a_g) in align_g.iter().enumerate() {
        let mut tr1 = servicefuns_l239_binarysearch1a(
            a_g.exons[0][EX_G],
            &transcriptome.tr_s,
            transcriptome.n_tr as i32,
        );
        if tr1 == -1 {
            continue;
        }

        let a_g_end = a_g.exons[a_g.n_exons as usize - 1][EX_G]
            + a_g.exons[a_g.n_exons as usize - 1][EX_L]
            - 1;

        tr1 += 1;
        loop {
            tr1 -= 1;
            let tr = tr1 as usize;
            let tr_stranded = if transcriptome.tr_str[tr] == 1 {
                a_g.str_
            } else {
                1 - a_g.str_
            };
            if a_g_end > transcriptome.tr_e[tr]
                || (p_solo.strand >= 0 && tr_stranded != p_solo.strand as u32)
            {
                if !(transcriptome.tr_e_max[tr] >= a_g_end && tr1 > 0) {
                    break;
                }
                continue;
            }

            let ex_i = transcriptome.tr_ex_i[tr] as usize;
            let mut dist_tr_ends = [0_u32; 2];
            let mut a_status = transcriptome_classifyalign_l8_aligntotranscript(
                a_g,
                transcriptome.tr_s[tr] as u64,
                transcriptome.tr_ex_n[tr],
                &transcriptome.ex_se[2 * ex_i..],
                &transcriptome.ex_len_cum[ex_i..],
                &mut dist_tr_ends,
            );

            if a_status == ALIGN_VS_TRANSCRIPT_CONCORDANT {
                let dist_tts = if transcriptome.tr_str[tr] == 1 {
                    dist_tr_ends[1]
                } else {
                    dist_tr_ends[0]
                };
                read_annot
                    .transcript_concordant
                    .push([tr1 as u32, dist_tts]);

                ann_feat.f_set.insert(transcriptome.tr_gene[tr]);
                ann_feat.f_align[iag].insert(transcriptome.tr_gene[tr]);
            }

            if velocyto_enabled {
                a_status = transcriptome_classifyalign_l93_aligntotranscriptminoverlap(
                    a_g,
                    transcriptome.tr_s[tr],
                    &transcriptome.ex_se[2 * ex_i..],
                    transcriptome.tr_ex_n[tr],
                    6,
                );

                if a_status >= 0 {
                    if re_ge != u32::MAX {
                        if re_ge == u32::MAX - 1 {
                            re_ge = transcriptome.tr_gene[tr];
                        }

                        if re_ge != transcriptome.tr_gene[tr] {
                            re_ge = u32::MAX;
                        } else if a_status != ALIGN_VS_TRANSCRIPT_EXON_INTRON_SPAN {
                            re_ann |= 1_u32 << ALIGN_VS_TRANSCRIPT_EXON_INTRON_SPAN;
                            re_ann |= 1_u32 << a_status;
                        }
                    }

                    let mut re_ann1 = 1_u8 << a_status;
                    if a_status == ALIGN_VS_TRANSCRIPT_EXON_INTRON_SPAN {
                        re_ann1 |= 1_u8 << ALIGN_VS_TRANSCRIPT_INTRON;
                        re_ann1 |= 1_u8 << ALIGN_VS_TRANSCRIPT_CONCORDANT;
                    }
                    read_annot
                        .tr_velocyto_type
                        .push(crate::solo_read_feature_record::TrTypeStruct {
                            tr: tr1 as u32,
                            type_: re_ann1,
                        });
                }
            }

            if !(transcriptome.tr_e_max[tr] >= a_g_end && tr1 > 0) {
                break;
            }
        }
    }

    if !ann_feat.f_set.is_empty() {
        ann_feat.ov_type = 1;
    }
    read_annot.gene_velocyto_simple[0] = if re_ge.wrapping_add(2) == 0 {
        u32::MAX
    } else {
        re_ge
    };
    read_annot.gene_velocyto_simple[1] = re_ann;
}
