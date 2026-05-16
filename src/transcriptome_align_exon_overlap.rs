#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original struct `GeneInfo1` at STAR/source/Transcriptome_alignExonOverlap.cpp:25."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GeneInfo1 {
    pub g: u32,
    pub ia: u32,
    pub ot: [bool; 6],
}

#[doc = "Original struct `GeneStrOverlapAlign` at STAR/source/Transcriptome_alignExonOverlap.cpp:12."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GeneStrOverlapAlign {
    pub g: u32,
    pub ov: i64,
    pub ia: u32,
    pub exl: u32,
    pub str: bool,
    pub sjc: bool,
}

#[doc = "Original `Transcriptome::alignExonOverlap` at STAR/source/Transcriptome_alignExonOverlap.cpp:10. Args: nA: uint, aAll: Transcript, strandType: int32, annFeat: ReadAnnotFeature"]
pub fn transcriptome_alignexonoverlap_l10_transcriptome_alignexonoverlap(
    transcriptome: &crate::transcriptome::Transcriptome,
    n_a: u32,
    a_all: &[crate::transcript::Transcript],
    strand_type: i32,
    ann_feat: &mut crate::read_annotations::ReadAnnotFeature,
) {
    let mut v_gene_info1 = Vec::<crate::transcriptome_align_exon_overlap::GeneInfo1>::new();
    let ot_as = [false, true, false, true, false, true];

    for iag in 0..n_a as usize {
        let a_g = &a_all[iag];
        let a_g_start = a_g.exons[0][EX_G] as u64;
        let last_exon = &a_g.exons[a_g.n_exons as usize - 1];
        let a_g_end = last_exon[EX_G] as u64 + last_exon[EX_L] as u64 - 1;

        let mut tr1 = servicefuns_l239_binarysearch1a(
            a_g_start as u32,
            &transcriptome.tr_s,
            transcriptome.n_tr as i32,
        );
        if tr1 == -1 {
            continue;
        }
        tr1 += 1;
        loop {
            tr1 -= 1;
            let tr = tr1 as usize;
            if a_g_end <= transcriptome.tr_e[tr] as u64 {
                let mut str1 = if strand_type == 0 {
                    a_g.str_ as i32
                } else {
                    1 - a_g.str_ as i32
                } == transcriptome.tr_str[tr] as i32 - 1;
                str1 = str1 || strand_type == -1;

                let ex_start = 2 * transcriptome.tr_ex_i[tr] as usize;
                let ex_len = 2 * transcriptome.tr_ex_n[tr] as usize;
                let (n_overlap, sj_concord) =
                    transcriptome_alignexonoverlap_l236_alignblocksoverlapexons(
                        a_g,
                        transcriptome.tr_ex_n[tr],
                        &transcriptome.ex_se[ex_start..ex_start + ex_len],
                        transcriptome.tr_s[tr] as u64,
                    );
                if n_overlap >= 0 {
                    let mut exl = 0_i32;
                    for iex in 0..a_g.n_exons as usize {
                        exl += a_g.exons[iex][EX_L] as i32;
                    }
                    v_gene_info1.push(crate::transcriptome_align_exon_overlap::GeneInfo1 {
                        g: transcriptome.tr_gene[tr],
                        ia: iag as u32,
                        ot: [
                            str1 && n_overlap == exl && sj_concord,
                            !str1 && n_overlap == exl && sj_concord,
                            str1 && n_overlap > exl / 2,
                            !str1 && n_overlap > exl / 2,
                            str1,
                            !str1,
                        ],
                    });
                }
            }

            if !(transcriptome.tr_e_max[tr] as u64 >= a_g_end && tr1 > 0) {
                break;
            }
        }
    }

    let mut ot_final = [false; 6];
    for v1 in &v_gene_info1 {
        for it in 0..ot_final.len() {
            if v1.ot[it] {
                ot_final[it] = true;
                break;
            }
        }
    }

    ann_feat.ov_type = if ot_final[0] {
        1
    } else if ot_final[1] {
        2
    } else if ot_final[2] {
        3
    } else if ot_final[3] {
        4
    } else if ot_final[4] {
        5
    } else if ot_final[5] {
        6
    } else {
        7
    };

    ann_feat.f_align.resize(n_a as usize, Default::default());
    for it in 0..ot_final.len() {
        if ot_final[it] {
            if ot_as[it] {
                return;
            }
            for v1 in &v_gene_info1 {
                if v1.ot[it] {
                    ann_feat.f_set.insert(v1.g);
                    ann_feat.f_align[v1.ia as usize].insert(v1.g);
                }
            }
            break;
        }
    }
}

#[doc = "Original `alignBlocksOverlapExons` at STAR/source/Transcriptome_alignExonOverlap.cpp:236. Args: aG: Transcript, exN1: uint16, exSE1: uint32, trStart1: uint64, sjConcord: bool"]
pub fn transcriptome_alignexonoverlap_l236_alignblocksoverlapexons(
    a_g: &crate::transcript::Transcript,
    ex_n1: u16,
    ex_se1: &[u32],
    tr_start1: u64,
) -> (i32, bool) {
    let mut i1 = 0_u64;
    let mut i2 = 0_u64;
    let mut n_overlap = 0_i32;
    let mut sj_concord = true;
    let tr_end1 = tr_start1 + ex_se1[2 * ex_n1 as usize - 1] as u64 + 1;

    while i1 < a_g.n_exons && i2 < ex_n1 as u64 {
        let rs1 = a_g.exons[i1 as usize][EX_G] as u64;
        let re1 = a_g.exons[i1 as usize][EX_G] as u64 + a_g.exons[i1 as usize][EX_L] as u64;
        let rs2 = tr_start1 + ex_se1[2 * i2 as usize] as u64;
        let re2 = tr_start1 + ex_se1[2 * i2 as usize + 1] as u64 + 1;

        if rs1 < tr_start1 || re1 > tr_end1 {
            return (-1, sj_concord);
        }

        if rs1 >= re2 {
            i2 += 1;
            if i1 > 0 && a_g.canon_sj[i1 as usize - 1] >= 0 {
                sj_concord = false;
            }
        } else if rs2 >= re1 {
            i1 += 1;
            sj_concord = false;
        } else {
            n_overlap += std::cmp::min(re1, re2) as i32 - std::cmp::max(rs1, rs2) as i32;

            if i1 > 0 && rs1 != rs2 && a_g.canon_sj[i1 as usize - 1] >= 0 {
                sj_concord = false;
            }
            if i1 < a_g.n_exons - 1 && re1 != re2 && a_g.canon_sj[i1 as usize] >= 0 {
                sj_concord = false;
            }

            if re1 >= re2 {
                i2 += 1;
            }
            if re2 >= re1 {
                i1 += 1;
            }
        }
    }
    (n_overlap, sj_concord)
}
