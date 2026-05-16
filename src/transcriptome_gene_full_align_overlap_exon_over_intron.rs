#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Transcriptome::geneFullAlignOverlap_ExonOverIntron` at STAR/source/Transcriptome_geneFullAlignOverlap_ExonOverIntron.cpp:5. Args: nA: uint, aAll: Transcript, strandType: int32, annFeat: ReadAnnotFeature, annFeatGeneConcordant: ReadAnnotFeature"]
pub fn transcriptome_genefullalignoverlap_exonoverintron_l5_transcriptome_genefullalignoverlap_exonoverintron(
    transcriptome: &crate::transcriptome::Transcriptome,
    n_a: u32,
    a_all: &[crate::transcript::Transcript],
    strand_type: i32,
    ann_feat: &mut crate::read_annotations::ReadAnnotFeature,
    ann_feat_gene_concordant: &crate::read_annotations::ReadAnnotFeature,
) {
    if !ann_feat_gene_concordant.f_set.is_empty() {
        *ann_feat = ann_feat_gene_concordant.clone();
        ann_feat.ov_type = 1;
        return;
    }

    ann_feat.f_align.resize(n_a as usize, Default::default());
    for i_a in 0..n_a as usize {
        let a = &a_all[i_a];
        let a_s = a.exons[0][EX_G] as u64;
        let last_exon = &a.exons[a.n_exons as usize - 1];
        let a_e = last_exon[EX_G] as u64 + last_exon[EX_L] as u64 - 1;

        let mut gi1 = servicefuns_l239_binarysearch1a(
            a_s,
            &transcriptome.gene_full.s,
            transcriptome.n_ge as i32,
        );

        while gi1 >= 0 && transcriptome.gene_full.e_max[gi1 as usize] >= a_e {
            let gi = gi1 as usize;
            if transcriptome.gene_full.e[gi] >= a_e {
                let str1 = if transcriptome.gene_full.str_[gi] == 1 {
                    a.str_ as i32
                } else {
                    1 - a.str_ as i32
                };
                if strand_type == -1 || strand_type == str1 {
                    ann_feat.f_set.insert(transcriptome.gene_full.g[gi]);
                    ann_feat.f_align[i_a].insert(transcriptome.gene_full.g[gi]);
                }
            }
            gi1 -= 1;
        }
    }
    if !ann_feat.f_set.is_empty() {
        ann_feat.ov_type = 5;
    }
}
