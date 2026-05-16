#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Transcriptome::geneFullAlignOverlap` at STAR/source/Transcriptome_geneFullAlignOverlap.cpp:5. Args: nA: uint, aAll: Transcript, strandType: int32, annFeat: ReadAnnotFeature"]
pub fn transcriptome_genefullalignoverlap_l5_transcriptome_genefullalignoverlap(
    transcriptome: &crate::transcriptome::Transcriptome,
    n_a: u64,
    a_all: &[crate::transcript::Transcript],
    strand_type: i32,
    ann_feat: &mut crate::read_annotations::ReadAnnotFeature,
) {
    ann_feat.f_align.resize(n_a as usize, Default::default());
    for i_a in 0..n_a as usize {
        let a = &a_all[i_a];
        let mut ib = a.n_exons as i64 - 1;
        while ib >= 0 {
            let exon = a.exons[ib as usize];
            let be1 = exon[EX_G] as u64 + exon[EX_L] as u64 - 1;
            let mut gi1 = servicefuns_l239_binarysearch1a(
                be1,
                &transcriptome.gene_full.s,
                transcriptome.n_ge as i32,
            );

            while gi1 >= 0 && transcriptome.gene_full.e_max[gi1 as usize] >= exon[EX_G] as u64 {
                let gi = gi1 as usize;
                if transcriptome.gene_full.e[gi] >= exon[EX_G] as u64 {
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
            ib -= 1;
        }
    }
}
