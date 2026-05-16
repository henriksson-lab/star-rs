#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Transcriptome::geneCountsAddAlign` at STAR/source/Transcriptome_geneCountsAddAlign.cpp:4. Args: nA: uint, aAll: Transcript, gene1: vector<int32>"]
pub fn transcriptome_genecountsaddalign_l4_transcriptome_genecountsaddalign(
    transcriptome: &mut crate::transcriptome::Transcriptome,
    n_a: u32,
    a_all: &[crate::transcript::Transcript],
    gene1: &mut Vec<i32>,
) {
    let n_type = transcriptome.quants.gene_counts.n_type as usize;
    gene1.clear();
    gene1.resize(n_type, -1);

    if n_a > 1 {
        transcriptome.quants.gene_counts.c_multi += 1;
        return;
    }

    let a = &a_all[0];
    let mut ib = a.n_exons as i64 - 1;
    while ib >= 0 {
        let exon = a.exons[ib as usize];
        let g1 = exon[EX_G] as u64 + exon[EX_L] as u64 - 1;
        let mut e1 = servicefuns_l239_binarysearch1a(
            g1,
            &transcriptome.ex_g.s,
            transcriptome.ex_g.n_ex as i32,
        );

        while e1 >= 0 && transcriptome.ex_g.e_max[e1 as usize] >= exon[EX_G] as u64 {
            let ie = e1 as usize;
            if transcriptome.ex_g.e[ie] >= exon[EX_G] as u64 {
                let str1 = (transcriptome.ex_g.str_[ie] as u64).wrapping_sub(1);
                for (itype, gene) in gene1.iter_mut().enumerate().take(n_type) {
                    if itype == 1 && a.str_ != str1 && str1 < 2 {
                        continue;
                    }
                    if itype == 2 && a.str_ == str1 && str1 < 2 {
                        continue;
                    }

                    let exon_gene = transcriptome.ex_g.g[ie] as i32;
                    if *gene == -1 {
                        *gene = exon_gene;
                    } else if *gene == -2 {
                        continue;
                    } else if *gene != exon_gene {
                        *gene = -2;
                    }
                }
            }
            e1 -= 1;
        }

        ib -= 1;
    }

    for (itype, gene) in gene1.iter().copied().enumerate().take(n_type) {
        if gene == -1 {
            transcriptome.quants.gene_counts.c_none[itype] += 1;
        } else if gene == -2 {
            transcriptome.quants.gene_counts.c_ambig[itype] += 1;
        } else {
            transcriptome.quants.gene_counts.g_count[itype][gene as usize] += 1;
        }
    }
}
