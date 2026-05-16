#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::quantTranscriptome` at STAR/source/ReadAlign_quantTranscriptome.cpp:7. Args: Tr: Transcriptome, nAlignG: uint, alignG: Transcript, alignT: Transcript"]
pub fn readalign_quanttranscriptome_l7_readalign_quanttranscriptome(
    read_align: &crate::read_align::ReadAlign,
    p: &crate::parameters_chimeric::Parameters,
    map_gen: &crate::genome::Genome,
    transcriptome: &crate::transcriptome::Transcriptome,
    align_g: &[crate::transcript::Transcript],
    primary_pick_fraction: f64,
) -> crate::quantifications::QuantTranscriptomeResult {
    let mut align_t = Vec::<crate::transcript::Transcript>::new();
    let mut n_align_t = 0_u32;

    for align1_in in align_g {
        if !p.quant_tr_sam_indel && (align1_in.n_del > 0 || align1_in.n_ins > 0) {
            continue;
        }
        if !p.quant_tr_sam_single_end
            && p.read_nmates == 2
            && align1_in.exons[0][EX_IFRAG]
                == align1_in.exons[align1_in.n_exons as usize - 1][EX_IFRAG]
        {
            continue;
        }

        let mut align2_storage;
        let align1 = if !p.quant_tr_sam_soft_clip {
            let mut n_mm1 = 0_u32;
            let r = if align1_in.ro_str == 0 {
                &read_align.read1[0]
            } else {
                &read_align.read1[2]
            };
            align2_storage = align1_in.clone();

            for iab in 0..align2_storage.n_exons as usize {
                let mut left1 = 0_u32;
                let mut right1 = 0_u32;
                if iab == 0 {
                    left1 = align2_storage.exons[iab][EX_R];
                } else if align2_storage.canon_sj[iab - 1] == -3 {
                    left1 = align2_storage.exons[iab][EX_R]
                        - read_align.read_length[align2_storage.exons[iab - 1][EX_IFRAG] as usize]
                        - 1;
                }

                if iab == align2_storage.n_exons as usize - 1 {
                    right1 = read_align.l_read
                        - align2_storage.exons[iab][EX_R]
                        - align2_storage.exons[iab][EX_L];
                } else if align2_storage.canon_sj[iab] == -3 {
                    right1 = read_align.read_length[align2_storage.exons[iab][EX_IFRAG] as usize]
                        - align2_storage.exons[iab][EX_R]
                        - align2_storage.exons[iab][EX_L];
                }

                for b in 1..=left1 {
                    let r1 = r[(align2_storage.exons[iab][EX_R] - b) as usize];
                    let g1 = map_gen.g[(align2_storage.exons[iab][EX_G] - b) as usize];
                    if r1 != g1 && r1 < 4 && g1 < 4 {
                        n_mm1 += 1;
                    }
                }
                for b in 0..right1 {
                    let r1 = r[(align2_storage.exons[iab][EX_R]
                        + align2_storage.exons[iab][EX_L]
                        + b) as usize];
                    let g1 = map_gen.g[(align2_storage.exons[iab][EX_G]
                        + align2_storage.exons[iab][EX_L]
                        + b) as usize];
                    if r1 != g1 && r1 < 4 && g1 < 4 {
                        n_mm1 += 1;
                    }
                }
                align2_storage.exons[iab][EX_R] -= left1;
                align2_storage.exons[iab][EX_G] -= left1;
                align2_storage.exons[iab][EX_L] += left1 + right1;
            }

            if align2_storage.n_mm + n_mm1
                > read_align
                    .out_filter_mismatch_nmax_total
                    .min((p.out_filter_mismatch_nover_lmax * (read_align.l_read - 1) as f64) as u32)
            {
                continue;
            }
            &mut align2_storage
        } else {
            align2_storage = align1_in.clone();
            &mut align2_storage
        };

        let added = transcriptome_quantalign_l91_transcriptome_quantalign(
            transcriptome,
            align1,
            &mut align_t,
        );
        n_align_t += added;
    }

    let mut bam_requests = Vec::new();
    if p.quant_tr_sam_bam_yes && n_align_t > 0 {
        let n_align_t_requested = n_align_t as usize;
        let n_align_t_available = align_t.len().min(n_align_t_requested);
        let mut primary_index = (primary_pick_fraction * n_align_t_available as f64) as usize;
        if primary_index >= n_align_t_available {
            primary_index = n_align_t_available - 1;
        }
        align_t[primary_index].primary_flag = true;

        for iatr in 0..n_align_t_available {
            bam_requests.push(crate::quantifications::QuantTranscriptomeBamRequest {
                transcript: align_t[iatr].clone(),
                n_align_t,
                i_align_t: iatr as u32,
            });
        }
    }

    crate::quantifications::QuantTranscriptomeResult {
        n_align_t,
        align_t,
        bam_requests,
    }
}
