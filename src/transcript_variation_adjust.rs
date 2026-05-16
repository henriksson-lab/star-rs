#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Transcript::variationAdjust` at STAR/source/Transcript_variationAdjust.cpp:4. Args: mapGen: Genome, R: char"]
pub fn transcript_variationadjust_l4_transcript_variationadjust(
    transcript: &mut crate::transcript::Transcript,
    map_gen: &crate::genome::Genome,
    r: &[u8],
) -> i32 {
    let var = &map_gen.var;
    if !var.yes {
        return 0;
    }

    let d_score = 0;

    for ie in 0..transcript.n_exons as usize {
        let exon = transcript.exons[ie];
        let mut isnp = servicefuns_l266_binarysearch1b(exon[EX_G], &var.snp.loci, var.snp.n as i32);
        if isnp >= 0 {
            while (isnp as u32) < var.snp.n && exon[EX_G] + exon[EX_L] > var.snp.loci[isnp as usize]
            {
                transcript.var_ind.push(isnp);
                transcript.var_gen_coord.push(
                    var.snp.loci[isnp as usize] - map_gen.chr_start[transcript.chr as usize] as u32,
                );

                let read_coord = exon[EX_R] + var.snp.loci[isnp as usize] - exon[EX_G];
                transcript.var_read_coord.push(read_coord);
                let nt_r = r[read_coord as usize];

                let mut igt = if nt_r > 3 { 4 } else { 1 };
                if nt_r <= 3 {
                    while igt < 3 {
                        if var.snp.nt[isnp as usize][igt as usize] == nt_r {
                            break;
                        }
                        igt += 1;
                    }
                }

                transcript.var_allele.push(igt);
                isnp += 1;
            }
        }
    }

    d_score
}
