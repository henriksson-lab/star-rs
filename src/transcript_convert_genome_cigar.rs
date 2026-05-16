#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Transcript::convertGenomeCigar` at STAR/source/Transcript_convertGenomeCigar.cpp:2. Args: genOut: Genome, A: Transcript"]
pub fn transcript_convertgenomecigar_l2_transcript_convertgenomecigar(
    source: &crate::transcript::Transcript,
    gen_out: &crate::genome::Genome,
    a: &mut crate::transcript::Transcript,
) -> bool {
    const BAM_CIGAR_M: u32 = 0;
    const BAM_CIGAR_I: u32 = 1;
    const BAM_CIGAR_D: u32 = 2;
    const BAM_CIGAR_N: u32 = 3;
    const BAM_CIGAR_S: u32 = 4;

    let co_bl = &gen_out.genome_out.conv_blocks;
    let mut ic_b = co_bl
        .partition_point(|block| block[0] <= source.g_start)
        .saturating_sub(1);

    a.g_start = source.g_start - co_bl[ic_b][0] + co_bl[ic_b][2];
    a.cigar.clear();
    a.cigar.reserve(source.cigar.len() + 100);

    let mut g_end = source.g_start;
    for &cc in &source.cigar {
        if cc[0] == BAM_CIGAR_M || cc[0] == BAM_CIGAR_D || cc[0] == BAM_CIGAR_N {
            let mut g_start1 = g_end;
            g_end += cc[1] as u64;
            while co_bl[ic_b][0] + co_bl[ic_b][1] <= g_end {
                let b_e = co_bl[ic_b][0] + co_bl[ic_b][1];
                let len1 = (b_e - g_start1) as u32;
                a.cigar.push([cc[0], len1]);

                let len1 = (co_bl[ic_b + 1][2] - co_bl[ic_b][2] - co_bl[ic_b][1]) as u32;
                a.cigar.push([BAM_CIGAR_N, len1]);

                g_start1 = b_e;
                ic_b += 1;
            }

            if g_end > g_start1 {
                a.cigar.push([cc[0], (g_end - g_start1) as u32]);
            }

            if g_end < co_bl[ic_b][0] {
                ic_b -= 1;
            }
        } else {
            a.cigar.push(cc);
        }
    }

    if !a.cigar.is_empty() {
        let mut c0 = 0usize;
        for c1 in 1..a.cigar.len() {
            if a.cigar[c1][0] == a.cigar[c0][0] {
                a.cigar[c0][1] += a.cigar[c1][1];
            } else if a.cigar[c1][0] == BAM_CIGAR_N
                && a.cigar[c0][0] == BAM_CIGAR_I
                && c0 > 0
                && a.cigar[c0 - 1][0] == BAM_CIGAR_N
            {
                a.cigar[c0 - 1][1] += a.cigar[c1][1];
            } else {
                c0 += 1;
                a.cigar[c0] = a.cigar[c1];
            }
        }
        a.cigar.truncate(c0 + 1);
    }

    if a.cigar.last().is_some_and(|cc| cc[0] == BAM_CIGAR_N) {
        a.cigar.pop();
    }
    if a.cigar.len() >= 2
        && a.cigar.last().is_some_and(|cc| cc[0] == BAM_CIGAR_S)
        && a.cigar[a.cigar.len() - 2][0] == BAM_CIGAR_N
    {
        let i = a.cigar.len() - 2;
        a.cigar.remove(i);
    }

    {
        let mut n_op = [0u32; 5];
        let mut c1 = 0usize;
        while c1 < a.cigar.len() {
            n_op[a.cigar[c1][0] as usize] += a.cigar[c1][1];
            if a.cigar[c1][0] == BAM_CIGAR_M {
                if n_op[BAM_CIGAR_M as usize] >= gen_out.align_sjdb_overhang_min {
                    break;
                }
            } else if a.cigar[c1][0] == BAM_CIGAR_N {
                c1 += 1;
                while c1 < a.cigar.len() && a.cigar[c1][0] != BAM_CIGAR_M {
                    n_op[a.cigar[c1][0] as usize] += a.cigar[c1][1];
                    c1 += 1;
                }
                if c1 < a.cigar.len() {
                    if c1 > 1 {
                        a.cigar.drain(0..c1 - 1);
                    }
                    a.cigar[0] = [
                        BAM_CIGAR_S,
                        n_op[BAM_CIGAR_M as usize]
                            + n_op[BAM_CIGAR_S as usize]
                            + n_op[BAM_CIGAR_I as usize],
                    ];
                    a.g_start += (n_op[BAM_CIGAR_M as usize]
                        + n_op[BAM_CIGAR_D as usize]
                        + n_op[BAM_CIGAR_N as usize]) as u64;
                }
                break;
            }
            c1 += 1;
        }
    }

    {
        let mut n_op = [0u32; 5];
        let mut c1 = a.cigar.len().saturating_sub(1);
        while c1 > 0 {
            n_op[a.cigar[c1][0] as usize] += a.cigar[c1][1];
            if a.cigar[c1][0] == BAM_CIGAR_M {
                if n_op[BAM_CIGAR_M as usize] >= gen_out.align_sjdb_overhang_min {
                    break;
                }
            } else if a.cigar[c1][0] == BAM_CIGAR_N {
                c1 -= 1;
                while a.cigar[c1][0] != BAM_CIGAR_M {
                    n_op[a.cigar[c1][0] as usize] += a.cigar[c1][1];
                    c1 -= 1;
                }
                a.cigar.truncate(c1 + 1);
                a.cigar.push([
                    BAM_CIGAR_S,
                    n_op[BAM_CIGAR_M as usize]
                        + n_op[BAM_CIGAR_S as usize]
                        + n_op[BAM_CIGAR_I as usize],
                ]);
                break;
            }
            c1 -= 1;
        }
    }

    a.g_length = 0;
    for cc in &a.cigar {
        if cc[0] == BAM_CIGAR_M || cc[0] == BAM_CIGAR_D || cc[0] == BAM_CIGAR_N {
            a.g_length += cc[1] as u64;
        }
    }

    a.str_ = source.str_;
    if a.g_start >= gen_out.genome_out.n_minus_strand_offset {
        a.str_ = 1 - a.str_;
        a.g_start =
            2 * gen_out.genome_out.n_minus_strand_offset - (a.g_start + a.g_length);
        a.cigar.reverse();
    }

    a.chr = gen_out.chr_bin[(a.g_start >> gen_out.p_ge.g_chr_bin_nbits) as usize] as u64;
    true
}
