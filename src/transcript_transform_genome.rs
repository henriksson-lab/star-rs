#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Transcript::transformGenome` at STAR/source/Transcript_transformGenome.cpp:4. Args: genOut: Genome, A: Transcript"]
pub fn transcript_transformgenome_l4_transcript_transformgenome(
    source: &crate::transcript::Transcript,
    gen_out: &crate::genome::Genome,
    a: &mut crate::transcript::Transcript,
) -> bool {
    let mut n_b = 0usize;
    let co_bl = &gen_out.genome_out.conv_blocks;

    for i_a in 0..source.n_exons as usize {
        let r1 = source.exons[i_a][EX_R] as u64;
        let len = source.exons[i_a][EX_L] as u64;
        let g1 = source.exons[i_a][EX_G] as u64;
        let g2 = g1 + len - 1;

        let mut c_bit = co_bl
            .partition_point(|block| block[0] <= g1)
            .saturating_sub(1);

        let b1 = co_bl[c_bit][0];
        let b2 = co_bl[c_bit][0] + co_bl[c_bit][1] - 1;
        let b1o = co_bl[c_bit][2];

        if g1 <= b2 {
            if a.exons.len() <= n_b {
                a.exons.resize(n_b + 1, [0; EX_SIZE]);
            }
            a.exons[n_b][EX_IFRAG] = source.exons[i_a][EX_IFRAG];
            a.exons[n_b][EX_R] = r1 as u32;
            a.exons[n_b][EX_G] = (b1o + g1 - b1) as u32;
            a.exons[n_b][EX_L] = if g2 <= b2 {
                len as u32
            } else {
                (b2 - g1 + 1) as u32
            };
            n_b += 1;
        }

        c_bit += 1;
        while g2 >= co_bl[c_bit][0] {
            if a.exons.len() <= n_b {
                a.exons.resize(n_b + 1, [0; EX_SIZE]);
            }
            a.exons[n_b][EX_IFRAG] = source.exons[i_a][EX_IFRAG];
            a.exons[n_b][EX_G] = co_bl[c_bit][2] as u32;
            a.exons[n_b][EX_R] = (r1 + co_bl[c_bit][0] - g1) as u32;
            a.exons[n_b][EX_L] = if g2 < co_bl[c_bit][0] + co_bl[c_bit][1] {
                (g2 - co_bl[c_bit][0] + 1) as u32
            } else {
                co_bl[c_bit][1] as u32
            };
            n_b += 1;
            c_bit += 1;
        }
    }

    if n_b == 0 {
        return false;
    }

    {
        let mut n_b1 = 1usize;
        for ib in 1..n_b {
            if n_b1 != ib {
                a.exons[n_b1] = a.exons[ib];
            }

            if a.exons[n_b1][EX_IFRAG] != a.exons[n_b1 - 1][EX_IFRAG] {
                n_b1 += 1;
                continue;
            }

            let gap_r = (a.exons[n_b1][EX_R] as u64)
                .wrapping_sub(a.exons[n_b1 - 1][EX_R] as u64)
                .wrapping_sub(a.exons[n_b1 - 1][EX_L] as u64);
            let gap_g = (a.exons[n_b1][EX_G] as u64)
                .wrapping_sub(a.exons[n_b1 - 1][EX_G] as u64)
                .wrapping_sub(a.exons[n_b1 - 1][EX_L] as u64);

            if gap_r == gap_g {
                a.exons[n_b1 - 1][EX_L] += a.exons[n_b1][EX_L] + gap_r as u32;
            } else {
                let min_gap = gap_r.min(gap_g);
                if min_gap > 0 {
                    a.exons[n_b1][EX_L] += min_gap as u32;
                    a.exons[n_b1][EX_G] -= min_gap as u32;
                    a.exons[n_b1][EX_R] -= min_gap as u32;
                }
                n_b1 += 1;
            }
        }
        n_b = n_b1;
    }

    a.n_exons = n_b as u32;
    a.exons.truncate(n_b);
    a.str_ = source.str_;
    a.chr = gen_out.chr_bin[(a.exons[0][EX_G] >> gen_out.p_ge.g_chr_bin_nbits) as usize];

    if a.canon_sj.len() < a.n_exons.saturating_sub(1) as usize {
        a.canon_sj.resize(a.n_exons.saturating_sub(1) as usize, 0);
    }
    if a.sj_annot.len() < a.n_exons.saturating_sub(1) as usize {
        a.sj_annot.resize(a.n_exons.saturating_sub(1) as usize, 0);
    }

    for ia in 0..a.n_exons.saturating_sub(1) as usize {
        a.sj_annot[ia] = 0;

        if a.exons[ia + 1][EX_IFRAG] != a.exons[ia][EX_IFRAG] {
            a.canon_sj[ia] = -3;
            continue;
        }

        let j_s = a.exons[ia][EX_G] + a.exons[ia][EX_L];
        let j_e = a.exons[ia + 1][EX_G] - 1;
        let sjdb_ind = binarysearch2_l3_binarysearch2(
            j_s,
            j_e,
            &gen_out.sjdb_start,
            &gen_out.sjdb_end,
            gen_out.sjdb_n as i32,
        );

        if sjdb_ind >= 0 {
            let sjdb_ind = sjdb_ind as usize;
            a.sj_annot[ia] = 1;
            a.canon_sj[ia] = gen_out.sjdb_motif[sjdb_ind] as i32;
            if gen_out.sjdb_motif[sjdb_ind] == 0 {
                let shift = gen_out.sjdb_shift_left[sjdb_ind] as u32;
                if a.exons[ia][EX_L] <= shift {
                    return false;
                }
                a.exons[ia][EX_L] -= shift;
                a.exons[ia + 1][EX_G] -= shift;
            }
        } else {
            let gap_g = j_e - j_s + 1;
            let gap_r = a.exons[ia + 1][EX_R] - a.exons[ia][EX_R] - a.exons[ia][EX_L];
            if gap_r > 0 {
                a.canon_sj[ia] = -2;
            } else if gap_g >= gen_out.align_intron_min {
                a.canon_sj[ia] = 0;
                let g = &gen_out.g;
                let js = j_s as usize;
                let je = j_e as usize;
                if g[js] == 2 && g[js + 1] == 3 && g[je - 1] == 0 && g[je] == 2 {
                    a.canon_sj[ia] = 1;
                } else if g[js] == 1 && g[js + 1] == 3 && g[je - 1] == 0 && g[je] == 1 {
                    a.canon_sj[ia] = 2;
                } else if g[js] == 2 && g[js + 1] == 1 && g[je - 1] == 0 && g[je] == 2 {
                    a.canon_sj[ia] = 3;
                } else if g[js] == 1 && g[js + 1] == 3 && g[je - 1] == 2 && g[je] == 1 {
                    a.canon_sj[ia] = 4;
                } else if g[js] == 0 && g[js + 1] == 3 && g[je - 1] == 0 && g[je] == 1 {
                    a.canon_sj[ia] = 5;
                } else if g[js] == 2 && g[js + 1] == 3 && g[je - 1] == 0 && g[je] == 3 {
                    a.canon_sj[ia] = 6;
                }
            } else {
                a.canon_sj[ia] = -1;
            }
        }
    }

    true
}
