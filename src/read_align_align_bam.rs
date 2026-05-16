#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::samAttrNM_MD` at STAR/source/ReadAlign_alignBAM.cpp:9. Args: trOut: Transcript, iEx1: uint, iEx2: uint, tagNM: uint, tagMD: string"]
pub fn readalign_alignbam_l9_readalign_samattrnm_md(
    read_align: &crate::read_align::ReadAlign,
    p: &crate::parameters_chimeric::Parameters,
    gen_out: &crate::genome::Genome,
    tr_out: &crate::transcript::Transcript,
    i_ex1: u64,
    i_ex2: u64,
) -> (u32, String) {
    let r = if tr_out.ro_str == 0 {
        &read_align.read1[0]
    } else {
        &read_align.read1[2]
    };
    let mut match_n = 0u64;
    let mut n_mm = 0u64;
    let mut n_i = 0u64;
    let mut n_d = 0u64;
    let mut tag_md = String::new();

    for iex in i_ex1 as usize..=i_ex2 as usize {
        for ii in 0..tr_out.exons[iex][EX_L] as usize {
            let r1 = r[ii + tr_out.exons[iex][EX_R] as usize];
            let g1 = gen_out.g[ii + tr_out.exons[iex][EX_G] as usize];
            if r1 != g1 || r1 == 4 || g1 == 4 {
                n_mm += 1;
                tag_md.push_str(&match_n.to_string());
                tag_md.push(p.genome_num_to_nt[g1 as usize] as char);
                match_n = 0;
            } else {
                match_n += 1;
            }
        }

        if iex < i_ex2 as usize {
            if tr_out.canon_sj[iex] < 0 {
                n_d += tr_out.exons[iex + 1][EX_G]
                    - (tr_out.exons[iex][EX_G] + tr_out.exons[iex][EX_L]);
            }
            n_i += tr_out.exons[iex + 1][EX_R] - tr_out.exons[iex][EX_R] - tr_out.exons[iex][EX_L];
            if tr_out.canon_sj[iex] == -1 {
                tag_md.push_str(&match_n.to_string());
                tag_md.push('^');
                let start = tr_out.exons[iex][EX_G] + tr_out.exons[iex][EX_L];
                let end = tr_out.exons[iex + 1][EX_G];
                for ii in start as usize..end as usize {
                    tag_md.push(p.genome_num_to_nt[gen_out.g[ii] as usize] as char);
                }
                match_n = 0;
            }
        }
    }

    tag_md.push_str(&match_n.to_string());
    ((n_mm + n_i + n_d) as u32, tag_md)
}

#[doc = "Original `ReadAlign::alignBAM` at STAR/source/ReadAlign_alignBAM.cpp:47. Args: trOut: Transcript, nTrOut: uint, iTrOut: uint, trChrStart: uint, mateChr: uint, mateStart: uint, mateStrand: char, alignType: int, mateMap: bool, outSAMattrOrder: vector<int>, outBAMarray: char, outBAMarrayN: uint"]
pub fn readalign_alignbam_l47_readalign_alignbam(
    read_align: &mut crate::read_align::ReadAlign,
    p: &crate::parameters_chimeric::Parameters,
    gen_out: &crate::genome::Genome,
    tr_out: &crate::transcript::Transcript,
    n_tr_out: u64,
    i_tr_out: u64,
    tr_chr_start: u64,
    mut mate_chr: u64,
    mut mate_start: u64,
    mut mate_strand: i8,
    align_type: i32,
    mate_map: Option<[bool; 2]>,
    out_sam_attr_order: &[i32],
) -> Result<crate::quantifications::AlignBamResult, String> {
    if p.out_sam_mode == "None" {
        return Ok(crate::quantifications::AlignBamResult::default());
    }

    let mut result = crate::quantifications::AlignBamResult::default();
    let mut i_ex_mate = 0_u64;
    let flag_paired = p.read_nmates == 2;
    let mut n_mates = 1_u64;
    if align_type < 0 {
        for ii in 0..tr_out.n_exons.saturating_sub(1) {
            if tr_out
                .canon_sj
                .get(ii as usize)
                .copied()
                .unwrap_or_default()
                == -3
            {
                n_mates = 2;
                i_ex_mate = ii;
                break;
            }
        }
    } else {
        n_mates = 0;
    }

    let mut t_len = 0_u64;
    let mut left_most_mate = 0_u64;
    if n_mates > 1 && p.out_sam_tlen == 2 {
        let left_end =
            tr_out.exons[i_ex_mate as usize][EX_G] + tr_out.exons[i_ex_mate as usize][EX_L];
        let right_end = tr_out.exons[tr_out.n_exons as usize - 1][EX_G]
            + tr_out.exons[tr_out.n_exons as usize - 1][EX_L];
        t_len = left_end.max(right_end)
            - tr_out.exons[0][EX_G].min(tr_out.exons[i_ex_mate as usize + 1][EX_G]);
        left_most_mate = if tr_out.exons[0][EX_G] <= tr_out.exons[i_ex_mate as usize + 1][EX_G] {
            0
        } else {
            1
        };
    }

    let left_mate = if flag_paired { tr_out.str_ } else { 0 };
    if p.out_sam_attr_present.mc {
        readalign_calccigar_l3_readalign_calccigar(
            read_align, tr_out, n_mates, i_ex_mate, left_mate,
        );
    }

    let mate_count = if align_type < 0 {
        n_mates
    } else {
        p.read_nmates as u64
    };
    let mate_map = mate_map.unwrap_or([false, false]);
    let default_read_bar = crate::solo_read_barcode::SoloReadBarcode::default();
    let read_bar = read_align
        .solo_read
        .read_bar
        .as_ref()
        .unwrap_or(&default_read_bar);

    for imate in 0..mate_count as usize {
        let mut i_ex1 = 0_u64;
        let mut i_ex2 = 0_u64;
        let mate: u64;
        let str_: u64;
        let mut packed_cigar: Vec<u32> = Vec::new();
        let mut mapq = 0_i32;
        let mut attr_out_array = vec![0_u8; 131_072];
        let mut attr_n = 0_usize;
        let mut trim_l1 = 0_u64;
        let mut trim_r1 = 0_u64;
        let mut sam_flag: u16;

        if align_type >= 0 {
            if mate_map.get(imate).copied().unwrap_or(false) {
                continue;
            }
            sam_flag = 0x4;
            if p.read_nmates == 2 {
                sam_flag |= 0x1 | if imate == 0 { 0x40 } else { 0x80 };
                if mate_map[1 - imate] {
                    if tr_out.str_ != (1 - imate) as u64 {
                        sam_flag |= 0x20;
                    }
                    mate_chr = tr_out.chr;
                    mate_start = tr_out.exons[0][EX_G] - tr_chr_start;
                    mate_strand = if tr_out.str_ == (1 - imate) as u64 {
                        0
                    } else {
                        1
                    };
                    if !tr_out.primary_flag && p.out_sam_unmapped_keep_pairs {
                        sam_flag |= 0x100;
                    }
                } else {
                    sam_flag |= 0x8;
                }
            }
            if read_align.read_filter == b'Y' as i32 {
                sam_flag |= 0x200;
            }
            mate = imate as u64;
            str_ = mate;

            attr_n += bamfunctions_l106_bamattrarraywrite(0, b"NH", &mut attr_out_array[attr_n..])
                as usize;
            attr_n += bamfunctions_l106_bamattrarraywrite(0, b"HI", &mut attr_out_array[attr_n..])
                as usize;
            attr_n += bamfunctions_l106_bamattrarraywrite(
                tr_out.max_score,
                b"AS",
                &mut attr_out_array[attr_n..],
            ) as usize;
            attr_n += bamfunctions_l106_bamattrarraywrite(
                tr_out.n_mm as i32,
                b"nM",
                &mut attr_out_array[attr_n..],
            ) as usize;
            attr_n += bamfunctions_l118_bamattrarraywrite(
                (align_type as u32).to_string().as_bytes()[0],
                b"uT",
                &mut attr_out_array[attr_n..],
            ) as usize;

            if p.out_sam_attr_present.rg && !p.out_sam_attr_rg.is_empty() {
                attr_n += bamfunctions_l124_bamattrarraywrite(
                    &p.out_sam_attr_rg[read_align.read_files_index as usize],
                    b"RG",
                    &mut attr_out_array[attr_n..],
                ) as usize;
            }
            if p.out_sam_attr_present.cr {
                attr_n += bamfunctions_l124_bamattrarraywrite(
                    &read_bar.cb_seq,
                    b"CR",
                    &mut attr_out_array[attr_n..],
                ) as usize;
            }
            if p.out_sam_attr_present.cy {
                attr_n += bamfunctions_l124_bamattrarraywrite(
                    &read_bar.cb_qual,
                    b"CY",
                    &mut attr_out_array[attr_n..],
                ) as usize;
            }
            if p.out_sam_attr_present.cb && !read_bar.cb_seq_corrected.is_empty() {
                attr_n += bamfunctions_l124_bamattrarraywrite(
                    &read_bar.cb_seq_corrected,
                    b"CB",
                    &mut attr_out_array[attr_n..],
                ) as usize;
            }
            if p.out_sam_attr_present.ur {
                attr_n += bamfunctions_l124_bamattrarraywrite(
                    &read_bar.umi_seq,
                    b"UR",
                    &mut attr_out_array[attr_n..],
                ) as usize;
            }
            if p.out_sam_attr_present.uy {
                attr_n += bamfunctions_l124_bamattrarraywrite(
                    &read_bar.umi_qual,
                    b"UY",
                    &mut attr_out_array[attr_n..],
                ) as usize;
            }
            if p.out_sam_attr_present.s_m {
                attr_n += bamfunctions_l106_bamattrarraywrite(
                    read_bar.cb_match,
                    b"sM",
                    &mut attr_out_array[attr_n..],
                ) as usize;
            }
            if p.out_sam_attr_present.s_s {
                attr_n += bamfunctions_l124_bamattrarraywrite(
                    &read_bar.b_seq,
                    b"sS",
                    &mut attr_out_array[attr_n..],
                ) as usize;
            }
            if p.out_sam_attr_present.s_q {
                attr_n += bamfunctions_l124_bamattrarraywrite(
                    &read_bar.b_qual,
                    b"sQ",
                    &mut attr_out_array[attr_n..],
                ) as usize;
            }
            if p.out_sam_attr_present.c_n {
                let clip = if read_align.clip_mates.len() > imate {
                    vec![
                        read_align.clip_mates[imate][0].clipped_n as i32,
                        read_align.clip_mates[imate][1].clipped_n as i32,
                    ]
                } else {
                    vec![0, 0]
                };
                attr_n += bamfunctions_l138_bamattrarraywrite(
                    &clip,
                    b"cN",
                    &mut attr_out_array[attr_n..],
                ) as usize;
            }
        } else {
            sam_flag = if flag_paired { 0x0001 } else { 0 };
            if flag_paired {
                if i_ex_mate == tr_out.n_exons - 1 {
                    if mate_chr > gen_out.n_chr_real as u64 {
                        sam_flag |= 0x0008;
                    }
                } else {
                    sam_flag |= 0x0002;
                }
            }
            if read_align.read_filter == b'Y' as i32 {
                sam_flag |= 0x200;
            }
            if align_type == -11 || align_type == -12 || align_type == -13 {
                sam_flag |= 0x800;
            }
            if !tr_out.primary_flag {
                sam_flag |= 0x100;
            }

            i_ex1 = if imate == 0 { 0 } else { i_ex_mate + 1 };
            i_ex2 = if imate == 0 {
                i_ex_mate
            } else {
                tr_out.n_exons - 1
            };
            mate = tr_out.exons[i_ex1 as usize][EX_IFRAG];
            str_ = tr_out.str_;

            if mate == 0 {
                sam_flag |= (str_ as u16) * 0x10;
                if n_mates == 2 {
                    sam_flag |= ((1 - str_) as u16) * 0x20;
                }
            } else {
                sam_flag |= ((1 - str_) as u16) * 0x10;
                if n_mates == 2 {
                    sam_flag |= (str_ as u16) * 0x20;
                }
            }
            if flag_paired {
                sam_flag |= if mate == 0 { 0x0040 } else { 0x0080 };
                if n_mates == 1 && mate_strand == 1 {
                    sam_flag |= 0x20;
                }
            }

            let trim_l = if str_ == 0 && mate == 0 {
                read_align.clip_mates[mate as usize][0].clipped_n
            } else if str_ == 0 && mate == 1 {
                read_align.clip_mates[mate as usize][1].clipped_n
            } else if str_ == 1 && mate == 0 {
                read_align.clip_mates[mate as usize][1].clipped_n
            } else {
                read_align.clip_mates[mate as usize][0].clipped_n
            };
            trim_l1 = trim_l as u64 + tr_out.exons[i_ex1 as usize][EX_R]
                - if tr_out.exons[i_ex1 as usize][EX_R] < read_align.read_length[left_mate as usize]
                {
                    0
                } else {
                    read_align.read_length[left_mate as usize] + 1
                };
            if trim_l1 > 0 {
                packed_cigar.push(
                    ((trim_l1 << 4) as u32)
                        | if align_type == -11 {
                            BAM_CIGAR_H
                        } else {
                            BAM_CIGAR_S
                        },
                );
            }

            let mut sj_intron = Vec::<i32>::new();
            let mut sj_motif = Vec::<u8>::new();
            for ii in i_ex1..=i_ex2 {
                let iiu = ii as usize;
                if ii > i_ex1 {
                    let prev = iiu - 1;
                    let gap_g = tr_out.exons[iiu][EX_G]
                        - (tr_out.exons[prev][EX_G] + tr_out.exons[prev][EX_L]);
                    let gap_r = tr_out.exons[iiu][EX_R]
                        - tr_out.exons[prev][EX_R]
                        - tr_out.exons[prev][EX_L];
                    if gap_r > 0 {
                        packed_cigar.push(((gap_r << 4) as u32) | BAM_CIGAR_I);
                    }
                    if tr_out.canon_sj[prev] >= 0 || tr_out.sj_annot[prev] == 1 {
                        packed_cigar.push(((gap_g << 4) as u32) | BAM_CIGAR_N);
                        sj_motif.push(
                            (tr_out.canon_sj[prev]
                                + if tr_out.sj_annot[prev] == 0 { 0 } else { 20 })
                                as u8,
                        );
                        sj_intron.push(
                            (tr_out.exons[prev][EX_G] + tr_out.exons[prev][EX_L] + 1 - tr_chr_start)
                                as i32,
                        );
                        sj_intron.push((tr_out.exons[iiu][EX_G] - tr_chr_start) as i32);
                    } else if gap_g > 0 {
                        packed_cigar.push(((gap_g << 4) as u32) | BAM_CIGAR_D);
                    }
                }
                if tr_out.exons[iiu][EX_L] > 0 {
                    packed_cigar.push(((tr_out.exons[iiu][EX_L] << 4) as u32) | BAM_CIGAR_M);
                }
            }
            if sj_motif.is_empty() {
                sj_motif.push((-1_i8) as u8);
                sj_intron.push(-1);
            }

            trim_r1 = (if tr_out.exons[i_ex1 as usize][EX_R]
                < read_align.read_length[left_mate as usize]
            {
                read_align.read_length_original[left_mate as usize]
            } else {
                read_align.read_length[left_mate as usize]
                    + 1
                    + read_align.read_length_original[mate as usize]
            }) - tr_out.exons[i_ex2 as usize][EX_R]
                - tr_out.exons[i_ex2 as usize][EX_L]
                - trim_l as u64;
            if trim_r1 > 0 {
                packed_cigar.push(
                    ((trim_r1 << 4) as u32)
                        | if align_type == -12 {
                            BAM_CIGAR_H
                        } else {
                            BAM_CIGAR_S
                        },
                );
            }

            mapq = p.out_sam_mapq_unique;
            if n_tr_out >= 5 {
                mapq = 0;
            } else if n_tr_out >= 3 {
                mapq = 1;
            } else if n_tr_out == 2 {
                mapq = 3;
            }

            let mut tag_nm: Option<u32> = None;
            let mut tag_md = String::new();
            for attr in out_sam_attr_order {
                match *attr {
                    ATTR_NH => {
                        attr_n += bamfunctions_l106_bamattrarraywrite(
                            n_tr_out as i32,
                            b"NH",
                            &mut attr_out_array[attr_n..],
                        ) as usize
                    }
                    ATTR_HI => {
                        attr_n += bamfunctions_l106_bamattrarraywrite(
                            (i_tr_out + p.out_sam_attr_ih_start as u64) as i32,
                            b"HI",
                            &mut attr_out_array[attr_n..],
                        ) as usize
                    }
                    ATTR_AS => {
                        attr_n += bamfunctions_l106_bamattrarraywrite(
                            tr_out.max_score,
                            b"AS",
                            &mut attr_out_array[attr_n..],
                        ) as usize
                    }
                    ATTR_NM_LOWER => {
                        attr_n += bamfunctions_l106_bamattrarraywrite(
                            tr_out.n_mm as i32,
                            b"nM",
                            &mut attr_out_array[attr_n..],
                        ) as usize
                    }
                    ATTR_JM => {
                        attr_n += bamfunctions_l130_bamattrarraywrite(
                            &sj_motif,
                            b"jM",
                            &mut attr_out_array[attr_n..],
                        ) as usize
                    }
                    ATTR_JI => {
                        attr_n += bamfunctions_l138_bamattrarraywrite(
                            &sj_intron,
                            b"jI",
                            &mut attr_out_array[attr_n..],
                        ) as usize
                    }
                    ATTR_XS => {
                        if tr_out.sj_motif_strand == 1 {
                            attr_n += bamfunctions_l118_bamattrarraywrite(
                                b'+',
                                b"XS",
                                &mut attr_out_array[attr_n..],
                            ) as usize;
                        } else if tr_out.sj_motif_strand == 2 {
                            attr_n += bamfunctions_l118_bamattrarraywrite(
                                b'-',
                                b"XS",
                                &mut attr_out_array[attr_n..],
                            ) as usize;
                        }
                    }
                    ATTR_NM => {
                        if tag_nm.is_none() {
                            let (nm, md) = readalign_alignbam_l9_readalign_samattrnm_md(
                                read_align, p, gen_out, tr_out, i_ex1, i_ex2,
                            );
                            tag_nm = Some(nm);
                            tag_md = md;
                        }
                        attr_n += bamfunctions_l106_bamattrarraywrite(
                            tag_nm.unwrap() as i32,
                            b"NM",
                            &mut attr_out_array[attr_n..],
                        ) as usize;
                    }
                    ATTR_MD => {
                        if tag_md.is_empty() {
                            let (nm, md) = readalign_alignbam_l9_readalign_samattrnm_md(
                                read_align, p, gen_out, tr_out, i_ex1, i_ex2,
                            );
                            tag_nm = Some(nm);
                            tag_md = md;
                        }
                        attr_n += bamfunctions_l124_bamattrarraywrite(
                            &tag_md,
                            b"MD",
                            &mut attr_out_array[attr_n..],
                        ) as usize;
                    }
                    ATTR_RG => {
                        attr_n += bamfunctions_l124_bamattrarraywrite(
                            &p.out_sam_attr_rg[read_align.read_files_index as usize],
                            b"RG",
                            &mut attr_out_array[attr_n..],
                        ) as usize
                    }
                    ATTR_RB => {
                        let mut rb = Vec::new();
                        for ii in i_ex1..=i_ex2 {
                            let ex = tr_out.exons[ii as usize];
                            rb.push((ex[EX_R] + 1) as i32);
                            rb.push((ex[EX_R] + ex[EX_L]) as i32);
                            rb.push(
                                (ex[EX_G] as u64 - gen_out.chr_start[tr_out.chr as usize] + 1)
                                    as i32,
                            );
                            rb.push(
                                (ex[EX_G] as u64 - gen_out.chr_start[tr_out.chr as usize]
                                    + ex[EX_L] as u64) as i32,
                            );
                        }
                        attr_n += bamfunctions_l138_bamattrarraywrite(
                            &rb,
                            b"rB",
                            &mut attr_out_array[attr_n..],
                        ) as usize;
                    }
                    ATTR_VG => {
                        if !tr_out.var_gen_coord.is_empty() {
                            let v: Vec<i32> =
                                tr_out.var_gen_coord.iter().map(|&x| x as i32).collect();
                            attr_n += bamfunctions_l138_bamattrarraywrite(
                                &v,
                                b"vG",
                                &mut attr_out_array[attr_n..],
                            ) as usize;
                        }
                    }
                    ATTR_VA => {
                        if !tr_out.var_allele.is_empty() {
                            attr_n += bamfunctions_l130_bamattrarraywrite(
                                &tr_out.var_allele,
                                b"vA",
                                &mut attr_out_array[attr_n..],
                            ) as usize;
                        }
                    }
                    ATTR_VW => {
                        if read_align.wasp_type != -1 {
                            attr_n += bamfunctions_l106_bamattrarraywrite(
                                read_align.wasp_type,
                                b"vW",
                                &mut attr_out_array[attr_n..],
                            ) as usize;
                        }
                    }
                    ATTR_CH => {
                        if align_type <= -10 {
                            attr_n += bamfunctions_l118_bamattrarraywrite(
                                b'1',
                                b"ch",
                                &mut attr_out_array[attr_n..],
                            ) as usize;
                        }
                    }
                    ATTR_MC => {
                        if n_mates > 1 {
                            attr_n += bamfunctions_l124_bamattrarraywrite(
                                &read_align.mates_cigar[1 - imate],
                                b"MC",
                                &mut attr_out_array[attr_n..],
                            ) as usize;
                        }
                    }
                    ATTR_CR => {
                        attr_n += bamfunctions_l124_bamattrarraywrite(
                            &read_bar.cb_seq,
                            b"CR",
                            &mut attr_out_array[attr_n..],
                        ) as usize;
                    }
                    ATTR_CY => {
                        attr_n += bamfunctions_l124_bamattrarraywrite(
                            &read_bar.cb_qual,
                            b"CY",
                            &mut attr_out_array[attr_n..],
                        ) as usize;
                    }
                    ATTR_UR => {
                        attr_n += bamfunctions_l124_bamattrarraywrite(
                            &read_bar.umi_seq,
                            b"UR",
                            &mut attr_out_array[attr_n..],
                        ) as usize;
                    }
                    ATTR_UY => {
                        attr_n += bamfunctions_l124_bamattrarraywrite(
                            &read_bar.umi_qual,
                            b"UY",
                            &mut attr_out_array[attr_n..],
                        ) as usize;
                    }
                    ATTR_SM => {
                        attr_n += bamfunctions_l106_bamattrarraywrite(
                            read_bar.cb_match,
                            b"sM",
                            &mut attr_out_array[attr_n..],
                        ) as usize;
                    }
                    ATTR_SS => {
                        attr_n += bamfunctions_l124_bamattrarraywrite(
                            &read_bar.b_seq,
                            b"sS",
                            &mut attr_out_array[attr_n..],
                        ) as usize;
                    }
                    ATTR_SQ => {
                        attr_n += bamfunctions_l124_bamattrarraywrite(
                            &read_bar.b_qual,
                            b"sQ",
                            &mut attr_out_array[attr_n..],
                        ) as usize;
                    }
                    ATTR_CB => {
                        if !read_bar.cb_seq_corrected.is_empty() {
                            attr_n += bamfunctions_l124_bamattrarraywrite(
                                &read_bar.cb_seq_corrected,
                                b"CB",
                                &mut attr_out_array[attr_n..],
                            ) as usize;
                        }
                    }
                    ATTR_UB | ATTR_GX | ATTR_GN => {}
                    _ => {
                        return Err(format!(
                            "EXITING because of FATAL BUG: unknown/unimplemented SAM/BAM atrribute (tag): {}\nSOLUTION: please make sure that SAM attributes in --outSAMattributes are compatible with the --outSAMtype option\n",
                            attr
                        ));
                    }
                }
            }
        }

        if p.read_files_type_n == 10 && !p.read_files_sam_attr_keep_none {
            let keep: Vec<u16> = p.read_files_sam_attr_keep.iter().copied().collect();
            let extra = read_align
                .read_name_extra
                .get(mate as usize)
                .map(String::as_str)
                .unwrap_or("");
            attr_n += bamfunctions_l147_bamattrarraywritesamtags(
                extra,
                &mut attr_out_array[attr_n..],
                p.read_files_sam_attr_keep_all,
                &keep,
            )? as usize;
        }

        let read0 = read_align
            .read0
            .get(mate as usize)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let qual0 = read_align
            .qual0
            .get(mate as usize)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let read_len_original = read_align
            .read_length_original
            .get(mate as usize)
            .copied()
            .unwrap_or(read0.len() as u64);
        let (mut seq_out, mut qual_out) = if mate == str_ {
            (
                read0[..read_len_original as usize].to_vec(),
                qual0[..read_len_original.min(qual0.len() as u64) as usize].to_vec(),
            )
        } else {
            let mut seq_rev = vec![0_u8; read_len_original as usize];
            sequencefuns_l16_revcomplementnucleotides(read0, &mut seq_rev, read_len_original);
            let mut qual_rev = qual0[..read_len_original.min(qual0.len() as u64) as usize].to_vec();
            qual_rev.reverse();
            (seq_rev, qual_rev)
        };

        let mut seq_mate_length = read_len_original;
        if align_type == -11 {
            seq_mate_length -= trim_l1;
            seq_out = seq_out[trim_l1 as usize..].to_vec();
            qual_out = qual_out[trim_l1 as usize..].to_vec();
        } else if align_type == -12 {
            seq_mate_length -= trim_r1;
            seq_out.truncate(seq_mate_length as usize);
            qual_out.truncate(seq_mate_length as usize);
        }

        let mut packed_seq = vec![0_u8; (seq_mate_length as usize + 1) / 2];
        sequencefuns_l122_nuclpackbam(&seq_out, &mut packed_seq, seq_mate_length);

        let mut core = [0_u32; 9];
        if align_type < 0 {
            core[1] = tr_out.chr as u32;
            core[2] = (tr_out.exons[i_ex1 as usize][EX_G] - tr_chr_start) as u32;
        } else {
            core[1] = u32::MAX;
            core[2] = u32::MAX;
        }
        let read_name_len = read_align.read_name.len() as u32;
        if align_type < 0 {
            let beg = (tr_out.exons[i_ex1 as usize][EX_G] - tr_chr_start) as u32;
            let end = (tr_out.exons[i_ex2 as usize][EX_G] + tr_out.exons[i_ex2 as usize][EX_L]
                - tr_chr_start) as u32;
            core[3] = ((bamfunctions_l95_reg2bin(beg as i32, end as i32) as u32) << 16)
                | ((mapq as u32) << 8)
                | read_name_len;
        } else {
            core[3] = ((bamfunctions_l95_reg2bin(-1, 0) as u32) << 16) | read_name_len;
        }
        let flag_out = (sam_flag & p.out_sam_flag_and) | p.out_sam_flag_or;
        core[4] = ((flag_out as u32) << 16) | packed_cigar.len() as u32;
        core[5] = seq_mate_length as u32;
        if n_mates > 1 {
            core[6] = tr_out.chr as u32;
            core[7] = (tr_out.exons[if imate == 0 {
                i_ex_mate as usize + 1
            } else {
                0
            }][EX_G]
                - tr_chr_start) as u32;
        } else if mate_chr < gen_out.n_chr_real as u64 {
            core[6] = mate_chr as u32;
            core[7] = mate_start as u32;
        } else {
            core[6] = u32::MAX;
            core[7] = u32::MAX;
        }
        if n_mates > 1 {
            if p.out_sam_tlen == 1 {
                let tl = (tr_out.exons[tr_out.n_exons as usize - 1][EX_G]
                    + tr_out.exons[tr_out.n_exons as usize - 1][EX_L]
                    - tr_out.exons[0][EX_G]) as i32;
                core[8] = if imate == 0 { tl } else { -tl } as u32;
            } else if p.out_sam_tlen == 2 {
                let tl = t_len as i32;
                core[8] = if imate as u64 == left_most_mate {
                    tl
                } else {
                    -tl
                } as u32;
            }
        } else {
            core[8] = 0;
        }

        let mut rec = vec![0_u8; 4];
        for value in core.iter().skip(1) {
            rec.extend_from_slice(&value.to_ne_bytes());
        }
        let qname = read_align.read_name.as_bytes();
        if qname.len() > 1 {
            rec.extend_from_slice(&qname[1..]);
        }
        rec.push(0);
        for cigar in &packed_cigar {
            rec.extend_from_slice(&cigar.to_ne_bytes());
        }
        rec.extend_from_slice(&packed_seq);
        if read_align.read_file_type == 2 && p.out_sam_mode != "NoQS" {
            for ii in 0..seq_mate_length as usize {
                rec.push(qual_out.get(ii).copied().unwrap_or(33).saturating_sub(33));
            }
        } else {
            rec.extend(std::iter::repeat_n(0xFF, seq_mate_length as usize));
        }
        rec.extend_from_slice(&attr_out_array[..attr_n]);
        let block_size = (rec.len() - 4) as u32;
        rec[0..4].copy_from_slice(&block_size.to_ne_bytes());

        result.record_sizes[imate] = rec.len() as u32;
        result.records[imate] = rec;
        result.sam_flags[imate] = flag_out;
        result.mapq[imate] = mapq;
        result.n_cigar[imate] = packed_cigar.len() as u32;
        if align_type < 0 {
            result.signal_records[imate] = Some(crate::parameters_chimeric::SignalFromBamRecord {
                tid: core[1] as i32,
                pos: core[2],
                flag: flag_out,
                cigar: packed_cigar,
                nh: Some(n_tr_out as u32),
            });
        }
    }

    result.n_lines = if result.record_sizes[1] == 0 { 1 } else { 2 };
    Ok(result)
}
