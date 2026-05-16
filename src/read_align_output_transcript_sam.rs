#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::outputTranscriptSAM` at STAR/source/ReadAlign_outputTranscriptSAM.cpp:5. Args: trOut: Transcript, nTrOut: uint, iTrOut: uint, mateChr: uint, mateStart: uint, mateStrand: char, unmapType: int, mateMap: bool, outStream: ostream"]
pub fn readalign_outputtranscriptsam_l5_readalign_outputtranscriptsam(
    tr_out: &crate::transcript::Transcript,
    n_tr_out: u32,
    i_tr_out: u32,
    mate_chr: u32,
    mate_start: u32,
    mate_strand: i8,
    unmap_type: i32,
    mate_map: Option<&[bool]>,
    out_stream: &mut String,
    read_filter: u8,
    read_name: &str,
    read0: &[String],
    read_file_type: i32,
    qual0: &[String],
    p: &crate::parameters_chimeric::Parameters,
    read_files_index: u32,
    read_name_extra: &[String],
    l_read: u32,
    read_length: &[u32],
    read_length_original: &[u32],
    clip_mates: &[Vec<crate::clip_mate::ClipMate>],
    read1: &[Vec<u8>; 3],
    gen_out: &crate::genome::Genome,
) -> Result<u64, String> {
    if p.out_sam_mode == "None" {
        return Ok(0);
    }

    let out_pos0 = out_stream.len();
    let read_name_out = read_name.get(1..).unwrap_or("");

    if unmap_type >= 0 {
        let mate_map = mate_map.unwrap_or(&[]);
        for imate in 0..p.read_nmates as usize {
            if !mate_map.get(imate).copied().unwrap_or(false) {
                let mut sam_flag = 0x4_u16;
                if p.read_nmates == 2 {
                    sam_flag |= 0x1 + if imate == 0 { 0x40 } else { 0x80 };
                    if mate_map.get(1 - imate).copied().unwrap_or(false) {
                        if tr_out.str_ != (1 - imate) as u32 {
                            sam_flag |= 0x20;
                        }
                    } else {
                        sam_flag |= 0x8;
                    }
                }

                if read_filter == b'Y' {
                    sam_flag |= 0x200;
                }

                if mate_map
                    .get(1usize.wrapping_sub(imate))
                    .copied()
                    .unwrap_or(false)
                    && !tr_out.primary_flag
                    && p.out_sam_unmapped_keep_pairs
                {
                    sam_flag |= 0x100;
                }

                out_stream.push_str(&format!("{}\t{}\t*\t0\t0\t*", read_name_out, sam_flag));

                if mate_map
                    .get(1usize.wrapping_sub(imate))
                    .copied()
                    .unwrap_or(false)
                {
                    out_stream.push_str(&format!(
                        "\t{}\t{}",
                        gen_out.chr_name[tr_out.chr as usize],
                        tr_out.exons[0][EX_G] as u64 + 1 - gen_out.chr_start[tr_out.chr as usize]
                    ));
                } else {
                    out_stream.push_str("\t*\t0");
                }

                out_stream.push_str(&format!(
                    "\t0\t{}\t{}\tNH:i:0\tHI:i:0\tAS:i:{}\tnM:i:{}\tuT:A:{}",
                    read0[imate],
                    if read_file_type == 2 {
                        qual0[imate].as_str()
                    } else {
                        "*"
                    },
                    tr_out.max_score,
                    tr_out.n_mm,
                    unmap_type
                ));
                if !p.out_sam_attr_rg.is_empty() {
                    out_stream.push_str(&format!(
                        "\tRG:Z:{}",
                        p.out_sam_attr_rg[read_files_index as usize]
                    ));
                }
                if p.read_files_type_n == 10 && !read_name_extra[imate].is_empty() {
                    out_stream.push('\t');
                    out_stream.push_str(&read_name_extra[imate]);
                }
                out_stream.push('\n');
            }
        }
        return Ok((out_stream.len() - out_pos0) as u64);
    }

    let flag_paired = p.read_nmates == 2;
    let mut i_ex_mate = 0_u32;
    let mut n_mates = 1_u32;
    while i_ex_mate < tr_out.n_exons - 1 {
        if tr_out.canon_sj[i_ex_mate as usize] == -3 {
            n_mates = 2;
            break;
        }
        i_ex_mate += 1;
    }

    let mut sam_flag_common = 0_u16;
    if flag_paired {
        sam_flag_common = 0x0001;
        if i_ex_mate == tr_out.n_exons - 1 {
            if mate_chr > gen_out.n_chr_real {
                sam_flag_common += 0x0008;
            }
        } else if p.align_ends_protrude_concordant_pair
            || (tr_out.exons[0][EX_G]
                <= tr_out.exons[i_ex_mate as usize + 1][EX_G] + tr_out.exons[0][EX_R]
                && tr_out.exons[i_ex_mate as usize][EX_G] + tr_out.exons[i_ex_mate as usize][EX_L]
                    <= tr_out.exons[tr_out.n_exons as usize - 1][EX_G] + l_read
                        - tr_out.exons[tr_out.n_exons as usize - 1][EX_R])
        {
            sam_flag_common += 0x0002;
        }
    }
    if read_filter == b'Y' {
        sam_flag_common += 0x200;
    }

    let str_ = tr_out.str_;
    let left_mate = if flag_paired { str_ } else { 0 };
    let mut mate_cigars = Vec::new();
    if p.out_sam_attr_present.mc {
        for imate in 0..n_mates {
            let i_ex1 = if imate == 0 { 0 } else { i_ex_mate + 1 };
            let i_ex2 = if imate == 0 {
                i_ex_mate
            } else {
                tr_out.n_exons - 1
            };
            let mate = tr_out.exons[i_ex1 as usize][EX_IFRAG] as usize;
            let trim_l = if str_ == 0 && mate == 0 {
                clip_mates[mate][0].clipped_n
            } else if str_ == 0 && mate == 1 {
                clip_mates[mate][1].clipped_n
            } else if str_ == 1 && mate == 0 {
                clip_mates[mate][1].clipped_n
            } else {
                clip_mates[mate][0].clipped_n
            };
            let mut cigar = String::new();
            let trim_l1 = trim_l + tr_out.exons[i_ex1 as usize][EX_R]
                - if tr_out.exons[i_ex1 as usize][EX_R] < read_length[left_mate as usize] {
                    0
                } else {
                    read_length[left_mate as usize] + 1
                };
            if trim_l1 > 0 {
                cigar.push_str(&format!("{}S", trim_l1));
            }
            for ii in i_ex1..=i_ex2 {
                let ii = ii as usize;
                if ii > i_ex1 as usize {
                    let gap_g = tr_out.exons[ii][EX_G]
                        - (tr_out.exons[ii - 1][EX_G] + tr_out.exons[ii - 1][EX_L]);
                    let gap_r = tr_out.exons[ii][EX_R]
                        - tr_out.exons[ii - 1][EX_R]
                        - tr_out.exons[ii - 1][EX_L];
                    if gap_r > 0 {
                        cigar.push_str(&format!("{}I", gap_r));
                    }
                    if tr_out.canon_sj[ii - 1] >= 0 || tr_out.sj_annot[ii - 1] == 1 {
                        cigar.push_str(&format!("{}N", gap_g));
                    } else if gap_g > 0 {
                        cigar.push_str(&format!("{}D", gap_g));
                    }
                }
                cigar.push_str(&format!("{}M", tr_out.exons[ii][EX_L]));
            }
            let trim_r1 = (if tr_out.exons[i_ex1 as usize][EX_R] < read_length[left_mate as usize] {
                read_length_original[left_mate as usize]
            } else {
                read_length[left_mate as usize] + 1 + read_length_original[mate]
            }) - tr_out.exons[i_ex2 as usize][EX_R]
                - tr_out.exons[i_ex2 as usize][EX_L]
                - trim_l;
            if trim_r1 > 0 {
                cigar.push_str(&format!("{}S", trim_r1));
            }
            mate_cigars.push(cigar);
        }
    }

    for imate in 0..n_mates {
        let mut sam_flag = sam_flag_common;
        let i_ex1 = if imate == 0 { 0 } else { i_ex_mate + 1 };
        let i_ex2 = if imate == 0 {
            i_ex_mate
        } else {
            tr_out.n_exons - 1
        };
        let mate = tr_out.exons[i_ex1 as usize][EX_IFRAG];

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

        if !tr_out.primary_flag {
            sam_flag |= 0x100;
        }

        let trim_l = if str_ == 0 && mate == 0 {
            clip_mates[mate as usize][0].clipped_n
        } else if str_ == 0 && mate == 1 {
            clip_mates[mate as usize][1].clipped_n
        } else if str_ == 1 && mate == 0 {
            clip_mates[mate as usize][1].clipped_n
        } else {
            clip_mates[mate as usize][0].clipped_n
        };

        let mut cigar = String::new();
        let trim_l1 = trim_l + tr_out.exons[i_ex1 as usize][EX_R]
            - if tr_out.exons[i_ex1 as usize][EX_R] < read_length[left_mate as usize] {
                0
            } else {
                read_length[left_mate as usize] + 1
            };
        if trim_l1 > 0 {
            cigar.push_str(&format!("{}S", trim_l1));
        }

        let mut sj_motif = String::new();
        let mut sj_intron = String::new();
        for ii in i_ex1..=i_ex2 {
            let ii = ii as usize;
            if ii > i_ex1 as usize {
                let gap_g = tr_out.exons[ii][EX_G]
                    - (tr_out.exons[ii - 1][EX_G] + tr_out.exons[ii - 1][EX_L]);
                let gap_r = tr_out.exons[ii][EX_R]
                    - tr_out.exons[ii - 1][EX_R]
                    - tr_out.exons[ii - 1][EX_L];
                if gap_r > 0 {
                    cigar.push_str(&format!("{}I", gap_r));
                }
                if tr_out.canon_sj[ii - 1] >= 0 || tr_out.sj_annot[ii - 1] == 1 {
                    cigar.push_str(&format!("{}N", gap_g));
                    sj_motif.push_str(&format!(
                        ",{}",
                        tr_out.canon_sj[ii - 1]
                            + if tr_out.sj_annot[ii - 1] == 0 {
                                0
                            } else {
                                SJ_SAM_ANNOTATED_MOTIF_SHIFT
                            }
                    ));
                    sj_intron.push_str(&format!(
                        ",{},{}",
                        tr_out.exons[ii - 1][EX_G] as u64 + tr_out.exons[ii - 1][EX_L] as u64 + 1
                            - gen_out.chr_start[tr_out.chr as usize],
                        tr_out.exons[ii][EX_G] as u64 - gen_out.chr_start[tr_out.chr as usize]
                    ));
                } else if gap_g > 0 {
                    cigar.push_str(&format!("{}D", gap_g));
                }
            }
            cigar.push_str(&format!("{}M", tr_out.exons[ii][EX_L]));
        }
        if sj_motif.is_empty() {
            sj_motif = ",-1".to_string();
            sj_intron = ",-1".to_string();
        }

        let trim_r1 = (if tr_out.exons[i_ex1 as usize][EX_R] < read_length[left_mate as usize] {
            read_length_original[left_mate as usize]
        } else {
            read_length[left_mate as usize] + 1 + read_length_original[mate as usize]
        }) - tr_out.exons[i_ex2 as usize][EX_R]
            - tr_out.exons[i_ex2 as usize][EX_L]
            - trim_l;
        if trim_r1 > 0 {
            cigar.push_str(&format!("{}S", trim_r1));
        }

        let (seq_out, qual_out) = if mate == str_ {
            (
                read0[mate as usize].clone(),
                qual0.get(mate as usize).cloned().unwrap_or_default(),
            )
        } else {
            let mut seq_mate = vec![0u8; read_length_original[mate as usize] as usize];
            sequencefuns_l16_revcomplementnucleotides(
                read0[mate as usize].as_bytes(),
                &mut seq_mate,
                read_length_original[mate as usize],
            );
            let qual_mate: String = qual0[mate as usize]
                .chars()
                .rev()
                .take(read_length_original[mate as usize] as usize)
                .collect();
            (String::from_utf8(seq_mate).unwrap(), qual_mate)
        };

        let mapq = if n_tr_out >= 5 {
            0
        } else if n_tr_out >= 3 {
            1
        } else if n_tr_out == 2 {
            3
        } else {
            p.out_sam_mapq_unique
        };

        out_stream.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            read_name_out,
            (sam_flag & p.out_sam_flag_and) | p.out_sam_flag_or,
            gen_out.chr_name[tr_out.chr as usize],
            tr_out.exons[i_ex1 as usize][EX_G] as u64 + 1 - gen_out.chr_start[tr_out.chr as usize],
            mapq,
            cigar
        ));

        if n_mates > 1 {
            out_stream.push_str(&format!(
                "\t=\t{}\t{}{}",
                tr_out.exons[if imate == 0 {
                    i_ex_mate as usize + 1
                } else {
                    0
                }][EX_G] as u64
                    + 1
                    - gen_out.chr_start[tr_out.chr as usize],
                if imate == 0 { "" } else { "-" },
                tr_out.exons[tr_out.n_exons as usize - 1][EX_G]
                    + tr_out.exons[tr_out.n_exons as usize - 1][EX_L]
                    - tr_out.exons[0][EX_G]
            ));
        } else if mate_chr < gen_out.n_chr_real {
            out_stream.push_str(&format!(
                "\t{}\t{}\t0",
                gen_out.chr_name[mate_chr as usize],
                mate_start as u64 + 1 - gen_out.chr_start[mate_chr as usize]
            ));
        } else {
            out_stream.push_str("\t*\t0\t0");
        }

        out_stream.push('\t');
        out_stream.push_str(&seq_out);
        out_stream.push('\t');
        if read_file_type == 2 && p.out_sam_mode != "NoQS" {
            out_stream.push_str(&qual_out);
        } else {
            out_stream.push('*');
        }

        let mut tag_nm = 0_u32;
        let mut tag_md = String::new();
        if p.out_sam_attr_present.nm || p.out_sam_attr_present.md {
            let r = if tr_out.ro_str == 0 {
                &read1[0]
            } else {
                &read1[2]
            };
            let mut match_n = 0_u32;
            for iex in i_ex1..=i_ex2 {
                let iex = iex as usize;
                for ii in 0..tr_out.exons[iex][EX_L] {
                    let r1 = r[(ii + tr_out.exons[iex][EX_R]) as usize];
                    let g1 = gen_out.g[(ii + tr_out.exons[iex][EX_G]) as usize];
                    if r1 != g1 || r1 == 4 || g1 == 4 {
                        tag_nm += 1;
                        tag_md.push_str(&match_n.to_string());
                        tag_md.push(p.genome_num_to_nt[g1 as usize] as char);
                        match_n = 0;
                    } else {
                        match_n += 1;
                    }
                }
                if iex < i_ex2 as usize {
                    if tr_out.canon_sj[iex] == -1 {
                        let gap = tr_out.exons[iex + 1][EX_G]
                            - (tr_out.exons[iex][EX_G] + tr_out.exons[iex][EX_L]);
                        tag_nm += gap;
                        tag_md.push_str(&format!("{}^", match_n));
                        for ii in tr_out.exons[iex][EX_G] + tr_out.exons[iex][EX_L]
                            ..tr_out.exons[iex + 1][EX_G]
                        {
                            tag_md
                                .push(p.genome_num_to_nt[gen_out.g[ii as usize] as usize] as char);
                        }
                        match_n = 0;
                    } else if tr_out.canon_sj[iex] == -2 {
                        tag_nm += tr_out.exons[iex + 1][EX_R]
                            - tr_out.exons[iex][EX_R]
                            - tr_out.exons[iex][EX_L];
                    }
                }
            }
            tag_md.push_str(&match_n.to_string());
        }

        for attr in &p.out_sam_attr_order {
            match *attr {
                ATTR_NH => out_stream.push_str(&format!("\tNH:i:{}", n_tr_out)),
                ATTR_HI => {
                    out_stream.push_str(&format!("\tHI:i:{}", i_tr_out + p.out_sam_attr_ih_start))
                }
                ATTR_AS => out_stream.push_str(&format!("\tAS:i:{}", tr_out.max_score)),
                ATTR_NM_LOWER => out_stream.push_str(&format!("\tnM:i:{}", tr_out.n_mm)),
                ATTR_JM => out_stream.push_str(&format!("\tjM:B:c{}", sj_motif)),
                ATTR_JI => out_stream.push_str(&format!("\tjI:B:i{}", sj_intron)),
                ATTR_XS => {
                    if tr_out.sj_motif_strand == 1 {
                        out_stream.push_str("\tXS:A:+");
                    } else if tr_out.sj_motif_strand == 2 {
                        out_stream.push_str("\tXS:A:-");
                    }
                }
                ATTR_NM => out_stream.push_str(&format!("\tNM:i:{}", tag_nm)),
                ATTR_MD => out_stream.push_str(&format!("\tMD:Z:{}", tag_md)),
                ATTR_RG => out_stream.push_str(&format!(
                    "\tRG:Z:{}",
                    p.out_sam_attr_rg[read_files_index as usize]
                )),
                ATTR_MC => {
                    if n_mates > 1 {
                        out_stream.push_str(&format!("\tMC:Z:{}", mate_cigars[1 - imate as usize]));
                    }
                }
                ATTR_HA => {
                    if gen_out.p_ge.transform.type_ == 2 {
                        out_stream.push_str(&format!("\tha:i:{}", tr_out.haplo_type));
                    }
                }
                ATTR_CH | ATTR_CR | ATTR_CY | ATTR_UR | ATTR_UY | ATTR_CB | ATTR_UB | ATTR_SM
                | ATTR_SS | ATTR_SQ | ATTR_RB | ATTR_VG | ATTR_VA | ATTR_VW | ATTR_GX | ATTR_GN => {
                }
                _ => {
                    return Err(format!(
                        "EXITING because of FATAL error: unknown/unimplemented SAM atrribute (tag): {}\nSOLUTION: contact Alex Dobin at dobin@cshl.edu\n",
                        attr
                    ));
                }
            }
        }

        if p.read_files_type_n == 10 && !read_name_extra[imate as usize].is_empty() {
            out_stream.push('\t');
            out_stream.push_str(&read_name_extra[imate as usize]);
        }
        out_stream.push('\n');
    }

    Ok((out_stream.len() - out_pos0) as u64)
}
