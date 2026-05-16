#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::outputSpliceGraphSAM` at STAR/source/ReadAlign_outputSpliceGraphSAM.cpp:5. Args: trOut: Transcript, nTrOut: uint, iTrOut: uint, outStream: ostream"]
pub fn readalign_outputsplicegraphsam_l5_readalign_outputsplicegraphsam(
    tr_out: &crate::transcript::Transcript,
    n_tr_out: u32,
    i_tr_out: u32,
    read_filter: u8,
    unmap_type: i32,
    read_name: &str,
    read0: &[String],
    read_file_type: i32,
    qual0: &[String],
    p: &crate::parameters_chimeric::Parameters,
    read_files_index: u32,
    read_name_extra: &[String],
    l_read: u32,
    gen_out: &crate::genome::Genome,
) -> Result<(String, u64), String> {
    let mut out_stream = String::new();
    let mut sam_flag: u16 = 0;
    if read_filter == b'Y' {
        sam_flag |= 0x200;
    }

    let read_name_out = read_name.get(1..).unwrap_or("");

    if unmap_type >= 0 {
        sam_flag |= 0x4;
        out_stream.push_str(&format!(
            "{}\t{}\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tNH:i:0\tHI:i:0\tAS:i:{}\tnM:i:{}\tuT:A:{}",
            read_name_out,
            sam_flag,
            read0[0],
            if read_file_type == 2 {
                qual0[0].as_str()
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
        if p.read_files_type_n == 10 && !read_name_extra[0].is_empty() {
            out_stream.push('\t');
            out_stream.push_str(&read_name_extra[0]);
        }
        out_stream.push('\n');
        let n = out_stream.len() as u64;
        return Ok((out_stream, n));
    }

    if !tr_out.primary_flag {
        sam_flag |= 0x100;
    }
    if tr_out.str_ == 1 {
        sam_flag |= 0x10;
    }

    let mut cigar = String::new();
    let cigar_chars = ['M', 'I', 'D', 'N', 'S'];
    for cc in &tr_out.cigar {
        cigar.push_str(&format!(
            "{}{}",
            cc[1],
            cigar_chars.get(cc[0] as usize).copied().unwrap_or('\0')
        ));
    }

    let (seq_out, qual_out) = if tr_out.str_ == 0 {
        (
            read0[0].chars().take(l_read as usize).collect::<String>(),
            qual0[0].chars().take(l_read as usize).collect::<String>(),
        )
    } else {
        let mut seq_rev = vec![0u8; l_read as usize];
        sequencefuns_l16_revcomplementnucleotides(read0[0].as_bytes(), &mut seq_rev, l_read);
        (
            String::from_utf8(seq_rev).unwrap(),
            qual0[0]
                .chars()
                .rev()
                .take(l_read as usize)
                .collect::<String>(),
        )
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

    let flag_out = (sam_flag & p.out_sam_flag_and) | p.out_sam_flag_or;
    out_stream.push_str(&format!(
        "{}\t{}\t{}\t{}\t{}\t{}",
        read_name_out,
        flag_out,
        gen_out.chr_name[tr_out.chr as usize],
        tr_out.g_start as u64 + 1 - gen_out.chr_start[tr_out.chr as usize],
        mapq,
        cigar
    ));
    out_stream.push_str("\t*\t0\t0");
    out_stream.push('\t');
    out_stream.push_str(&seq_out);
    out_stream.push('\t');
    if read_file_type == 2 && p.out_sam_mode != "NoQS" {
        out_stream.push_str(&qual_out);
    } else {
        out_stream.push('*');
    }

    let tag_nm = tr_out.n_mm + tr_out.l_ins + tr_out.l_del;
    for attr in &p.out_sam_attr_order {
        match *attr {
            ATTR_NH => out_stream.push_str(&format!("\tNH:i:{}", n_tr_out)),
            ATTR_HI => {
                out_stream.push_str(&format!("\tHI:i:{}", i_tr_out + p.out_sam_attr_ih_start))
            }
            ATTR_AS => out_stream.push_str(&format!("\tAS:i:{}", tr_out.max_score)),
            ATTR_NM_LOWER => out_stream.push_str(&format!("\tnM:i:{}", tr_out.n_mm)),
            ATTR_NM => out_stream.push_str(&format!("\tNM:i:{}", tag_nm)),
            ATTR_RG => out_stream.push_str(&format!(
                "\tRG:Z:{}",
                p.out_sam_attr_rg[read_files_index as usize]
            )),
            ATTR_JM | ATTR_JI | ATTR_XS | ATTR_MD | ATTR_CH | ATTR_CR | ATTR_CY | ATTR_UR
            | ATTR_UY | ATTR_CB | ATTR_UB | ATTR_SM | ATTR_SS | ATTR_SQ | ATTR_RB | ATTR_VG
            | ATTR_VA | ATTR_VW => {}
            _ => {
                return Err(format!(
                    "EXITING because of FATAL error: unknown/unimplemented SAM atrribute (tag): {}\nSOLUTION: contact Alex Dobin at dobin@cshl.edu\n",
                    attr
                ));
            }
        }
    }

    if p.read_files_type_n == 10 && !read_name_extra[0].is_empty() {
        out_stream.push('\t');
        out_stream.push_str(&read_name_extra[0]);
    }
    out_stream.push('\n');
    let n = out_stream.len() as u64;
    Ok((out_stream, n))
}
