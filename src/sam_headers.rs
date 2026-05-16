#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `samHeaders` at STAR/source/samHeaders.cpp:5. Args: P: Parameters, genomeOut: Genome, transcriptomeMain: Transcriptome"]
pub fn samheaders_l5_samheaders(
    p: &mut crate::parameters_chimeric::Parameters,
    genome_out: &mut crate::genome::Genome,
    transcriptome_main: &crate::transcriptome::Transcriptome,
    extra_references_txt: &str,
    comment_file_txt: &str,
) {
    const STAR_VERSION: &str = "2.7.11b";

    if p.quant_tr_sam_bam_yes {
        let mut sam_header_stream = String::new();
        let mut tr_length = Vec::new();
        for ii in 0..transcriptome_main.tr_id.len() {
            let iex1 = transcriptome_main.tr_ex_i[ii] as usize
                + transcriptome_main.tr_ex_n[ii] as usize
                - 1;
            let length = transcriptome_main.ex_len_cum[iex1]
                + transcriptome_main.ex_se[2 * iex1 + 1]
                - transcriptome_main.ex_se[2 * iex1]
                + 1;
            tr_length.push(length);
            sam_header_stream.push_str("@SQ\tSN:");
            sam_header_stream.push_str(&transcriptome_main.tr_id[ii]);
            sam_header_stream.push_str("\tLN:");
            sam_header_stream.push_str(&length.to_string());
            sam_header_stream.push('\n');
        }
        for rg in &p.out_sam_attr_rgline_split {
            sam_header_stream.push_str("@RG\t");
            sam_header_stream.push_str(rg);
            sam_header_stream.push('\n');
        }
        p.out_quant_bam_header = bamfunctions_l77_outbamwriteheader(
            &sam_header_stream,
            &transcriptome_main.tr_id,
            &tr_length,
        );
    }

    if p.out_sam_mode == "None" || p.out_sam_type[0] == "None" {
        return;
    }

    let mut sam_header_stream = String::new();

    for ii in 0..genome_out.n_chr_real as usize {
        sam_header_stream.push_str("@SQ\tSN:");
        sam_header_stream.push_str(&genome_out.chr_name[ii]);
        sam_header_stream.push_str("\tLN:");
        sam_header_stream.push_str(&genome_out.chr_length[ii].to_string());
        sam_header_stream.push('\n');
    }

    genome_out.chr_name_all = genome_out.chr_name.clone();
    genome_out.chr_length_all = genome_out
        .chr_length
        .iter()
        .map(|length| *length as u32)
        .collect();

    for line1 in extra_references_txt.lines() {
        let mut stream1 = line1.split_whitespace();
        let field1 = stream1.next().unwrap_or("");
        if field1 != "" {
            sam_header_stream.push_str(line1);
            sam_header_stream.push('\n');
            let name_field = stream1.next().unwrap_or("");
            genome_out
                .chr_name_all
                .push(name_field.get(3..).unwrap_or("").to_string());
            let length_field = stream1.next().unwrap_or("");
            let length = length_field
                .get(3..)
                .unwrap_or("")
                .parse::<u64>()
                .unwrap_or(0) as u32;
            genome_out.chr_length_all.push(length);
        }
    }

    if p.out_sam_header_pg.first().map(String::as_str) != Some("-") {
        sam_header_stream.push_str(
            p.out_sam_header_pg
                .first()
                .map(String::as_str)
                .unwrap_or(""),
        );
        for ii in 1..p.out_sam_header_pg.len() {
            sam_header_stream.push('\t');
            sam_header_stream.push_str(&p.out_sam_header_pg[ii]);
        }
        if !p.out_sam_header_pg.is_empty() {
            sam_header_stream.push('\n');
        }
    }

    sam_header_stream.push_str("@PG\tID:STAR\tPN:STAR\tVN:");
    sam_header_stream.push_str(STAR_VERSION);
    sam_header_stream.push_str("\tCL:");
    sam_header_stream.push_str(&p.command_line_full);
    sam_header_stream.push('\n');

    if p.out_sam_header_comment_file != "-" {
        for line1 in comment_file_txt.lines() {
            if line1
                .chars()
                .any(|c| !" \t\n\u{000B}\u{000C}\r".contains(c))
            {
                sam_header_stream.push_str(line1);
                sam_header_stream.push('\n');
            }
        }
    }

    for rg in &p.out_sam_attr_rgline_split {
        sam_header_stream.push_str("@RG\t");
        sam_header_stream.push_str(rg);
        sam_header_stream.push('\n');
    }

    sam_header_stream.push_str("@CO\tuser command line: ");
    sam_header_stream.push_str(&p.command_line);
    sam_header_stream.push('\n');

    sam_header_stream.push_str(&p.sam_header_extra);

    if p.out_sam_header_hd.first().map(String::as_str) != Some("-") {
        p.sam_header_hd = p
            .out_sam_header_hd
            .first()
            .cloned()
            .unwrap_or_else(|| "@HD\tVN:1.4".to_string());
        for ii in 1..p.out_sam_header_hd.len() {
            p.sam_header_hd.push('\t');
            p.sam_header_hd.push_str(&p.out_sam_header_hd[ii]);
        }
    } else {
        p.sam_header_hd = "@HD\tVN:1.4".to_string();
    }

    p.sam_header = p.sam_header_hd.clone() + "\n" + &sam_header_stream;
    p.sam_header_sorted_coord = p.sam_header_hd.clone()
        + if p.out_sam_header_hd.is_empty() {
            ""
        } else {
            "\tSO:coordinate"
        }
        + "\n"
        + &sam_header_stream;

    if p.out_sam_bool {
        p.out_sam_contents.push_str(&p.sam_header);
    }
    if p.out_bam_unsorted {
        p.out_bam_unsorted_header = bamfunctions_l77_outbamwriteheader(
            &p.sam_header,
            &genome_out.chr_name_all,
            &genome_out.chr_length_all,
        );
    }
}
