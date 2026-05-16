#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ChimericAlign::chimericJunctionOutput` at STAR/source/ChimericAlign_chimericJunctionOutput.cpp:4. Args: outStream: fstream, chimN: uint, maxNonChimAlignScore: int, PEmerged_flag: bool, chimScoreBest: int, maxPossibleAlignScore: int"]
pub fn chimericalign_chimericjunctionoutput_l4_chimericalign_chimericjunctionoutput(
    chim: &crate::chimeric_align::ChimericAlign,
    map_gen: &crate::genome::Genome,
    p: &crate::parameters_chimeric::Parameters,
    read_name: &str,
    read_files_index: u32,
    chim_n: u64,
    max_non_chim_align_score: i32,
    pe_merged_flag: bool,
    chim_score_best: i32,
    max_possible_align_score: i32,
    solo_bar: Option<&crate::solo_read_barcode::SoloReadBarcode>,
) -> String {
    let read_name_out = read_name.get(1..).unwrap_or("");
    let mut out_stream = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        map_gen.chr_name[chim.al1.chr as usize],
        chim.chim_j1 - map_gen.chr_start[chim.al1.chr as usize] + 1,
        if chim.al1.str_ == 0 { "+" } else { "-" },
        map_gen.chr_name[chim.al2.chr as usize],
        chim.chim_j2 - map_gen.chr_start[chim.al2.chr as usize] + 1,
        if chim.al2.str_ == 0 { "+" } else { "-" },
        chim.chim_motif,
        chim.chim_repeat1,
        chim.chim_repeat2,
        read_name_out,
        chim.al1.exons[0][EX_G] as u64 - map_gen.chr_start[chim.al1.chr as usize] + 1,
        transcript_generatecigarp_l4_transcript_generatecigarp(&chim.al1),
        chim.al2.exons[0][EX_G] as u64 - map_gen.chr_start[chim.al2.chr as usize] + 1,
        transcript_generatecigarp_l4_transcript_generatecigarp(&chim.al2),
        chim_n,
        max_possible_align_score,
        max_non_chim_align_score,
        chim.chim_score,
        chim_score_best,
    );
    out_stream.push_str(&format!("\t{}", if pe_merged_flag { 1 } else { 0 }));

    if p.out_sam_attr_present.rg {
        out_stream.push_str(&format!(
            "\t{}",
            p.out_sam_attr_rg[read_files_index as usize]
        ));
    }
    if p.p_solo.solo_type > 0 {
        if let Some(solo_bar) = solo_bar {
            out_stream.push_str(&format!("\t{}\t{}", solo_bar.cb_seq, solo_bar.umi_seq));
        }
    }
    out_stream.push('\n');
    out_stream
}
