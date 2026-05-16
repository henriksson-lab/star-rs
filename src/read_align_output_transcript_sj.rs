#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::outputTranscriptSJ` at STAR/source/ReadAlign_outputTranscriptSJ.cpp:4. Args: trOut: Transcript, nTrOut: uint, chunkOutSJ: OutSJ, sjReadStartN: uint"]
pub fn readalign_outputtranscriptsj_l4_readalign_outputtranscriptsj(
    tr_out: &crate::transcript::Transcript,
    n_tr_out: u32,
    chunk_out_sj: &mut crate::out_sj::OutSJ,
    sj_read_start_n: u64,
) {
    if tr_out.n_exons == 0 {
        return;
    }

    for iex in 0..tr_out.n_exons.saturating_sub(1) as usize {
        if tr_out.canon_sj[iex] >= 0 {
            let start = tr_out.exons[iex][EX_G] + tr_out.exons[iex][EX_L];
            let gap = tr_out.exons[iex + 1][EX_G] - start;
            let overhang =
                std::cmp::min(tr_out.exons[iex][EX_L], tr_out.exons[iex + 1][EX_L]) as u16;

            let mut duplicate_sj = false;
            for ii in sj_read_start_n as usize..chunk_out_sj.n as usize {
                if chunk_out_sj.junctions[ii].start == start
                    && chunk_out_sj.junctions[ii].gap == gap
                {
                    duplicate_sj = true;
                    if chunk_out_sj.junctions[ii].overhang_left < overhang {
                        chunk_out_sj.junctions[ii].overhang_left = overhang;
                        chunk_out_sj.junctions[ii].overhang_right = overhang;
                    }
                    break;
                }
            }
            if duplicate_sj {
                continue;
            }

            let motif = tr_out.canon_sj[iex];
            chunk_out_sj
                .junctions
                .push(crate::out_sj::JunctionRecord {
                    start,
                    gap,
                    overhang_left: overhang,
                    overhang_right: overhang,
                    motif,
                    strand: if motif == 0 {
                        0
                    } else {
                        ((motif + 1) % 2 + 1) as i8
                    },
                    annot: tr_out.sj_annot[iex],
                    count_unique: if n_tr_out == 1 { 1 } else { 0 },
                    count_multiple: if n_tr_out == 1 { 0 } else { 1 },
                });
            chunk_out_sj.n += 1;
        }
    }
}
