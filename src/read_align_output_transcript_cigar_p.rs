#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::outputTranscriptCIGARp` at STAR/source/ReadAlign_outputTranscriptCIGARp.cpp:4. Args: trOut: Transcript"]
pub fn readalign_outputtranscriptcigarp_l4_readalign_outputtranscriptcigarp(
    read_align: &crate::read_align::ReadAlign,
    tr_out: &crate::transcript::Transcript,
    read_files_in_n: usize,
) -> String {
    let mut tr = tr_out.clone();
    tr.read_nmates = if read_files_in_n > 1 { 2 } else { 1 };
    tr.read_length_original = read_align.read_length_original.clone();
    tr.read_length_pair_original = read_align.read_length_pair_original;
    transcript_generatecigarp_l4_transcript_generatecigarp(&tr)
}
