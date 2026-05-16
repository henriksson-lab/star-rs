#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ChimericDetection` at STAR/source/ChimericDetection.h:12."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericDetection {
    pub p: Parameters,
    pub ra: Option<ReadAlign>,
    pub tr_all: Vec<Vec<Transcript>>,
    pub n_w: u32,
    pub n_win_tr: Vec<u32>,
    pub read1: [Vec<u8>; 2],
    pub out_gen: Genome,
    pub chim_aligns: Vec<ChimericAlign>,
    pub ostream_chim_junction_attached: bool,
}

#[doc = "Original `ChimericDetection::ChimericDetection` at STAR/source/ChimericDetection.cpp:3. Args: Pin: Parameters, trAll: Transcript, nWinTr: uint, Read1in: char, mapGenIn: Genome, ostreamChimJunctionIn: fstream, RAin: ReadAlign"]
pub fn chimericdetection_l3_chimericdetection_chimericdetection(
    p: crate::parameters_chimeric::Parameters,
    tr_all: Vec<Vec<crate::transcript::Transcript>>,
    n_win_tr: Vec<u32>,
    read1: [Vec<u8>; 2],
    out_gen: crate::genome::Genome,
    ostream_chim_junction_attached: bool,
    ra: crate::read_align::ReadAlign,
) -> crate::chimeric_detection::ChimericDetection {
    crate::chimeric_detection::ChimericDetection {
        p,
        ra: Some(ra),
        n_w: n_win_tr.len() as u32,
        tr_all,
        n_win_tr,
        read1,
        out_gen,
        ostream_chim_junction_attached,
        chim_aligns: Vec::new(),
    }
}
