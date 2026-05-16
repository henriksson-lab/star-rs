#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;
use std::borrow::Cow;

#[doc = "Original class `ChimericDetection` at STAR/source/ChimericDetection.h:12."]
#[derive(Clone, Debug, PartialEq)]
pub struct ChimericDetection<'a> {
    pub p: Cow<'a, Parameters>,
    pub ra: Option<Cow<'a, ReadAlign>>,
    pub tr_all: Cow<'a, [Vec<Transcript>]>,
    pub n_w: u64,
    pub n_win_tr: Vec<u32>,
    pub read1: [Cow<'a, [u8]>; 2],
    pub out_gen: Cow<'a, Genome>,
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
) -> crate::chimeric_detection::ChimericDetection<'static> {
    let [read1_0, read1_1] = read1;
    crate::chimeric_detection::ChimericDetection {
        p: Cow::Owned(p),
        ra: Some(Cow::Owned(ra)),
        n_w: n_win_tr.len() as u64,
        tr_all: Cow::Owned(tr_all),
        n_win_tr,
        read1: [Cow::Owned(read1_0), Cow::Owned(read1_1)],
        out_gen: Cow::Owned(out_gen),
        ostream_chim_junction_attached,
        chim_aligns: Vec::new(),
    }
}

pub fn chimericdetection_borrowed<'a>(
    p: &'a crate::parameters_chimeric::Parameters,
    tr_all: &'a [Vec<crate::transcript::Transcript>],
    n_win_tr: Vec<u32>,
    read1: [&'a [u8]; 2],
    out_gen: &'a crate::genome::Genome,
    ostream_chim_junction_attached: bool,
    ra: &'a crate::read_align::ReadAlign,
) -> crate::chimeric_detection::ChimericDetection<'a> {
    crate::chimeric_detection::ChimericDetection {
        p: Cow::Borrowed(p),
        ra: Some(Cow::Borrowed(ra)),
        n_w: n_win_tr.len() as u64,
        tr_all: Cow::Borrowed(tr_all),
        n_win_tr,
        read1: [Cow::Borrowed(read1[0]), Cow::Borrowed(read1[1])],
        out_gen: Cow::Borrowed(out_gen),
        ostream_chim_junction_attached,
        chim_aligns: Vec::new(),
    }
}
