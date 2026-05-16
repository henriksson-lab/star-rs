#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ClipMate` at STAR/source/ClipMate.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ClipMate {
    pub type_: i32,
    pub n: u32,
    pub n_after_ad: u32,
    pub ad_seq: String,
    pub ad_seq_num: Vec<u8>,
    pub ad_mmp: f64,
    pub clipped_info: u8,
    pub clipped_ad_n: u32,
    pub clipped_ad_mm: u32,
    pub clipped_n: u32,
    pub cr4: Option<ClipCR4>,
}
