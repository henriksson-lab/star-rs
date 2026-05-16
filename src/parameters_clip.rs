#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ParametersClip` at STAR/source/ParametersClip.h:21."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParametersClip {
    pub adapter_type: Vec<String>,
    pub in_: [ReadClipInput; 2],
    pub read_nmates: u32,
    pub read_nends: u32,
}

#[doc = "Original class `ReadClipInput` at STAR/source/ParametersClip.h:12."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadClipInput {
    pub n: Vec<u32>,
    pub n_after_ad: Vec<u32>,
    pub ad_seq: Vec<String>,
    pub ad_mmp: Vec<f64>,
}
