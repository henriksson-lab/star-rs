#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ReadAnnotFeature` at STAR/source/ReadAnnotations.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAnnotFeature {
    pub f_set: std::collections::BTreeSet<u32>,
    pub f_align: Vec<std::collections::BTreeSet<u32>>,
    pub ov_type: u32,
}

#[doc = "Original class `ReadAnnotations` at STAR/source/ReadAnnotations.h:20."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAnnotations {
    pub annot_features: Vec<ReadAnnotFeature>,
    pub gene_exon_overlap: Vec<i32>,
    pub transcript_concordant: Vec<[u32; 2]>,
    pub gene_velocyto_simple: [u32; 2],
    pub tr_velocyto_type: Vec<TrTypeStruct>,
}
