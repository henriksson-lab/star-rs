#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::outputVariation` at STAR/source/ReadAlign_outputVariation.cpp:3. Args: Var: Variation, Tr: Transcript, iTr: uint, nTr: uint"]
pub fn readalign_outputvariation_l3_readalign_outputvariation(
    var: &crate::variation::Variation,
    _tr: &crate::transcript::Transcript,
    _i_tr: u32,
    _n_tr: u32,
) {
    if !var.yes {
        return;
    }
}
