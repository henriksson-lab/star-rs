#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `sjAlignSplit` at STAR/source/sjAlignSplit.cpp:3. Args: a1: uint, aLength: uint, mapGen: Genome, a1D: uint, aLengthD: uint, a1A: uint, aLengthA: uint, isj: uint"]
pub fn sjalignsplit_l3_sjalignsplit(
    a1: u64,
    a_length: u64,
    map_gen: &crate::genome::Genome,
    a1_d: &mut u64,
    a_length_d: &mut u64,
    a1_a: &mut u64,
    a_length_a: &mut u64,
    isj: &mut u64,
) -> bool {
    let sj1 = (a1 - map_gen.sj_gstart) % map_gen.sjdb_length as u64;
    if sj1 < map_gen.sjdb_overhang as u64 && sj1 + a_length > map_gen.sjdb_overhang as u64 {
        *isj = (a1 - map_gen.sj_gstart) / map_gen.sjdb_length as u64;
        *a_length_d = map_gen.sjdb_overhang as u64 - sj1;
        *a_length_a = a_length - *a_length_d;
        *a1_d = map_gen.sj_dstart[*isj as usize] + sj1;
        *a1_a = map_gen.sj_astart[*isj as usize];
        true
    } else {
        false
    }
}
