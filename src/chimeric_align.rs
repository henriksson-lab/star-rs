#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ChimericAlign` at STAR/source/ChimericAlign.h:14."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericAlign {
    pub seg1: ChimericSegment,
    pub seg2: ChimericSegment,
    pub chim_score: i32,
    pub chim_j1: u32,
    pub chim_j2: u32,
    pub chim_repeat1: u32,
    pub chim_repeat2: u32,
    pub chim_motif: i32,
    pub chim_str: i32,
    pub stitching_done: bool,
    pub al1: Transcript,
    pub al2: Transcript,
    pub ex1: u32,
    pub ex2: u32,
    pub junction_overhang_min: u32,
}

#[doc = "Original `ChimericAlign::ChimericAlign` at STAR/source/ChimericAlign.cpp:3. Args: seg1in: ChimericSegment, seg2in: ChimericSegment, chimScoreIn: int, genomeIn: Genome, RAin: ReadAlign"]
pub fn chimericalign_l3_chimericalign_chimericalign(
    seg1: crate::chimeric_segment::ChimericSegment,
    seg2: crate::chimeric_segment::ChimericSegment,
    chim_score: i32,
    junction_overhang_min: u32,
) -> crate::chimeric_align::ChimericAlign {
    let mut al1 = seg1.align.clone();
    let mut al2 = seg2.align.clone();
    if al1.ro_start > al2.ro_start {
        std::mem::swap(&mut al1, &mut al2);
    }
    let ex1 = if al1.str_ == 1 {
        0
    } else {
        al1.n_exons.saturating_sub(1)
    };
    let ex2 = if al2.str_ == 0 {
        0
    } else {
        al2.n_exons.saturating_sub(1)
    };
    crate::chimeric_align::ChimericAlign {
        seg1,
        seg2,
        chim_score,
        chim_j1: 0,
        chim_j2: 0,
        chim_repeat1: 0,
        chim_repeat2: 0,
        chim_motif: 0,
        chim_str: 0,
        stitching_done: false,
        al1,
        al2,
        ex1,
        ex2,
        junction_overhang_min,
    }
}

#[doc = "Original `ChimericAlign::chimericCheck` at STAR/source/ChimericAlign.cpp:17. Args: "]
pub fn chimericalign_l17_chimericalign_chimericcheck(
    chim: &crate::chimeric_align::ChimericAlign,
) -> bool {
    let ex1 = chim.ex1 as usize;
    let ex2 = chim.ex2 as usize;
    let mut chim_good = true;
    chim_good = chim_good && chim.al1.exons[ex1][EX_IFRAG] <= chim.al2.exons[ex2][EX_IFRAG];
    chim_good = chim_good
        && (chim.al1.exons[ex1][EX_IFRAG] < chim.al2.exons[ex2][EX_IFRAG]
            || (chim.al1.exons[ex1][EX_L] >= chim.junction_overhang_min
                && chim.al2.exons[ex2][EX_L] >= chim.junction_overhang_min));
    chim_good
}
