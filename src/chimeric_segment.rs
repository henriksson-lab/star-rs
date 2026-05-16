#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ChimericSegment` at STAR/source/ChimericSegment.h:9."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericSegment {
    pub align: Transcript,
    pub ro_s: u32,
    pub ro_e: u32,
    pub str_: u32,
    pub segment_min: u32,
    pub segment_read_gap_max: u32,
}

#[doc = "Original `ChimericSegment::ChimericSegment` at STAR/source/ChimericSegment.cpp:3. Args: Pin: Parameters, alignIn: Transcript"]
pub fn chimericsegment_l3_chimericsegment_chimericsegment(
    p_ch: &crate::parameters_chimeric::ParametersChimeric,
    align_in: &crate::transcript::Transcript,
) -> crate::chimeric_segment::ChimericSegment {
    let mut seg = crate::chimeric_segment::ChimericSegment {
        align: align_in.clone(),
        segment_min: p_ch.segment_min,
        segment_read_gap_max: p_ch.segment_read_gap_max,
        ..Default::default()
    };

    if (seg.align.intron_motifs[1] == 0 && seg.align.intron_motifs[2] == 0)
        || (seg.align.intron_motifs[1] > 0 && seg.align.intron_motifs[2] > 0)
    {
        seg.str_ = 0;
    } else if (seg.align.str_ == 0) == (seg.align.intron_motifs[1] > 0) {
        seg.str_ = 1;
    } else {
        seg.str_ = 2;
    }

    seg.ro_s = if seg.align.str_ == 0 {
        seg.align.exons[0][EX_R]
    } else {
        seg.align.l_read
            - seg.align.exons[seg.align.n_exons as usize - 1][EX_R]
            - seg.align.exons[seg.align.n_exons as usize - 1][EX_L]
    };
    seg.ro_e = if seg.align.str_ == 0 {
        seg.align.exons[seg.align.n_exons as usize - 1][EX_R]
            + seg.align.exons[seg.align.n_exons as usize - 1][EX_L]
            - 1
    } else {
        seg.align.l_read - seg.align.exons[0][EX_R] - 1
    };
    if seg.ro_s > seg.align.read_length[0] {
        seg.ro_s -= 1;
    }
    if seg.ro_e > seg.align.read_length[0] {
        seg.ro_e -= 1;
    }
    seg
}

#[doc = "Original `ChimericSegment::segmentCheck` at STAR/source/ChimericSegment.cpp:19. Args: "]
pub fn chimericsegment_l19_chimericsegment_segmentcheck(
    seg: &crate::chimeric_segment::ChimericSegment,
) -> bool {
    let mut seg_good = true;
    seg_good = seg_good && seg.align.r_length >= seg.segment_min;
    seg_good = seg_good && seg.align.intron_motifs[0] == 0;
    seg_good
}
