#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::mappedFilter` at STAR/source/ReadAlign_mappedFilter.cpp:3. Args: "]
pub fn readalign_mappedfilter_l3_readalign_mappedfilter(
    read_align: &mut crate::read_align::ReadAlign,
    out_filter_score_min: i32,
    out_filter_score_min_over_lread: f64,
    out_filter_match_nmin: u32,
    out_filter_match_nmin_over_lread: f64,
    out_filter_mismatch_nover_lmax: f64,
    out_filter_multimap_nmax: u64,
) {
    read_align.unmap_type = -1;
    if read_align.n_w == 0 {
        read_align.stats_ra.unmapped_other += 1;
        read_align.unmap_type = 0;
    } else if read_align.tr_best.max_score < out_filter_score_min
        || read_align.tr_best.max_score
            < (out_filter_score_min_over_lread * (read_align.l_read - 1) as f64) as i32
        || read_align.tr_best.n_match < out_filter_match_nmin
        || read_align.tr_best.n_match
            < (out_filter_match_nmin_over_lread * (read_align.l_read - 1) as f64) as u32
    {
        read_align.stats_ra.unmapped_short += 1;
        read_align.unmap_type = 1;
    } else if read_align.tr_best.n_mm > read_align.out_filter_mismatch_nmax_total
        || read_align.tr_best.n_mm as f64 / read_align.tr_best.r_length as f64
            > out_filter_mismatch_nover_lmax
    {
        read_align.stats_ra.unmapped_mismatch += 1;
        read_align.unmap_type = 2;
    } else if read_align.n_tr > out_filter_multimap_nmax {
        read_align.stats_ra.unmapped_multi += 1;
        read_align.unmap_type = 3;
    }
}
