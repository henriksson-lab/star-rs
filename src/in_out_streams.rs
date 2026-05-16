#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `InOutStreams` at STAR/source/InOutStreams.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InOutStreams {
    pub log_stdout_attached: bool,
    pub out_sam_attached: bool,
    pub log_stdout_flushed: bool,
    pub out_sam_flushed: bool,
    pub log_stdout_file_flushed: bool,
    pub out_sam_file_flushed: bool,
    pub out_chim_sam_flushed: bool,
    pub out_chim_junction_flushed: bool,
    pub log_progress_flushed: bool,
    pub log_main_flushed: bool,
    pub log_final_flushed: bool,
    pub out_local_chains_flushed: bool,
    pub out_sam_file_closed: bool,
    pub out_chim_sam_closed: bool,
    pub out_chim_junction_closed: bool,
    pub log_progress_closed: bool,
    pub log_final_closed: bool,
    pub out_local_chains_closed: bool,
    pub out_unmapped_reads_open: [bool; 2],
    pub out_unmapped_reads_flushed: [bool; 2],
    pub out_unmapped_reads_closed: [bool; 2],
}

#[doc = "Original `InOutStreams::InOutStreams` at STAR/source/InOutStreams.cpp:4. Args: "]
pub fn inoutstreams_l4_inoutstreams_inoutstreams() -> crate::in_out_streams::InOutStreams {
    crate::in_out_streams::InOutStreams {
        log_stdout_attached: false,
        out_sam_attached: false,
        ..Default::default()
    }
}

#[doc = "Original `InOutStreams::~InOutStreams` at STAR/source/InOutStreams.cpp:11. Args: "]
pub fn inoutstreams_l11_inoutstreams_inoutstreams(
    streams: &mut crate::in_out_streams::InOutStreams,
) {
    if streams.log_stdout_attached {
        streams.log_stdout_flushed = true;
    }
    if streams.out_sam_attached {
        streams.out_sam_flushed = true;
    }

    streams.log_stdout_file_flushed = true;
    streams.out_sam_file_flushed = true;
    streams.out_chim_sam_flushed = true;
    streams.out_chim_junction_flushed = true;
    streams.log_progress_flushed = true;
    streams.log_main_flushed = true;
    streams.log_final_flushed = true;
    streams.out_local_chains_flushed = true;

    streams.out_sam_file_closed = true;
    streams.out_chim_sam_closed = true;
    streams.out_chim_junction_closed = true;
    streams.log_progress_closed = true;
    streams.log_final_closed = true;
    streams.out_local_chains_closed = true;

    for ii in 0..2 {
        if streams.out_unmapped_reads_open[ii] {
            streams.out_unmapped_reads_flushed[ii] = true;
            streams.out_unmapped_reads_closed[ii] = true;
            streams.out_unmapped_reads_open[ii] = false;
        }
    }
}
