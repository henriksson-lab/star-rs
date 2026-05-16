#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `exitWithError` at STAR/source/ErrorWarning.cpp:8. Args: messageOut: string, streamOut1: ostream, streamOut2: ostream, errorInt: int, P: Parameters"]
pub fn errorwarning_l8_exitwitherror(
    message_out: &str,
    stream_out1_good: bool,
    stream_out2_good: bool,
    error_int: i32,
    p: &crate::parameters_chimeric::Parameters,
    raw_time: libc::time_t,
) -> crate::out_sj::ExitWithErrorResult {
    let fatal_line = format!(
        "\n{}\n{} ...... FATAL ERROR, exiting\n",
        message_out,
        timefunctions_l14_timemonthdaytime(raw_time)
    );
    crate::out_sj::ExitWithErrorResult {
        stream_out1: if stream_out1_good {
            fatal_line.clone()
        } else {
            String::new()
        },
        stream_out2: if stream_out2_good {
            fatal_line
        } else {
            String::new()
        },
        error_int,
        thread_mutex_locked: p.run_thread_n > 1,
        in_out_deleted: true,
    }
}

#[doc = "Original `warningMessage` at STAR/source/ErrorWarning.cpp:25. Args: messageOut: string, streamOut1: ostream, streamOut2: ostream, P: Parameters"]
pub fn errorwarning_l25_warningmessage(message_out: &str) -> String {
    format!("!!!!! WARNING: {}\n", message_out)
}
