#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `timeMonthDayTime` at STAR/source/TimeFunctions.cpp:4. Args: "]
pub fn timefunctions_l4_timemonthdaytime() -> String {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or_default();
    format_local_time_month_day_time(now)
}

#[doc = "Original `timeMonthDayTime` at STAR/source/TimeFunctions.cpp:14. Args: rawTime: time_t"]
pub fn timefunctions_l14_timemonthdaytime(raw_time: libc::time_t) -> String {
    format_local_time_month_day_time(raw_time)
}

#[cfg(not(windows))]
pub fn format_local_time_month_day_time_impl(raw_time: libc::time_t) -> String {
    let mut tm = std::mem::MaybeUninit::<libc::tm>::uninit();
    let tm = unsafe {
        if libc::localtime_r(&raw_time, tm.as_mut_ptr()).is_null() {
            return String::new();
        }
        tm.assume_init()
    };

    let mut buf = [0 as libc::c_char; 32];
    let len = unsafe {
        libc::strftime(
            buf.as_mut_ptr(),
            buf.len(),
            c"%b %d %H:%M:%S".as_ptr(),
            &tm,
        )
    };
    if len == 0 {
        return String::new();
    }

    unsafe { std::ffi::CStr::from_ptr(buf.as_ptr()) }
        .to_string_lossy()
        .into_owned()
}

#[cfg(windows)]
pub fn format_local_time_month_day_time_impl(raw_time: libc::time_t) -> String {
    use chrono::{Local, TimeZone};

    let timestamp = raw_time as i64;
    match Local.timestamp_opt(timestamp, 0).single() {
        Some(local_time) => local_time.format("%b %d %H:%M:%S").to_string(),
        None => String::new(),
    }
}
