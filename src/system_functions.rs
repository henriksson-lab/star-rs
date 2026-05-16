#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `linuxProcMemory` at STAR/source/systemFunctions.cpp:6. Args: "]
pub fn systemfunctions_l6_linuxprocmemory() -> String {
    let input = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    let mut out_string = String::new();
    for str1 in input.lines() {
        if str1.starts_with("VmPeak")
            || str1.starts_with("VmSize")
            || str1.starts_with("VmHWM")
            || str1.starts_with("VmRSS")
        {
            out_string += str1;
            out_string += "; ";
        }
    }
    out_string += "\n";
    out_string
}
