#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ThreadControl` at STAR/source/ThreadControl.h:9."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ThreadControl {
    pub chunk_in_n: u32,
    pub chunk_out_n: u32,
}

#[doc = "Original `ThreadControl::ThreadControl` at STAR/source/ThreadControl.cpp:3. Args: "]
pub fn threadcontrol_l3_threadcontrol_threadcontrol() -> crate::thread_control::ThreadControl {
    crate::thread_control::ThreadControl {
        chunk_in_n: 0,
        chunk_out_n: 0,
    }
}
