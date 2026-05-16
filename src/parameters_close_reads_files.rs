#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Parameters::closeReadsFiles` at STAR/source/Parameters_closeReadsFiles.cpp:5. Args: "]
pub fn parameters_closereadsfiles_l5_parameters_closereadsfiles(
    p: &mut crate::parameters_chimeric::Parameters,
) -> Vec<i32> {
    let mut killed_pids = Vec::new();
    for imate in 0..p.read_files_in.len() {
        if p.read_in_open.get(imate).copied().unwrap_or(false) {
            p.read_in_open[imate] = false;
        }
        if p.read_files_command_pid.get(imate).copied().unwrap_or(0) > 0 {
            killed_pids.push(p.read_files_command_pid[imate]);
        }
    }
    killed_pids
}
