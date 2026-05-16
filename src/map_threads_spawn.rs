#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `mapThreadsSpawn` at STAR/source/mapThreadsSpawn.cpp:6. Args: P: Parameters, RAchunk: ReadAlignChunk"]
pub fn mapthreadsspawn_l6_mapthreadsspawn<F>(
    run_thread_n: i32,
    create_status: &[i32],
    join_status: &[i32],
    mut process_chunks_main: F,
) -> Result<String, String>
where
    F: FnMut() -> Result<String, String>,
{
    let mut log_main = String::new();

    for ithread in 1..run_thread_n {
        let thread_status = create_status
            .get(ithread as usize)
            .copied()
            .unwrap_or_default();
        if thread_status > 0 {
            return Err(format!(
                "EXITING because of FATAL ERROR: phtread error while creating thread # {}, error code: {}",
                ithread, thread_status
            ));
        }
        log_main.push_str(&format!("Created thread # {}\n", ithread));
    }

    log_main.push_str(&process_chunks_main()?);

    for ithread in 1..run_thread_n {
        let thread_status = join_status
            .get(ithread as usize)
            .copied()
            .unwrap_or_default();
        if thread_status > 0 {
            return Err(format!(
                "EXITING because of FATAL ERROR: phtread error while joining thread # {}, error code: {}",
                ithread, thread_status
            ));
        }
        log_main.push_str(&format!("Joined thread # {}\n", ithread));
    }

    Ok(log_main)
}
