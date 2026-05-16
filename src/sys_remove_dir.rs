#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `removeFileOrDir` at STAR/source/sysRemoveDir.cpp:7. Args: fpath: char, sb: struct stat, typeflag: int, ftwbuf: struct FTW"]
pub fn sysremovedir_l7_removefileordir(fpath: &std::path::Path, typeflag: i32) -> i32 {
    const FTW_F: i32 = 0;
    const FTW_DP: i32 = 5;
    if typeflag == FTW_F {
        let _ = std::fs::remove_file(fpath);
    } else if typeflag == FTW_DP {
        let _ = std::fs::remove_dir(fpath);
    } else {
        return -1;
    }
    0
}

#[doc = "Original `sysRemoveDir` at STAR/source/sysRemoveDir.cpp:25. Args: dirName: std::string"]
pub fn sysremovedir_l25_sysremovedir(dir_name: &std::path::Path) -> std::io::Result<()> {
    std::fs::remove_dir_all(dir_name)
}
