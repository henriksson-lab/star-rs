#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `funCompareUintAndSuffixesMemcmp` at STAR/source/funCompareUintAndSuffixesMemcmp.cpp:7. Args: a: void, b: void"]
pub fn funcompareuintandsuffixesmemcmp_l7_funcompareuintandsuffixesmemcmp(
    a: &[u64],
    b: &[u64],
    g: &[u8],
    l: usize,
) -> i32 {
    if a[0] > b[0] {
        1
    } else if a[0] < b[0] {
        -1
    } else {
        let ia = a[1] as usize;
        let ib = b[1] as usize;
        if l == 0 {
            return if a[1] > b[1] { 1 } else { -1 };
        }
        let comp = unsafe {
            libc::memcmp(
                g.as_ptr().add(ia) as *const libc::c_void,
                g.as_ptr().add(ib) as *const libc::c_void,
                l,
            )
        };
        if comp != 0 {
            comp
        } else if a[1] > b[1] {
            1
        } else {
            -1
        }
    }
}
