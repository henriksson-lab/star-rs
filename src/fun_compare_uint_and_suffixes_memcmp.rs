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
        for ii in 0..l {
            if g[ia + ii] != g[ib + ii] {
                return g[ia + ii] as i32 - g[ib + ii] as i32;
            }
        }
        if a[1] > b[1] { 1 } else { -1 }
    }
}
