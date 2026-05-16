#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `funCompareUintAndSuffixes` at STAR/source/funCompareUintAndSuffixes.cpp:6. Args: a: void, b: void"]
pub fn funcompareuintandsuffixes_l6_funcompareuintandsuffixes(
    a: &[u64],
    b: &[u64],
    g: &[u8],
) -> i32 {
    if a[0] > b[0] {
        1
    } else if a[0] < b[0] {
        -1
    } else {
        let ia = a[1] as usize;
        let ib = b[1] as usize;
        let mut ig = 0usize;
        loop {
            if g[ia + ig] > g[ib + ig] {
                return 1;
            } else if g[ia + ig] < g[ib + ig] {
                return -1;
            } else if g[ia + ig] == 5 {
                if a[1] > b[1] {
                    return 1;
                } else {
                    return -1;
                }
            } else {
                ig += 1;
            }
        }
    }
}
