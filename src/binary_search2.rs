#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `binarySearch2` at STAR/source/binarySearch2.cpp:3. Args: x: uint, y: uint, X: uint, Y: uint, N: int"]
pub fn binarysearch2_l3_binarysearch2(x: u32, y: u32, x_arr: &[u32], y_arr: &[u32], n: i32) -> i32 {
    let n = if n <= 0 {
        0
    } else {
        std::cmp::min(n as usize, std::cmp::min(x_arr.len(), y_arr.len())) as i32
    };
    if n == 0 || x > x_arr[n as usize - 1] || x < x_arr[0] {
        return -1;
    }

    let mut i1 = 0;
    let mut i2 = n - 1;
    while i2 > i1 + 1 {
        let i3 = (i1 + i2) / 2;
        if x_arr[i3 as usize] > x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }

    let i3;
    if x == x_arr[i1 as usize] {
        i3 = i1;
    } else if x == x_arr[i2 as usize] {
        i3 = i2;
    } else {
        return -1;
    }

    for jj in (0..=i3).rev() {
        if x != x_arr[jj as usize] {
            break;
        } else if y == y_arr[jj as usize] {
            return jj;
        }
    }

    for jj in i3..n {
        if x != x_arr[jj as usize] {
            return -1;
        } else if y == y_arr[jj as usize] {
            return jj;
        }
    }

    -2
}
