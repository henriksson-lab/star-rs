#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `sum1D` at STAR/source/serviceFuns.cpp:7. Args: a: T, N: uint"]
pub fn servicefuns_l7_sum1d<T>(a: &[T], n: u32) -> T
where
    T: Copy + Default + std::ops::AddAssign,
{
    let mut s = T::default();
    for ii in 0..std::cmp::min(n as usize, a.len()) {
        s += a[ii];
    }
    s
}

#[doc = "Original `funCompareNumbers` at STAR/source/serviceFuns.cpp:13. Args: a: void, b: void"]
pub fn servicefuns_l13_funcomparenumbers<T: Copy + Ord>(a: &T, b: &T) -> i32 {
    let va = *a;
    let vb = *b;
    if va > vb {
        1
    } else if va == vb {
        0
    } else {
        -1
    }
}

#[doc = "Original `funCompareNumbersReverse` at STAR/source/serviceFuns.cpp:26. Args: a: void, b: void"]
pub fn servicefuns_l26_funcomparenumbersreverse<T: Copy + Ord>(a: &T, b: &T) -> i32 {
    let va = *a;
    let vb = *b;
    if va > vb {
        -1
    } else if va == vb {
        0
    } else {
        1
    }
}

#[doc = "Original `funCompareNumbersReverseShift` at STAR/source/serviceFuns.cpp:39. Args: a: void, b: void"]
pub fn servicefuns_l39_funcomparenumbersreverseshift<T: Copy + Ord, const SHIFT: usize>(
    a: &[T],
    b: &[T],
) -> i32 {
    match (a.get(SHIFT), b.get(SHIFT)) {
        (Some(&va), Some(&vb)) => {
            if va > vb {
                -1
            } else if va == vb {
                0
            } else {
                1
            }
        }
        _ => 0,
    }
}

#[doc = "Original `funCompareUint1` at STAR/source/serviceFuns.cpp:53. Args: a: void, b: void"]
pub fn servicefuns_l53_funcompareuint1(a: &u32, b: &u32) -> i32 {
    let va = *a;
    let vb = *b;
    if va > vb {
        1
    } else if va == vb {
        0
    } else {
        -1
    }
}

#[doc = "Original `funCompareUint2` at STAR/source/serviceFuns.cpp:66. Args: a: void, b: void"]
pub fn servicefuns_l66_funcompareuint2(a: &[u32], b: &[u32]) -> i32 {
    match (a.get(0), b.get(0)) {
        (Some(&va), Some(&vb)) if va > vb => 1,
        (Some(&va), Some(&vb)) if va < vb => -1,
        (Some(_), Some(_)) => match (a.get(1), b.get(1)) {
            (Some(&va1), Some(&vb1)) if va1 > vb1 => 1,
            (Some(&va1), Some(&vb1)) if va1 < vb1 => -1,
            _ => 0,
        },
        _ => 0,
    }
}

#[doc = "Original `funCompareArrays` at STAR/source/serviceFuns.cpp:84. Args: a: void, b: void"]
pub fn servicefuns_l84_funcomparearrays<T: Copy + Ord, const ARRAY_SIZE: usize>(
    a: &[T],
    b: &[T],
) -> i32 {
    for ii in 0..std::cmp::min(ARRAY_SIZE, std::cmp::min(a.len(), b.len())) {
        if a[ii] > b[ii] {
            return 1;
        } else if a[ii] < b[ii] {
            return -1;
        }
    }
    0
}

#[doc = "Original `funCompareArraysReverse` at STAR/source/serviceFuns.cpp:101. Args: a: void, b: void"]
pub fn servicefuns_l101_funcomparearraysreverse<T: Copy + Ord, const ARRAY_SIZE: usize>(
    a: &[T],
    b: &[T],
) -> i32 {
    for ii in 0..std::cmp::min(ARRAY_SIZE, std::cmp::min(a.len(), b.len())) {
        if a[ii] > b[ii] {
            return -1;
        } else if a[ii] < b[ii] {
            return 1;
        }
    }
    0
}

#[doc = "Original `funCompareArraysShift` at STAR/source/serviceFuns.cpp:118. Args: a: void, b: void"]
pub fn servicefuns_l118_funcomparearraysshift<
    T: Copy + Ord,
    const ARRAY_SIZE: usize,
    const SHIFT: usize,
>(
    a: &[T],
    b: &[T],
) -> i32 {
    let available = a.len().min(b.len()).saturating_sub(SHIFT);
    for ii in 0..std::cmp::min(ARRAY_SIZE, available) {
        if a[ii + SHIFT] > b[ii + SHIFT] {
            return 1;
        } else if a[ii + SHIFT] < b[ii + SHIFT] {
            return -1;
        }
    }
    0
}

#[doc = "Original `funCompareTypeSecondFirst` at STAR/source/serviceFuns.cpp:135. Args: a: void, b: void"]
pub fn servicefuns_l135_funcomparetypesecondfirst<T: Copy + Ord>(a: &[T], b: &[T]) -> i32 {
    match (a.get(1), b.get(1)) {
        (Some(&va), Some(&vb)) if va > vb => 1,
        (Some(&va), Some(&vb)) if va < vb => -1,
        (Some(_), Some(_)) => match (a.get(0), b.get(0)) {
            (Some(&va1), Some(&vb1)) if va1 > vb1 => 1,
            (Some(&va1), Some(&vb1)) if va1 < vb1 => -1,
            _ => 0,
        },
        _ => 0,
    }
}

#[doc = "Original `funCompareTypeShift` at STAR/source/serviceFuns.cpp:153. Args: a: void, b: void"]
pub fn servicefuns_l153_funcomparetypeshift<T: Copy + Ord, const SHIFT: usize>(
    a: &[T],
    b: &[T],
) -> i32 {
    match (a.get(SHIFT), b.get(SHIFT)) {
        (Some(&va), Some(&vb)) => {
            if va > vb {
                1
            } else if va == vb {
                0
            } else {
                -1
            }
        }
        _ => 0,
    }
}

#[doc = "Original `splitString` at STAR/source/serviceFuns.cpp:167. Args: s: std::string, delim: char, elems: std::vector<std::string>"]
pub fn servicefuns_l167_splitstring(s: &str, delim: char, elems: &mut Vec<String>) -> i32 {
    let mut max_l = 0;
    elems.clear();
    for item in s.split_terminator(delim) {
        max_l = max_l.max(item.len() as i32);
        elems.push(item.to_string());
    }
    max_l
}

#[doc = "Original `binarySearch1` at STAR/source/serviceFuns.cpp:192. Args: x: argType, X: argType, N: uint32"]
pub fn servicefuns_l192_binarysearch1<T: Copy + Ord>(x: T, x_arr: &[T], n: u32) -> u32 {
    let n = std::cmp::min(n as usize, x_arr.len());
    if n == 0 {
        return u32::MAX;
    }
    if x > x_arr[n as usize - 1] || x < x_arr[0] {
        return u32::MAX;
    }
    let mut i1 = 0;
    let mut i2 = n as u32 - 1;
    while i2 > i1 + 1 {
        let i3 = (i1 + i2) / 2;
        if x_arr[i3 as usize] > x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    while i1 < n as u32 - 1 && x == x_arr[i1 as usize + 1] {
        i1 += 1;
    }
    i1
}

#[doc = "Original `binarySearch_leLeft` at STAR/source/serviceFuns.cpp:212. Args: x: argType, X: argType, N: uint32, i1: uint32"]
pub fn servicefuns_l212_binarysearch_leleft<T: Copy + Ord>(
    x: T,
    x_arr: &[T],
    n: u32,
    i1_out: &mut u32,
) -> bool {
    let n = std::cmp::min(n as usize, x_arr.len());
    if n == 0 {
        return false;
    }
    if x > x_arr[n as usize - 1] || x < x_arr[0] {
        return false;
    }
    let mut i1 = 0;
    let mut i2 = n as u32 - 1;
    while i2 > i1 + 1 {
        let i3 = (i1 + i2) / 2;
        if x_arr[i3 as usize] > x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    while i1 > 0 && x == x_arr[i1 as usize - 1] {
        i1 -= 1;
    }
    *i1_out = i1;
    true
}

#[doc = "Original `binarySearch1a` at STAR/source/serviceFuns.cpp:239. Args: x: argType, X: argType, N: int32"]
pub fn servicefuns_l239_binarysearch1a<T: Copy + Ord>(x: T, x_arr: &[T], n: i32) -> i32 {
    let n = if n <= 0 {
        0
    } else {
        std::cmp::min(n as usize, x_arr.len()) as i32
    };
    if n == 0 {
        return -1;
    }
    if x > x_arr[n as usize - 1] {
        return n - 1;
    } else if x < x_arr[0] {
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
    while i1 < n - 1 && x == x_arr[i1 as usize + 1] {
        i1 += 1;
    }
    i1
}

#[doc = "Original `binarySearch1b` at STAR/source/serviceFuns.cpp:266. Args: x: argType, X: argType, N: int32"]
pub fn servicefuns_l266_binarysearch1b<T: Copy + Ord>(x: T, x_arr: &[T], n: i32) -> i32 {
    let n = if n <= 0 {
        0
    } else {
        std::cmp::min(n as usize, x_arr.len()) as i32
    };
    if n == 0 {
        return -1;
    }
    if x > x_arr[n as usize - 1] {
        return -1;
    } else if x <= x_arr[0] {
        return 0;
    }
    let mut i1 = 0;
    let mut i2 = n - 1;
    while i2 > i1 + 1 {
        let i3 = (i1 + i2) / 2;
        if x_arr[i3 as usize] >= x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    i2
}

#[doc = "Original `binarySearchExact` at STAR/source/serviceFuns.cpp:294. Args: x: argType, X: argType, N: uint64"]
pub fn servicefuns_l294_binarysearchexact<T: Copy + Ord>(x: T, x_arr: &[T], n: u64) -> i64 {
    let n = std::cmp::min(n as usize, x_arr.len());
    if n == 0 {
        return -1;
    }
    if x > x_arr[n as usize - 1] || x < x_arr[0] {
        return -1;
    }
    let mut i1 = 0i32;
    let mut i2 = n as i32 - 1;
    while i2 > i1 + 1 {
        let i3 = (i1 + i2) / 2;
        if x_arr[i3 as usize] >= x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    if x == x_arr[i2 as usize] {
        i2 as i64
    } else if x == x_arr[i1 as usize] {
        i1 as i64
    } else {
        -1
    }
}
