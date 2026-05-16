#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `medianUint2` at STAR/source/SuffixArrayFuns.cpp:4. Args: a: uint, b: uint"]
pub fn suffixarrayfuns_l4_medianuint2(a: u64, b: u64) -> u64 {
    a / 2 + b / 2 + (a % 2 + b % 2) / 2
}

#[doc = "Original `compareSeqToGenome` at STAR/source/SuffixArrayFuns.cpp:10. Args: mapGen: Genome, s2: char, S: uint, N: uint, L: uint, iSA: uint, dirR: bool, compRes: bool"]
pub fn suffixarrayfuns_l10_compareseqtogenome(
    map_gen: &crate::genome::Genome,
    s2: [&[u8]; 2],
    s: u64,
    n: u64,
    l: u64,
    i_sa: u64,
    dir_r: bool,
    comp_res: &mut bool,
) -> u64 {
    if l >= n {
        return n;
    }
    let mut sastr = genome_sa_index_value(map_gen, i_sa);
    let dir_g = (sastr >> map_gen.gstrand_bit) == 0;
    sastr &= map_gen.gstrand_mask as u64;
    let s_usize = s as usize;
    let l_usize = l as usize;
    let n_minus_l = n - l;

    if dir_r && dir_g {
        let s_start = s_usize + l_usize;
        let g_start = sastr as usize + l_usize;
        let s_ptr = unsafe { s2[0].as_ptr().add(s_start) };
        let g_ptr = unsafe { map_gen.g.as_ptr().add(g_start) };
        for ii in 0..n_minus_l as usize {
            let sb = unsafe { *s_ptr.add(ii) };
            let gb = unsafe { *g_ptr.add(ii) };
            if sb != gb {
                *comp_res = sb > gb;
                return ii as u64 + l;
            }
        }
        n
    } else if dir_r && !dir_g {
        let s_start = s_usize + l_usize;
        let Some(g_start) = map_gen
            .n_genome
            .checked_sub(1)
            .and_then(|v| v.checked_sub(sastr))
            .and_then(|v| v.checked_sub(l))
            .map(|v| v as usize)
        else {
            *comp_res = false;
            return l;
        };
        let s_ptr = unsafe { s2[1].as_ptr().add(s_start) };
        let g_ptr = unsafe { map_gen.g.as_ptr().add(g_start) };
        for ii in 0..n_minus_l as usize {
            let sb = unsafe { *s_ptr.add(ii) };
            let gb = unsafe { *g_ptr.sub(ii) };
            if sb != gb {
                *comp_res = !(sb > gb || gb > 3);
                return ii as u64 + l;
            }
        }
        n
    } else if !dir_r && dir_g {
        let Some(s_start) = s_usize.checked_sub(l_usize) else {
            *comp_res = false;
            return l;
        };
        let g_start = sastr as usize + l_usize;
        let s_ptr = unsafe { s2[1].as_ptr().add(s_start) };
        let g_ptr = unsafe { map_gen.g.as_ptr().add(g_start) };
        for ii in 0..n_minus_l as usize {
            let sb = unsafe { *s_ptr.sub(ii) };
            let gb = unsafe { *g_ptr.add(ii) };
            if sb != gb {
                *comp_res = sb > gb;
                return ii as u64 + l;
            }
        }
        n
    } else {
        let Some(s_start) = s_usize.checked_sub(l_usize) else {
            *comp_res = false;
            return l;
        };
        let Some(g_start) = map_gen
            .n_genome
            .checked_sub(1)
            .and_then(|v| v.checked_sub(sastr))
            .and_then(|v| v.checked_sub(l))
            .map(|v| v as usize)
        else {
            *comp_res = false;
            return l;
        };
        let s_ptr = unsafe { s2[0].as_ptr().add(s_start) };
        let g_ptr = unsafe { map_gen.g.as_ptr().add(g_start) };
        for ii in 0..n_minus_l as usize {
            let sb = unsafe { *s_ptr.sub(ii) };
            let gb = unsafe { *g_ptr.sub(ii) };
            if sb != gb {
                *comp_res = !(sb > gb || gb > 3);
                return ii as u64 + l;
            }
        }
        n
    }
}

#[doc = "Original `findMultRange` at STAR/source/SuffixArrayFuns.cpp:106. Args: mapGen: Genome, i3: uint, L3: uint, i1: uint, L1: uint, i1a: uint, L1a: uint, i1b: uint, L1b: uint, s: char, dirR: bool, S: uint"]
pub fn suffixarrayfuns_l106_findmultrange(
    map_gen: &crate::genome::Genome,
    i3: u64,
    l3: u64,
    i1: u64,
    l1: u64,
    mut i1a: u64,
    l1a: u64,
    mut i1b: u64,
    mut l1b: u64,
    s: [&[u8]; 2],
    dir_r: bool,
    s_offset: u64,
) -> u64 {
    if l1 < l3 {
        l1b = l1;
        i1b = i1;
        i1a = i3;
    } else if l1a < l1 {
        l1b = l1a;
        i1b = i1a;
        i1a = i1;
    }

    let mut comp_res = false;
    while i1b + 1 < i1a || i1b > i1a + 1 {
        let i1c = suffixarrayfuns_l4_medianuint2(i1a, i1b);
        let l1c = suffixarrayfuns_l10_compareseqtogenome(
            map_gen,
            s,
            s_offset,
            l3,
            l1b,
            i1c,
            dir_r,
            &mut comp_res,
        );
        if l1c == l3 {
            i1a = i1c;
        } else {
            i1b = i1c;
            l1b = l1c;
        }
    }
    i1a
}

#[doc = "Original `maxMappableLength` at STAR/source/SuffixArrayFuns.cpp:133. Args: mapGen: Genome, s: char, S: uint, N: uint, i1: uint, i2: uint, dirR: bool, L: uint, indStartEnd: uint"]
pub fn suffixarrayfuns_l133_maxmappablelength(
    map_gen: &crate::genome::Genome,
    s: [&[u8]; 2],
    s_offset: u64,
    n: u64,
    mut i1: u64,
    mut i2: u64,
    dir_r: bool,
    l: &mut u64,
    ind_start_end: &mut [u64; 2],
) -> u64 {
    let mut comp_res = false;
    let initial_l = *l;
    let original_i1 = i1;
    let original_i2 = i2;

    let mut l1 = suffixarrayfuns_l10_compareseqtogenome(
        map_gen,
        s,
        s_offset,
        n,
        *l,
        i1,
        dir_r,
        &mut comp_res,
    );
    let mut l2 = suffixarrayfuns_l10_compareseqtogenome(
        map_gen,
        s,
        s_offset,
        n,
        *l,
        i2,
        dir_r,
        &mut comp_res,
    );

    *l = std::cmp::min(l1, l2);

    let mut l1a = l1;
    let mut l1b = l1;
    let mut i1a = i1;
    let mut i1b = i1;
    let mut l2a = l2;
    let mut l2b = l2;
    let mut i2a = i2;
    let mut i2b = i2;

    let mut i3 = i1;
    let mut l3 = l1;
    while i1 + 1 < i2 {
        i3 = suffixarrayfuns_l4_medianuint2(i1, i2);
        l3 = suffixarrayfuns_l10_compareseqtogenome(
            map_gen,
            s,
            s_offset,
            n,
            *l,
            i3,
            dir_r,
            &mut comp_res,
        );

        if l3 == n {
            break;
        }

        if comp_res {
            if l3 > l1 {
                l1b = l1a;
                l1a = l1;
                i1b = i1a;
                i1a = i1;
            }
            i1 = i3;
            l1 = l3;
        } else {
            if l3 > l2 {
                l2b = l2a;
                l2a = l2;
                i2b = i2a;
                i2a = i2;
            }
            i2 = i3;
            l2 = l3;
        }
        *l = std::cmp::min(l1, l2);
    }

    if l3 < n {
        if l1 > l2 {
            i3 = i1;
            l3 = l1;
        } else {
            i3 = i2;
            l3 = l2;
        }
    }

    i1 = suffixarrayfuns_l106_findmultrange(
        map_gen, i3, l3, i1, l1, i1a, l1a, i1b, l1b, s, dir_r, s_offset,
    );
    i2 = suffixarrayfuns_l106_findmultrange(
        map_gen, i3, l3, i2, l2, i2a, l2a, i2b, l2b, s, dir_r, s_offset,
    );

    *l = l3;
    ind_start_end[0] = i1;
    ind_start_end[1] = i2;

    let n_rep = i2 - i1 + 1;
    if initial_l <= 4
        && original_i2 >= original_i1
        && original_i2 - original_i1 <= 10_000
        && *l == initial_l
        && n_rep > 8
    {
        let mut best_l = initial_l;
        let mut best_start = original_i1;
        let mut best_end = original_i1;
        let mut seen_best = false;
        for i in original_i1..=original_i2 {
            let li = suffixarrayfuns_l10_compareseqtogenome(
                map_gen,
                s,
                s_offset,
                n,
                initial_l,
                i,
                dir_r,
                &mut comp_res,
            );
            if li > best_l {
                best_l = li;
                best_start = i;
                best_end = i;
                seen_best = true;
            } else if li == best_l {
                if !seen_best {
                    best_start = i;
                    seen_best = true;
                }
                best_end = i;
            }
        }
        if seen_best && best_l > initial_l {
            *l = best_l;
            ind_start_end[0] = best_start;
            ind_start_end[1] = best_end;
            return best_end - best_start + 1;
        }
    }

    n_rep
}

#[doc = "Original `compareRefEnds` at STAR/source/SuffixArrayFuns.cpp:210. Args: mapGen: Genome, SAstr: uint64, gInsert: uint64, strG: bool, strR: bool"]
pub fn suffixarrayfuns_l210_comparerefends(
    n_genome: u64,
    sastr: u64,
    g_insert: u64,
    str_g: bool,
    str_r: bool,
) -> i32 {
    if str_g {
        if str_r {
            if sastr < g_insert { 1 } else { -1 }
        } else {
            1
        }
    } else if str_r {
        -1
    } else if g_insert == u64::MAX {
        -1
    } else if sastr < n_genome - g_insert {
        1
    } else {
        -1
    }
}

#[inline(always)]
pub fn genome_sa_index_value(map_gen: &crate::genome::Genome, i_sa: u64) -> u64 {
    if let Some(value) = map_gen.sa.get(i_sa as usize) {
        *value
    } else {
        packedarray_h18_packedarray_index(&map_gen.sa_packed, i_sa)
    }
}

#[doc = "Original `compareSeqToGenome1` at STAR/source/SuffixArrayFuns.cpp:221. Args: mapGen: Genome, s2: char, S: uint, N: uint, L: uint, iSA: uint, dirR: bool, gInsert: uint64, compRes: int"]
pub fn suffixarrayfuns_l221_compareseqtogenome1(
    map_gen: &crate::genome::Genome,
    s2: [&[u8]; 2],
    s: u64,
    n: u64,
    l: u64,
    i_sa: u64,
    dir_r: bool,
    g_insert: u64,
    comp_res: &mut i32,
) -> u64 {
    if l >= n {
        return n;
    }
    let mut sastr = genome_sa_index_value(map_gen, i_sa);
    let dir_g = (sastr >> map_gen.gstrand_bit) == 0;
    sastr &= map_gen.gstrand_mask as u64;
    let s_start = s as usize + l as usize;

    if dir_g {
        let g_start = sastr as usize + l as usize;
        for ii in 0..(n - l) as usize {
            let sb = s2[0][s_start + ii];
            let gb = map_gen.g[g_start + ii];
            if sb != gb {
                *comp_res = if sb > gb { 1 } else { -1 };
                return ii as u64 + l;
            } else if sb == GENOME_SPACING_CHAR {
                *comp_res = suffixarrayfuns_l210_comparerefends(
                    map_gen.n_genome,
                    sastr,
                    g_insert,
                    dir_g,
                    dir_r,
                );
                return ii as u64 + l;
            }
        }
        n
    } else {
        let Some(g_start) = map_gen
            .n_genome
            .checked_sub(1)
            .and_then(|v| v.checked_sub(sastr))
            .and_then(|v| v.checked_sub(l))
            .map(|v| v as usize)
        else {
            *comp_res = -1;
            return l;
        };
        for ii in 0..(n - l) as usize {
            let sb = s2[1][s_start + ii];
            let Some(gb) = g_start
                .checked_sub(ii)
                .and_then(|pos| map_gen.g.get(pos))
                .copied()
            else {
                *comp_res = -1;
                return ii as u64 + l;
            };
            if sb != gb {
                let mut s1 = sb;
                let mut g1 = gb;
                if s1 < 4 {
                    s1 = 3 - s1;
                }
                if g1 < 4 {
                    g1 = 3 - g1;
                }
                *comp_res = if s1 > g1 { 1 } else { -1 };
                return ii as u64 + l;
            } else if sb == GENOME_SPACING_CHAR {
                *comp_res = suffixarrayfuns_l210_comparerefends(
                    map_gen.n_genome,
                    sastr,
                    g_insert,
                    dir_g,
                    dir_r,
                );
                return ii as u64 + l;
            }
        }
        n
    }
}

#[doc = "Original `suffixArraySearch1` at STAR/source/SuffixArrayFuns.cpp:297. Args: mapGen: Genome, s: char, S: uint, N: uint, gInsert: uint64, strR: bool, i1: uint, i2: uint, L: uint"]
pub fn suffixarrayfuns_l297_suffixarraysearch1(
    map_gen: &crate::genome::Genome,
    s: [&[u8]; 2],
    s_offset: u64,
    n: u64,
    g_insert: u64,
    str_r: bool,
    mut i1: u64,
    mut i2: u64,
    l: &mut u64,
) -> u64 {
    let mut comp_res = 0i32;

    let mut l1 = suffixarrayfuns_l221_compareseqtogenome1(
        map_gen,
        s,
        s_offset,
        n,
        *l,
        i1,
        str_r,
        g_insert,
        &mut comp_res,
    );
    if comp_res < 0 {
        *l = l1;
        return 0;
    }

    let mut l2 = suffixarrayfuns_l221_compareseqtogenome1(
        map_gen,
        s,
        s_offset,
        n,
        *l,
        i2,
        str_r,
        g_insert,
        &mut comp_res,
    );
    if comp_res > 0 {
        *l = l2;
        return u64::MAX - 1;
    }

    *l = std::cmp::min(l1, l2);
    while i1 + 1 < i2 {
        let i3 = suffixarrayfuns_l4_medianuint2(i1, i2);
        let l3 = suffixarrayfuns_l221_compareseqtogenome1(
            map_gen,
            s,
            s_offset,
            n,
            *l,
            i3,
            str_r,
            g_insert,
            &mut comp_res,
        );
        if l3 == n {
            *l = n;
            return i3;
        }

        if comp_res > 0 {
            i1 = i3;
            l1 = l3;
        } else if comp_res < 0 {
            i2 = i3;
            l2 = l3;
        }
        *l = std::cmp::min(l1, l2);
    }
    i2
}

#[doc = "Original `funCalcSAiFromSA` at STAR/source/SuffixArrayFuns.cpp:353. Args: gSeq: char, gSA: PackedArray, mapGen: Genome, iSA: uint, L: int, iL4: int"]
pub fn suffixarrayfuns_l353_funcalcsaifromsa(
    g_seq: &[u8],
    g_sa: &[u64],
    n_genome: u64,
    gstrand_bit: u8,
    gstrand_mask: u64,
    i_sa: u64,
    l: i32,
    i_l4: &mut i32,
) -> u64 {
    let mut sastr = g_sa[i_sa as usize];
    let dir_g = (sastr >> gstrand_bit) == 0;
    sastr &= gstrand_mask;
    *i_l4 = -1;
    let mut saind = 0u64;

    if dir_g {
        for ii in 0..l {
            let g2 = g_seq[sastr as usize + ii as usize];
            if g2 > 3 {
                *i_l4 = ii;
                saind <<= 2 * (l - ii) as u32;
                return saind;
            }
            saind <<= 2;
            saind += g2 as u64;
        }
        saind
    } else {
        let start = n_genome as i128 - sastr as i128 - 16;
        for ii in 0..l {
            let pos = start + (15 - ii) as i128;
            let g2 = if pos >= 0 {
                g_seq
                    .get(pos as usize)
                    .copied()
                    .unwrap_or(GENOME_SPACING_CHAR)
            } else {
                GENOME_SPACING_CHAR
            };
            if g2 > 3 {
                *i_l4 = ii;
                saind <<= 2 * (l - ii) as u32;
                return saind;
            }
            saind <<= 2;
            saind += (3 - g2) as u64;
        }
        saind
    }
}

pub fn suffixarrayfuns_l353_funcalcsaifrompacked(
    g_seq: &[u8],
    g_sa: &crate::packed_array::PackedArray,
    n_genome: u64,
    gstrand_bit: u8,
    gstrand_mask: u64,
    i_sa: u64,
    l: i32,
    i_l4: &mut i32,
) -> u64 {
    let mut sastr = packedarray_h18_packedarray_index(g_sa, i_sa);
    let dir_g = (sastr >> gstrand_bit) == 0;
    sastr &= gstrand_mask;
    *i_l4 = -1;
    let mut saind = 0u64;

    if dir_g {
        for ii in 0..l {
            let g2 = g_seq[sastr as usize + ii as usize];
            if g2 > 3 {
                *i_l4 = ii;
                saind <<= 2 * (l - ii) as u32;
                return saind;
            }
            saind <<= 2;
            saind += g2 as u64;
        }
        saind
    } else {
        let start = n_genome as i128 - sastr as i128 - 16;
        for ii in 0..l {
            let pos = start + (15 - ii) as i128;
            let g2 = if pos >= 0 {
                g_seq
                    .get(pos as usize)
                    .copied()
                    .unwrap_or(GENOME_SPACING_CHAR)
            } else {
                GENOME_SPACING_CHAR
            };
            if g2 > 3 {
                *i_l4 = ii;
                saind <<= 2 * (l - ii) as u32;
                return saind;
            }
            saind <<= 2;
            saind += (3 - g2) as u64;
        }
        saind
    }
}

#[doc = "Original `funCalcSAi` at STAR/source/SuffixArrayFuns.cpp:397. Args: G: char, iL: uint"]
pub fn suffixarrayfuns_l397_funcalcsai(g: &[u8], il: u32) -> i64 {
    let mut ind1 = 0i64;
    for il1 in 0..=il as usize {
        let base = g[il1] as i64;
        if base > 3 {
            return -ind1;
        } else {
            ind1 <<= 2;
            ind1 += base;
        }
    }
    ind1
}
