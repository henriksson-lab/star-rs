#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `complementSeqNumbers` at STAR/source/SequenceFuns.cpp:4. Args: ReadsIn: char, ReadsOut: char, Lread: uint"]
pub fn sequencefuns_l4_complementseqnumbers(reads_in: &[u8], reads_out: &mut [u8], lread: u32) {
    let lread = std::cmp::min(
        lread as usize,
        std::cmp::min(reads_in.len(), reads_out.len()),
    );
    for jj in 0..lread {
        reads_out[jj] = match reads_in[jj] {
            3 => 0,
            2 => 1,
            1 => 2,
            0 => 3,
            value => value,
        };
    }
}

#[doc = "Original `revComplementNucleotides` at STAR/source/SequenceFuns.cpp:16. Args: ReadsIn: char, ReadsOut: char, Lread: uint"]
pub fn sequencefuns_l16_revcomplementnucleotides(
    reads_in: &[u8],
    reads_out: &mut [u8],
    lread: u32,
) {
    let lread = std::cmp::min(
        lread as usize,
        std::cmp::min(reads_in.len(), reads_out.len()),
    );
    for jj in 0..lread {
        reads_out[jj] = match reads_in[lread - 1 - jj] {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'S' => b'S',
            b'W' => b'W',
            b'B' => b'V',
            b'D' => b'H',
            b'V' => b'B',
            b'H' => b'D',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            b'n' => b'n',
            b'r' => b'y',
            b'y' => b'r',
            b'k' => b'm',
            b'm' => b'k',
            b's' => b's',
            b'w' => b'w',
            b'b' => b'v',
            b'd' => b'h',
            b'v' => b'b',
            b'h' => b'd',
            value => value,
        };
    }
}

#[doc = "Original `revComplementNucleotides` at STAR/source/SequenceFuns.cpp:56. Args: seq: string"]
pub fn sequencefuns_l56_revcomplementnucleotides(seq: &mut String) {
    let seq1 = seq.as_bytes().to_vec();
    let mut out = vec![0; seq1.len()];
    for jj in 0..seq1.len() {
        out[jj] = match seq1[seq1.len() - 1 - jj] {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'S' => b'S',
            b'W' => b'W',
            b'B' => b'V',
            b'D' => b'H',
            b'V' => b'B',
            b'H' => b'D',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            b'n' => b'n',
            b'r' => b'y',
            b'y' => b'r',
            b'k' => b'm',
            b'm' => b'k',
            b's' => b's',
            b'w' => b'w',
            b'b' => b'v',
            b'd' => b'h',
            b'v' => b'b',
            b'h' => b'd',
            value => value,
        };
    }
    *seq = String::from_utf8(out).unwrap();
}

#[doc = "Original `nuclToNumBAM` at STAR/source/SequenceFuns.cpp:99. Args: cc: char"]
pub fn sequencefuns_l99_nucltonumbam(cc: u8) -> u8 {
    match cc {
        b'=' => 0,
        b'A' | b'a' => 1,
        b'C' | b'c' => 2,
        b'M' | b'm' => 3,
        b'G' | b'g' => 4,
        b'R' | b'r' => 5,
        b'S' | b's' => 6,
        b'V' | b'v' => 7,
        b'T' | b't' => 8,
        b'W' | b'w' => 9,
        b'Y' | b'y' => 10,
        b'H' | b'h' => 11,
        b'K' | b'k' => 12,
        b'D' | b'd' => 13,
        b'B' | b'b' => 14,
        b'N' | b'n' => 15,
        _ => 15,
    }
}

#[doc = "Original `nuclPackBAM` at STAR/source/SequenceFuns.cpp:122. Args: ReadsIn: char, ReadsOut: char, Lread: uint"]
pub fn sequencefuns_l122_nuclpackbam(reads_in: &[u8], reads_out: &mut [u8], lread: u32) {
    let lread = std::cmp::min(lread as usize, reads_in.len());
    let packed_len = lread.div_ceil(2);
    let lread = std::cmp::min(lread, reads_out.len().saturating_mul(2));
    for jj in 0..std::cmp::min(lread / 2, packed_len) {
        reads_out[jj] = sequencefuns_l99_nucltonumbam(reads_in[2 * jj]) << 4
            | sequencefuns_l99_nucltonumbam(reads_in[2 * jj + 1]);
    }
    if lread % 2 == 1 && lread / 2 < reads_out.len() {
        reads_out[lread / 2] = sequencefuns_l99_nucltonumbam(reads_in[lread - 1]) << 4;
    }
}

#[doc = "Original `convertNucleotidesToNumbers` at STAR/source/SequenceFuns.cpp:131. Args: R0: char, R1: char, Lread: uint"]
pub fn sequencefuns_l131_convertnucleotidestonumbers(r0: &[u8], r1: &mut [u8], lread: u32) {
    let lread = std::cmp::min(lread as usize, std::cmp::min(r0.len(), r1.len()));
    for jj in 0..lread {
        r1[jj] = match r0[jj] {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4,
        };
    }
}

#[doc = "Original `convertCapitalBasesToNum` at STAR/source/SequenceFuns.cpp:148. Args: rS: uint8_t, N: uint64_t"]
pub fn sequencefuns_l148_convertcapitalbasestonum(rs: &mut [u8], n: u64) {
    let n = std::cmp::min(n as usize, rs.len());
    for ib in 0..n {
        rs[ib] = match rs[ib] {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        };
    }
}

#[doc = "Original `convertNucleotidesToNumbersRemoveControls` at STAR/source/SequenceFuns.cpp:170. Args: R0: char, R1: char, Lread: uint"]
pub fn sequencefuns_l170_convertnucleotidestonumbersremovecontrols(
    r0: &[u8],
    r1: &mut [u8],
    lread: u32,
) -> u32 {
    let mut ir1 = 0;
    let lread = std::cmp::min(lread as usize, std::cmp::min(r0.len(), r1.len()));
    for jj in 0..lread {
        match r0[jj] {
            b'A' | b'a' => r1[jj] = 0,
            b'C' | b'c' => r1[jj] = 1,
            b'G' | b'g' => r1[jj] = 2,
            b'T' | b't' => r1[jj] = 3,
            value => {
                if value < 32 {
                    continue;
                } else {
                    r1[jj] = 4;
                }
            }
        }
        ir1 += 1;
    }
    ir1
}

#[doc = "Original `convertNt01234` at STAR/source/SequenceFuns.cpp:195. Args: R0: char"]
pub fn sequencefuns_l195_convertnt01234(r0: u8) -> u8 {
    match r0 {
        b'a' | b'A' => 0,
        b'c' | b'C' => 1,
        b'g' | b'G' => 2,
        b't' | b'T' => 3,
        _ => 4,
    }
}

#[doc = "Original `convertNuclStrToInt32` at STAR/source/SequenceFuns.cpp:219. Args: S: string, intOut: uint32"]
pub fn sequencefuns_l219_convertnuclstrtoint32(s: &str, int_out: &mut u32) -> i32 {
    *int_out = 0;
    let mut pos_n = -1;
    for (ii, value) in s.bytes().enumerate() {
        let mut nt = sequencefuns_l195_convertnt01234(value) as u32;
        if nt > 3 {
            if pos_n >= 0 {
                return -2;
            }
            pos_n = ii as i32;
            nt = 0;
        }
        *int_out = *int_out << 2;
        *int_out += nt;
    }
    pos_n
}

#[doc = "Original `convertNuclInt32toString` at STAR/source/SequenceFuns.cpp:237. Args: nuclNum: uint32, L: uint32"]
pub fn sequencefuns_l237_convertnuclint32tostring(mut nucl_num: u32, l: u32) -> String {
    let mut nucl_out = vec![b'N'; l as usize];
    let nucl_char = b"ACGT";
    for ii in 1..=l as usize {
        nucl_out[l as usize - ii] = nucl_char[(nucl_num & 3) as usize];
        nucl_num >>= 2;
    }
    String::from_utf8(nucl_out).unwrap()
}

#[doc = "Original `convertNuclStrToInt64` at STAR/source/SequenceFuns.cpp:249. Args: S: string, intOut: uint64"]
pub fn sequencefuns_l249_convertnuclstrtoint64(s: &str, int_out: &mut u64) -> i64 {
    *int_out = 0;
    let mut pos_n = -1;
    for (ii, value) in s.bytes().enumerate() {
        let mut nt = sequencefuns_l195_convertnt01234(value) as u64;
        if nt > 3 {
            if pos_n >= 0 {
                return -2;
            }
            pos_n = ii as i64;
            nt = 0;
        }
        *int_out = *int_out << 2;
        *int_out += nt;
    }
    pos_n
}

#[doc = "Original `convertNuclInt64toString` at STAR/source/SequenceFuns.cpp:267. Args: nuclNum: uint64, L: uint32"]
pub fn sequencefuns_l267_convertnuclint64tostring(mut nucl_num: u64, l: u32) -> String {
    let mut nucl_out = vec![b'N'; l as usize];
    let nucl_char = b"ACGT";
    for ii in 1..=l as usize {
        nucl_out[l as usize - ii] = nucl_char[(nucl_num & 3) as usize];
        nucl_num >>= 2;
    }
    String::from_utf8(nucl_out).unwrap()
}

#[doc = "Original `chrFind` at STAR/source/SequenceFuns.cpp:280. Args: Start: uint, i2: uint, chrStart: uint"]
pub fn sequencefuns_l280_chrfind(start: u32, mut i2: u32, chr_start: &[u32]) -> u32 {
    let mut i1 = 0;
    while i1 + 1 < i2 {
        let i3 = (i1 + i2) / 2;
        if chr_start[i3 as usize] > start {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    i1
}

#[doc = "Original `localSearch` at STAR/source/SequenceFuns.cpp:293. Args: x: char, nx: uint, y: char, ny: uint, pMM: double"]
pub fn sequencefuns_l293_localsearch(x: &[u8], nx: u32, y: &[u8], ny: u32, p_mm: f64) -> u32 {
    let mut n_match_best = 0;
    let mut n_mm_best = 0;
    let nx = std::cmp::min(nx as usize, x.len());
    let ny = std::cmp::min(ny as usize, y.len());
    let mut ix_best = nx;
    for ix in 0..nx {
        let mut n_match = 0;
        let mut n_mm = 0;
        for iy in 0..ny.min(nx - ix) {
            if x[ix + iy] > 3 {
                continue;
            }
            if x[ix + iy] == y[iy] {
                n_match += 1;
            } else {
                n_mm += 1;
            }
        }
        if (n_match > n_match_best || (n_match == n_match_best && n_mm < n_mm_best))
            && (n_mm as f64) / (n_match as f64) <= p_mm
        {
            ix_best = ix;
            n_match_best = n_match;
            n_mm_best = n_mm;
        }
    }
    ix_best as u32
}

#[doc = "Original `localSearchNisMM` at STAR/source/SequenceFuns.cpp:317. Args: x: char, nx: uint, y: char, ny: uint, pMM: double"]
pub fn sequencefuns_l317_localsearchnismm(x: &[u8], nx: u32, y: &[u8], ny: u32, p_mm: f64) -> u32 {
    let mut n_match_best = 0;
    let mut n_mm_best = 0;
    let nx = std::cmp::min(nx as usize, x.len());
    let ny = std::cmp::min(ny as usize, y.len());
    let mut ix_best = nx;
    for ix in 0..nx {
        let mut n_match = 0;
        let mut n_mm = 0;
        for iy in 0..ny.min(nx - ix) {
            if x[ix + iy] == y[iy] && y[iy] < 4 {
                n_match += 1;
            } else {
                n_mm += 1;
            }
        }
        if (n_match > n_match_best || (n_match == n_match_best && n_mm < n_mm_best))
            && (n_mm as f64) / (n_match as f64) <= p_mm
        {
            ix_best = ix;
            n_match_best = n_match;
            n_mm_best = n_mm;
        }
    }
    ix_best as u32
}

#[doc = "Original `localAlignHammingDist` at STAR/source/SequenceFuns.cpp:341. Args: text: string, query: string, pos: uint32"]
pub fn sequencefuns_l341_localalignhammingdist(text: &str, query: &str, pos: &mut u32) -> u32 {
    let text_b = text.as_bytes();
    let query_b = query.as_bytes();
    let mut dist_best = query_b.len() as u32;
    if text_b.len() < query_b.len() {
        return text_b.len() as u32 + 1;
    }
    for ii in 0..=text_b.len() - query_b.len() {
        let mut dist1 = 0;
        for jj in 0..query_b.len() {
            if query_b[jj] != b'N' && text_b[jj + ii] != query_b[jj] {
                dist1 += 1;
            }
        }
        if dist1 < dist_best {
            dist_best = dist1;
            *pos = ii as u32;
        }
    }
    dist_best
}

#[doc = "Original `qualitySplit` at STAR/source/SequenceFuns.cpp:411. Args: r: char, L: uint, maxNsplit: uint, minLsplit: uint, splitR: uint"]
pub fn sequencefuns_l411_qualitysplit(
    r: &[u8],
    l: u32,
    max_nsplit: u32,
    min_lsplit: u32,
    split_r: &mut [Vec<u32>; 3],
) -> u32 {
    let mut ir = 0;
    let mut is = 0;
    let mut lgood_min = 0;
    let mut ifrag = 0;
    let l = std::cmp::min(l, r.len() as u32);
    let max_nsplit = std::cmp::min(
        max_nsplit,
        split_r.iter().map(|values| values.len()).min().unwrap_or(0) as u32,
    );
    while (ir < l) & (is < max_nsplit) {
        while ir < l && r[ir as usize] > 3 {
            if r[ir as usize] == MARK_FRAG_SPACER_BASE {
                ifrag += 1;
            }
            ir += 1;
        }
        if ir == l {
            break;
        }
        let ir1 = ir;
        while ir < l && r[ir as usize] <= 3 {
            ir += 1;
        }
        if ir - ir1 > lgood_min {
            lgood_min = ir - ir1;
        }
        if ir - ir1 < min_lsplit {
            continue;
        }
        split_r[0][is as usize] = ir1;
        split_r[1][is as usize] = ir - ir1;
        split_r[2][is as usize] = ifrag;
        is += 1;
    }
    if is == 0 {
        split_r[1][0] = lgood_min;
    }
    is
}
