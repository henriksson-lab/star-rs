#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `funCompareNames` at STAR/source/bamRemoveDuplicates.cpp:13. Args: a: void, b: void"]
pub fn bamremoveduplicates_l13_funcomparenames(a: &[u32], b: &[u32]) -> i32 {
    let la = (a[3] << 24) >> 24;
    let lb = (b[3] << 24) >> 24;
    if la > lb {
        return 1;
    } else if la < lb {
        return -1;
    } else {
        let ca = &a[9..];
        let cb = &b[9..];
        for ii in 0..la as usize {
            let byte_a = ca[ii / 4].to_ne_bytes()[ii % 4];
            let byte_b = cb[ii / 4].to_ne_bytes()[ii % 4];
            if byte_a > byte_b {
                return 1;
            } else if byte_a < byte_b {
                return -1;
            }
        }
        let fa = a[4] >> 16;
        let fb = b[4] >> 16;
        if (fa & 0x80) > (fb & 0x80) {
            return 1;
        } else if (fa & 0x80) < (fb & 0x80) {
            return -1;
        }
        0
    }
}

#[doc = "Original `funStartExtendS` at STAR/source/bamRemoveDuplicates.cpp:34. Args: p: uint32"]
pub fn bamremoveduplicates_l34_funstartextends(p: &[u32]) -> u32 {
    let name_len = ((p[3] << 24) >> 24) as usize;
    let cig_byte = 9 * 4 + name_len;
    let cig0 =
        unsafe { std::ptr::read_unaligned((p.as_ptr() as *const u8).add(cig_byte) as *const u32) };
    if ((cig0 << 28) >> 28) == 4 {
        p[2] - (cig0 >> 4)
    } else {
        p[2]
    }
}

#[doc = "Original `funCigarExtendS` at STAR/source/bamRemoveDuplicates.cpp:43. Args: p: uint32, cout: uint32"]
pub fn bamremoveduplicates_l43_funcigarextends(p: &[u32], cout: &mut [u32]) -> u32 {
    let name_len = ((p[3] << 24) >> 24) as usize;
    let cig_byte = 9 * 4 + name_len;
    let n = (p[4] << 16) >> 16;
    let mut n1 = n;
    let cig0 =
        unsafe { std::ptr::read_unaligned((p.as_ptr() as *const u8).add(cig_byte) as *const u32) };
    if ((cig0 << 28) >> 28) == 4 {
        n1 -= 1;
        for i in 0..n1 as usize {
            cout[i] = unsafe {
                std::ptr::read_unaligned(
                    (p.as_ptr() as *const u8).add(cig_byte + (i + 1) * 4) as *const u32
                )
            };
        }
        cout[0] += (cig0 >> 4) << 4;
    } else {
        for i in 0..n as usize {
            cout[i] = unsafe {
                std::ptr::read_unaligned(
                    (p.as_ptr() as *const u8).add(cig_byte + i * 4) as *const u32
                )
            };
        }
    }
    let cign_1 = unsafe {
        std::ptr::read_unaligned(
            (p.as_ptr() as *const u8).add(cig_byte + (n as usize - 1) * 4) as *const u32,
        )
    };
    if ((cign_1 << 28) >> 28) == 4 {
        n1 -= 1;
        cout[n1 as usize - 1] += (cign_1 >> 4) << 4;
    }
    n1
}

#[doc = "Original `funCompareCigarsExtendS` at STAR/source/bamRemoveDuplicates.cpp:61. Args: pa: uint32, pb: uint32"]
pub fn bamremoveduplicates_l61_funcomparecigarsextends(pa: &[u32], pb: &[u32]) -> i32 {
    let mut ca = [0u32; 100];
    let mut cb = [0u32; 100];
    let na = bamremoveduplicates_l43_funcigarextends(pa, &mut ca);
    let nb = bamremoveduplicates_l43_funcigarextends(pb, &mut cb);
    if na > nb {
        return 1;
    } else if na < nb {
        return -1;
    }
    for ii in 0..na as usize {
        if ca[ii] > cb[ii] {
            return 1;
        } else if ca[ii] < cb[ii] {
            return -1;
        }
    }
    0
}

#[doc = "Original `funCompareCoordFlagCigarSeq` at STAR/source/bamRemoveDuplicates.cpp:72. Args: a: void, b: void"]
pub fn bamremoveduplicates_l72_funcomparecoordflagcigarseq(
    pa1: &[u32],
    pa2: &[u32],
    pb1: &[u32],
    pb2: &[u32],
    mate2_bases_n: u32,
) -> i32 {
    let start_a1 = bamremoveduplicates_l34_funstartextends(pa1);
    let start_b1 = bamremoveduplicates_l34_funstartextends(pb1);
    if start_a1 > start_b1 {
        return 1;
    } else if start_a1 < start_b1 {
        return -1;
    }

    let start_a2 = bamremoveduplicates_l34_funstartextends(pa2);
    let start_b2 = bamremoveduplicates_l34_funstartextends(pb2);
    if start_a2 > start_b2 {
        return 1;
    } else if start_a2 < start_b2 {
        return -1;
    }

    let flag_a1 = pa1[4] >> 16;
    let flag_b1 = pb1[4] >> 16;
    if flag_a1 > flag_b1 {
        return 1;
    } else if flag_a1 < flag_b1 {
        return -1;
    }

    let flag_a2 = pa2[4] >> 16;
    let flag_b2 = pb2[4] >> 16;
    if flag_a2 > flag_b2 {
        return 1;
    } else if flag_a2 < flag_b2 {
        return -1;
    }

    let mut ret1 = bamremoveduplicates_l61_funcomparecigarsextends(pa1, pb1);
    if ret1 != 0 {
        return ret1;
    }
    ret1 = bamremoveduplicates_l61_funcomparecigarsextends(pa2, pb2);
    if ret1 != 0 {
        return ret1;
    }

    let pa2_name_len = ((pa2[3] << 24) >> 24) as usize;
    let pb2_name_len = ((pb2[3] << 24) >> 24) as usize;
    let pa2_n_cigar = ((pa2[4] << 16) >> 16) as usize;
    let pb2_n_cigar = ((pb2[4] << 16) >> 16) as usize;
    let sa_byte = 9 * 4 + pa2_name_len + pa2_n_cigar * 4;
    let sb_byte = 9 * 4 + pb2_name_len + pb2_n_cigar * 4;
    let sa = pa2.as_ptr() as *const u8;
    let sb = pb2.as_ptr() as *const u8;

    if ((pa2[4] >> 16) & 0x10) == 0 {
        let mut ii = 1;
        while ii < mate2_bases_n {
            let a_byte = unsafe { *sa.add(sa_byte + (ii / 2) as usize) };
            let b_byte = unsafe { *sb.add(sb_byte + (ii / 2) as usize) };
            if a_byte > b_byte {
                return 1;
            } else if a_byte < b_byte {
                return -1;
            }
            ii += 2;
        }
        if mate2_bases_n % 2 > 0 {
            let a_byte = unsafe { *sa.add(sa_byte + (ii / 2) as usize) } >> 4;
            let b_byte = unsafe { *sb.add(sb_byte + (ii / 2) as usize) } >> 4;
            if a_byte > b_byte {
                return 1;
            } else if a_byte < b_byte {
                return -1;
            }
        }
    } else {
        let mut ii = pa2[5] - mate2_bases_n;
        if ii % 2 > 0 {
            let a_byte = unsafe { *sa.add(sa_byte + (ii / 2) as usize) } & 15;
            let b_byte = unsafe { *sb.add(sb_byte + (ii / 2) as usize) } & 15;
            if a_byte > b_byte {
                return 1;
            } else if a_byte < b_byte {
                return -1;
            }
            ii += 1;
        }
        while ii < pa2[5] {
            let a_byte = unsafe { *sa.add(sa_byte + (ii / 2) as usize) };
            let b_byte = unsafe { *sb.add(sb_byte + (ii / 2) as usize) };
            if a_byte > b_byte {
                return 1;
            } else if a_byte < b_byte {
                return -1;
            }
            ii += 2;
        }
    }

    0
}

#[doc = "Original `bamRemoveDuplicates` at STAR/source/bamRemoveDuplicates.cpp:114. Args: bamFileName: string, bamFileNameOut: string, P: Parameters"]
pub fn bamremoveduplicates_l114_bamremoveduplicates(
    bam_records: &mut [Vec<u32>],
    mate2_bases_n: u32,
    mark_multi: bool,
) -> Result<(), String> {
    let validate_record = |record: &[u32]| -> Result<(), String> {
        if record.len() < 9 {
            return Err(
                "EXITING because of fatal ERROR: malformed BAM record for deduplication"
                    .to_string(),
            );
        }
        let name_len = ((record[3] << 24) >> 24) as usize;
        let n_cigar = ((record[4] << 16) >> 16) as usize;
        let seq_len = record[5] as usize;
        if n_cigar == 0 || n_cigar > 100 {
            return Err(
                "EXITING because of fatal ERROR: malformed BAM CIGAR for deduplication".to_string(),
            );
        }
        if seq_len < mate2_bases_n as usize {
            return Err(
                "EXITING because of fatal ERROR: BAM read sequence is shorter than bamRemoveDuplicatesMate2basesN"
                    .to_string(),
            );
        }
        let record_bytes = record.len().checked_mul(4).ok_or_else(|| {
            "EXITING because of fatal ERROR: malformed BAM record for deduplication".to_string()
        })?;
        let aux_start = 9usize
            .checked_mul(4)
            .and_then(|v| v.checked_add(name_len))
            .and_then(|v| v.checked_add(n_cigar.checked_mul(4)?))
            .and_then(|v| v.checked_add(seq_len.div_ceil(2)))
            .and_then(|v| v.checked_add(seq_len))
            .ok_or_else(|| {
                "EXITING because of fatal ERROR: malformed BAM record for deduplication".to_string()
            })?;
        if aux_start > record_bytes {
            return Err(
                "EXITING because of fatal ERROR: malformed BAM record for deduplication"
                    .to_string(),
            );
        }
        Ok(())
    };

    for record in bam_records.iter() {
        validate_record(record)?;
    }

    let record_bytes = |record: &[u32]| -> Vec<u8> {
        let mut bytes = Vec::with_capacity(record.len() * 4);
        for word in record {
            bytes.extend_from_slice(&word.to_ne_bytes());
        }
        bytes
    };

    let aux_start = |record: &[u32]| -> usize {
        let name_len = ((record[3] << 24) >> 24) as usize;
        let n_cigar = ((record[4] << 16) >> 16) as usize;
        let seq_len = record[5] as usize;
        9 * 4 + name_len + n_cigar * 4 + seq_len.div_ceil(2) + seq_len
    };

    let bam_aux_i = |record: &[u32], tag: &[u8; 2]| -> Option<i32> {
        let bytes = record_bytes(record);
        let mut pos = aux_start(record);
        while pos + 3 <= bytes.len() {
            if bytes[pos] == 0 && bytes[pos + 1] == 0 {
                break;
            }
            let matches = bytes[pos] == tag[0] && bytes[pos + 1] == tag[1];
            let type_ = bytes[pos + 2];
            pos += 3;
            match type_ {
                b'A' | b'c' => {
                    if pos + 1 > bytes.len() {
                        return None;
                    }
                    let value = bytes[pos] as i8 as i32;
                    if matches {
                        return Some(value);
                    }
                    pos += 1;
                }
                b'C' => {
                    if pos + 1 > bytes.len() {
                        return None;
                    }
                    let value = bytes[pos] as i32;
                    if matches {
                        return Some(value);
                    }
                    pos += 1;
                }
                b's' => {
                    if pos + 2 > bytes.len() {
                        return None;
                    }
                    let value = i16::from_ne_bytes([bytes[pos], bytes[pos + 1]]) as i32;
                    if matches {
                        return Some(value);
                    }
                    pos += 2;
                }
                b'S' => {
                    if pos + 2 > bytes.len() {
                        return None;
                    }
                    let value = u16::from_ne_bytes([bytes[pos], bytes[pos + 1]]) as i32;
                    if matches {
                        return Some(value);
                    }
                    pos += 2;
                }
                b'i' => {
                    if pos + 4 > bytes.len() {
                        return None;
                    }
                    let value = i32::from_ne_bytes([
                        bytes[pos],
                        bytes[pos + 1],
                        bytes[pos + 2],
                        bytes[pos + 3],
                    ]);
                    if matches {
                        return Some(value);
                    }
                    pos += 4;
                }
                b'I' => {
                    if pos + 4 > bytes.len() {
                        return None;
                    }
                    let value = u32::from_ne_bytes([
                        bytes[pos],
                        bytes[pos + 1],
                        bytes[pos + 2],
                        bytes[pos + 3],
                    ]) as i32;
                    if matches {
                        return Some(value);
                    }
                    pos += 4;
                }
                b'f' => {
                    if pos + 4 > bytes.len() {
                        return None;
                    }
                    if matches {
                        return Some(f32::from_ne_bytes([
                            bytes[pos],
                            bytes[pos + 1],
                            bytes[pos + 2],
                            bytes[pos + 3],
                        ]) as i32);
                    }
                    pos += 4;
                }
                b'Z' | b'H' => {
                    while pos < bytes.len() && bytes[pos] != 0 {
                        pos += 1;
                    }
                    if pos >= bytes.len() {
                        return None;
                    }
                    pos += 1;
                    if matches {
                        return None;
                    }
                }
                _ => return None,
            }
        }
        None
    };

    const DUP_FLAG_WORD: u32 = 0x400 << 16;
    let mut unique = Vec::new();
    for (idx, record) in bam_records.iter_mut().enumerate() {
        let n_mult = bam_aux_i(record, b"NH").ok_or_else(|| {
            "EXITING because of fatal ERROR: SAM tag NH is missing from a read, but it's required for deduplication. \nSOLUTION: re-generate BAM file with NH and AS tags.".to_string()
        })?;
        if n_mult == 1 || (n_mult > 1 && mark_multi) {
            record[4] |= DUP_FLAG_WORD;
        }
        if n_mult == 1 {
            unique.push(idx);
        }
    }

    unique.sort_by(|&a, &b| {
        match bamremoveduplicates_l13_funcomparenames(&bam_records[a], &bam_records[b]) {
            x if x < 0 => std::cmp::Ordering::Less,
            x if x > 0 => std::cmp::Ordering::Greater,
            _ => std::cmp::Ordering::Equal,
        }
    });

    let pair_n = unique.len() / 2;
    let mut pairs: Vec<[usize; 2]> = (0..pair_n)
        .map(|ii| [unique[ii * 2], unique[ii * 2 + 1]])
        .collect();
    pairs.sort_by(|a, b| {
        match bamremoveduplicates_l72_funcomparecoordflagcigarseq(
            &bam_records[a[0]],
            &bam_records[a[1]],
            &bam_records[b[0]],
            &bam_records[b[1]],
            mate2_bases_n,
        ) {
            x if x < 0 => std::cmp::Ordering::Less,
            x if x > 0 => std::cmp::Ordering::Greater,
            _ => std::cmp::Ordering::Equal,
        }
    });

    let mut best_score = -999;
    let mut best_pair = 0usize;
    for pp in 0..pairs.len() {
        let score = bam_aux_i(&bam_records[pairs[pp][0]], b"AS").ok_or_else(|| {
            "EXITING because of fatal ERROR: SAM tag AS is missing from a read, but it's required for deduplication. \nSOLUTION: re-generate BAM file with NH and AS tags.".to_string()
        })?;
        if score > best_score {
            best_score = score;
            best_pair = pp;
        }
        let last_or_next_different = pp == pairs.len() - 1
            || bamremoveduplicates_l72_funcomparecoordflagcigarseq(
                &bam_records[pairs[pp][0]],
                &bam_records[pairs[pp][1]],
                &bam_records[pairs[pp + 1][0]],
                &bam_records[pairs[pp + 1][1]],
                mate2_bases_n,
            ) != 0;
        if last_or_next_different {
            bam_records[pairs[best_pair][1]][4] ^= DUP_FLAG_WORD;
            bam_records[pairs[best_pair][0]][4] ^= DUP_FLAG_WORD;
            best_score = -999;
        }
    }

    Ok(())
}
