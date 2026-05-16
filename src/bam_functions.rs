#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `bam_cigarString` at STAR/source/BAMfunctions.cpp:5. Args: b: bam1_t"]
pub fn bamfunctions_l5_bam_cigarstring(cigar: &[u32], n_cigar: i32) -> String {
    let mut cigar_string = String::new();
    if n_cigar > 0 {
        for cigar_op in cigar.iter().take(n_cigar as usize) {
            cigar_string.push_str(&(cigar_op >> 4).to_string());
            cigar_string.push(match cigar_op & 0xf {
                0 => 'M',
                1 => 'I',
                2 => 'D',
                3 => 'N',
                4 => 'S',
                5 => 'H',
                6 => 'P',
                7 => '=',
                8 => 'X',
                9 => 'B',
                _ => '?',
            });
        }
    }
    cigar_string
}

#[doc = "Original `bam_read1_fromArray` at STAR/source/BAMfunctions.cpp:30. Args: bamChar: char, b: bam1_t"]
pub fn bamfunctions_l30_bam_read1_fromarray(
    bam_char: &[u8],
    b: &mut crate::bam_output::Bam1,
) -> i32 {
    if bam_char.len() < 36 {
        return -4;
    }
    let block_len = i32::from_ne_bytes(bam_char[0..4].try_into().unwrap());
    if block_len < 32 {
        return -4;
    }
    let record_len = match (block_len as usize).checked_add(4) {
        Some(len) => len,
        None => return -4,
    };
    if record_len > bam_char.len() {
        return -4;
    }
    let mut x = [0u32; 8];
    for ii in 0..8 {
        let start = 4 + ii * 4;
        x[ii] = u32::from_ne_bytes(bam_char[start..start + 4].try_into().unwrap());
    }

    b.core.tid = x[0] as i32;
    b.core.pos = x[1] as i32;
    b.core.bin = x[2] >> 16;
    b.core.qual = (x[2] >> 8) & 0xff;
    b.core.l_qname = x[2] & 0xff;
    b.core.flag = x[3] >> 16;
    b.core.n_cigar = x[3] & 0xffff;
    b.core.l_qseq = x[4] as i32;
    b.core.mtid = x[5] as i32;
    b.core.mpos = x[6] as i32;
    b.core.isize = x[7] as i32;
    b.l_data = block_len - 32;
    if b.l_data < 0 || b.core.l_qseq < 0 {
        return -4;
    }
    if 36usize
        .checked_add(b.l_data as usize)
        .is_none_or(|len| len > bam_char.len())
    {
        return -4;
    }

    if b.m_data < b.l_data {
        b.m_data = b.l_data;
        if b.m_data > 0 {
            b.m_data = (b.m_data as u32).next_power_of_two() as i32;
        }
    }
    b.data_offset = 4 * 9;

    4 + block_len
}

#[doc = "Original `outBAMwriteHeader` at STAR/source/BAMfunctions.cpp:77. Args: fp: BGZF, samh: string, chrn: vector <string>, chrl: vector <uint>"]
pub fn bamfunctions_l77_outbamwriteheader(samh: &str, chrn: &[String], chrl: &[u32]) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"BAM\x01");
    let hlen = samh.len() as i32;
    out.extend_from_slice(&hlen.to_ne_bytes());
    out.extend_from_slice(samh.as_bytes());
    let nchr = chrn.len() as i32;
    out.extend_from_slice(&nchr.to_ne_bytes());
    for ii in 0..chrn.len() {
        let rlen = chrn[ii].len() as i32 + 1;
        let slen = chrl[ii] as i32;
        out.extend_from_slice(&rlen.to_ne_bytes());
        out.extend_from_slice(chrn[ii].as_bytes());
        out.push(0);
        out.extend_from_slice(&slen.to_ne_bytes());
    }
    out
}

#[doc = "Original `reg2bin` at STAR/source/BAMfunctions.cpp:95. Args: beg: int, end: int"]
pub fn bamfunctions_l95_reg2bin(beg: i32, mut end: i32) -> i32 {
    end -= 1;
    if beg >> 14 == end >> 14 {
        return ((1 << 15) - 1) / 7 + (beg >> 14);
    }
    if beg >> 17 == end >> 17 {
        return ((1 << 12) - 1) / 7 + (beg >> 17);
    }
    if beg >> 20 == end >> 20 {
        return ((1 << 9) - 1) / 7 + (beg >> 20);
    }
    if beg >> 23 == end >> 23 {
        return ((1 << 6) - 1) / 7 + (beg >> 23);
    }
    if beg >> 26 == end >> 26 {
        return 1 + (beg >> 26);
    }
    0
}

#[doc = "Original `bamAttrArrayWrite` at STAR/source/BAMfunctions.cpp:106. Args: attr: int32, tagName: char, attrArray: char"]
pub fn bamfunctions_l106_bamattrarraywrite(
    attr: i32,
    tag_name: &[u8],
    attr_array: &mut [u8],
) -> i32 {
    attr_array[0] = tag_name[0];
    attr_array[1] = tag_name[1];
    attr_array[2] = b'i';
    attr_array[3..7].copy_from_slice(&attr.to_ne_bytes());
    7
}

#[doc = "Original `bamAttrArrayWrite` at STAR/source/BAMfunctions.cpp:112. Args: attr: float, tagName: char, attrArray: char"]
pub fn bamfunctions_l112_bamattrarraywrite(
    attr: f32,
    tag_name: &[u8],
    attr_array: &mut [u8],
) -> i32 {
    attr_array[0] = tag_name[0];
    attr_array[1] = tag_name[1];
    attr_array[2] = b'f';
    attr_array[3..7].copy_from_slice(&attr.to_ne_bytes());
    7
}

#[doc = "Original `bamAttrArrayWrite` at STAR/source/BAMfunctions.cpp:118. Args: attr: char, tagName: char, attrArray: char"]
pub fn bamfunctions_l118_bamattrarraywrite(
    attr: u8,
    tag_name: &[u8],
    attr_array: &mut [u8],
) -> i32 {
    attr_array[0] = tag_name[0];
    attr_array[1] = tag_name[1];
    attr_array[2] = b'A';
    attr_array[3] = attr;
    4
}

#[doc = "Original `bamAttrArrayWrite` at STAR/source/BAMfunctions.cpp:124. Args: attr: string, tagName: char, attrArray: char"]
pub fn bamfunctions_l124_bamattrarraywrite(
    attr: &str,
    tag_name: &[u8],
    attr_array: &mut [u8],
) -> i32 {
    attr_array[0] = tag_name[0];
    attr_array[1] = tag_name[1];
    attr_array[2] = b'Z';
    attr_array[3..3 + attr.len()].copy_from_slice(attr.as_bytes());
    attr_array[3 + attr.len()] = 0;
    (3 + attr.len() + 1) as i32
}

#[doc = "Original `bamAttrArrayWrite` at STAR/source/BAMfunctions.cpp:130. Args: attr: vector<char>, tagName: char, attrArray: char"]
pub fn bamfunctions_l130_bamattrarraywrite(
    attr: &[u8],
    tag_name: &[u8],
    attr_array: &mut [u8],
) -> i32 {
    attr_array[0] = tag_name[0];
    attr_array[1] = tag_name[1];
    attr_array[2] = b'B';
    attr_array[3] = b'c';
    attr_array[4..8].copy_from_slice(&(attr.len() as i32).to_ne_bytes());
    attr_array[8..8 + attr.len()].copy_from_slice(attr);
    (8 + attr.len()) as i32
}

#[doc = "Original `bamAttrArrayWrite` at STAR/source/BAMfunctions.cpp:138. Args: attr: vector<int32>, tagName: char, attrArray: char"]
pub fn bamfunctions_l138_bamattrarraywrite(
    attr: &[i32],
    tag_name: &[u8],
    attr_array: &mut [u8],
) -> i32 {
    attr_array[0] = tag_name[0];
    attr_array[1] = tag_name[1];
    attr_array[2] = b'B';
    attr_array[3] = b'i';
    attr_array[4..8].copy_from_slice(&(attr.len() as i32).to_ne_bytes());
    for (i, value) in attr.iter().enumerate() {
        attr_array[8 + i * 4..12 + i * 4].copy_from_slice(&value.to_ne_bytes());
    }
    (8 + 4 * attr.len()) as i32
}

#[doc = "Original `bamAttrArrayWriteSAMtags` at STAR/source/BAMfunctions.cpp:147. Args: attrStr: string, attrArray: char, P: Parameters"]
pub fn bamfunctions_l147_bamattrarraywritesamtags(
    attr_str: &str,
    attr_array: &mut [u8],
    sam_attr_keep_all: bool,
    sam_attr_keep: &[u16],
) -> Result<i32, String> {
    let mut nattr = 0usize;
    for attr1 in attr_str.split('\t') {
        if attr1.is_empty() {
            continue;
        }

        let attr_bytes = attr1.as_bytes();
        if attr_bytes.len() < 5 || attr_bytes[2] != b':' || attr_bytes[4] != b':' {
            return Err(format!(
                "EXITING because of FATAL ERROR: malformed SAM attribute: {}\n",
                attr1
            ));
        }
        let tagn = u16::from_ne_bytes([attr_bytes[0], attr_bytes[1]]);
        if !sam_attr_keep_all && !sam_attr_keep.contains(&tagn) {
            continue;
        }

        match attr_bytes[3] {
            b'i' => {
                let a1 = attr1[5..].parse::<i32>().map_err(|_| {
                    format!(
                        "EXITING because of FATAL ERROR: malformed SAM integer attribute: {}\n",
                        attr1
                    )
                })?;
                nattr +=
                    bamfunctions_l106_bamattrarraywrite(a1, attr_bytes, &mut attr_array[nattr..])
                        as usize;
            }
            b'A' => {
                if attr_bytes.len() < 6 {
                    return Err(format!(
                        "EXITING because of FATAL ERROR: malformed SAM character attribute: {}\n",
                        attr1
                    ));
                }
                let a1 = attr_bytes[5];
                nattr +=
                    bamfunctions_l118_bamattrarraywrite(a1, attr_bytes, &mut attr_array[nattr..])
                        as usize;
            }
            b'Z' => {
                nattr += bamfunctions_l124_bamattrarraywrite(
                    &attr1[5..],
                    attr_bytes,
                    &mut attr_array[nattr..],
                ) as usize;
            }
            b'f' => {
                let a1 = attr1[5..].parse::<f32>().map_err(|_| {
                    format!(
                        "EXITING because of FATAL ERROR: malformed SAM float attribute: {}\n",
                        attr1
                    )
                })?;
                nattr +=
                    bamfunctions_l112_bamattrarraywrite(a1, attr_bytes, &mut attr_array[nattr..])
                        as usize;
            }
            _ => {}
        }
    }

    Ok(nattr as i32)
}
