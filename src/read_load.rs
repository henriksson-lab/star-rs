#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `readLoad` at STAR/source/readLoad.cpp:4. Args: readInStream: istream, P: Parameters, Lread: uint, LreadOriginal: uint, readName: char, Seq: char, SeqNum: char, Qual: char, clipOneMate: vector<ClipMate>, iReadAll: uint, readFilesIndex: uint32, readFilter: char, readNameExtra: string"]
pub fn readload_l4_readload<R: std::io::BufRead + ?Sized>(
    read_in_stream: &mut R,
    p: &crate::parameters_chimeric::Parameters,
    l_read: &mut u64,
    l_read_original: &mut u64,
    read_name: &mut String,
    seq: &mut String,
    seq_num: &mut Vec<u8>,
    qual: &mut String,
    clip_one_mate: &mut [crate::clip_mate::ClipMate],
    i_read_all: &mut u64,
    read_files_index: &mut u32,
    read_filter: &mut u8,
    read_name_extra: &mut String,
) -> Result<i32, String> {
    let first = match read_in_stream.fill_buf().map_err(|e| e.to_string())? {
        [] => return Ok(-1),
        buf => buf[0],
    };
    if first != b'@' && first != b'>' {
        return Ok(-1);
    }

    let mut line1 = String::new();
    read_in_stream
        .read_line(&mut line1)
        .map_err(|e| e.to_string())?;
    if line1.ends_with('\n') {
        line1.pop();
    }

    let mut line1_parts = line1.split_whitespace();
    read_name.clear();
    if let Some(name) = line1_parts.next() {
        read_name.push_str(name);
    }
    if read_name.len() >= DEF_READ_NAME_LENGTH_MAX - 1 {
        return Err(format!(
            "EXITING because of FATAL ERROR in reads input: read name is too long:0\nRead Name={}\nDEF_readNameLengthMax={}\nSOLUTION: increase DEF_readNameLengthMax in IncludeDefine.h and re-compile STAR\n",
            read_name, DEF_READ_NAME_LENGTH_MAX
        ));
    }
    if let Some(value) = line1_parts.next() {
        if let Ok(parsed) = value.parse::<u64>() {
            *i_read_all = parsed;
        }
    }
    if let Some(value) = line1_parts.next() {
        if let Some(byte) = value.as_bytes().first() {
            *read_filter = *byte;
        }
    }
    if let Some(value) = line1_parts.next() {
        if let Ok(parsed) = value.parse::<u32>() {
            *read_files_index = parsed;
        }
    }
    read_name_extra.clear();
    let extra = line1_parts.collect::<Vec<_>>().join(" ");
    if !extra.is_empty() {
        read_name_extra.push_str(&extra);
    }

    seq.clear();
    read_in_stream.read_line(seq).map_err(|e| e.to_string())?;
    if seq.ends_with('\n') {
        seq.pop();
    }
    *l_read = seq.len() as u64;

    if *l_read < 1 {
        return Err(format!(
            "EXITING because of FATAL ERROR in reads input: short read sequence line: {}\nRead Name={}\nRead Sequence=\"{}\"\nDEF_readNameLengthMax={}\nDEF_readSeqLengthMax={}\n",
            *l_read, read_name, seq, DEF_READ_NAME_LENGTH_MAX, DEF_READ_SEQ_LENGTH_MAX
        ));
    }
    if *l_read as usize > DEF_READ_SEQ_LENGTH_MAX {
        return Err(format!(
            "EXITING because of FATAL ERROR in reads input: Lread>={}   while DEF_readSeqLengthMax={}\nRead Name={}\nSOLUTION: increase DEF_readSeqLengthMax in IncludeDefine.h and re-compile STAR\n",
            *l_read, DEF_READ_SEQ_LENGTH_MAX, read_name
        ));
    }

    *l_read_original = *l_read;
    seq_num.clear();
    seq_num.resize(*l_read as usize, 0);
    sequencefuns_l131_convertnucleotidestonumbers(seq.as_bytes(), seq_num, *l_read as u64);

    let mut read_file_type = 0;
    if read_name.as_bytes().first() == Some(&b'@') {
        read_file_type = 2;
        let mut plus_line = String::new();
        read_in_stream
            .read_line(&mut plus_line)
            .map_err(|e| e.to_string())?;
        if let Some(byte) = plus_line.as_bytes().first() {
            clip_one_mate[0].clipped_info = *byte;
        }
    }

    clipmate_clip_l5_clipmate_clip(&mut clip_one_mate[0], l_read, seq_num);
    clipmate_clip_l5_clipmate_clip(&mut clip_one_mate[1], l_read, seq_num);

    if read_name.as_bytes().first() == Some(&b'@') {
        qual.clear();
        read_in_stream.read_line(qual).map_err(|e| e.to_string())?;
        if qual.ends_with('\n') {
            qual.pop();
        }
        if qual.len() as u64 != *l_read_original {
            return Err(format!(
                "EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length\n{}\n{}\n{}\nSOLUTION: fix your fastq file\n",
                read_name, seq, qual
            ));
        }
        if p.out_qs_conversion_add != 0 {
            let mut qual_bytes = qual.as_bytes().to_vec();
            for q in qual_bytes.iter_mut().take(*l_read_original as usize) {
                let mut qs = *q as i32 + p.out_qs_conversion_add;
                qs = qs.clamp(33, 126);
                *q = qs as u8;
            }
            *qual = String::from_utf8(qual_bytes).unwrap();
        }
    } else if read_name.as_bytes().first() == Some(&b'>') {
        read_file_type = 1;
        qual.clear();
        qual.extend(std::iter::repeat('A').take(*l_read_original as usize));
    } else {
        return Err(format!(
            "Unknown reads file format: header line does not start with @ or > : {}\n",
            read_name
        ));
    }

    for sep in &p.read_name_separator_char {
        if let Some(pos) = read_name.find(*sep) {
            read_name.truncate(pos);
        }
    }

    Ok(read_file_type)
}
