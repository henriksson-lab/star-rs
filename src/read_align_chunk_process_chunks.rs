#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlignChunk::processChunks` at STAR/source/ReadAlignChunk_processChunks.cpp:11. Args: "]
pub fn readalignchunk_processchunks_l11_readalignchunk_processchunks<F>(
    chunk: &mut crate::read_align_chunk::ReadAlignChunk,
    p: &mut crate::parameters_chimeric::Parameters,
    thread_chunks: &mut crate::thread_control::ThreadControl,
    stats_all: &mut crate::stats::Stats,
    time_current: libc::time_t,
    input_mates: &[String],
    mut one_read: F,
    mut real_context: Option<(
        &crate::genome::Genome,
        &mut crate::transcriptome::Transcriptome,
    )>,
) -> Result<crate::read_align_chunk::ReadAlignChunkProcessChunksResult, String>
where
    F: FnMut(&mut crate::read_align::ReadAlign) -> i32,
{
    use std::io::BufRead;

    let mut result = crate::read_align_chunk::ReadAlignChunkProcessChunksResult::default();
    chunk.no_reads_left = false;
    let mut new_file = false;
    let mut read_in: Vec<std::io::Cursor<&[u8]>> = input_mates
        .iter()
        .map(|s| std::io::Cursor::new(s.as_bytes()))
        .collect();
    if read_in.len() < p.read_nends as usize {
        read_in.resize_with(p.read_nends as usize, || std::io::Cursor::new(&[][..]));
    }
    if chunk.chunk_in.len() < p.read_nends as usize {
        chunk.chunk_in.resize(p.read_nends as usize, Vec::new());
    }
    if chunk.chunk_in_size_bytes_total.len() < p.read_nends as usize {
        chunk
            .chunk_in_size_bytes_total
            .resize(p.read_nends as usize, 0);
    }

    let retain_chunk_inputs = real_context.is_none();
    while !chunk.no_reads_left {
        if p.out_filter_by_sjout_stage < 2 {
            for imate in 0..p.read_nends as usize {
                chunk.chunk_in[imate].clear();
                chunk.chunk_in_size_bytes_total[imate] = 0;
            }

            loop {
                let next_char = match read_in[0].fill_buf().map_err(|e| e.to_string())? {
                    [] => break,
                    buf => buf[0],
                };
                let chunk_limit = p.chunk_in_size_bytes.max(1);
                if chunk.chunk_in_size_bytes_total[0] >= chunk_limit
                    || chunk.chunk_in_size_bytes_total
                        [1usize.min(chunk.chunk_in_size_bytes_total.len() - 1)]
                        >= chunk_limit
                {
                    break;
                }
                if p.read_map_number != 0 && p.i_read_all == p.read_map_number {
                    break;
                } else if p.read_files_type_n == 10 && p.out_filter_by_sjout_stage != 2 {
                    if next_char == b'@' {
                        let mut discard = String::new();
                        read_in[0]
                            .read_line(&mut discard)
                            .map_err(|e| e.to_string())?;
                        continue;
                    }

                    let mut lines = Vec::new();
                    let mut line1 = String::new();
                    if read_in[0]
                        .read_line(&mut line1)
                        .map_err(|e| e.to_string())?
                        == 0
                    {
                        break;
                    }
                    if line1.trim_end() == "FILE" {
                        new_file = true;
                    } else {
                        lines.push(line1);
                        p.i_read_all += 1;
                        for _ in 1..p.read_nmates {
                            let mut line = String::new();
                            read_in[0].read_line(&mut line).map_err(|e| e.to_string())?;
                            lines.push(line);
                        }

                        let mut first_name = String::new();
                        let mut first_flag = 0u64;
                        let mut imate1 = 0usize;
                        for (imate, line) in lines.iter().enumerate() {
                            let fields: Vec<&str> = line.split_whitespace().collect();
                            if fields.len() < 11 {
                                return Err(format!(
                                    "EXITING because of FATAL ERROR in input reads: wrong SAM line format for read # {}\n",
                                    p.i_read_all
                                ));
                            }
                            let str1 = fields[0].to_string();
                            let flag1 = fields[1].parse::<u64>().map_err(|_| {
                                format!(
                                    "EXITING because of FATAL ERROR in input BAM file: could not parse FLAG for read {}\n",
                                    str1
                                )
                            })?;
                            if imate == 0 {
                                first_name = str1.clone();
                                first_flag = flag1;
                            } else {
                                if first_name != str1 {
                                    return Err(format!(
                                        "EXITING because of FATAL ERROR in input BAM file: the consecutive lines in paired-end BAM have different read IDs:\n{}   vs   {}\n\n SOLUTION: fix BAM file formatting. Paired-end reads should be always consecutive lines, with exactly 2 lines per paired-end read",
                                        first_name, str1
                                    ));
                                }
                                if !(((first_flag & 0x40) != 0 && (flag1 & 0x80) != 0)
                                    || ((flag1 & 0x40) != 0 && (first_flag & 0x80) != 0))
                                {
                                    return Err(format!(
                                        "EXITING because of FATAL ERROR in input BAM file: the consecutive lines in paired-end BAM have wrong mate FLAG bits:\n{}   {}   vs   {}   {}\n\n SOLUTION: fix BAM file formatting. Paired-end reads should be always consecutive lines, with exactly 2 lines per paired-end read. Mate1 should have 0x40 bit set in the FLAG, Mate2 should have 0x80 bit set in the FLAG",
                                        first_name, first_flag, str1, flag1
                                    ));
                                }
                            }
                            let pass_filter_illumina = if (flag1 & 0x800) != 0 { 'Y' } else { 'N' };
                            if imate == 1 {
                                imate1 = 1 - imate1;
                            } else if p.read_nmates == 2 && (flag1 & 0x80) != 0 {
                                imate1 = 1;
                            } else {
                                imate1 = 0;
                            }
                            let read_id = if p.out_sam_read_id == "Number" {
                                format!("@{}", p.i_read_all)
                            } else {
                                format!("@{}", str1)
                            };
                            chunk.chunk_in[imate1].extend_from_slice(
                                format!(
                                    "{} {} {} {}",
                                    read_id, p.i_read_all, pass_filter_illumina, p.read_files_index
                                )
                                .as_bytes(),
                            );
                            let mut seq1 = fields[9].to_string();
                            let mut qual1 = fields[10].to_string();
                            if (flag1 & 0x10) != 0 {
                                sequencefuns_l56_revcomplementnucleotides(&mut seq1);
                                qual1 = qual1.chars().rev().collect();
                            }
                            let attrs = if fields.len() > 11 {
                                format!(" {}", fields[11..].join(" "))
                            } else {
                                String::new()
                            };
                            chunk.chunk_in[imate1].extend_from_slice(
                                format!("{}\n{}\n+\n{}\n", attrs, seq1, qual1).as_bytes(),
                            );
                            chunk.chunk_in_size_bytes_total[imate1] =
                                chunk.chunk_in[imate1].len() as u64;
                        }
                    }
                } else if next_char == b'@' {
                    p.i_read_all += 1;
                    if p.out_filter_by_sjout_stage != 2 {
                        let mut name_line = String::new();
                        read_in[0]
                            .read_line(&mut name_line)
                            .map_err(|e| e.to_string())?;
                        let mut fields = name_line.split_whitespace();
                        let mut read_id = fields.next().unwrap_or_default().to_string();
                        readalignchunk_processchunks_l298_removestringendcontrol(&mut read_id);
                        if p.out_sam_read_id_number {
                            read_id = format!("@{}", p.i_read_all);
                        }
                        let mut pass_filter_illumina = 'N';
                        if let Some(field2) = fields.next() {
                            let bytes = field2.as_bytes();
                            if bytes.len() >= 4
                                && bytes[1] == b':'
                                && bytes[2] == b'Y'
                                && bytes[3] == b':'
                            {
                                pass_filter_illumina = 'Y';
                            }
                        }
                        read_id.push_str(&format!(
                            " {} {} {}",
                            p.i_read_all, pass_filter_illumina, p.read_files_index
                        ));
                        for imate in 1..p.read_nends as usize {
                            let mut discard = String::new();
                            read_in[imate]
                                .read_line(&mut discard)
                                .map_err(|e| e.to_string())?;
                        }
                        for imate in 0..p.read_nends as usize {
                            chunk.chunk_in[imate].extend_from_slice(read_id.as_bytes());
                            chunk.chunk_in[imate].push(b'\n');
                        }
                    }
                    for imate in 0..p.read_nends as usize {
                        if p.out_filter_by_sjout_stage == 2 {
                            let mut line = Vec::new();
                            readalignchunk_processchunks_l284_fastqreadoneline(
                                &mut read_in[imate],
                                &mut line,
                            )
                            .map_err(|e| e.to_string())?;
                            chunk.chunk_in[imate].extend_from_slice(&line);
                        }
                        let mut line = Vec::new();
                        readalignchunk_processchunks_l284_fastqreadoneline(
                            &mut read_in[imate],
                            &mut line,
                        )
                        .map_err(|e| e.to_string())?;
                        chunk.chunk_in[imate].extend_from_slice(&line);
                        let mut discard = String::new();
                        read_in[imate]
                            .read_line(&mut discard)
                            .map_err(|e| e.to_string())?;
                        chunk.chunk_in[imate].extend_from_slice(b"+\n");
                        line.clear();
                        readalignchunk_processchunks_l284_fastqreadoneline(
                            &mut read_in[imate],
                            &mut line,
                        )
                        .map_err(|e| e.to_string())?;
                        chunk.chunk_in[imate].extend_from_slice(&line);
                        chunk.chunk_in_size_bytes_total[imate] = chunk.chunk_in[imate].len() as u64;
                    }
                } else if next_char == b'>' {
                    p.i_read_all += 1;
                    for imate in 0..p.read_nends as usize {
                        if p.out_filter_by_sjout_stage != 2 {
                            let mut header = String::new();
                            read_in[imate]
                                .read_line(&mut header)
                                .map_err(|e| e.to_string())?;
                            let head_word = header.split_whitespace().next().unwrap_or_default();
                            let read_id = if p.out_sam_read_id == "Number" {
                                format!(">{}", p.i_read_all)
                            } else {
                                head_word.to_string()
                            };
                            chunk.chunk_in[imate].extend_from_slice(
                                format!(
                                    "{} {} {} {} \n",
                                    read_id, p.i_read_all, 'N', p.read_files_index
                                )
                                .as_bytes(),
                            );
                        }
                        loop {
                            let fasta_next =
                                match read_in[imate].fill_buf().map_err(|e| e.to_string())? {
                                    [] => break,
                                    buf => buf[0],
                                };
                            if matches!(fasta_next, b'@' | b'>' | b' ' | b'\n') {
                                break;
                            }
                            let mut seq_line = Vec::new();
                            read_in[imate]
                                .read_until(b'\n', &mut seq_line)
                                .map_err(|e| e.to_string())?;
                            while seq_line
                                .last()
                                .is_some_and(|&byte| byte == b'\n' || byte < 33)
                            {
                                seq_line.pop();
                            }
                            chunk.chunk_in[imate].extend_from_slice(&seq_line);
                        }
                        chunk.chunk_in[imate].push(b'\n');
                        chunk.chunk_in_size_bytes_total[imate] = chunk.chunk_in[imate].len() as u64;
                    }
                } else if next_char == b' ' || next_char == b'\n' {
                    result.log_main.push_str(&format!(
                        "Thread #{} end of input stream, nextChar={}\n",
                        chunk.i_thread, next_char
                    ));
                    break;
                } else {
                    let mut marker_line = String::new();
                    read_in[0]
                        .read_line(&mut marker_line)
                        .map_err(|e| e.to_string())?;
                    let mut parts = marker_line.splitn(2, char::is_whitespace);
                    let word1 = parts.next().unwrap_or_default();
                    if word1 == "FILE" {
                        new_file = true;
                        p.read_files_index = parts
                            .next()
                            .and_then(|v| v.split_whitespace().next())
                            .and_then(|v| v.parse::<u32>().ok())
                            .unwrap_or(p.read_files_index);
                    } else {
                        return Err(format!(
                            "EXITING because of FATAL ERROR in input reads: wrong read ID line format: the read ID lines should start with @ or > \nOffending line for read # {}\n{} {}\nSOLUTION: verify and correct the input read files\n",
                            p.i_read_all + 1,
                            word1,
                            parts.next().unwrap_or_default()
                        ));
                    }
                }

                if new_file {
                    result
                        .log_main
                        .push_str(&format!("Starting to map file # {}\n", p.read_files_index));
                    for imate in 0..p.read_files_names.len() {
                        let file_name = p.read_files_names[imate]
                            .get(p.read_files_index as usize)
                            .cloned()
                            .unwrap_or_default();
                        result
                            .log_main
                            .push_str(&format!("mate {}:   {}\n", imate + 1, file_name));
                        if imate > 0 && imate < read_in.len() {
                            let mut marker_line = String::new();
                            read_in[imate]
                                .read_line(&mut marker_line)
                                .map_err(|e| e.to_string())?;
                            let mut parts = marker_line.split_whitespace();
                            if parts.next() == Some("FILE") {
                                p.read_files_index = parts
                                    .next()
                                    .and_then(|v| v.parse::<u32>().ok())
                                    .unwrap_or(p.read_files_index);
                            }
                        }
                    }
                    new_file = false;
                }
            }

            if chunk
                .chunk_in_size_bytes_total
                .first()
                .copied()
                .unwrap_or_default()
                == 0
            {
                chunk.no_reads_left = true;
                chunk.i_chunk_in = thread_chunks.chunk_in_n;
                thread_chunks.chunk_in_n += 1;
            } else {
                chunk.no_reads_left = false;
                chunk.i_chunk_in = thread_chunks.chunk_in_n;
                thread_chunks.chunk_in_n += 1;
                result.chunks_read += 1;
            }
            for imate in 0..p.read_nends as usize {
                chunk.chunk_in[imate].push(b'\n');
                chunk.chunk_in_size_bytes_total[imate] = chunk.chunk_in[imate].len() as u64 - 1;
            }
        } else {
            chunk.no_reads_left = true;
            for imate in 0..p.read_nends as usize {
                let bytes = chunk
                    .chunk_out_filter_by_sjout_files
                    .get(imate)
                    .map(|s| s.as_bytes().to_vec())
                    .unwrap_or_default();
                chunk.chunk_in[imate] = bytes;
                chunk.chunk_in_size_bytes_total[imate] = chunk.chunk_in[imate].len() as u64;
            }
        }

        if retain_chunk_inputs {
            result.chunk_inputs.push(chunk.chunk_in.clone());
        }
        let map_result = readalignchunk_mapchunk_l7_readalignchunk_mapchunk(
            chunk,
            p,
            stats_all,
            time_current,
            &mut one_read,
            match real_context.as_mut() {
                Some((map_gen, transcriptome)) => Some((*map_gen, &mut **transcriptome)),
                None => None,
            },
        )?;
        result.log_main.push_str(&map_result.log_main);
        result.map_chunks.push(map_result);

        if chunk.i_thread == 0 && p.run_thread_n > 1 && p.out_sam_order == "PairedKeepInputOrder" {
            result
                .paired_keep_input_order_cat_after_chunks
                .push(thread_chunks.chunk_out_n);
        }
    }

    if p.out_filter_by_sjout_stage != 1 && chunk.ra.i_read > 0 {
        if p.out_bam_unsorted {
            result.flushed_bam_unsorted = true;
        }
        if p.out_bam_coord {
            result.flushed_bam_coord = true;
        }
        if chunk.chunk_out_bam_quant.is_some() {
            result.flushed_bam_quant = true;
        }
        if p.p_ch.segment_min > 0 {
            result.chim_sam_cat_path = chunk.chunk_out_chim_sam_path.clone();
            result.chim_junction_cat_path = chunk.chunk_out_chim_junction_path.clone();
        }
        if p.out_reads_unmapped == "Fastx" {
            result.unmapped_fastx_cat_paths = chunk.chunk_out_unmapped_reads_paths.clone();
        }
    }
    result
        .log_main
        .push_str(&format!("Completed: thread #{}\n", chunk.i_thread));
    chunk.log_main.push_str(&result.log_main);
    Ok(result)
}

#[doc = "Original `fastqReadOneLine` at STAR/source/ReadAlignChunk_processChunks.cpp:284. Args: streamIn: ifstream, arrIn: char"]
pub fn readalignchunk_processchunks_l284_fastqreadoneline<R: std::io::BufRead>(
    stream_in: &mut R,
    arr_in: &mut Vec<u8>,
) -> std::io::Result<u64> {
    arr_in.clear();
    let mut line = Vec::new();
    if stream_in.read_until(b'\n', &mut line)? == 0 {
        return Ok(0);
    }

    if line.last() == Some(&b'\n') {
        line.pop();
    }
    if line.last().is_some_and(|&byte| byte < 33) {
        line.pop();
    }
    line.push(b'\n');

    arr_in.extend_from_slice(&line);
    Ok(arr_in.len() as u64)
}

#[doc = "Original `removeStringEndControl` at STAR/source/ReadAlignChunk_processChunks.cpp:298. Args: str: string"]
pub fn readalignchunk_processchunks_l298_removestringendcontrol(str_: &mut String) {
    if str_.as_bytes().last().is_some_and(|&byte| byte < 33) {
        str_.pop();
    }
}
