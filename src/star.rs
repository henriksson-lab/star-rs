#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `usage` at STAR/source/STAR.cpp:36. Args: usageType: int"]
pub fn star_l36_usage(
    usage_type: i32,
    star_version: &str,
    compilation_time_place: &str,
    parameters_default: &[u8],
) -> String {
    let mut out = String::new();
    out.push_str("Usage: STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq\n");
    out.push_str("Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2022\n\n");
    out.push_str("STAR version=");
    out.push_str(star_version);
    out.push('\n');
    out.push_str("STAR compilation time,server,dir=");
    out.push_str(compilation_time_place);
    out.push('\n');
    out.push_str("For more details see:\n");
    out.push_str("<https://github.com/alexdobin/STAR>\n");
    out.push_str("<https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>\n");

    if usage_type == 0 {
        out.push_str("\nTo list all parameters, run STAR --help\n");
    } else if usage_type == 1 {
        out.push_str(&String::from_utf8_lossy(parameters_default));
    }

    out
}

#[doc = "Original `main` at STAR/source/STAR.cpp:58. Args: argInN: int, argIn: char"]
pub fn star_l58_main(
    arg_in: &[String],
    mut p: crate::parameters_chimeric::Parameters,
    parameters_default: &[u8],
    genome_main_in: Option<crate::genome::Genome>,
    transcriptome_main_in: Option<crate::transcriptome::Transcriptome>,
    pass1_chunks_in: Option<Vec<crate::read_align_chunk::ReadAlignChunk>>,
    map_chunks_in: Option<Vec<crate::read_align_chunk::ReadAlignChunk>>,
    existing_read_files: &std::collections::BTreeSet<String>,
    temp_files_by_bam_bin: &[Vec<Vec<u8>>],
    signal_records: &[crate::parameters_chimeric::SignalFromBamRecord],
) -> Result<crate::parameters_chimeric::StarMainResult, String> {
    let mut result = crate::parameters_chimeric::StarMainResult {
        parameters: p.clone(),
        ..Default::default()
    };

    if arg_in.len() <= 1 {
        result.usage = star_l36_usage(0, "2.7.11b", "COMPILATION_TIME_PLACE", parameters_default);
        result.exit_code = 0;
        return Ok(result);
    }
    if arg_in.len() == 2 && (arg_in[1] == "-h" || arg_in[1] == "--help") {
        result.usage = star_l36_usage(1, "2.7.11b", "COMPILATION_TIME_PLACE", parameters_default);
        result.exit_code = 0;
        return Ok(result);
    }

    let raw_time_start = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    let run_mode = p
        .run_mode_in
        .first()
        .map(|s| s.as_str())
        .unwrap_or("alignReads")
        .to_string();

    result
        .log_stdout
        .push_str(&format!("\t{}\n", p.command_line));
    result.log_stdout.push_str(&format!(
        "\tSTAR version: {}   compiled: {}\n",
        "2.7.11b", "COMPILATION_TIME_PLACE"
    ));
    result.log_stdout.push_str(&format!(
        "{} ..... started STAR run\n",
        timefunctions_l14_timemonthdaytime(raw_time_start)
    ));

    if run_mode == "genomeGenerate" {
        let gtf_contents = if p.p_ge.sjdb_gtf_file != "-" && p.p_ge.sjdb_overhang > 0 {
            Some(
                crate::io_utils::read_to_string_auto_gzip(&p.p_ge.sjdb_gtf_file).map_err(|_| {
                    format!(
                        "FATAL error, could not open file pGe.sjdbGTFfile={}\n",
                        p.p_ge.sjdb_gtf_file
                    )
                })?,
            )
        } else {
            None
        };
        let mut genome_main = genome_l15_genome_genome(p.p_ge.clone());
        let gen_out = genome_genomegenerate_l98_genome_genomegenerate(
            &mut genome_main,
            &mut p,
            gtf_contents.as_deref(),
        )?;
        result.log_main.push_str(&gen_out.log_main);
        result.log_stdout.push_str(&gen_out.log_stdout);
        result.genome_generate.push(gen_out);

        if p.p_ge.transform.type_ > 0 {
            p.p_ge.transform.type_ = 0;
            p.p_ge.transform.type_string = "None".to_string();
            p.p_ge.transform.vcf_file = "-".to_string();
            p.p_ge.g_dir.push_str("/OriginalGenome/");
            let mut genome_orig = genome_l15_genome_genome(p.p_ge.clone());
            let gen_orig = genome_genomegenerate_l98_genome_genomegenerate(
                &mut genome_orig,
                &mut p,
                gtf_contents.as_deref(),
            )?;
            result.log_main.push_str(&gen_orig.log_main);
            result.log_stdout.push_str(&gen_orig.log_stdout);
            result.genome_generate.push(gen_orig);
        }

        if !p.out_file_tmp.is_empty() {
            if sysremovedir_l25_sysremovedir(std::path::Path::new(&p.out_file_tmp)).is_ok() {
                result.removed_tmp = true;
            }
        }
        result
            .log_main
            .push_str("DONE: Genome generation, EXITING\n");
        result.exit_code = 0;
        result.parameters = p;
        result.genome = Some(genome_main);
        return Ok(result);
    } else if run_mode == "liftOver" {
        for (ii, chain_file) in p.p_ge.g_chain_files.iter().enumerate() {
            let chain = chain_l5_chain_chain(chain_file)?;
            let out_file_name = format!("{}GTFliftOver_{}.gtf", p.out_file_name_prefix, ii + 1);
            if let Some(parent) = std::path::Path::new(&out_file_name).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|err| err.to_string())?;
            }
            chain_l58_chain_liftovergtf(&chain, &p.p_ge.sjdb_gtf_file, &out_file_name)?;
        }
        result
            .log_main
            .push_str("DONE: lift-over of GTF file, EXITING\n");
        result.exit_code = 0;
        result.parameters = p;
        return Ok(result);
    } else if run_mode == "inputAlignmentsFromBAM" {
        let bam_file_name = p.input_bam_file.clone();
        let input_bam_bytes = std::fs::read(&bam_file_name).map_err(|_| {
            format!(
                "EXITING because of fatal ERROR: could not open --inputBAMfile {}\n",
                bam_file_name
            )
        })?;
        let bam_bytes = if input_bam_bytes.starts_with(b"BAM\x01") {
            input_bam_bytes
        } else if input_bam_bytes.starts_with(&[0x1f, 0x8b]) {
            let mut bam_bytes = Vec::<u8>::new();
            let mut block_start = 0usize;
            while block_start < input_bam_bytes.len() {
                if block_start + 18 > input_bam_bytes.len()
                    || input_bam_bytes[block_start] != 0x1f
                    || input_bam_bytes[block_start + 1] != 0x8b
                    || input_bam_bytes[block_start + 2] != 8
                    || input_bam_bytes[block_start + 3] & 4 == 0
                {
                    return Err(format!(
                        "EXITING because of fatal ERROR: could not read BAM header from --inputBAMfile {}\n",
                        bam_file_name
                    ));
                }
                let xlen = u16::from_le_bytes(
                    input_bam_bytes[block_start + 10..block_start + 12]
                        .try_into()
                        .unwrap(),
                ) as usize;
                let extra_start = block_start + 12;
                let extra_end = extra_start.checked_add(xlen).ok_or_else(|| {
                    format!(
                        "EXITING because of fatal ERROR: malformed BAM header in --inputBAMfile {}\n",
                        bam_file_name
                    )
                })?;
                if extra_end > input_bam_bytes.len() {
                    return Err(format!(
                        "EXITING because of fatal ERROR: truncated BAM header in --inputBAMfile {}\n",
                        bam_file_name
                    ));
                }
                let mut bsize = None;
                let mut extra_pos = extra_start;
                while extra_pos + 4 <= extra_end {
                    let slen = u16::from_le_bytes(
                        input_bam_bytes[extra_pos + 2..extra_pos + 4]
                            .try_into()
                            .unwrap(),
                    ) as usize;
                    extra_pos += 4;
                    if extra_pos + slen > extra_end {
                        return Err(format!(
                            "EXITING because of fatal ERROR: truncated BAM header in --inputBAMfile {}\n",
                            bam_file_name
                        ));
                    }
                    if input_bam_bytes[extra_pos - 4] == b'B'
                        && input_bam_bytes[extra_pos - 3] == b'C'
                        && slen == 2
                    {
                        bsize = Some(u16::from_le_bytes(
                            input_bam_bytes[extra_pos..extra_pos + 2]
                                .try_into()
                                .unwrap(),
                        ) as usize);
                    }
                    extra_pos += slen;
                }
                let block_size = bsize
                    .ok_or_else(|| {
                        format!(
                            "EXITING because of fatal ERROR: could not read BAM header from --inputBAMfile {}\n",
                            bam_file_name
                        )
                    })?
                    .checked_add(1)
                    .ok_or_else(|| {
                        format!(
                            "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                            bam_file_name
                        )
                    })?;
                let block_end = block_start.checked_add(block_size).ok_or_else(|| {
                    format!(
                        "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                        bam_file_name
                    )
                })?;
                if block_end > input_bam_bytes.len() || block_size < 26 {
                    return Err(format!(
                        "EXITING because of fatal ERROR: truncated BAM record in --inputBAMfile {}\n",
                        bam_file_name
                    ));
                }
                let mut compressed_start = extra_end;
                if input_bam_bytes[block_start + 3] & 8 != 0 {
                    while compressed_start < block_end && input_bam_bytes[compressed_start] != 0 {
                        compressed_start += 1;
                    }
                    compressed_start = compressed_start
                        .checked_add(usize::from(compressed_start < block_end))
                        .ok_or_else(|| {
                            format!(
                                "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                                bam_file_name
                            )
                        })?;
                }
                if input_bam_bytes[block_start + 3] & 16 != 0 {
                    while compressed_start < block_end && input_bam_bytes[compressed_start] != 0 {
                        compressed_start += 1;
                    }
                    compressed_start = compressed_start
                        .checked_add(usize::from(compressed_start < block_end))
                        .ok_or_else(|| {
                            format!(
                                "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                                bam_file_name
                            )
                        })?;
                }
                if input_bam_bytes[block_start + 3] & 2 != 0 {
                    compressed_start = compressed_start.checked_add(2).ok_or_else(|| {
                        format!(
                            "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                            bam_file_name
                        )
                    })?;
                }
                let compressed_end = block_end - 8;
                if compressed_start > compressed_end {
                    return Err(format!(
                        "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                        bam_file_name
                    ));
                }
                let mut decoder = flate2::read::DeflateDecoder::new(
                    &input_bam_bytes[compressed_start..compressed_end],
                );
                std::io::Read::read_to_end(&mut decoder, &mut bam_bytes).map_err(|_| {
                    format!(
                        "EXITING because of fatal ERROR: could not read BAM header from --inputBAMfile {}\n",
                        bam_file_name
                    )
                })?;
                block_start = block_end;
            }
            bam_bytes
        } else {
            input_bam_bytes
        };
        if bam_bytes.len() < 12 || &bam_bytes[0..4] != b"BAM\x01" {
            return Err(format!(
                "EXITING because of fatal ERROR: could not read BAM header from --inputBAMfile {}\n",
                bam_file_name
            ));
        }

        let mut bam_pos = 4usize;
        let header_len_i32 =
            i32::from_ne_bytes(bam_bytes[bam_pos..bam_pos + 4].try_into().unwrap());
        if header_len_i32 < 0 {
            return Err(format!(
                "EXITING because of fatal ERROR: malformed BAM header in --inputBAMfile {}\n",
                bam_file_name
            ));
        }
        let header_len = header_len_i32 as usize;
        bam_pos += 4;
        if bam_pos
            .checked_add(header_len)
            .and_then(|pos| pos.checked_add(4))
            .is_none_or(|end| end > bam_bytes.len())
        {
            return Err(format!(
                "EXITING because of fatal ERROR: truncated BAM header in --inputBAMfile {}\n",
                bam_file_name
            ));
        }
        bam_pos += header_len;
        let ref_n = i32::from_ne_bytes(bam_bytes[bam_pos..bam_pos + 4].try_into().unwrap());
        if ref_n < 0 {
            return Err(format!(
                "EXITING because of fatal ERROR: malformed BAM reference header in --inputBAMfile {}\n",
                bam_file_name
            ));
        }
        bam_pos += 4;
        let mut target_names = Vec::<String>::new();
        let mut target_lens = Vec::<u32>::new();
        for _ in 0..ref_n {
            if bam_pos + 4 > bam_bytes.len() {
                return Err(format!(
                    "EXITING because of fatal ERROR: truncated BAM reference header in --inputBAMfile {}\n",
                    bam_file_name
                ));
            }
            let name_len_i32 =
                i32::from_ne_bytes(bam_bytes[bam_pos..bam_pos + 4].try_into().unwrap());
            if name_len_i32 <= 0 {
                return Err(format!(
                    "EXITING because of fatal ERROR: malformed BAM reference header in --inputBAMfile {}\n",
                    bam_file_name
                ));
            }
            let name_len = name_len_i32 as usize;
            bam_pos += 4;
            if bam_pos
                .checked_add(name_len)
                .and_then(|pos| pos.checked_add(4))
                .is_none_or(|end| end > bam_bytes.len())
            {
                return Err(format!(
                    "EXITING because of fatal ERROR: truncated BAM reference header in --inputBAMfile {}\n",
                    bam_file_name
                ));
            }
            let name_bytes = &bam_bytes[bam_pos..bam_pos + name_len];
            let name_end = name_bytes
                .iter()
                .position(|byte| *byte == 0)
                .unwrap_or(name_bytes.len());
            target_names.push(String::from_utf8_lossy(&name_bytes[..name_end]).to_string());
            bam_pos += name_len;
            target_lens.push(u32::from_ne_bytes(
                bam_bytes[bam_pos..bam_pos + 4].try_into().unwrap(),
            ));
            bam_pos += 4;
        }
        let records_start = bam_pos;
        let mut signal_records = Vec::<crate::parameters_chimeric::SignalFromBamRecord>::new();
        let mut duplicate_records = Vec::<Vec<u32>>::new();
        let mut duplicate_record_lens = Vec::<usize>::new();
        while bam_pos < bam_bytes.len() {
            if bam_pos + 4 > bam_bytes.len() {
                return Err(format!(
                    "EXITING because of fatal ERROR: truncated BAM record in --inputBAMfile {}\n",
                    bam_file_name
                ));
            }
            let block_len = i32::from_ne_bytes(bam_bytes[bam_pos..bam_pos + 4].try_into().unwrap());
            if block_len < 32 {
                return Err(format!(
                    "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                    bam_file_name
                ));
            }
            let record_len = 4usize.checked_add(block_len as usize).ok_or_else(|| {
                format!(
                    "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                    bam_file_name
                )
            })?;
            if bam_pos
                .checked_add(record_len)
                .is_none_or(|end| end > bam_bytes.len())
            {
                return Err(format!(
                    "EXITING because of fatal ERROR: truncated BAM record in --inputBAMfile {}\n",
                    bam_file_name
                ));
            }
            let record = &bam_bytes[bam_pos..bam_pos + record_len];
            let mut bam1 = crate::bam_output::Bam1::default();
            if bamfunctions_l30_bam_read1_fromarray(record, &mut bam1) < 0 {
                return Err(format!(
                    "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                    bam_file_name
                ));
            }
            let cigar_start = 4usize
                .checked_add(32)
                .and_then(|v| v.checked_add(bam1.core.l_qname as usize))
                .ok_or_else(|| {
                    format!(
                        "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                        bam_file_name
                    )
                })?;
            let cigar_bytes = (bam1.core.n_cigar as usize).checked_mul(4).ok_or_else(|| {
                format!(
                    "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                    bam_file_name
                )
            })?;
            let seq_start = cigar_start.checked_add(cigar_bytes).ok_or_else(|| {
                format!(
                    "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                    bam_file_name
                )
            })?;
            let seq_bytes = (bam1.core.l_qseq as usize).div_ceil(2);
            let qual_start = seq_start.checked_add(seq_bytes).ok_or_else(|| {
                format!(
                    "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                    bam_file_name
                )
            })?;
            let aux_start = qual_start
                .checked_add(bam1.core.l_qseq.max(0) as usize)
                .ok_or_else(|| {
                    format!(
                        "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                        bam_file_name
                    )
                })?;
            if aux_start > record.len() {
                return Err(format!(
                    "EXITING because of fatal ERROR: malformed BAM record in --inputBAMfile {}\n",
                    bam_file_name
                ));
            }
            let mut cigar = Vec::<u32>::new();
            for ic in 0..bam1.core.n_cigar as usize {
                let start = cigar_start + ic * 4;
                if start + 4 <= record.len() {
                    cigar.push(u32::from_ne_bytes(
                        record[start..start + 4].try_into().unwrap(),
                    ));
                }
            }
            let mut nh = None;
            let mut aux_pos = aux_start;
            while aux_pos + 3 <= record.len() {
                let tag0 = record[aux_pos];
                let tag1 = record[aux_pos + 1];
                let type_ = record[aux_pos + 2];
                aux_pos += 3;
                match type_ {
                    b'c' | b'C' | b'A' => {
                        if aux_pos + 1 > record.len() {
                            break;
                        }
                        let value = record[aux_pos] as u32;
                        if tag0 == b'N' && tag1 == b'H' {
                            nh = Some(value);
                        }
                        aux_pos += 1;
                    }
                    b's' | b'S' => {
                        if aux_pos + 2 > record.len() {
                            break;
                        }
                        let value =
                            u16::from_ne_bytes(record[aux_pos..aux_pos + 2].try_into().unwrap())
                                as u32;
                        if tag0 == b'N' && tag1 == b'H' {
                            nh = Some(value);
                        }
                        aux_pos += 2;
                    }
                    b'i' | b'I' => {
                        if aux_pos + 4 > record.len() {
                            break;
                        }
                        let value =
                            u32::from_ne_bytes(record[aux_pos..aux_pos + 4].try_into().unwrap());
                        if tag0 == b'N' && tag1 == b'H' {
                            nh = Some(value);
                        }
                        aux_pos += 4;
                    }
                    b'f' => {
                        if aux_pos + 4 > record.len() {
                            break;
                        }
                        aux_pos += 4;
                    }
                    b'Z' | b'H' => {
                        while aux_pos < record.len() && record[aux_pos] != 0 {
                            aux_pos += 1;
                        }
                        aux_pos += usize::from(aux_pos < record.len());
                    }
                    _ => break,
                }
            }
            signal_records.push(crate::parameters_chimeric::SignalFromBamRecord {
                tid: bam1.core.tid,
                pos: bam1.core.pos.max(0) as u32,
                flag: bam1.core.flag as u16,
                cigar,
                nh,
            });
            let words = record_len.div_ceil(4);
            let mut record_words = vec![0u32; words];
            for iw in 0..words {
                let start = iw * 4;
                let end = (start + 4).min(record_len);
                let mut word = [0u8; 4];
                word[..end - start].copy_from_slice(&record[start..end]);
                record_words[iw] = u32::from_ne_bytes(word);
            }
            duplicate_records.push(record_words);
            duplicate_record_lens.push(record_len);
            bam_pos += record_len;
        }

        if p.out_wig_flags.yes {
            result.log_stdout.push_str(&format!(
                "{} ..... reading from BAM, output wiggle\n",
                timefunctions_l4_timemonthdaytime()
            ));
            result.log_main.push_str(&format!(
                "{} ..... reading from BAM, output wiggle\n",
                timefunctions_l4_timemonthdaytime()
            ));
            result.signal = Some(signalfrombam_l5_signalfrombam(
                &format!("{}Signal", p.out_file_name_prefix),
                &p,
                &target_names,
                &target_lens,
                &signal_records,
            )?);
            result.log_stdout.push_str(&format!(
                "{} ..... done\n",
                timefunctions_l4_timemonthdaytime()
            ));
            result.log_main.push_str(&format!(
                "{} ..... done\n",
                timefunctions_l4_timemonthdaytime()
            ));
        } else if p.bam_remove_duplicates_type != "-" {
            result.log_stdout.push_str(&format!(
                "{} ..... reading from BAM, remove duplicates, output BAM\n",
                timefunctions_l4_timemonthdaytime()
            ));
            result.log_main.push_str(&format!(
                "{} ..... reading from BAM, remove duplicates, output BAM\n",
                timefunctions_l4_timemonthdaytime()
            ));
            bamremoveduplicates_l114_bamremoveduplicates(
                &mut duplicate_records,
                p.bam_remove_duplicates_mate2bases_n,
                p.bam_remove_duplicates_mark_multi,
            )?;
            result
                .processed_bam_output
                .extend_from_slice(&bam_bytes[..records_start]);
            for (record, record_len) in duplicate_records.iter().zip(duplicate_record_lens.iter()) {
                let mut record_bytes = Vec::with_capacity(record.len() * 4);
                for word in record {
                    record_bytes.extend_from_slice(&word.to_ne_bytes());
                }
                result
                    .processed_bam_output
                    .extend_from_slice(&record_bytes[..*record_len]);
            }
            result.log_stdout.push_str(&format!(
                "{} ..... done\n",
                timefunctions_l4_timemonthdaytime()
            ));
            result.log_main.push_str(&format!(
                "{} ..... done\n",
                timefunctions_l4_timemonthdaytime()
            ));
        } else {
            return Err("EXITING because of fatal INPUT ERROR: at the moment --runMode inputFromBAM only works with --outWigType bedGraph OR --bamRemoveDuplicatesType Identical\n".to_string());
        }

        if !p.out_file_tmp.is_empty() {
            if sysremovedir_l25_sysremovedir(std::path::Path::new(&p.out_file_tmp)).is_ok() {
                result.removed_tmp = true;
            }
        }
        result.exit_code = 0;
        result.parameters = p;
        return Ok(result);
    } else if run_mode == "soloCellFiltering" {
        let input_prefix = p.run_mode_in.get(1).cloned().unwrap_or_default();
        let input_prefix = format!("{}/", input_prefix);
        let matrix_file = format!(
            "{}{}",
            input_prefix,
            p.p_solo
                .out_file_names
                .get(3)
                .cloned()
                .unwrap_or_else(|| "matrix.mtx".to_string())
        );
        let barcodes_file = format!(
            "{}{}",
            input_prefix,
            p.p_solo
                .out_file_names
                .get(2)
                .cloned()
                .unwrap_or_else(|| "barcodes.tsv".to_string())
        );
        let features_file = format!(
            "{}{}",
            input_prefix,
            p.p_solo
                .out_file_names
                .get(1)
                .cloned()
                .unwrap_or_else(|| "features.tsv".to_string())
        );
        let matrix_contents = std::fs::read_to_string(&matrix_file).map_err(|err| {
            format!(
                "EXITING because of fatal INPUT FILE error: could not open {}\n{}",
                matrix_file, err
            )
        })?;
        let barcodes_contents = std::fs::read_to_string(&barcodes_file).map_err(|err| {
            format!(
                "EXITING because of fatal INPUT FILE error: could not open {}\n{}",
                barcodes_file, err
            )
        })?;
        let features_contents = std::fs::read_to_string(&features_file).map_err(|err| {
            format!(
                "EXITING because of fatal INPUT FILE error: could not open {}\n{}",
                features_file, err
            )
        })?;
        let transcriptome_main = transcriptome_main_in.unwrap_or_default();
        let current_dir = std::env::current_dir()
            .map(|path| path.to_string_lossy().to_string())
            .unwrap_or_else(|_| ".".to_string());
        let solo_result = solo_l23_solo_solo(
            &p,
            &transcriptome_main,
            &matrix_contents,
            &barcodes_contents,
            &features_contents,
            &current_dir,
        )?;
        result.log_stdout.push_str(&solo_result.log_stdout);
        result.log_main.push_str(&solo_result.log_main);
        result.solo_cell_filtering = Some(solo_result);
        result.exit_code = 0;
        result.parameters = p;
        result.transcriptome = Some(transcriptome_main);
        return Ok(result);
    } else if run_mode != "alignReads" && run_mode != "soloCellFiltering" {
        result.log_main.push_str(&format!(
            "EXITING because of INPUT ERROR: unknown value of input parameter runMode={}\n",
            run_mode
        ));
        result.exit_code = 1;
        result.parameters = p;
        return Ok(result);
    }

    let mut transcriptome_main = transcriptome_main_in.unwrap_or_default();
    let mut genome_main = if let Some(genome_main) = genome_main_in {
        genome_main
    } else {
        crate::cli::load_genome_from_parameters(&mut p)?
    };

    let solo_cell_filter = solo_l23_solo_solo(&p, &transcriptome_main, "", "", "", ".")?;
    result.log_stdout.push_str(&solo_cell_filter.log_stdout);
    result.log_main.push_str(&solo_cell_filter.log_main);
    if solo_cell_filter.exited {
        result.exit_code = 0;
        result.parameters = p;
        result.genome = Some(genome_main);
        return Ok(result);
    }

    let mut sjdb_loci = crate::sjdb_class::SjdbClass::default();
    if p.sjdb_insert_yes {
        let genome_main1 = genome_sjdb_insert_snapshot(&genome_main);
        let sjdb = sjdbinsertjunctions_l11_sjdbinsertjunctions(
            &mut p,
            &mut genome_main,
            &genome_main1,
            &mut sjdb_loci,
        )?;
        result.log_main.push_str(&sjdb.log_main);
        result.sjdb_insert = Some(sjdb);
    }

    result
        .log_progress
        .push_str(&stats_l62_stats_progressreportheader());

    let two_pass = twopassrunpass1_l9_twopassrunpass1(
        &mut p,
        &mut genome_main,
        Some(&transcriptome_main),
        &mut sjdb_loci,
        pass1_chunks_in,
        existing_read_files,
    )?;
    result.log_progress.push_str(&two_pass.log_progress);
    result.log_stdout.push_str(&two_pass.log_stdout);
    result.log_main.push_str(&two_pass.log_main);
    if two_pass.sjdb_insert.is_some() || !two_pass.log_progress.is_empty() {
        result.two_pass = Some(two_pass);
    }

    let mut stats_all = crate::stats::Stats::default();
    stats_l4_stats_resetn(&mut stats_all);
    let raw_time_map = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    stats_all.time_start = raw_time_start;
    stats_all.time_start_map = raw_time_map;
    stats_all.time_last_report = raw_time_map;
    result.log_stdout.push_str(&format!(
        "{} ..... started mapping\n",
        timefunctions_l14_timemonthdaytime(raw_time_map)
    ));

    if p.quant_tr_sam_yes && transcriptome_main.tr_id.is_empty() {
        transcriptome_main = crate::transcriptome::Transcriptome::default();
    }
    if p.out_sam_type.is_empty() {
        p.out_sam_type.push("None".to_string());
    }
    let sam_header_comment_txt = if p.out_sam_header_comment_file != "-" {
        std::fs::read_to_string(&p.out_sam_header_comment_file).unwrap_or_default()
    } else {
        String::new()
    };
    samheaders_l5_samheaders(
        &mut p,
        &mut genome_main,
        &transcriptome_main,
        "",
        &sam_header_comment_txt,
    );
    let sam_header_for_chim = p.sam_header.clone();
    let mut p_ch = std::mem::take(&mut p.p_ch);
    parameterschimeric_initialize_l6_parameterschimeric_initialize(
        &mut p_ch,
        &mut p,
        &sam_header_for_chim,
    )?;
    p.p_ch = p_ch;

    let map_chunks_injected = map_chunks_in.is_some();
    let mut ra_chunks = if let Some(chunks) = map_chunks_in {
        chunks
    } else {
        let mut chunks = Vec::with_capacity(p.run_thread_n.max(0) as usize);
        for ii in 0..p.run_thread_n {
            chunks.push(readalignchunk_l5_readalignchunk_readalignchunk(
                &p,
                &genome_main,
                Some(&transcriptome_main),
                ii,
            )?);
        }
        chunks
    };

    let mut collected_signal_records = Vec::<crate::parameters_chimeric::SignalFromBamRecord>::new();
    let mut thread_chunks = crate::thread_control::ThreadControl::default();

    if p.run_restart_type != 1 {
        let input_mates = if map_chunks_injected {
            Vec::new()
        } else if !p.read_files_command_string.is_empty() {
            let mut input_mates = vec![String::new(); p.read_nends as usize];
            for imate in 0..p.read_nends as usize {
                for ifile in 0..p.read_files_n as usize {
                    let file_name = p
                        .read_files_names
                        .get(imate)
                        .and_then(|v| v.get(ifile))
                        .ok_or_else(|| {
                            "EXITING: because of fatal INPUT file error: could not open read file: \nSOLUTION: check that this file exists and has read permision.\n".to_string()
                        })?;
                    input_mates[imate].push_str(&format!("FILE {}\n", ifile));
                    let contents = if p.read_files_command.first().map(|s| s.as_str()) == Some("-")
                        && p.read_files_command_string.trim() == "cat"
                    {
                        crate::io_utils::read_to_string_auto_gzip(file_name).map_err(|_| {
                            format!(
                                "EXITING: because of fatal INPUT file error: could not open read file: {}\nSOLUTION: check that this file exists and has read permision.\n",
                                file_name
                            )
                        })?
                    } else {
                        let command_words = if p.read_files_command.first().map(|s| s.as_str())
                            == Some("-")
                            && p.read_files_command_string.trim() == "cat"
                        {
                            vec!["cat".to_string()]
                        } else {
                            p.read_files_command.clone()
                        };
                        if command_words.is_empty() {
                            crate::io_utils::read_to_string_auto_gzip(file_name).map_err(|_| {
                                format!(
                                    "EXITING because of fatal input ERROR: could not open readFilesIn={}\n",
                                    file_name
                                )
                            })?
                        } else {
                            let output = std::process::Command::new(&command_words[0])
                                .args(&command_words[1..])
                                .arg(file_name)
                                .output()
                                .map_err(|e| {
                                    format!(
                                        "EXITING: because of fatal EXECUTION error: failed to execute readFilesCommand {} for {}\n{}\n",
                                        command_words.join(" "),
                                        file_name,
                                        e
                                    )
                                })?;
                            if !output.status.success() {
                                return Err(format!(
                                    "EXITING: because of fatal EXECUTION error: readFilesCommand {} failed for {}\n",
                                    command_words.join(" "),
                                    file_name
                                ));
                            }
                            String::from_utf8_lossy(&output.stdout).into_owned()
                        }
                    };
                    input_mates[imate].push_str(&contents);
                    if !contents.ends_with('\n') {
                        input_mates[imate].push('\n');
                    }
                }
            }
            input_mates
        } else {
            let mut input_mates = Vec::new();
            for imate in 0..p.read_nends as usize {
                let rf_name = p
                    .read_files_names
                    .get(imate)
                    .and_then(|v| v.first())
                    .cloned()
                    .unwrap_or_else(|| {
                        format!(
                            "{}{}",
                            p.read_files_prefix_final,
                            p.read_files_in.get(imate).cloned().unwrap_or_default()
                        )
                    });
                input_mates.push(crate::io_utils::read_to_string_auto_gzip(&rf_name).map_err(
                    |_| {
                        format!(
                            "EXITING because of fatal input ERROR: could not open readFilesIn={}\n",
                            rf_name
                        )
                    },
                )?);
            }
            input_mates
        };
        let mut process_chunks_result =
            Option::<crate::read_align_chunk::ReadAlignChunkProcessChunksResult>::None;
        result
            .log_main
            .push_str(&mapthreadsspawn_l6_mapthreadsspawn(
                p.run_thread_n,
                &vec![0; p.run_thread_n.max(0) as usize],
                &vec![0; p.run_thread_n.max(0) as usize],
                || {
                    if map_chunks_injected || ra_chunks.is_empty() {
                        return Ok(String::new());
                    }
                    let process = readalignchunk_processchunks_l11_readalignchunk_processchunks(
                        &mut ra_chunks[0],
                        &mut p,
                        &mut thread_chunks,
                        &mut stats_all,
                        raw_time_map,
                        &input_mates,
                        |_ra| -1,
                        Some((&genome_main, &mut transcriptome_main)),
                    )?;
                    let log_main = process.log_main.clone();
                    process_chunks_result = Some(process);
                    Ok(log_main)
                },
            )?);
        if let Some(process) = process_chunks_result {
            if p.out_sam_bool {
                for map_chunk in &process.map_chunks {
                    p.out_sam_contents
                        .push_str(&String::from_utf8_lossy(&map_chunk.direct_sam_output));
                }
            }
            for map_chunk in &process.map_chunks {
                collected_signal_records.extend(map_chunk.signal_records.iter().cloned());
            }
            if let Some(chunk) = ra_chunks.get(0) {
                stats_all.read_n += chunk.ra.stats_ra.read_n;
                stats_all.read_bases += chunk.ra.stats_ra.read_bases;
                stats_all.mapped_reads_u += chunk.ra.stats_ra.mapped_reads_u;
                stats_all.mapped_reads_m += chunk.ra.stats_ra.mapped_reads_m;
                stats_all.unmapped_mismatch += chunk.ra.stats_ra.unmapped_mismatch;
                stats_all.unmapped_short += chunk.ra.stats_ra.unmapped_short;
                stats_all.unmapped_other += chunk.ra.stats_ra.unmapped_other;
            }
            result.process_chunks.push(process);
        }
    }

    if p.out_filter_by_sjout_stage == 1 {
        result
            .log_main
            .push_str("Completed stage 1 mapping of outFilterBySJout mapping\n");
        let sj1 = outputsj_l20_outputsj(&ra_chunks, &mut p, &genome_main)?;
        result.log_main.push_str(&sj1.log_main);
        p.read_files_index = u32::MAX;
        p.out_filter_by_sjout_stage = 2;
        if p.out_bam_coord {
            for chunk in &mut ra_chunks {
                bamoutput_l179_bamoutput_coordunmappedpreparebysjout(
                    &mut chunk.chunk_out_bam_coord,
                    &p,
                );
            }
        }
        result
            .log_main
            .push_str(&mapthreadsspawn_l6_mapthreadsspawn(
                p.run_thread_n,
                &vec![0; p.run_thread_n.max(0) as usize],
                &vec![0; p.run_thread_n.max(0) as usize],
                || {
                    if map_chunks_injected || ra_chunks.is_empty() {
                        return Ok(String::new());
                    }
                    let process = readalignchunk_processchunks_l11_readalignchunk_processchunks(
                        &mut ra_chunks[0],
                        &mut p,
                        &mut thread_chunks,
                        &mut stats_all,
                        raw_time_map,
                        &[],
                        |_ra| -1,
                        Some((&genome_main, &mut transcriptome_main)),
                    )?;
                    let log_main = process.log_main.clone();
                    if p.out_sam_bool {
                        for map_chunk in &process.map_chunks {
                            p.out_sam_contents
                                .push_str(&String::from_utf8_lossy(&map_chunk.direct_sam_output));
                        }
                    }
                    for map_chunk in &process.map_chunks {
                        collected_signal_records.extend(map_chunk.signal_records.iter().cloned());
                    }
                    result.process_chunks.push(process);
                    Ok(log_main)
                },
            )?);
    }

    if p.out_bam_coord && p.limit_bam_sort_ram == 0 {
        p.limit_bam_sort_ram = genome_main.n_genome
            + genome_main.sa_packed.length_byte
            + genome_main.sai_packed.length_byte;
    }
    let raw_time_finish_map = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    stats_all.time_finish_map = raw_time_finish_map;
    result.log_stdout.push_str(&format!(
        "{} ..... finished mapping\n",
        timefunctions_l14_timemonthdaytime(raw_time_finish_map)
    ));
    result.log_main.push_str(&format!(
        "{} ..... finished mapping\nRAM after mapping:\n{}\n",
        timefunctions_l14_timemonthdaytime(raw_time_finish_map),
        systemfunctions_l6_linuxprocmemory()
    ));

    genome_l33_genome_freememory(&mut genome_main);

    if p.run_restart_type != 1 && p.out_sj {
        let sj = outputsj_l20_outputsj(&ra_chunks, &mut p, &genome_main)?;
        result.log_main.push_str(&sj.log_main);
        result.output_sj = Some(sj);
    }

    let mut solo_main = solo_l5_solo_solo(&p, &transcriptome_main);
    if solo_main.p_solo.solo_type != SOLO_TYPE_NONE {
        let count_p = p.clone();
        let count_p_solo = p.p_solo.clone();
        let quant_p_solo = p.p_solo.clone();
        let quant_transcriptome = transcriptome_main.clone();
        let quant_cluster_contents = if p.p_solo.cluster_cb_file != "-" {
            std::fs::read_to_string(&p.p_solo.cluster_cb_file).unwrap_or_default()
        } else {
            String::new()
        };
        let raw_time_solo = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|duration| duration.as_secs() as libc::time_t)
            .unwrap_or(raw_time_finish_map);
        let time_counting_start = timefunctions_l4_timemonthdaytime();
        let time_counting_finish = timefunctions_l4_timemonthdaytime();
        let time_process_start = timefunctions_l4_timemonthdaytime();
        let time_writing_raw_matrix = timefunctions_l4_timemonthdaytime();
        let time_cell_filtering = timefunctions_l4_timemonthdaytime();
        let time_finished_redistribution = timefunctions_l4_timemonthdaytime();
        let time_finished_collapsing = timefunctions_l4_timemonthdaytime();
        let time_finished_counting = timefunctions_l4_timemonthdaytime();
        let linux_proc_memory = systemfunctions_l6_linuxprocmemory();
        let solo_process = solo_l48_solo_processandoutput(
            &mut solo_main,
            &mut p,
            &transcriptome_main,
            &mut ra_chunks,
            &stats_all,
            ".",
            &[],
            &time_counting_start,
            &time_counting_finish,
            &time_process_start,
            &time_writing_raw_matrix,
            &time_cell_filtering,
            &time_finished_redistribution,
            &time_finished_collapsing,
            &time_finished_counting,
            &linux_proc_memory,
            |_, solo_feature| {
                solofeature_countcbgeneumi_l7_solofeature_countcbgeneumi(
                    solo_feature,
                    &count_p,
                    &count_p_solo,
                    raw_time_solo,
                )
                .unwrap_or_else(|err| err)
            },
            |_, solo_feature| {
                let quant = solofeature_quanttranscript_l12_solofeature_quanttranscript(
                    solo_feature,
                    &quant_p_solo,
                    &quant_transcriptome,
                    count_p.run_thread_n,
                    &quant_cluster_contents,
                    &time_process_start,
                    &time_cell_filtering,
                    &time_finished_counting,
                );
                quant.log_main
            },
        )?;
        result.log_stdout.push_str(&solo_process.log_stdout);
        result.log_main.push_str(&solo_process.log_main);
        result.solo_process_and_output = Some(solo_process);
    }

    if p.quant_ge_count_yes && map_chunks_injected {
        let mut summed_transcriptome = None;
        if !ra_chunks.is_empty() {
            let (first_chunk, rest_chunks) = ra_chunks.split_at_mut(1);
            if let Some(first_chunk_tr) = first_chunk[0].chunk_tr.as_mut() {
                for chunk in rest_chunks
                    .iter()
                    .take(p.run_thread_n.max(0).saturating_sub(1) as usize)
                {
                    if let Some(chunk_tr) = &chunk.chunk_tr {
                        quantifications_l25_quantifications_addquants(
                            &mut first_chunk_tr.quants,
                            &chunk_tr.quants,
                        );
                    }
                }
                summed_transcriptome = Some(first_chunk_tr.clone());
            }
        }
        if let Some(first_chunk_tr) = summed_transcriptome {
            transcriptome_main = first_chunk_tr;
        }
    }

    if p.run_thread_n > 1 && p.out_sam_order == "PairedKeepInputOrder" {
        let mut chunk_cat = Vec::new();
        let mut i_c = 0_u32;
        readalignchunk_l151_readalignchunk_chunkfilescat(
            &mut chunk_cat,
            &format!("{}/Aligned.out.sam.chunk", p.out_file_tmp),
            &mut i_c,
        )
        .map_err(|err| err.to_string())?;
        p.out_sam_contents
            .push_str(&String::from_utf8_lossy(&chunk_cat));
    }

    let mut chunk_bam_bin_files = Vec::<Vec<Vec<u8>>>::new();
    let bam_temp_files = if temp_files_by_bam_bin.is_empty() && p.out_bam_coord {
        let n_bins = p.out_bam_coord_nbins as usize;
        chunk_bam_bin_files.resize_with(n_bins, Vec::new);
        for ibin in 0..n_bins {
            for chunk in ra_chunks.iter().take(p.run_thread_n as usize) {
                let stream = chunk
                    .chunk_out_bam_coord
                    .bin_streams
                    .get(ibin)
                    .cloned()
                    .unwrap_or_default();
                chunk_bam_bin_files[ibin].push(stream);
                if ibin == n_bins - 1 {
                    chunk_bam_bin_files[ibin].push(Vec::new());
                }
            }
        }
        chunk_bam_bin_files.as_slice()
    } else {
        temp_files_by_bam_bin
    };
    let bam_sort = bamsortbycoordinate_l8_bamsortbycoordinate(
        &mut p,
        &ra_chunks,
        &genome_main,
        bam_temp_files,
    )?;
    if p.out_bam_coord {
        result.bam_sort = Some(bam_sort);
    }

    if p.out_wig_flags.yes {
        result.log_stdout.push_str(&format!(
            "{} ..... started wiggle output\n",
            timefunctions_l4_timemonthdaytime()
        ));
        result.log_main.push_str(&format!(
            "{} ..... started wiggle output\n",
            timefunctions_l4_timemonthdaytime()
        ));
        let signal = signalfrombam_l5_signalfrombam(
            &format!("{}Signal", p.out_file_name_prefix),
            &p,
            &genome_main.chr_name_all,
            &genome_main.chr_length_all,
            if signal_records.is_empty() {
                &collected_signal_records
            } else {
                signal_records
            },
        )?;
        result.signal = Some(signal);
    }

    if p.p_ch.out_chim_sam_opened {
        for process_chunk in &result.process_chunks {
            for map_chunk in &process_chunk.map_chunks {
                p.p_ch
                    .out_chim_sam_contents
                    .push_str(&map_chunk.chimeric_sam_output);
            }
        }
    }

    if p.p_ch.out_chim_junction_opened {
        for process_chunk in &result.process_chunks {
            for map_chunk in &process_chunk.map_chunks {
                p.p_ch
                    .out_chim_junction_contents
                    .push_str(&map_chunk.chimeric_junction_output);
            }
        }
        p.p_ch
            .out_chim_junction_contents
            .push_str(&stats_l147_stats_writelines(
                &stats_all,
                &p.p_ch.out_junction_format,
                "#",
                &format!("2.7.11b   {}", p.command_line),
            ));
    }

    if let Some(progress) = stats_l73_stats_progressreport(&mut stats_all, raw_time_finish_map + 60)
    {
        result.log_progress.push_str(&progress);
    }
    result.log_progress.push_str("ALL DONE!\n");
    result.log_final_out = stats_l99_stats_reportfinal(&mut stats_all, raw_time_finish_map);
    result.log_stdout.push_str(&format!(
        "{} ..... finished successfully\n",
        timefunctions_l14_timemonthdaytime(stats_all.time_finish)
    ));
    result.log_main.push_str("ALL DONE!\n");

    if p.out_tmp_keep == "None" && !p.out_file_tmp.is_empty() {
        if sysremovedir_l25_sysremovedir(std::path::Path::new(&p.out_file_tmp)).is_ok() {
            result.removed_tmp = true;
        }
    }
    result.killed_read_command_pids =
        parameters_closereadsfiles_l5_parameters_closereadsfiles(&mut p);
    result.exit_code = 0;
    result.parameters = p;
    result.genome = Some(genome_main);
    result.transcriptome = Some(transcriptome_main);
    result.stats_all = stats_all;
    result.read_chunks = ra_chunks;
    Ok(result)
}
