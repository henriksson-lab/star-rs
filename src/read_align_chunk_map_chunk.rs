#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlignChunk::mapChunk` at STAR/source/ReadAlignChunk_mapChunk.cpp:7. Args: "]
pub fn readalignchunk_mapchunk_l7_readalignchunk_mapchunk<F>(
    chunk: &mut crate::read_align_chunk::ReadAlignChunk,
    p: &crate::parameters_chimeric::Parameters,
    stats_all: &mut crate::stats::Stats,
    time_current: libc::time_t,
    one_read: F,
    real_context: Option<(
        &crate::genome::Genome,
        &mut crate::transcriptome::Transcriptome,
    )>,
) -> Result<crate::read_align_chunk::ReadAlignChunkMapChunkResult, String>
where
    F: FnMut(&mut crate::read_align::ReadAlign) -> i32,
{
    if !chunk.ra.clip_mates.is_empty()
        && !chunk.ra.clip_mates[0].is_empty()
        && !chunk.chunk_in.is_empty()
    {
        let ch_size = chunk
            .chunk_in_size_bytes_total
            .first()
            .copied()
            .unwrap_or(chunk.chunk_in[0].len() as u64);
        clipmate_clipchunk_l7_clipmate_clipchunk(
            &mut chunk.ra.clip_mates[0][0],
            &mut chunk.chunk_in[0],
            ch_size,
        )?;
    }

    let chunk_inputs: Vec<&[u8]> = chunk
        .chunk_in
        .iter()
        .map(|input| {
            let ptr = input.as_ptr();
            let len = input.len();
            // SAFETY: the cursors below only read from chunk_in while mapChunk
            // mutates disjoint output/alignment fields. The existing STAR chunk
            // input buffers are not resized or written after this point.
            unsafe { std::slice::from_raw_parts(ptr, len) }
        })
        .collect();
    let mut real_read_inputs: Vec<std::io::Cursor<&[u8]>> = chunk_inputs
        .iter()
        .map(|input| std::io::Cursor::new(*input))
        .collect();
    let mut read_streams: Vec<&mut dyn std::io::BufRead> = real_read_inputs
        .iter_mut()
        .map(|input| input as &mut dyn std::io::BufRead)
        .collect();
    readalignchunk_mapchunk_l7_readalignchunk_mapchunk_with_read_streams(
        chunk,
        p,
        stats_all,
        time_current,
        one_read,
        real_context,
        &mut read_streams,
    )
}

#[doc = "STAR direct API helper: map a chunk using caller-supplied FASTQ-like read streams instead of chunk_in buffers."]
pub fn readalignchunk_mapchunk_l7_readalignchunk_mapchunk_with_read_streams<F>(
    chunk: &mut crate::read_align_chunk::ReadAlignChunk,
    p: &crate::parameters_chimeric::Parameters,
    stats_all: &mut crate::stats::Stats,
    time_current: libc::time_t,
    one_read: F,
    real_context: Option<(
        &crate::genome::Genome,
        &mut crate::transcriptome::Transcriptome,
    )>,
    read_streams: &mut [&mut dyn std::io::BufRead],
) -> Result<crate::read_align_chunk::ReadAlignChunkMapChunkResult, String>
where
    F: FnMut(&mut crate::read_align::ReadAlign) -> i32,
{
    readalignchunk_mapchunk_l7_readalignchunk_mapchunk_with_next_read(
        chunk,
        p,
        stats_all,
        time_current,
        one_read,
        real_context,
        |ra,
         p_one_read,
         map_gen,
         transcriptome,
         pe_merge_ra,
         wasp_ra,
         chunk_out_sj,
         chunk_out_sj1,
         chunk_out_filter_by_sjout_files,
         chunk_out_unmapped_reads_stream,
         out_sam_stream| {
            readalign_oneread_l8_readalign_oneread(
                ra,
                read_streams,
                p_one_read,
                map_gen,
                transcriptome,
                None,
                None,
                None,
                pe_merge_ra,
                wasp_ra,
                None,
                &[],
                chunk_out_sj,
                chunk_out_sj1,
                chunk_out_filter_by_sjout_files,
                chunk_out_unmapped_reads_stream,
                out_sam_stream,
                0.0,
                None,
            )
        },
    )
}

pub fn readalignchunk_mapchunk_l7_readalignchunk_mapchunk_with_next_read<F, N>(
    chunk: &mut crate::read_align_chunk::ReadAlignChunk,
    p: &crate::parameters_chimeric::Parameters,
    stats_all: &mut crate::stats::Stats,
    time_current: libc::time_t,
    mut one_read: F,
    mut real_context: Option<(
        &crate::genome::Genome,
        &mut crate::transcriptome::Transcriptome,
    )>,
    mut next_one_read: N,
) -> Result<crate::read_align_chunk::ReadAlignChunkMapChunkResult, String>
where
    F: FnMut(&mut crate::read_align::ReadAlign) -> i32,
    N: FnMut(
        &mut crate::read_align::ReadAlign,
        &mut crate::parameters_chimeric::Parameters,
        &crate::genome::Genome,
        &mut crate::transcriptome::Transcriptome,
        &mut crate::read_align::ReadAlign,
        &mut crate::read_align::ReadAlign,
        &mut crate::out_sj::OutSJ,
        &mut crate::out_sj::OutSJ,
        &mut [String],
        &mut [String],
        &mut String,
    ) -> Result<crate::quantifications::ReadAlignOneReadResult, String>,
{
    stats_l4_stats_resetn(&mut chunk.ra.stats_ra);

    let paired_keep_input_order = p.out_sam_order == "PairedKeepInputOrder" && p.run_thread_n > 1;
    let mut result = crate::read_align_chunk::ReadAlignChunkMapChunkResult::default();
    if paired_keep_input_order {
        chunk.chunk_out_bam_file_name = format!(
            "{}/Aligned.tmp.sam.chunk{}",
            p.out_file_tmp, chunk.i_chunk_in
        );
        result.paired_keep_input_order_tmp_name = Some(chunk.chunk_out_bam_file_name.clone());
    }

    let mut read_status = 0;
    let mut pe_merge_ra = crate::read_align::ReadAlign::default();
    let mut wasp_ra = crate::read_align::ReadAlign::default();
    let mut chunk_out_unmapped_reads_stream = vec![String::new(); p.read_nends as usize];
    let mut p_one_read = p.clone();
    while read_status == 0 {
        chunk.ra.out_bam_bytes = 0;
        if let Some((map_gen, transcriptome)) = real_context.as_mut() {
            if chunk.chunk_out_filter_by_sjout_files.len() < p.read_nends as usize {
                chunk
                    .chunk_out_filter_by_sjout_files
                    .resize(p.read_nends as usize, String::new());
            }
            let mut out_sam_stream = String::new();
            let one_read_result = next_one_read(
                &mut chunk.ra,
                &mut p_one_read,
                *map_gen,
                &mut **transcriptome,
                &mut pe_merge_ra,
                &mut wasp_ra,
                &mut chunk.chunk_out_sj,
                &mut chunk.chunk_out_sj1,
                &mut chunk.chunk_out_filter_by_sjout_files,
                &mut chunk_out_unmapped_reads_stream,
                &mut out_sam_stream,
            )?;
            if one_read_result.status == 0 {
                if let Some(chimeric_detection) = &one_read_result.pe_overlap.chimeric_detection {
                    if let Some(old_output) = &chimeric_detection.old_output {
                        result.chimeric_sam_output.push_str(&old_output.chim_sam);
                        result
                            .chimeric_junction_output
                            .push_str(&old_output.chim_junction);
                        for request in &old_output.bam_requests {
                            let chim_bam =
                                chimericalign_chimericbamoutput_l7_chimericalign_chimericbamoutput(
                                    &request.al1,
                                    &request.al2,
                                    &chunk.ra,
                                    map_gen,
                                    request.i_tr,
                                    request.chim_n,
                                    request.is_best_chim_align,
                                    &p_one_read,
                                );
                            for bam_request in &chim_bam.bam_requests {
                                let bam = readalign_alignbam_l47_readalign_alignbam(
                                    &mut chunk.ra,
                                    &p_one_read,
                                    map_gen,
                                    &bam_request.transcript,
                                    bam_request.n_tr_out,
                                    bam_request.i_tr_out,
                                    bam_request.tr_chr_start,
                                    bam_request.mate_chr,
                                    bam_request.mate_start,
                                    bam_request.mate_strand,
                                    bam_request.align_type,
                                    bam_request.mate_map,
                                    &p_one_read.out_sam_attr_order,
                                )?;
                                for iline in 0..bam.n_lines as usize {
                                    let bam_size = bam.record_sizes[iline] as u64;
                                    if bam_size == 0 {
                                        continue;
                                    }
                                    if let Some(out_bam_unsorted) =
                                        chunk.chunk_out_bam_unsorted.as_mut()
                                    {
                                        bamoutput_l52_bamoutput_unsortedonealign(
                                            out_bam_unsorted,
                                            &bam.records[iline],
                                            bam_size,
                                            bam_size,
                                        )?;
                                    }
                                    if p_one_read.out_bam_coord {
                                        bamoutput_l77_bamoutput_coordonealign(
                                            &mut chunk.chunk_out_bam_coord,
                                            &mut p_one_read,
                                            &bam.records[iline],
                                            bam_size,
                                            chunk.ra.i_read_all,
                                        )?;
                                    }
                                }
                            }
                        }
                    }
                    if let Some(mult_output) = &chimeric_detection.mult_output {
                        result
                            .chimeric_junction_output
                            .push_str(&mult_output.chim_junction);
                        for chim_bam in &mult_output.bam_outputs {
                            for bam_request in &chim_bam.bam_requests {
                                let bam = readalign_alignbam_l47_readalign_alignbam(
                                    &mut chunk.ra,
                                    &p_one_read,
                                    map_gen,
                                    &bam_request.transcript,
                                    bam_request.n_tr_out,
                                    bam_request.i_tr_out,
                                    bam_request.tr_chr_start,
                                    bam_request.mate_chr,
                                    bam_request.mate_start,
                                    bam_request.mate_strand,
                                    bam_request.align_type,
                                    bam_request.mate_map,
                                    &p_one_read.out_sam_attr_order,
                                )?;
                                for iline in 0..bam.n_lines as usize {
                                    let bam_size = bam.record_sizes[iline] as u64;
                                    if bam_size == 0 {
                                        continue;
                                    }
                                    if let Some(out_bam_unsorted) =
                                        chunk.chunk_out_bam_unsorted.as_mut()
                                    {
                                        bamoutput_l52_bamoutput_unsortedonealign(
                                            out_bam_unsorted,
                                            &bam.records[iline],
                                            bam_size,
                                            bam_size,
                                        )?;
                                    }
                                    if p_one_read.out_bam_coord {
                                        bamoutput_l77_bamoutput_coordonealign(
                                            &mut chunk.chunk_out_bam_coord,
                                            &mut p_one_read,
                                            &bam.records[iline],
                                            bam_size,
                                            chunk.ra.i_read_all,
                                        )?;
                                    }
                                }
                            }
                        }
                    }
                }
                if let Some(chimeric_detection) = &one_read_result.chimeric_detection {
                    if let Some(old_output) = &chimeric_detection.old_output {
                        result.chimeric_sam_output.push_str(&old_output.chim_sam);
                        result
                            .chimeric_junction_output
                            .push_str(&old_output.chim_junction);
                        for request in &old_output.bam_requests {
                            let chim_bam =
                                chimericalign_chimericbamoutput_l7_chimericalign_chimericbamoutput(
                                    &request.al1,
                                    &request.al2,
                                    &chunk.ra,
                                    map_gen,
                                    request.i_tr,
                                    request.chim_n,
                                    request.is_best_chim_align,
                                    &p_one_read,
                                );
                            for bam_request in &chim_bam.bam_requests {
                                let bam = readalign_alignbam_l47_readalign_alignbam(
                                    &mut chunk.ra,
                                    &p_one_read,
                                    map_gen,
                                    &bam_request.transcript,
                                    bam_request.n_tr_out,
                                    bam_request.i_tr_out,
                                    bam_request.tr_chr_start,
                                    bam_request.mate_chr,
                                    bam_request.mate_start,
                                    bam_request.mate_strand,
                                    bam_request.align_type,
                                    bam_request.mate_map,
                                    &p_one_read.out_sam_attr_order,
                                )?;
                                for iline in 0..bam.n_lines as usize {
                                    let bam_size = bam.record_sizes[iline] as u64;
                                    if bam_size == 0 {
                                        continue;
                                    }
                                    if let Some(out_bam_unsorted) =
                                        chunk.chunk_out_bam_unsorted.as_mut()
                                    {
                                        bamoutput_l52_bamoutput_unsortedonealign(
                                            out_bam_unsorted,
                                            &bam.records[iline],
                                            bam_size,
                                            bam_size,
                                        )?;
                                    }
                                    if p_one_read.out_bam_coord {
                                        bamoutput_l77_bamoutput_coordonealign(
                                            &mut chunk.chunk_out_bam_coord,
                                            &mut p_one_read,
                                            &bam.records[iline],
                                            bam_size,
                                            chunk.ra.i_read_all,
                                        )?;
                                    }
                                }
                            }
                        }
                    }
                    if let Some(mult_output) = &chimeric_detection.mult_output {
                        result
                            .chimeric_junction_output
                            .push_str(&mult_output.chim_junction);
                        for chim_bam in &mult_output.bam_outputs {
                            for bam_request in &chim_bam.bam_requests {
                                let bam = readalign_alignbam_l47_readalign_alignbam(
                                    &mut chunk.ra,
                                    &p_one_read,
                                    map_gen,
                                    &bam_request.transcript,
                                    bam_request.n_tr_out,
                                    bam_request.i_tr_out,
                                    bam_request.tr_chr_start,
                                    bam_request.mate_chr,
                                    bam_request.mate_start,
                                    bam_request.mate_strand,
                                    bam_request.align_type,
                                    bam_request.mate_map,
                                    &p_one_read.out_sam_attr_order,
                                )?;
                                for iline in 0..bam.n_lines as usize {
                                    let bam_size = bam.record_sizes[iline] as u64;
                                    if bam_size == 0 {
                                        continue;
                                    }
                                    if let Some(out_bam_unsorted) =
                                        chunk.chunk_out_bam_unsorted.as_mut()
                                    {
                                        bamoutput_l52_bamoutput_unsortedonealign(
                                            out_bam_unsorted,
                                            &bam.records[iline],
                                            bam_size,
                                            bam_size,
                                        )?;
                                    }
                                    if p_one_read.out_bam_coord {
                                        bamoutput_l77_bamoutput_coordonealign(
                                            &mut chunk.chunk_out_bam_coord,
                                            &mut p_one_read,
                                            &bam.records[iline],
                                            bam_size,
                                            chunk.ra.i_read_all,
                                        )?;
                                    }
                                }
                            }
                        }
                    }
                }
                if let Some(output_alignments) = &one_read_result.output_alignments {
                    for request in &output_alignments.write_sam.bam_requests {
                        let bam = readalign_alignbam_l47_readalign_alignbam(
                            &mut chunk.ra,
                            &p_one_read,
                            map_gen,
                            &request.transcript,
                            request.n_tr_out,
                            request.i_tr_out,
                            request.tr_chr_start,
                            request.mate_chr,
                            request.mate_start,
                            request.mate_strand,
                            request.align_type,
                            request.mate_map,
                            &p_one_read.out_sam_attr_order,
                        )?;
                        for iline in 0..bam.n_lines as usize {
                            let bam_size = bam.record_sizes[iline] as u64;
                            if bam_size == 0 {
                                continue;
                            }
                            if let Some(out_bam_unsorted) = chunk.chunk_out_bam_unsorted.as_mut() {
                                bamoutput_l52_bamoutput_unsortedonealign(
                                    out_bam_unsorted,
                                    &bam.records[iline],
                                    bam_size,
                                    bam_size,
                                )?;
                            }
                            if p_one_read.out_bam_coord {
                                bamoutput_l77_bamoutput_coordonealign(
                                    &mut chunk.chunk_out_bam_coord,
                                    &mut p_one_read,
                                    &bam.records[iline],
                                    bam_size,
                                    chunk.ra.i_read_all,
                                )?;
                            }
                            if p_one_read.out_wig_flags.yes
                                && let Some(signal_record) = bam.signal_records[iline].clone()
                            {
                                result.signal_records.push(signal_record);
                            }
                        }
                    }
                    if let Some(quant_transcriptome) = &output_alignments.quant_transcriptome {
                        for request in &quant_transcriptome.bam_requests {
                            let bam = readalign_alignbam_l47_readalign_alignbam(
                                &mut chunk.ra,
                                &p_one_read,
                                map_gen,
                                &request.transcript,
                                request.n_align_t,
                                request.i_align_t,
                                0,
                                u64::MAX,
                                u64::MAX,
                                0,
                                -1,
                                None,
                                &p_one_read.out_sam_attr_order_quant,
                            )?;
                            for iline in 0..bam.n_lines as usize {
                                let bam_size = bam.record_sizes[iline] as u64;
                                if bam_size == 0 {
                                    continue;
                                }
                                if let Some(out_bam_quant) = chunk.chunk_out_bam_quant.as_mut() {
                                    bamoutput_l52_bamoutput_unsortedonealign(
                                        out_bam_quant,
                                        &bam.records[iline],
                                        bam_size,
                                        bam_size,
                                    )?;
                                }
                            }
                        }
                    }
                }
            }
            if !out_sam_stream.is_empty() {
                let bytes = out_sam_stream.as_bytes();
                let start = chunk.chunk_out_bam_total as usize;
                let end = start + bytes.len();
                if chunk.chunk_out_bam.len() < end {
                    chunk.chunk_out_bam.resize(end, 0);
                }
                chunk.chunk_out_bam[start..end].copy_from_slice(bytes);
                chunk.ra.out_bam_bytes = bytes.len() as u64;
            }
            read_status = one_read_result.status;
        } else {
            read_status = one_read(&mut chunk.ra);
        }

        if read_status == 0 {
            chunk.ra.i_read += 1;
            result.reads_processed += 1;
            chunk.chunk_out_bam_total += chunk.ra.out_bam_bytes;
        }

        if p.out_sam_bool {
            if chunk.chunk_out_bam_total > p.chunk_out_bam_size_bytes {
                return Err("EXITING because of fatal error: buffer size for SAM/BAM output is too small\nSolution: increase input parameter --limitOutSAMoneReadBytes\n".to_string());
            } else if chunk.chunk_out_bam_total + p.limit_out_sam_one_read_bytes
                > p.chunk_out_bam_size_bytes
                || (read_status == -1 && chunk.no_reads_left)
            {
                let bytes = chunk.chunk_out_bam_total as usize;
                let out = &chunk.chunk_out_bam[..bytes.min(chunk.chunk_out_bam.len())];
                if paired_keep_input_order {
                    result.paired_keep_input_order_tmp.extend_from_slice(out);
                } else {
                    result.direct_sam_output.extend_from_slice(out);
                }
                chunk.chunk_out_bam_total = 0;
            }
        }

        if p.out_sj {
            if chunk.chunk_out_sj.n > chunk.chunk_out_sj.n_store {
                return Err("EXITING because of fatal error: buffer size for SJ output is too small\nSolution: increase input parameter --limitOutSJoneRead\n".to_string());
            } else if chunk.chunk_out_sj.n + p.limit_out_sj_one_read > chunk.chunk_out_sj.n_store
                || (read_status == -1 && chunk.no_reads_left)
            {
                outsj_l36_outsj_collapsesj(&mut chunk.chunk_out_sj)?;
                if chunk.chunk_out_sj.n + 2 * p.limit_out_sj_one_read > chunk.chunk_out_sj.n_store {
                    outsj_l62_outsj_datasizeincrease(&mut chunk.chunk_out_sj);
                    result.log_main.push_str(&format!(
                        "Increased the size of chunkOutSJ to {}\n",
                        chunk.chunk_out_sj.n_store
                    ));
                }
            }
        }

        if p.out_filter_by_sjout_stage == 1 {
            if chunk.chunk_out_sj1.n > chunk.chunk_out_sj.n_store {
                return Err("EXITING because of fatal error: buffer size for SJ output is too small\nSolution: increase input parameter --limitOutSJoneRead\n".to_string());
            } else if chunk.chunk_out_sj1.n + p.limit_out_sj_one_read > chunk.chunk_out_sj.n_store
                || (read_status == -1 && chunk.no_reads_left)
            {
                outsj_l36_outsj_collapsesj(&mut chunk.chunk_out_sj1)?;
                if chunk.chunk_out_sj1.n + 2 * p.limit_out_sj_one_read > chunk.chunk_out_sj.n_store
                {
                    outsj_l62_outsj_datasizeincrease(&mut chunk.chunk_out_sj);
                    result.log_main.push_str(&format!(
                        "Increased the size of chunkOutSJ to {}\n",
                        chunk.chunk_out_sj.n_store
                    ));
                }
            }
        }
    }

    if p.out_sam_bool && paired_keep_input_order {
        let bytes = chunk.chunk_out_bam_total as usize;
        let out = &chunk.chunk_out_bam[..bytes.min(chunk.chunk_out_bam.len())];
        result.paired_keep_input_order_tmp.extend_from_slice(out);
        chunk.chunk_out_bam_total = 0;
        result.paired_keep_input_order_final_name = Some(format!(
            "{}/Aligned.out.sam.chunk{}",
            p.out_file_tmp, chunk.i_chunk_in
        ));
        if let Some(parent) = std::path::Path::new(&chunk.chunk_out_bam_file_name).parent()
            && !parent.as_os_str().is_empty()
        {
            std::fs::create_dir_all(parent).map_err(|err| err.to_string())?;
        }
        std::fs::write(
            &chunk.chunk_out_bam_file_name,
            &result.paired_keep_input_order_tmp,
        )
        .map_err(|err| err.to_string())?;
        if let Some(final_name) = &result.paired_keep_input_order_final_name {
            std::fs::rename(&chunk.chunk_out_bam_file_name, final_name)
                .map_err(|err| err.to_string())?;
        }
    }
    if let Some(out_bam_unsorted) = chunk.chunk_out_bam_unsorted.as_mut() {
        bamoutput_l70_bamoutput_unsortedflush(out_bam_unsorted);
    }
    if let Some(out_bam_quant) = chunk.chunk_out_bam_quant.as_mut() {
        bamoutput_l70_bamoutput_unsortedflush(out_bam_quant);
    }
    if p.out_bam_coord {
        // p_one_read carries the bin-sorting mutations accumulated during the
        // per-record coordOneAlign calls above; flush against that state, not
        // a fresh clone of p (which would discard them).
        bamoutput_l168_bamoutput_coordflush(&mut chunk.chunk_out_bam_coord, &mut p_one_read)?;
    }
    if p.out_reads_unmapped == "Fastx" {
        result.unmapped_fastx_outputs = chunk_out_unmapped_reads_stream;
    }
    if let Some(out_bam_quant) = &chunk.chunk_out_bam_quant {
        result.quant_bam_output = out_bam_quant.bgzf_bam.clone();
    }

    stats_l21_stats_addstats(stats_all, &chunk.ra.stats_ra);
    result.progress_report = stats_l73_stats_progressreport(stats_all, time_current);
    Ok(result)
}
