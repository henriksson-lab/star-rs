use std::collections::BTreeSet;
use std::fs;

use crate::cli::{PARAMETERS_DEFAULT, load_genome_from_parameters, parameter_files_from_args};
use crate::generated::functions::{
    parameters_l310_parameters_inputparameters,
    parameterschimeric_initialize_l6_parameterschimeric_initialize,
    readalignchunk_l5_readalignchunk_readalignchunk,
    readalignchunk_mapchunk_l7_readalignchunk_mapchunk, samheaders_l5_samheaders,
    sjdbinsertjunctions_l11_sjdbinsertjunctions, solo_l23_solo_solo, stats_l4_stats_resetn,
    stats_l62_stats_progressreportheader, timefunctions_l14_timemonthdaytime,
    twopassrunpass1_l9_twopassrunpass1,
};
use crate::generated::structs::{
    Genome, Parameters, ReadAlignChunk, ReadAlignChunkProcessChunksResult, StarMainResult, Stats,
    Transcriptome,
};

/// Intentional helper API for direct access to the translated aligner core.
///
/// Used by Bascet.
///
/// The main translation keeps STAR's original file-oriented control flow for
/// parity and auditability. This type is deliberately outside that one-function
/// per original-function scaffold and is a user-approved deviation from CCC's
/// original-function mapping: it lets callers feed borrowed read buffers into
/// the existing `ReadAlignChunk` machinery without staging FASTQ files.
pub struct DirectReadPair<'a> {
    pub name: &'a str,
    pub r1: &'a [u8],
    pub q1: &'a [u8],
    pub r2: &'a [u8],
    pub q2: &'a [u8],
}

/// Stateful direct aligner wrapper around STAR's translated core structures.
///
/// Used by Bascet.
///
/// This is an intentional helper layer for zero-copy direct access from other
/// Rust code and a user-approved CCC deviation rather than a missing original
/// STAR function. It initializes the same translated parameters, genome,
/// headers, stats, and `ReadAlignChunk` state that the CLI path uses, then
/// exposes a small in-memory chunk API so external callers can drive the
/// aligner without shelling out or writing temporary read files.
pub struct DirectStarRun {
    parameters: Parameters,
    genome: Genome,
    transcriptome: Transcriptome,
    stats_all: Stats,
    read_chunks: Vec<ReadAlignChunk>,
    process: ReadAlignChunkProcessChunksResult,
    raw_time_start: libc::time_t,
    raw_time_map: libc::time_t,
    chunk_index: u32,
    result: StarMainResult,
}

impl DirectStarRun {
    /// Used by Bascet.
    pub fn new(args: &[String]) -> Result<Self, String> {
        let mut parameters = Parameters::default();
        let parameter_files = parameter_files_from_args(args)?;
        let _scan_state = parameters_l310_parameters_inputparameters(
            &mut parameters,
            args,
            PARAMETERS_DEFAULT,
            &parameter_files,
            None,
            &[],
        )?;

        if parameters.two_pass_yes {
            return Err(
                "STAR two-pass mode is not supported by the direct in-memory STAR integration"
                    .to_string(),
            );
        }

        if !parameters.out_file_tmp.is_empty() {
            fs::create_dir_all(&parameters.out_file_tmp).map_err(|e| e.to_string())?;
        }

        let raw_time_start = current_unix_time();
        let mut result = StarMainResult {
            parameters: parameters.clone(),
            ..Default::default()
        };
        result
            .log_stdout
            .push_str(&format!("\t{}\n", parameters.command_line));
        result.log_stdout.push_str(&format!(
            "{} ..... started STAR run\n",
            timefunctions_l14_timemonthdaytime(raw_time_start)
        ));

        let mut genome = load_genome_from_parameters(&mut parameters)?;
        let mut transcriptome = Transcriptome::default();

        let solo_cell_filter = solo_l23_solo_solo(&parameters, &transcriptome, "", "", "", ".")?;
        result.log_stdout.push_str(&solo_cell_filter.log_stdout);
        result.log_main.push_str(&solo_cell_filter.log_main);
        if solo_cell_filter.exited {
            result.exit_code = 0;
            result.parameters = parameters.clone();
            result.genome = Some(genome.clone());
        }

        let mut sjdb_loci = crate::generated::structs::SjdbClass::default();
        if parameters.sjdb_insert_yes {
            let genome_main1 = genome.clone();
            let sjdb = sjdbinsertjunctions_l11_sjdbinsertjunctions(
                &mut parameters,
                &mut genome,
                &genome_main1,
                &mut sjdb_loci,
            )?;
            result.log_main.push_str(&sjdb.log_main);
            result.sjdb_insert = Some(sjdb);
        }

        result
            .log_progress
            .push_str(&stats_l62_stats_progressreportheader());
        let existing_read_files = BTreeSet::new();
        let two_pass = twopassrunpass1_l9_twopassrunpass1(
            &mut parameters,
            &mut genome,
            Some(&transcriptome),
            &mut sjdb_loci,
            None,
            &existing_read_files,
        )?;
        result.log_progress.push_str(&two_pass.log_progress);
        result.log_stdout.push_str(&two_pass.log_stdout);
        result.log_main.push_str(&two_pass.log_main);

        let mut stats_all = Stats::default();
        stats_l4_stats_resetn(&mut stats_all);
        let raw_time_map = current_unix_time();
        stats_all.time_start = raw_time_start;
        stats_all.time_start_map = raw_time_map;
        stats_all.time_last_report = raw_time_map;
        result.log_stdout.push_str(&format!(
            "{} ..... started mapping\n",
            timefunctions_l14_timemonthdaytime(raw_time_map)
        ));

        if parameters.quant_tr_sam_yes && transcriptome.tr_id.is_empty() {
            transcriptome = Transcriptome::default();
        }
        if parameters.out_sam_type.is_empty() {
            parameters.out_sam_type.push("None".to_string());
        }
        samheaders_l5_samheaders(&mut parameters, &mut genome, &transcriptome, "", "");
        let sam_header_for_chim = parameters.sam_header.clone();
        let mut p_ch = std::mem::take(&mut parameters.p_ch);
        parameterschimeric_initialize_l6_parameterschimeric_initialize(
            &mut p_ch,
            &mut parameters,
            &sam_header_for_chim,
        )?;
        parameters.p_ch = p_ch;

        let mut read_chunks = Vec::with_capacity(parameters.run_thread_n.max(0) as usize);
        for ii in 0..parameters.run_thread_n {
            read_chunks.push(readalignchunk_l5_readalignchunk_readalignchunk(
                &parameters,
                &genome,
                Some(&transcriptome),
                ii,
            )?);
        }
        if read_chunks.is_empty() {
            return Err("STAR direct run created no read chunks".to_string());
        }

        Ok(Self {
            parameters,
            genome,
            transcriptome,
            stats_all,
            read_chunks,
            process: ReadAlignChunkProcessChunksResult::default(),
            raw_time_start,
            raw_time_map,
            chunk_index: 0,
            result,
        })
    }

    /// Clear the current in-memory FASTQ chunk while preserving allocated
    /// buffers for repeated direct aligner calls.
    ///
    /// Used by Bascet.
    pub fn clear_chunk_input(&mut self) {
        let read_nends = self.parameters.read_nends as usize;
        let chunk = &mut self.read_chunks[0];
        if chunk.chunk_in.len() < read_nends {
            chunk.chunk_in.resize(read_nends, Vec::new());
        }
        if chunk.chunk_in_size_bytes_total.len() < read_nends {
            chunk.chunk_in_size_bytes_total.resize(read_nends, 0);
        }
        for imate in 0..read_nends {
            chunk.chunk_in[imate].clear();
            chunk.chunk_in_size_bytes_total[imate] = 0;
        }
    }

    /// Check whether appending this borrowed read pair would reach STAR's
    /// configured chunk byte limit.
    ///
    /// Used by Bascet.
    pub fn read_pair_would_exceed_chunk(&self, read: &DirectReadPair<'_>) -> bool {
        let chunk = &self.read_chunks[0];
        let limit = self.parameters.chunk_in_size_bytes.max(1);
        let r1_size = chunk.chunk_in[0].len() as u64 + read.r1.len() as u64 + read.q1.len() as u64;
        let r2_size = if self.parameters.read_nends > 1 {
            chunk.chunk_in[1].len() as u64 + read.r2.len() as u64 + read.q2.len() as u64
        } else {
            0
        };
        r1_size >= limit || r2_size >= limit
    }

    /// Append one borrowed read pair to the current STAR chunk.
    ///
    /// The chunk bytes are shaped exactly like STAR's internal FASTQ buffers so
    /// the translated core aligner can consume them unchanged.
    ///
    /// Used by Bascet.
    pub fn append_read_pair(&mut self, read: &DirectReadPair<'_>) {
        self.parameters.i_read_all += 1;
        let read_id = if self.parameters.out_sam_read_id_number {
            format!(
                "@{} {} N {}",
                self.parameters.i_read_all,
                self.parameters.i_read_all,
                self.parameters.read_files_index
            )
        } else {
            format!(
                "@{} {} N {}",
                read.name, self.parameters.i_read_all, self.parameters.read_files_index
            )
        };

        let chunk = &mut self.read_chunks[0];
        write_star_mate_to_chunk(&mut chunk.chunk_in[0], &read_id, read.r1, read.q1);
        if self.parameters.read_nends > 1 {
            write_star_mate_to_chunk(&mut chunk.chunk_in[1], &read_id, read.r2, read.q2);
        }
        for imate in 0..self.parameters.read_nends as usize {
            chunk.chunk_in_size_bytes_total[imate] = chunk.chunk_in[imate].len() as u64;
        }
    }

    /// Report whether the current direct chunk has reached STAR's byte limit.
    ///
    /// Used by Bascet.
    pub fn chunk_reached_limit(&self) -> bool {
        let chunk = &self.read_chunks[0];
        let limit = self.parameters.chunk_in_size_bytes.max(1);
        chunk.chunk_in_size_bytes_total[0] >= limit
            || (self.parameters.read_nends > 1 && chunk.chunk_in_size_bytes_total[1] >= limit)
    }

    /// Finalize the current direct chunk and pass it to STAR's translated
    /// `ReadAlignChunk::mapChunk` logic.
    ///
    /// Used by Bascet.
    pub fn finalize_and_map_chunk(&mut self) -> Result<(), String> {
        let read_nends = self.parameters.read_nends as usize;
        let chunk = &mut self.read_chunks[0];
        for imate in 0..read_nends {
            chunk.chunk_in[imate].push(b'\n');
            chunk.chunk_in_size_bytes_total[imate] = chunk.chunk_in[imate].len() as u64 - 1;
        }
        chunk.no_reads_left = false;
        chunk.i_chunk_in = self.chunk_index;
        self.process.chunks_read += 1;

        let map_result = readalignchunk_mapchunk_l7_readalignchunk_mapchunk(
            chunk,
            &self.parameters,
            &mut self.stats_all,
            self.raw_time_map,
            |_ra| -1,
            Some((&self.genome, &mut self.transcriptome)),
        )?;
        self.process.log_main.push_str(&map_result.log_main);
        self.process.map_chunks.push(map_result);
        self.chunk_index += 1;
        Ok(())
    }

    /// Finish the direct run and return the same aggregate result shape used by
    /// the CLI-facing STAR entry point.
    ///
    /// Used by Bascet.
    pub fn finish(mut self) -> StarMainResult {
        self.read_chunks[0].no_reads_left = true;
        if self.parameters.out_sam_bool
            && self.parameters.out_sam_order != "PairedKeepInputOrder"
            && self.read_chunks[0].chunk_out_bam_total > 0
            && let Some(last_map_chunk) = self.process.map_chunks.last_mut()
        {
            let bytes = self.read_chunks[0].chunk_out_bam_total as usize;
            let out = &self.read_chunks[0].chunk_out_bam
                [..bytes.min(self.read_chunks[0].chunk_out_bam.len())];
            last_map_chunk.direct_sam_output.extend_from_slice(out);
            self.read_chunks[0].chunk_out_bam_total = 0;
        }
        self.result.log_main.push_str(&self.process.log_main);
        self.result.process_chunks.push(self.process);
        self.result.exit_code = 0;
        self.result.parameters = self.parameters;
        self.result.genome = Some(self.genome);
        self.result.transcriptome = Some(self.transcriptome);
        self.result.stats_all = self.stats_all;
        self.result.read_chunks = self.read_chunks;
        self.result
    }

    pub fn raw_time_start(&self) -> libc::time_t {
        self.raw_time_start
    }
}

// Used by Bascet.
//
// User-approved CCC deviation for direct access: this encodes borrowed
// read/quality slices into STAR's in-memory chunk format without requiring
// caller-side FASTQ files.
fn write_star_mate_to_chunk(chunk_in: &mut Vec<u8>, read_id: &str, seq: &[u8], qual: &[u8]) {
    chunk_in.extend_from_slice(read_id.as_bytes());
    chunk_in.push(b'\n');
    chunk_in.extend_from_slice(seq);
    chunk_in.extend_from_slice(b"\n+\n");
    chunk_in.extend_from_slice(qual);
    chunk_in.push(b'\n');
}

// Used by Bascet.
//
// User-approved CCC deviation matching the timestamp type expected by the
// translated STAR stats/logging code.
fn current_unix_time() -> libc::time_t {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0)
}
