use std::collections::BTreeSet;
use std::fs;
use std::sync::Arc;

use crate::cli::{PARAMETERS_DEFAULT, load_genome_from_parameters, parameter_files_from_args};
use crate::{
    Genome, Parameters, ReadAlignChunk, ReadAlignChunkProcessChunksResult, StarMainResult, Stats,
    Transcriptome,
};
use crate::{
    clipmate_clip_l5_clipmate_clip, parameters_l310_parameters_inputparameters,
    parameterschimeric_initialize_l6_parameterschimeric_initialize,
    readalign_oneread_l8_readalign_oneread_loaded, readalignchunk_l5_readalignchunk_readalignchunk,
    readalignchunk_mapchunk_l7_readalignchunk_mapchunk,
    readalignchunk_mapchunk_l7_readalignchunk_mapchunk_with_next_read, samheaders_l5_samheaders,
    sequencefuns_l131_convertnucleotidestonumbers, sjdbinsertjunctions_l11_sjdbinsertjunctions,
    solo_l23_solo_solo, stats_l4_stats_resetn, stats_l21_stats_addstats,
    stats_l62_stats_progressreportheader, timefunctions_l14_timemonthdaytime,
    twopassrunpass1_l9_twopassrunpass1,
};

// STAR-DEVIATION:
// The types below are not translations of original STAR C++ classes. They are
// an embedded/library API layer around the translated core so callers can own
// read parsing, chunk scheduling, and worker reuse without duplicating the
// genome index or using STAR's file-oriented thread scheduler.
#[derive(Clone, Copy, Debug)]
pub struct StarReadMate<'a> {
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>,
}

#[derive(Clone, Copy, Debug)]
pub struct StarReadPair<'a> {
    pub name: &'a str,
    pub mate1: StarReadMate<'a>,
    pub mate2: Option<StarReadMate<'a>>,
    pub read_number: u64,
    pub read_files_index: u32,
    pub filter: u8,
    pub extra: &'a str,
}

// STAR-DEVIATION:
// Trait boundary for already-parsed input chunks. The generic associated
// iterator keeps the direct path monomorphized for library callers while the
// normal STAR CLI path continues to use the original translated chunk reader.
pub trait StarReadChunk {
    type Iter<'a>: Iterator<Item = StarReadPair<'a>>
    where
        Self: 'a;

    fn chunk_index(&self) -> u32;
    fn reads(&self) -> Self::Iter<'_>;

    fn estimated_input_bytes(&self) -> usize {
        0
    }
}

// STAR-DEVIATION:
// Convenience owned implementation of `StarReadChunk` for embedders. Serious
// callers can provide their own offset-backed chunk storage to avoid per-read
// Vec/String ownership while using the same trait.
#[derive(Clone, Debug, Default)]
pub struct OwnedStarReadPair {
    pub name: String,
    pub r1: Vec<u8>,
    pub q1: Vec<u8>,
    pub r2: Vec<u8>,
    pub q2: Vec<u8>,
    pub read_number: u64,
    pub read_files_index: u32,
    pub filter: u8,
    pub extra: String,
}

#[derive(Clone, Debug, Default)]
pub struct OwnedStarReadChunk {
    pub chunk_index: u32,
    pub reads: Vec<OwnedStarReadPair>,
}

pub struct OwnedStarReadChunkIter<'a> {
    inner: std::slice::Iter<'a, OwnedStarReadPair>,
}

impl<'a> Iterator for OwnedStarReadChunkIter<'a> {
    type Item = StarReadPair<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let read = self.inner.next()?;
        Some(StarReadPair {
            name: &read.name,
            mate1: StarReadMate {
                seq: &read.r1,
                qual: Some(&read.q1),
            },
            mate2: if read.r2.is_empty() {
                None
            } else {
                Some(StarReadMate {
                    seq: &read.r2,
                    qual: Some(&read.q2),
                })
            },
            read_number: read.read_number,
            read_files_index: read.read_files_index,
            filter: read.filter,
            extra: &read.extra,
        })
    }
}

impl StarReadChunk for OwnedStarReadChunk {
    type Iter<'a> = OwnedStarReadChunkIter<'a>;

    fn chunk_index(&self) -> u32 {
        self.chunk_index
    }

    fn reads(&self) -> Self::Iter<'_> {
        OwnedStarReadChunkIter {
            inner: self.reads.iter(),
        }
    }

    fn estimated_input_bytes(&self) -> usize {
        self.reads
            .iter()
            .map(|read| {
                read.name.len()
                    + read.r1.len()
                    + read.q1.len()
                    + read.r2.len()
                    + read.q2.len()
                    + read.extra.len()
                    + 32
            })
            .sum()
    }
}

// STAR-DEVIATION:
// Shared immutable STAR state for embedded callers. Original STAR stores this
// behind process-wide/file-oriented control flow; this wrapper makes the split
// explicit so an external scheduler can share reference state across workers.
pub struct DirectStarContext {
    parameters: Arc<Parameters>,
    genome: Arc<Genome>,
    transcriptome_template: Arc<Transcriptome>,
    raw_time_start: libc::time_t,
    raw_time_map: libc::time_t,
    result_template: StarMainResult,
}

// STAR-DEVIATION:
// Per-worker mutable STAR state. This deliberately mirrors STAR's real unit of
// parallelism (`ReadAlignChunk` plus per-worker transcriptome/stat scratch)
// instead of creating one full STAR run per caller thread.
pub struct DirectStarWorker {
    worker_id: i32,
    chunk: ReadAlignChunk,
    transcriptome: Transcriptome,
}

// STAR-DEVIATION:
// Result envelope that returns both the mapped output and the reusable worker
// to the embedding scheduler. Original STAR writes/merges chunks internally;
// embedders can instead order and consume chunks themselves.
pub struct DirectStarMappedChunk<C> {
    pub chunk_index: u32,
    pub worker_id: i32,
    pub input: C,
    pub worker: DirectStarWorker,
    pub map_result: crate::read_align_chunk::ReadAlignChunkMapChunkResult,
    pub stats: Stats,
}

impl DirectStarContext {
    pub fn new(args: &[String]) -> Result<Self, String> {
        let run = DirectStarRun::new(args)?;
        Ok(Self {
            parameters: Arc::new(run.parameters),
            genome: Arc::new(run.genome),
            transcriptome_template: Arc::new(run.transcriptome),
            raw_time_start: run.raw_time_start,
            raw_time_map: run.raw_time_map,
            result_template: run.result,
        })
    }

    pub fn run_thread_n(&self) -> usize {
        self.parameters.run_thread_n.max(0) as usize
    }

    pub fn read_nends(&self) -> usize {
        self.parameters.read_nends as usize
    }

    pub fn chunk_input_limit_bytes(&self) -> u64 {
        self.parameters.chunk_in_size_bytes
    }

    pub fn sam_header(&self) -> &str {
        &self.parameters.sam_header
    }

    pub fn make_worker(&self, worker_id: i32) -> Result<DirectStarWorker, String> {
        let chunk = readalignchunk_l5_readalignchunk_readalignchunk(
            self.parameters.as_ref(),
            self.genome.as_ref(),
            Some(self.transcriptome_template.as_ref()),
            worker_id,
        )?;
        Ok(DirectStarWorker {
            worker_id,
            chunk,
            transcriptome: (*self.transcriptome_template).clone(),
        })
    }

    pub fn map_read_chunk<C: StarReadChunk>(
        &self,
        mut worker: DirectStarWorker,
        input: C,
    ) -> Result<DirectStarMappedChunk<C>, String> {
        worker.chunk.no_reads_left = true;
        worker.chunk.i_chunk_in = input.chunk_index();

        let mut stats = Stats::default();
        stats_l4_stats_resetn(&mut stats);
        stats.time_start = self.raw_time_start;
        stats.time_start_map = self.raw_time_map;
        stats.time_last_report = self.raw_time_map;

        let map_result = if can_stream_star_read_chunk(&worker.chunk) {
            map_star_read_chunk_direct(
                &mut worker.chunk,
                self.parameters.as_ref(),
                &input,
                &mut stats,
                self.raw_time_map,
                self.genome.as_ref(),
                &mut worker.transcriptome,
            )?
        } else {
            // STAR-DEVIATION:
            // Fallback for rare configurations that need STAR's chunk-level input
            // mutation before readLoad, e.g. clipChunk. The common direct path
            // below streams generated FASTQ records without materializing chunk_in.
            pack_star_read_chunk(&mut worker.chunk, &self.parameters, &input)?;
            readalignchunk_mapchunk_l7_readalignchunk_mapchunk(
                &mut worker.chunk,
                self.parameters.as_ref(),
                &mut stats,
                self.raw_time_map,
                |_ra| -1,
                Some((self.genome.as_ref(), &mut worker.transcriptome)),
            )?
        };

        Ok(DirectStarMappedChunk {
            chunk_index: input.chunk_index(),
            worker_id: worker.worker_id,
            input,
            worker,
            map_result,
            stats,
        })
    }

    pub fn finish_mapped_chunks<C>(
        &self,
        mapped_chunks_in_order: impl IntoIterator<Item = DirectStarMappedChunk<C>>,
    ) -> StarMainResult {
        // STAR-DEVIATION:
        // Caller supplies mapped chunks in deterministic output order. This
        // replaces STAR's internal thread/chunk concatenation for embedded use.
        let mut result = self.result_template.clone();
        let mut process = ReadAlignChunkProcessChunksResult::default();
        let mut stats_all = Stats::default();
        stats_l4_stats_resetn(&mut stats_all);
        stats_all.time_start = self.raw_time_start;
        stats_all.time_start_map = self.raw_time_map;
        stats_all.time_last_report = self.raw_time_map;

        let mut read_chunks = Vec::new();
        for mapped in mapped_chunks_in_order {
            process.log_main.push_str(&mapped.map_result.log_main);
            process.map_chunks.push(mapped.map_result);
            process.chunks_read += 1;
            stats_l21_stats_addstats(&mut stats_all, &mapped.stats);
            read_chunks.push(mapped.worker.chunk);
        }

        result.process_chunks.push(process);
        result.exit_code = 0;
        result.parameters = self.parameters.as_ref().clone();
        result.genome = Some(self.genome.as_ref().clone());
        result.transcriptome = Some(self.transcriptome_template.as_ref().clone());
        result.stats_all = stats_all;
        result.read_chunks = read_chunks;
        result
    }
}

impl DirectStarWorker {
    pub fn worker_id(&self) -> i32 {
        self.worker_id
    }

    pub fn chunk(&self) -> &ReadAlignChunk {
        &self.chunk
    }

    pub fn chunk_mut(&mut self) -> &mut ReadAlignChunk {
        &mut self.chunk
    }
}

fn can_stream_star_read_chunk(chunk: &ReadAlignChunk) -> bool {
    chunk.ra.clip_mates.is_empty() || chunk.ra.clip_mates[0].is_empty() || chunk.chunk_in.is_empty()
}

fn map_star_read_chunk_direct<C: StarReadChunk>(
    chunk: &mut ReadAlignChunk,
    parameters: &Parameters,
    input: &C,
    stats: &mut Stats,
    raw_time_map: libc::time_t,
    genome: &Genome,
    transcriptome: &mut Transcriptome,
) -> Result<crate::read_align_chunk::ReadAlignChunkMapChunkResult, String> {
    let start_read_number = chunk.ra.i_read_all + 1;
    let mut reads = input.reads();
    let mut next_read_number = start_read_number;
    readalignchunk_mapchunk_l7_readalignchunk_mapchunk_with_next_read(
        chunk,
        parameters,
        stats,
        raw_time_map,
        |_ra| -1,
        Some((genome, transcriptome)),
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
            let Some(loaded) =
                load_next_star_read_pair(ra, p_one_read, &mut reads, &mut next_read_number)?
            else {
                let mut result = crate::quantifications::ReadAlignOneReadResult::default();
                result.status = -1;
                return Ok(result);
            };
            readalign_oneread_l8_readalign_oneread_loaded(
                ra,
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
                loaded.read0,
                loaded.qual0,
                loaded.read_name_mates,
                loaded.read_status,
            )
        },
    )
}

struct LoadedStarReadPair {
    read0: Vec<String>,
    qual0: Vec<String>,
    read_name_mates: Vec<String>,
    read_status: i32,
}

fn load_next_star_read_pair<'a, I>(
    read_align: &mut crate::read_align::ReadAlign,
    parameters: &crate::parameters_chimeric::Parameters,
    reads: &mut I,
    next_read_number: &mut u64,
) -> Result<Option<LoadedStarReadPair>, String>
where
    I: Iterator<Item = StarReadPair<'a>>,
{
    let Some(read) = reads.next() else {
        return Ok(None);
    };
    let read_nends = parameters.read_nends as usize;
    prepare_direct_read_align_scratch(read_align, parameters);

    validate_star_read_mate(read.name, 0, read.mate1).map_err(|err| err.to_string())?;
    if read_nends > 1 {
        let mate2 = read
            .mate2
            .ok_or_else(|| format!("STAR direct read {} is missing mate 2", read.name))?;
        validate_star_read_mate(read.name, 1, mate2).map_err(|err| err.to_string())?;
    }

    let read_number = if read.read_number == 0 {
        let value = *next_read_number;
        *next_read_number += 1;
        value
    } else {
        read.read_number
    };
    let read_files_index = if read.read_files_index == 0 {
        parameters.read_files_index
    } else {
        read.read_files_index
    };
    let read_filter = if read.filter == 0 { b'N' } else { read.filter };
    read_align.i_read_all = read_number;
    read_align.read_filter = read_filter as i32;
    read_align.read_files_index = read_files_index;
    if read_align.read_name_extra.len() < read_nends {
        read_align.read_name_extra.resize(read_nends, String::new());
    }

    let mut read0 = std::mem::take(&mut read_align.read0_text);
    if read0.len() < read_nends {
        read0.resize(read_nends, String::new());
    }
    let mut qual0 = std::mem::take(&mut read_align.qual0_text);
    if qual0.len() < read_nends {
        qual0.resize(read_nends, String::new());
    }
    let mut read_name_mates = std::mem::take(&mut read_align.read_name_mates_text);
    if read_name_mates.len() < read_nends {
        read_name_mates.resize(read_nends, String::new());
    }

    let mut star_name = if parameters.out_sam_read_id_number {
        read_number.to_string()
    } else {
        read.name.strip_prefix('@').unwrap_or(read.name).to_string()
    };
    trim_star_read_name(&mut star_name, &parameters.read_name_separator_char);
    for mate_index in 0..read_nends {
        let mate = if mate_index == 0 {
            read.mate1
        } else {
            read.mate2.expect("mate2 presence was validated")
        };
        load_star_mate_into_scratch(
            read_align,
            parameters,
            mate_index,
            mate,
            read.extra,
            &star_name,
            &mut read0[mate_index],
            &mut qual0[mate_index],
            &mut read_name_mates[mate_index],
        )?;
    }

    Ok(Some(LoadedStarReadPair {
        read0,
        qual0,
        read_name_mates,
        read_status: 2,
    }))
}

fn trim_star_read_name(read_name: &mut String, separators: &[char]) {
    for separator in separators {
        if let Some(pos) = read_name.find(*separator) {
            read_name.truncate(pos);
        }
    }
}

fn prepare_direct_read_align_scratch(
    read_align: &mut crate::read_align::ReadAlign,
    parameters: &crate::parameters_chimeric::Parameters,
) {
    let read_nends = parameters.read_nends as usize;
    if read_align.read_length.len() < read_nends.max(2) {
        read_align.read_length.resize(read_nends.max(2), 0);
    }
    if read_align.read_length_original.len() < read_nends.max(2) {
        read_align.read_length_original.resize(read_nends.max(2), 0);
    }
    if read_align.read0.len() < read_nends {
        read_align.read0.resize(read_nends, Vec::new());
    }
    if read_align.qual0.len() < read_nends {
        read_align.qual0.resize(read_nends, Vec::new());
    }
    if read_align.read_name_mates.len() < read_nends {
        read_align
            .read_name_mates
            .resize(read_nends, vec![0; crate::DEF_READ_NAME_LENGTH_MAX]);
    }
    if read_align.clip_mates.len() < read_nends {
        read_align
            .clip_mates
            .resize(read_nends, vec![crate::clip_mate::ClipMate::default(); 2]);
    }
    if read_align.qual_hist.len() < parameters.read_nmates as usize {
        read_align
            .qual_hist
            .resize(parameters.read_nmates as usize, vec![0; 256]);
    }
    for mate_index in 0..read_nends {
        if read_align.clip_mates[mate_index].len() < 2 {
            read_align.clip_mates[mate_index].resize(2, crate::clip_mate::ClipMate::default());
        }
    }
}

fn load_star_mate_into_scratch(
    read_align: &mut crate::read_align::ReadAlign,
    parameters: &crate::parameters_chimeric::Parameters,
    mate_index: usize,
    mate: StarReadMate<'_>,
    extra: &str,
    star_name: &str,
    read0: &mut String,
    qual0: &mut String,
    read_name_mates: &mut String,
) -> Result<(), String> {
    read_name_mates.clear();
    read_name_mates.push('@');
    read_name_mates.push_str(star_name);
    read0.clear();
    read0.push_str(std::str::from_utf8(mate.seq).map_err(|err| err.to_string())?);
    qual0.clear();
    if let Some(qual) = mate.qual {
        if parameters.out_qs_conversion_add == 0 {
            qual0.push_str(std::str::from_utf8(qual).map_err(|err| err.to_string())?);
        } else {
            for byte in qual {
                let qs = (*byte as i32 + parameters.out_qs_conversion_add).clamp(33, 126) as u8;
                qual0.push(qs as char);
            }
        }
    } else {
        qual0.extend(std::iter::repeat('A').take(mate.seq.len()));
    }
    read_align.read_name_extra[mate_index].clear();
    read_align.read_name_extra[mate_index].push_str(extra);
    read_align.read_length[mate_index] = mate.seq.len() as u64;
    read_align.read_length_original[mate_index] = mate.seq.len() as u64;
    read_align.read1[mate_index].clear();
    read_align.read1[mate_index].resize(mate.seq.len(), 0);
    sequencefuns_l131_convertnucleotidestonumbers(
        mate.seq,
        &mut read_align.read1[mate_index],
        mate.seq.len() as u64,
    );
    clipmate_clip_l5_clipmate_clip(
        &mut read_align.clip_mates[mate_index][0],
        &mut read_align.read_length[mate_index],
        &mut read_align.read1[mate_index],
    );
    clipmate_clip_l5_clipmate_clip(
        &mut read_align.clip_mates[mate_index][1],
        &mut read_align.read_length[mate_index],
        &mut read_align.read1[mate_index],
    );
    Ok(())
}

fn validate_star_read_mate(
    read_name: &str,
    mate_index: usize,
    mate: StarReadMate<'_>,
) -> std::io::Result<()> {
    if mate.seq.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "STAR direct read {} has an empty mate {} sequence",
                read_name,
                mate_index + 1
            ),
        ));
    }
    if let Some(qual) = mate.qual
        && qual.len() != mate.seq.len()
    {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "STAR direct read {} mate {} quality length differs from sequence length",
                read_name,
                mate_index + 1
            ),
        ));
    }
    Ok(())
}

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

        let mut sjdb_loci = crate::SjdbClass::default();
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

// STAR-DEVIATION:
// Adapter from embedded parsed-read chunks into STAR's current in-memory FASTQ
// chunk representation. This is intentionally isolated so replacing it with a
// lower-copy direct `ReadAlign` loader will not affect the public trait API or
// the normal STAR CLI path.
fn pack_star_read_chunk<C: StarReadChunk>(
    chunk: &mut ReadAlignChunk,
    parameters: &Parameters,
    input: &C,
) -> Result<(), String> {
    let read_nends = parameters.read_nends as usize;
    if chunk.chunk_in.len() < read_nends {
        chunk.chunk_in.resize(read_nends, Vec::new());
    }
    if chunk.chunk_in_size_bytes_total.len() < read_nends {
        chunk.chunk_in_size_bytes_total.resize(read_nends, 0);
    }
    for imate in 0..read_nends {
        chunk.chunk_in[imate].clear();
        chunk.chunk_in[imate].reserve(input.estimated_input_bytes() / read_nends.max(1));
        chunk.chunk_in_size_bytes_total[imate] = 0;
    }

    let mut next_read_number = chunk.ra.i_read_all + 1;
    for read in input.reads() {
        if read.mate1.seq.is_empty() {
            return Err(format!(
                "STAR direct read {} has an empty mate 1 sequence",
                read.name
            ));
        }
        if let Some(qual) = read.mate1.qual
            && qual.len() != read.mate1.seq.len()
        {
            return Err(format!(
                "STAR direct read {} mate 1 quality length differs from sequence length",
                read.name
            ));
        }
        if read_nends > 1 {
            let Some(mate2) = read.mate2 else {
                return Err(format!("STAR direct read {} is missing mate 2", read.name));
            };
            if mate2.seq.is_empty() {
                return Err(format!(
                    "STAR direct read {} has an empty mate 2 sequence",
                    read.name
                ));
            }
            if let Some(qual) = mate2.qual
                && qual.len() != mate2.seq.len()
            {
                return Err(format!(
                    "STAR direct read {} mate 2 quality length differs from sequence length",
                    read.name
                ));
            }
        }

        let read_name = read.name.strip_prefix('@').unwrap_or(read.name);
        let read_filter = if read.filter == 0 { b'N' } else { read.filter };
        let read_number = if read.read_number == 0 {
            let value = next_read_number;
            next_read_number += 1;
            value
        } else {
            read.read_number
        };
        let read_files_index = if read.read_files_index == 0 {
            parameters.read_files_index
        } else {
            read.read_files_index
        };
        let mut read_id = if parameters.out_sam_read_id_number {
            format!(
                "@{} {} {} {}",
                read_number, read_number, read_filter as char, read_files_index
            )
        } else {
            format!(
                "@{} {} {} {}",
                read_name, read_number, read_filter as char, read_files_index
            )
        };
        if !read.extra.is_empty() {
            read_id.push(' ');
            read_id.push_str(read.extra);
        }

        write_star_mate_to_chunk_optional(
            &mut chunk.chunk_in[0],
            &read_id,
            read.mate1.seq,
            read.mate1.qual,
        );
        if read_nends > 1 {
            let mate2 = read.mate2.expect("mate2 presence was validated");
            write_star_mate_to_chunk_optional(
                &mut chunk.chunk_in[1],
                &read_id,
                mate2.seq,
                mate2.qual,
            );
        }
    }

    for imate in 0..read_nends {
        chunk.chunk_in[imate].push(b'\n');
        chunk.chunk_in_size_bytes_total[imate] = chunk.chunk_in[imate].len() as u64 - 1;
    }
    Ok(())
}

// STAR-DEVIATION:
// Direct API helper that synthesizes a FASTQ record from parsed mate slices.
// Missing qualities are represented as STAR-compatible dummy 'A' qualities,
// matching the existing FASTA fallback behavior in the translated loader.
fn write_star_mate_to_chunk_optional(
    chunk_in: &mut Vec<u8>,
    read_id: &str,
    seq: &[u8],
    qual: Option<&[u8]>,
) {
    chunk_in.extend_from_slice(read_id.as_bytes());
    chunk_in.push(b'\n');
    chunk_in.extend_from_slice(seq);
    chunk_in.extend_from_slice(b"\n+\n");
    if let Some(qual) = qual {
        chunk_in.extend_from_slice(qual);
    } else {
        chunk_in.extend(std::iter::repeat(b'A').take(seq.len()));
    }
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
