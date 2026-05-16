#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ReadAlignChunk` at STAR/source/ReadAlignChunk.h:12."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignChunk {
    pub i_thread: i32,
    pub i_chunk_in: u32,
    pub no_reads_left: bool,
    pub chunk_tr: Option<Transcriptome>,
    pub ra: ReadAlign,
    pub chunk_in: Vec<Vec<u8>>,
    pub chunk_in_size_bytes_total: Vec<u64>,
    pub read_in_stream_n: usize,
    pub chunk_out_bam: Vec<u8>,
    pub chunk_out_bam_total: u64,
    pub chunk_out_bam_file_name: String,
    pub chunk_out_bam_unsorted: Option<BAMoutput>,
    pub chunk_out_bam_coord: BAMoutput,
    pub chunk_out_bam_quant: Option<BAMoutput>,
    pub chunk_out_sj: OutSJ,
    pub chunk_out_sj1: OutSJ,
    pub chunk_out_chim_sam_path: Option<String>,
    pub chunk_out_chim_junction_path: Option<String>,
    pub chunk_out_unmapped_reads_paths: Vec<String>,
    pub chunk_out_filter_by_sjout_files: Vec<String>,
    pub wasp_ra_present: bool,
    pub pe_merge_ra_present: bool,
    pub log_main: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignChunkMapChunkResult {
    pub direct_sam_output: Vec<u8>,
    pub paired_keep_input_order_tmp: Vec<u8>,
    pub chimeric_sam_output: String,
    pub chimeric_junction_output: String,
    pub unmapped_fastx_outputs: Vec<String>,
    pub signal_records: Vec<SignalFromBamRecord>,
    pub quant_bam_output: Vec<u8>,
    pub paired_keep_input_order_tmp_name: Option<String>,
    pub paired_keep_input_order_final_name: Option<String>,
    pub progress_report: Option<String>,
    pub log_main: String,
    pub reads_processed: u64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignChunkProcessChunksResult {
    pub map_chunks: Vec<ReadAlignChunkMapChunkResult>,
    pub chunk_inputs: Vec<Vec<Vec<u8>>>,
    pub log_main: String,
    pub chunks_read: u32,
    pub paired_keep_input_order_cat_after_chunks: Vec<u32>,
    pub flushed_bam_unsorted: bool,
    pub flushed_bam_coord: bool,
    pub flushed_bam_quant: bool,
    pub chim_sam_cat_path: Option<String>,
    pub chim_junction_cat_path: Option<String>,
    pub unmapped_fastx_cat_paths: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BamSortByCoordinateResult {
    pub output_bam: Vec<u8>,
    pub bin_outputs: Vec<Vec<u8>>,
    pub bin_names: Vec<String>,
    pub removed_files: Vec<String>,
    pub max_mem: u64,
    pub unmapped_reads_n: u64,
}

#[doc = "Original `ReadAlignChunk::ReadAlignChunk` at STAR/source/ReadAlignChunk.cpp:5. Args: Pin: Parameters, genomeIn: Genome, TrIn: Transcriptome, iChunk: int"]
pub fn readalignchunk_l5_readalignchunk_readalignchunk(
    p: &crate::parameters_chimeric::Parameters,
    genome: &crate::genome::Genome,
    tr_in: Option<&crate::transcriptome::Transcriptome>,
    i_chunk: i32,
) -> Result<crate::read_align_chunk::ReadAlignChunk, String> {
    let mut chunk_tr = None;
    if p.quant_yes {
        let mut tr = tr_in.cloned().unwrap_or_default();
        transcriptome_l150_transcriptome_quantsallocate(&mut tr, p.quant_ge_count_yes);
        chunk_tr = Some(tr);
    }

    let mut ra = readalign_l6_readalign_readalign(p, genome, chunk_tr.as_ref().or(tr_in), i_chunk);
    ra.i_read = 0;

    // STAR allocates each chunkIn[ii] as a buffer of chunkInSizeBytesArray bytes
    // and memset's it to '\n' so that empty/incomplete reads stream as blank lines.
    let mut chunk_in = Vec::with_capacity(p.read_nends as usize);
    for _ in 0..p.read_nends {
        chunk_in.push(vec![b'\n'; p.chunk_in_size_bytes_array as usize]);
    }

    let mut chunk_out_bam = Vec::new();
    let mut chunk_out_bam_total = 0_u64;
    if p.out_sam_bool {
        chunk_out_bam = vec![0_u8; p.chunk_out_bam_size_bytes as usize];
        chunk_out_bam_total = 0;
    }

    let chunk_out_bam_unsorted = if p.out_bam_unsorted {
        Some(bamoutput_l36_bamoutput_bamoutput(Vec::new(), p))
    } else {
        None
    };

    let chunk_out_bam_coord = if p.out_bam_coord {
        bamoutput_l9_bamoutput_bamoutput(i_chunk, &p.out_bam_sort_tmp_dir, p)
    } else {
        crate::bam_output::BAMoutput::default()
    };

    let chunk_out_bam_quant = if p.quant_tr_sam_bam_yes {
        Some(bamoutput_l36_bamoutput_bamoutput(Vec::new(), p))
    } else {
        None
    };

    let chunk_out_sj = if p.out_sj {
        outsj_l4_outsj_outsj(p.limit_out_sj_collapsed)
    } else {
        crate::out_sj::OutSJ::default()
    };
    let chunk_out_sj1 = if p.out_filter_by_sjout_stage == 1 {
        outsj_l4_outsj_outsj(p.limit_out_sj_collapsed)
    } else {
        crate::out_sj::OutSJ::default()
    };

    let mut log_main = String::new();
    let mut chunk_out_chim_sam_path = None;
    let mut chunk_out_chim_junction_path = None;
    if p.p_ch.segment_min > 0 {
        if p.p_ch.out_sam_old {
            let prefix = format!("{}/Chimeric.out.sam.thread", p.out_file_tmp);
            let (_file, file_name, log) =
                readalignchunk_l116_readalignchunk_chunkfstreamopen(&prefix, i_chunk)?;
            log_main.push_str(&log);
            chunk_out_chim_sam_path = Some(file_name);
        }
        if p.p_ch.out_junctions {
            let prefix = format!("{}/Chimeric.out.junction.thread", p.out_file_tmp);
            let (_file, file_name, log) =
                readalignchunk_l116_readalignchunk_chunkfstreamopen(&prefix, i_chunk)?;
            log_main.push_str(&log);
            chunk_out_chim_junction_path = Some(file_name);
        }
    }

    let mut chunk_out_unmapped_reads_paths = Vec::new();
    if p.out_reads_unmapped == "Fastx" {
        for imate in 0..p.read_nends {
            let prefix = format!("{}/Unmapped.out.mate{}.thread", p.out_file_tmp, imate);
            let (_file, file_name, log) =
                readalignchunk_l116_readalignchunk_chunkfstreamopen(&prefix, i_chunk)?;
            log_main.push_str(&log);
            chunk_out_unmapped_reads_paths.push(file_name);
        }
    }

    let mut chunk_out_filter_by_sjout_files = Vec::new();
    if p.out_filter_type == "BySJout" {
        let prefix = format!("{}/FilterBySJoutFiles.mate1.thread", p.out_file_tmp);
        let (_file, file_name, log) =
            readalignchunk_l116_readalignchunk_chunkfstreamopen(&prefix, i_chunk)?;
        log_main.push_str(&log);
        chunk_out_filter_by_sjout_files.push(file_name);
        if p.read_nends == 2 {
            let prefix = format!("{}/FilterBySJoutFiles.mate2.thread", p.out_file_tmp);
            let (_file, file_name, log) =
                readalignchunk_l116_readalignchunk_chunkfstreamopen(&prefix, i_chunk)?;
            log_main.push_str(&log);
            chunk_out_filter_by_sjout_files.push(file_name);
        }
    }

    Ok(crate::read_align_chunk::ReadAlignChunk {
        i_thread: i_chunk,
        i_chunk_in: 0,
        no_reads_left: false,
        chunk_tr,
        ra,
        chunk_in,
        chunk_in_size_bytes_total: vec![0; p.read_nends as usize],
        read_in_stream_n: p.read_nends as usize,
        chunk_out_bam,
        chunk_out_bam_total,
        chunk_out_bam_file_name: String::new(),
        chunk_out_bam_unsorted,
        chunk_out_bam_coord,
        chunk_out_bam_quant,
        chunk_out_sj,
        chunk_out_sj1,
        chunk_out_chim_sam_path,
        chunk_out_chim_junction_path,
        chunk_out_unmapped_reads_paths,
        chunk_out_filter_by_sjout_files,
        wasp_ra_present: p.wasp_yes,
        pe_merge_ra_present: p.pe_overlap_nbases_min > 0,
        log_main,
    })
}

#[doc = "Original `ReadAlignChunk::chunkFstreamOpen` at STAR/source/ReadAlignChunk.cpp:116. Args: filePrefix: string, iChunk: int, fstreamOut: fstream"]
pub fn readalignchunk_l116_readalignchunk_chunkfstreamopen(
    file_prefix: &str,
    i_chunk: i32,
) -> Result<(std::fs::File, String, String), String> {
    let f_name1 = format!("{}{}", file_prefix, i_chunk);
    let mut log_main = format!("Opening the file: {} ... ", f_name1);

    match std::fs::remove_file(&f_name1) {
        Ok(()) => {}
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {}
        Err(_) => {}
    }

    match std::fs::OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(&f_name1)
    {
        Ok(file) => drop(file),
        Err(_) => {
            log_main.push_str("failed!\n");
            return Err(format!(
                "EXITING because of FATAL ERROR: could not create output file {}\nSolution: check that you have permission to write this file\n",
                f_name1
            ));
        }
    }

    match std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .open(&f_name1)
    {
        Ok(file) => {
            log_main.push_str("ok\n");
            Ok((file, f_name1, log_main))
        }
        Err(_) => {
            log_main.push_str("failed!\n");
            Err(format!(
                "EXITING because of FATAL ERROR: could not create output file {}\nSolution: check that you have permission to write this file\n",
                f_name1
            ))
        }
    }
}

#[doc = "Original `ReadAlignChunk::chunkFstreamCat` at STAR/source/ReadAlignChunk.cpp:137. Args: chunkOut: fstream, allOut: ofstream, mutexFlag: bool, mutexVal: pthread_mutex_t"]
pub fn readalignchunk_l137_readalignchunk_chunkfstreamcat<C, A>(
    chunk_out: &mut C,
    all_out: &mut A,
    _mutex_flag: bool,
) -> std::io::Result<()>
where
    C: std::io::Read + std::io::Write + std::io::Seek,
    A: std::io::Write,
{
    chunk_out.flush()?;
    chunk_out.seek(std::io::SeekFrom::Start(0))?;
    let mut buf = Vec::new();
    chunk_out.read_to_end(&mut buf)?;
    all_out.write_all(&buf)?;
    all_out.flush()?;
    chunk_out.seek(std::io::SeekFrom::Start(0))?;
    Ok(())
}

#[doc = "Original `ReadAlignChunk::chunkFilesCat` at STAR/source/ReadAlignChunk.cpp:151. Args: allOut: ostream, filePrefix: string, iC: uint"]
pub fn readalignchunk_l151_readalignchunk_chunkfilescat<W: std::io::Write>(
    all_out: &mut W,
    file_prefix: &str,
    i_c: &mut u32,
) -> std::io::Result<()> {
    loop {
        let name1 = format!("{}{}", file_prefix, *i_c);
        let mut file_chunk_in = match std::fs::File::open(&name1) {
            Ok(file) => file,
            Err(_) => break,
        };
        std::io::copy(&mut file_chunk_in, all_out)?;
        all_out.flush()?;
        drop(file_chunk_in);
        std::fs::remove_file(&name1)?;
        *i_c += 1;
    }
    Ok(())
}
