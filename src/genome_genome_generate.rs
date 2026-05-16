#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `funCompareSuffixes` at STAR/source/Genome_genomeGenerate.cpp:29. Args: a: void, b: void"]
pub fn genome_genomegenerate_l29_funcomparesuffixes(
    a: u64,
    b: u64,
    global_g: &[u8],
    global_l: u64,
) -> i32 {
    let mut jj = 0u64;
    while jj < global_l {
        let start_a = (a - 7 - jj * 8) as usize;
        let start_b = (b - 7 - jj * 8) as usize;
        let va = u64::from_ne_bytes(global_g[start_a..start_a + 8].try_into().unwrap());
        let vb = u64::from_ne_bytes(global_g[start_b..start_b + 8].try_into().unwrap());
        let has5_va = ((va ^ 0x0505_0505_0505_0505u64).wrapping_sub(0x0101_0101_0101_0101u64)
            & !(va ^ 0x0505_0505_0505_0505u64)
            & 0x8080_8080_8080_8080u64)
            != 0;
        let has5_vb = ((vb ^ 0x0505_0505_0505_0505u64).wrapping_sub(0x0101_0101_0101_0101u64)
            & !(vb ^ 0x0505_0505_0505_0505u64)
            & 0x8080_8080_8080_8080u64)
            != 0;

        if has5_va && has5_vb {
            let va1 = va.to_ne_bytes();
            let vb1 = vb.to_ne_bytes();
            for ii in (0..8).rev() {
                if va1[ii] > vb1[ii] {
                    return 1;
                } else if va1[ii] < vb1[ii] {
                    return -1;
                } else if va1[ii] == 5 {
                    if a > b {
                        return -1;
                    } else {
                        return 1;
                    }
                }
            }
        } else if va > vb {
            return 1;
        } else if va < vb {
            return -1;
        }
        jj += 1;
    }

    if a > b { -1 } else { 1 }
}

#[doc = "Original `funG2strLocus` at STAR/source/Genome_genomeGenerate.cpp:91. Args: SAstr: uint, N: uint, GstrandBit: char, GstrandMask: uint"]
pub fn genome_genomegenerate_l91_fung2strlocus(
    mut sastr: u32,
    n: u32,
    gstrand_bit: u8,
    gstrand_mask: u32,
) -> u32 {
    let strand_g = (sastr >> gstrand_bit) == 0;
    sastr &= gstrand_mask;
    if !strand_g {
        sastr += n;
    }
    sastr
}

pub fn compare_suffix_positions_star_wordwise(
    g: &[u8],
    a: usize,
    b: usize,
    global_l: usize,
) -> std::cmp::Ordering {
    if a == b {
        return std::cmp::Ordering::Equal;
    }

    let ptr = g.as_ptr();
    for jj in 0..global_l {
        let offset = jj * std::mem::size_of::<u64>();
        let Some(start_a) = a.checked_sub(7 + offset) else {
            return compare_suffix_positions_star_bytewise(g, a, b, offset);
        };
        let Some(start_b) = b.checked_sub(7 + offset) else {
            return compare_suffix_positions_star_bytewise(g, a, b, offset);
        };
        if start_a + 8 > g.len() || start_b + 8 > g.len() {
            return compare_suffix_positions_star_bytewise(g, a, b, offset);
        }

        // STAR casts `(globalG - 7) + suffix_index - 8*jj` to an
        // unsigned-long-long pointer. Suffix positions are byte offsets, so
        // unaligned reads are expected.
        let va = unsafe { std::ptr::read_unaligned(ptr.add(start_a) as *const u64) };
        let vb = unsafe { std::ptr::read_unaligned(ptr.add(start_b) as *const u64) };
        if has_genome_spacing_byte(va) && has_genome_spacing_byte(vb) {
            let va_bytes = va.to_ne_bytes();
            let vb_bytes = vb.to_ne_bytes();
            for ii in (0..8).rev() {
                if va_bytes[ii] != vb_bytes[ii] {
                    return va_bytes[ii].cmp(&vb_bytes[ii]);
                }
                if va_bytes[ii] == GENOME_SPACING_CHAR {
                    return b.cmp(&a);
                }
            }
        } else if va != vb {
            return va.cmp(&vb);
        }
    }

    b.cmp(&a)
}

pub fn compare_suffix_positions_star_bytewise(
    g: &[u8],
    a: usize,
    b: usize,
    start_offset: usize,
) -> std::cmp::Ordering {
    let ptr = g.as_ptr();
    let mut offset = start_offset;
    loop {
        let ba = if let Some(pos) = a.checked_sub(offset) {
            if pos < g.len() {
                unsafe { *ptr.add(pos) }
            } else {
                GENOME_SPACING_CHAR
            }
        } else {
            GENOME_SPACING_CHAR
        };
        let bb = if let Some(pos) = b.checked_sub(offset) {
            if pos < g.len() {
                unsafe { *ptr.add(pos) }
            } else {
                GENOME_SPACING_CHAR
            }
        } else {
            GENOME_SPACING_CHAR
        };
        if ba != bb {
            return ba.cmp(&bb);
        }
        if ba == GENOME_SPACING_CHAR {
            return b.cmp(&a);
        }
        offset += 1;
    }
}

pub fn has_genome_spacing_byte(v: u64) -> bool {
    let x = v ^ 0x0505_0505_0505_0505u64;
    ((x.wrapping_sub(0x0101_0101_0101_0101u64)) & !x & 0x8080_8080_8080_8080u64) != 0
}

#[doc = "Original `Genome::genomeGenerate` at STAR/source/Genome_genomeGenerate.cpp:98. Args: "]
pub fn genome_genomegenerate_l98_genome_genomegenerate(
    genome: &mut crate::genome::Genome,
    p: &mut crate::parameters_chimeric::Parameters,
    gtf_contents: Option<&str>,
) -> Result<crate::parameters_chimeric::GenomeGenerateResult, String> {
    use std::io::Write;

    let mut result = crate::parameters_chimeric::GenomeGenerateResult::default();
    std::fs::create_dir_all(&genome.p_ge.g_dir).map_err(|e| e.to_string())?;

    if genome.sjdb_overhang == 0
        && (genome
            .p_ge
            .sjdb_file_chr_start_end
            .first()
            .map(|s| s.as_str())
            != Some("-")
            || genome.p_ge.sjdb_gtf_file != "-")
    {
        return Err("EXITING because of FATAL INPUT PARAMETER ERROR: for generating genome with annotations (--sjdbFileChrStartEnd or --sjdbGTFfile options)\nyou need to specify >0 --sjdbOverhang\nSOLUTION: re-run genome generation specifying non-zero --sjdbOverhang, which ideally should be equal to OneMateLength-1, or could be chosen generically as ~100\n".to_string());
    }
    if genome
        .p_ge
        .sjdb_file_chr_start_end
        .first()
        .map(|s| s.as_str())
        == Some("-")
        && genome.p_ge.sjdb_gtf_file == "-"
    {
        genome.sjdb_overhang = 0;
        genome.sjdb_length = 0;
    }

    let raw_time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    let start_msg = format!(
        "{} ... starting to generate Genome files\n",
        timefunctions_l14_timemonthdaytime(raw_time)
    );
    result.log_main.push_str(&start_msg);
    result.log_stdout.push_str(&start_msg);

    genome.genome_chr_bin_nbases = 1_u32.checked_shl(genome.p_ge.g_chr_bin_nbits).unwrap_or(0);
    if genome.genome_chr_bin_nbases == 0 {
        genome.genome_chr_bin_nbases = 1;
    }

    genome.chr_name.clear();
    genome.chr_start.clear();
    genome.chr_length.clear();
    genome.chr_name_index.clear();
    let n_genome_scan = genomescanfastafiles_l5_genomescanfastafiles(p, &mut [], false, genome)?;
    let raw_time_alloc = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    let (_n_g1_alloc, mut g_alloc, _g_offset): (u64, Vec<u8>, usize) =
        genome_l219_genome_genomesequenceallocate(n_genome_scan as u64, p, raw_time_alloc)
            .map_err(|e| e.stream_out2)?;
    genomescanfastafiles_l5_genomescanfastafiles(p, g_alloc.as_mut_slice(), true, genome)?;
    genome.n_genome = n_genome_scan as u64;
    genome.g = g_alloc;

    result.n_genome_true = genome.chr_length.iter().sum();
    result.log_main.push_str(&format!(
        "Genome sequence total length = {}\n",
        result.n_genome_true
    ));
    result
        .log_main
        .push_str(&format!("Genome size with padding = {}\n", genome.n_genome));

    let mut sjdb_loci = crate::sjdb_class::SjdbClass::default();
    let genome_dir = genome.p_ge.g_dir.clone();
    let (mut main_gtf, gtf_log) = gtf_l7_gtf_gtf(genome, p, &genome_dir, gtf_contents)?;
    result.log_main.push_str(&gtf_log);

    let transform_out =
        genome_transformgenome_l33_genome_transformgenome(genome, &mut main_gtf, "");
    result.log_main.push_str(&transform_out.log_main);

    let super_out = gtf_supertranscript_l9_gtf_supertranscript(&mut main_gtf, genome, p);
    if !super_out.log_main.is_empty()
        || !super_out.transcript_sequences_fasta.is_empty()
        || !super_out.super_transcript_sequences_fasta.is_empty()
    {
        result.log_main.push_str(&super_out.log_main);
        result.super_transcriptome = Some(super_out);
    }

    genome_l209_genome_chrbinfill(genome);
    let gtf_tg_out = gtf_transcriptgenesj_l23_gtf_transcriptgenesj(
        &mut main_gtf,
        genome,
        &mut sjdb_loci,
        &genome_dir,
        &mut result.log_main,
    );
    if main_gtf.gtf_yes {
        std::fs::create_dir_all(&genome_dir).map_err(|e| e.to_string())?;
        for (name, contents) in [
            ("geneInfo.tab", gtf_tg_out.gene_info_tab.as_str()),
            (
                "transcriptInfo.tab",
                gtf_tg_out.transcript_info_tab.as_str(),
            ),
            ("exonInfo.tab", gtf_tg_out.exon_info_tab.as_str()),
            ("exonGeTrInfo.tab", gtf_tg_out.exon_ge_tr_info_tab.as_str()),
            (
                "sjdbList.fromGTF.out.tab",
                gtf_tg_out.sjdb_list_from_gtf_out_tab.as_str(),
            ),
        ] {
            let path = format!("{}/{}", genome_dir, name);
            std::fs::write(&path, contents).map_err(|e| e.to_string())?;
            result.files_written.push(path);
        }
    }
    if main_gtf.gtf_yes || !gtf_tg_out.sjdb_list_from_gtf_out_tab.is_empty() {
        result.gtf = Some(gtf_tg_out);
    }
    result
        .log_main
        .push_str(&sjdbloadfromfiles_l6_sjdbloadfromfiles(p, &mut sjdb_loci)?);

    if result.n_genome_true > 0
        && genome.p_ge.g_saindex_nbases as f64 > (result.n_genome_true as f64).log2() / 2.0 - 1.0
    {
        result.log_main.push_str(&format!(
            "WARNING: --genomeSAindexNbases {} is too large for the genome size={}\n",
            genome.p_ge.g_saindex_nbases, result.n_genome_true
        ));
    }

    genome_genomegenerate_l417_genome_writechrinfo(genome, &genome.p_ge.g_dir)?;
    result.files_written.extend(
        [
            "chrName.txt",
            "chrStart.txt",
            "chrLength.txt",
            "chrNameLength.txt",
        ]
        .iter()
        .map(|name| format!("{}/{}", genome.p_ge.g_dir, name)),
    );

    let n_genome = genome.n_genome as usize;
    if genome.g.len() < 2 * n_genome {
        genome.g.resize(2 * n_genome, GENOME_SPACING_CHAR);
    }
    for ii in 0..n_genome {
        genome.g[2 * n_genome - 1 - ii] = if genome.g[ii] < 4 {
            3 - genome.g[ii]
        } else {
            genome.g[ii]
        };
    }

    let sparse_d = std::cmp::max(1, genome.p_ge.g_sasparse_d) as usize;
    let g_len = 2 * n_genome;
    for ii in 0..n_genome {
        genome.g.swap(ii, g_len - 1 - ii);
    }
    // STAR's funCompareSuffixes packs the 4-byte prefix as nibbles:
    // p1 = (G[ii]<<12) + (G[ii-1]<<8) + (G[ii-2]<<4) + G[ii-3].
    // Base-5 packing (an earlier translation bug) made chunk order disagree
    // with the suffix comparator, which then trips the monotonicity check in
    // genomeSAindexChunk.
    let sa_prefix = |ii: usize| -> usize {
        let mut prefix = 0usize;
        for offset in 0..4 {
            let b = ii
                .checked_sub(offset)
                .and_then(|pos| genome.g.get(pos))
                .copied()
                .unwrap_or(GENOME_SPACING_CHAR);
            prefix = (prefix << 4) | (b as usize & 0xF);
        }
        prefix
    };
    let ind_pref_n = 1usize << 16;
    let mut ind_pref_count = vec![0usize; ind_pref_n];
    let mut n_sa_usize = 0usize;
    for ii in (0..g_len).step_by(sparse_d) {
        if genome.g[ii] < 4 {
            ind_pref_count[sa_prefix(ii)] += 1;
            n_sa_usize += 1;
        }
    }
    genome.n_sa = n_sa_usize as u64;

    let total_with_sj =
        genome.n_genome + p.limit_sjdb_insert_nsj as u64 * genome.sjdb_length as u64;
    let mut gstrand_bit = if total_with_sj > 0 {
        (total_with_sj as f64).log2().floor() as u32 + 1
    } else {
        1
    };
    if gstrand_bit < 32 {
        gstrand_bit = 32;
    }
    genome.gstrand_bit = gstrand_bit;
    genome.gstrand_mask = if genome.gstrand_bit >= 64 {
        u64::MAX
    } else {
        !(1_u64 << genome.gstrand_bit)
    };
    result.log_main.push_str(&format!(
        "Estimated genome size with padding and SJs: total=genome+SJ={} = {} + {}\n",
        total_with_sj,
        genome.n_genome,
        p.limit_sjdb_insert_nsj as u64 * genome.sjdb_length as u64
    ));
    result
        .log_main
        .push_str(&format!("GstrandBit={}\n", genome.gstrand_bit));

    packedarray_l8_packedarray_definebits(
        &mut genome.sa_packed,
        genome.gstrand_bit as u64 + 1,
        genome.n_sa as u64,
    );
    let sa_pass1_len = if p.sjdb_insert_yes {
        genome.n_sa
            + 2 * genome.sjdb_length as u64
                * std::cmp::min(sjdb_loci.chr.len() as u64, p.limit_sjdb_insert_nsj as u64)
    } else {
        genome.n_sa
    };
    packedarray_l8_packedarray_definebits(
        &mut genome.sa_insert,
        genome.gstrand_bit as u64 + 1,
        sa_pass1_len,
    );
    result
        .log_main
        .push_str(&format!("Number of SA indices: {}\n", genome.n_sa));

    let raw_time_sort = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    let sort_msg = format!(
        "{} ... starting to sort Suffix Array. This may take a long time...\n",
        timefunctions_l14_timemonthdaytime(raw_time_sort)
    );
    result.log_main.push_str(&sort_msg);
    result.log_stdout.push_str(&sort_msg);

    let sa_chunk_size_from_ram = p
        .limit_genome_generate_ram
        .saturating_sub(genome.g.len() as u64)
        / 8
        / std::cmp::max(1, p.run_thread_n as u64);
    let mut sa_chunk_size = (sa_chunk_size_from_ram * 6 / 10) as usize;
    if p.run_thread_n > 1 && genome.n_sa > 0 {
        sa_chunk_size = std::cmp::min(
            sa_chunk_size,
            genome.n_sa as usize / (p.run_thread_n as usize - 1),
        );
    }
    if sa_chunk_size == 0 {
        sa_chunk_size = std::cmp::max(1, genome.n_sa as usize);
    }

    let mut ind_pref_start = Vec::<usize>::new();
    let mut ind_pref_chunk_count = Vec::<usize>::new();
    ind_pref_start.push(0);
    let mut chunk_size1 = ind_pref_count[0];
    for (ii, count) in ind_pref_count.iter().enumerate().skip(1) {
        chunk_size1 += *count;
        if chunk_size1 > sa_chunk_size {
            ind_pref_start.push(ii);
            ind_pref_chunk_count.push(chunk_size1 - *count);
            chunk_size1 = *count;
        }
    }
    ind_pref_start.push(ind_pref_n);
    ind_pref_chunk_count.push(chunk_size1);
    let sa_chunk_n = ind_pref_chunk_count.len();
    result.log_main.push_str(&format!(
        "Number of chunks: {};   chunks size limit: {} bytes\n",
        sa_chunk_n,
        sa_chunk_size * 8
    ));

    let raw_time_chunks = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    let chunks_msg = format!(
        "{} ... sorting Suffix Array chunks and saving them to disk...\n",
        timefunctions_l14_timemonthdaytime(raw_time_chunks)
    );
    result.log_main.push_str(&chunks_msg);
    result.log_stdout.push_str(&chunks_msg);

    let global_l = g_len / std::mem::size_of::<u64>() + 2;

    // Parallelize SA chunk collection + sort + write. Each chunk is
    // independent: it reads the genome immutably and writes to its own
    // SA_<i_chunk> file. STAR's C++ uses OpenMP across chunks for the same
    // reason. We honor --runThreadN by building a scoped Rayon pool.
    let g_bytes: &[u8] = &genome.g;
    let n_genome_val = genome.n_genome;
    let g_dir = genome.p_ge.g_dir.clone();
    let ind_pref_start_ref = &ind_pref_start;
    let ind_pref_chunk_count_ref = &ind_pref_chunk_count;
    let n_threads = std::cmp::max(1, p.run_thread_n as usize);

    let sort_one_chunk = |i_chunk: usize| -> Result<(), String> {
        let pref_start = ind_pref_start_ref[i_chunk];
        let pref_end = ind_pref_start_ref[i_chunk + 1];
        let mut sa_chunk = Vec::<u64>::with_capacity(ind_pref_chunk_count_ref[i_chunk]);
        for ii in (0..g_len).step_by(sparse_d) {
            if g_bytes[ii] < 4 {
                let mut prefix = 0usize;
                for offset in 0..4 {
                    let b = ii
                        .checked_sub(offset)
                        .and_then(|pos| g_bytes.get(pos))
                        .copied()
                        .unwrap_or(GENOME_SPACING_CHAR);
                    prefix = (prefix << 4) | (b as usize & 0xF);
                }
                if prefix >= pref_start && prefix < pref_end {
                    sa_chunk.push(ii as u64);
                }
            }
        }
        if sa_chunk.len() != ind_pref_chunk_count_ref[i_chunk] {
            return Err(format!(
                "EXITING because of FATAL problem while generating the suffix array\nChunk {i_chunk} contains {} indices, expected {}\n",
                sa_chunk.len(),
                ind_pref_chunk_count_ref[i_chunk]
            ));
        }
        sa_chunk.sort_by(|a, b| {
            compare_suffix_positions_star_wordwise(g_bytes, *a as usize, *b as usize, global_l)
        });
        let chunk_file_name = format!("{}/SA_{}", g_dir, i_chunk);
        let mut chunk_file = streamfuns_l91_ofstropen(&chunk_file_name, "ERROR_OUT")?;
        let mut chunk_bytes = Vec::<u8>::with_capacity(
            std::cmp::min(sa_chunk.len(), 1 << 20) * std::mem::size_of::<u64>(),
        );
        for sa_pos in sa_chunk {
            let sa_in = 2 * n_genome_val - 1 - sa_pos;
            chunk_bytes.extend_from_slice(&sa_in.to_ne_bytes());
            if chunk_bytes.len() >= (1 << 20) * std::mem::size_of::<u64>() {
                streamfuns_l51_fstreamwritebig(
                    &mut chunk_file,
                    &chunk_bytes,
                    chunk_bytes.len() as u64,
                )
                .map_err(|e| e.to_string())?;
                chunk_bytes.clear();
            }
        }
        if !chunk_bytes.is_empty() {
            streamfuns_l51_fstreamwritebig(&mut chunk_file, &chunk_bytes, chunk_bytes.len() as u64)
                .map_err(|e| e.to_string())?;
        }
        Ok(())
    };

    if n_threads <= 1 || sa_chunk_n <= 1 {
        for i_chunk in 0..sa_chunk_n {
            sort_one_chunk(i_chunk)?;
        }
    } else {
        use rayon::iter::IntoParallelIterator;
        use rayon::iter::ParallelIterator;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build()
            .map_err(|e| e.to_string())?;
        pool.install(|| {
            (0..sa_chunk_n)
                .into_par_iter()
                .try_for_each(|i_chunk| sort_one_chunk(i_chunk))
        })?;
    }

    let raw_time_pack = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    let pack_msg = format!(
        "{} ... loading chunks from disk, packing SA...\n",
        timefunctions_l14_timemonthdaytime(raw_time_pack)
    );
    result.log_main.push_str(&pack_msg);
    result.log_stdout.push_str(&pack_msg);

    packedarray_l31_packedarray_allocatearray(&mut genome.sa_packed);
    let n2bit = 1_u64 << genome.gstrand_bit;
    genome.sa.clear();
    let mut packed_ind = 0u64;
    let mut read_buffer = vec![0u8; (1 << 20) * std::mem::size_of::<u64>()];
    for i_chunk in 0..sa_chunk_n {
        let chunk_file_name = format!("{}/SA_{}", genome.p_ge.g_dir, i_chunk);
        let mut chunk_file = std::fs::File::open(&chunk_file_name).map_err(|e| e.to_string())?;
        let mut leftover = 0usize;
        loop {
            let n_more = std::io::Read::read(&mut chunk_file, &mut read_buffer[leftover..])
                .map_err(|e| e.to_string())?;
            if n_more == 0 {
                if leftover != 0 {
                    return Err(format!(
                        "EXITING because of FATAL problem while generating the suffix array\nChunk {i_chunk} has a partial suffix-array record\n"
                    ));
                }
                break;
            }
            let avail = leftover + n_more;
            let usable = avail - (avail % std::mem::size_of::<u64>());
            for bytes in read_buffer[..usable].chunks_exact(std::mem::size_of::<u64>()) {
                let sa_in = u64::from_ne_bytes(bytes.try_into().unwrap());
                let packed = if sa_in < genome.n_genome {
                    sa_in
                } else {
                    (sa_in - genome.n_genome) | n2bit
                };
                packedarray_l17_packedarray_writepacked(&mut genome.sa_packed, packed_ind, packed);
                packed_ind += 1;
            }
            leftover = avail - usable;
            if leftover != 0 {
                read_buffer.copy_within(usable..avail, 0);
            }
        }
        let _ = std::fs::remove_file(&chunk_file_name);
    }
    if packed_ind != genome.n_sa as u64 {
        return Err(format!(
            "EXITING because of FATAL problem while generating the suffix array\nThe number of indices read from chunks = {packed_ind} is not equal to expected nSA={}\nSOLUTION: try to re-run suffix array generation, if it still does not work, report this problem to the author\n",
            genome.n_sa
        ));
    }
    genome.n_sa_byte = genome.sa_packed.length_byte;
    for ii in 0..n_genome {
        genome.g.swap(ii, g_len - 1 - ii);
    }

    let raw_time_done = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    let done_msg = format!(
        "{} ... finished generating suffix array\n",
        timefunctions_l14_timemonthdaytime(raw_time_done)
    );
    result.log_main.push_str(&done_msg);
    result.log_stdout.push_str(&done_msg);

    let sai = genomesaindex_from_packed(genome)?;
    if !sai.is_empty() {
        packedarray_l8_packedarray_definebits(
            &mut genome.sai_packed,
            genome.gstrand_bit as u64 + 3,
            sai.len() as u64,
        );
        packedarray_l31_packedarray_allocatearray(&mut genome.sai_packed);
        for (ii, value) in sai.iter().enumerate() {
            packedarray_l17_packedarray_writepacked(
                &mut genome.sai_packed,
                ii as u64,
                *value as u64,
            );
        }
    }

    genome.sjdb_n = 0;
    if p.sjdb_insert_yes {
        p.sjdb_insert_out_dir = genome.p_ge.g_dir.clone();
        p.two_pass_pass2 = false;
        let genome1 = genome_sjdb_insert_snapshot(genome);
        let sjdb_out =
            sjdbinsertjunctions_l11_sjdbinsertjunctions(p, genome, &genome1, &mut sjdb_loci)?;
        result.log_main.push_str(&sjdb_out.log_main);
        result.sjdb_insert = Some(sjdb_out);
    }

    genome.p_ge.g_file_sizes.clear();
    genome.p_ge.g_file_sizes.push(genome.n_genome);
    genome.p_ge.g_file_sizes.push(genome.sa_packed.length_byte);

    let genome_parameters = format!("{}/genomeParameters.txt", genome.p_ge.g_dir);
    genomeparameterswrite_l4_genomeparameterswrite(&genome_parameters, p, "ERROR_OUT", genome)?;
    result.files_written.push(genome_parameters);

    let raw_time_write = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    let write_msg = format!(
        "{} ... writing Genome to disk ...\n",
        timefunctions_l14_timemonthdaytime(raw_time_write)
    );
    result.log_main.push_str(&write_msg);
    result.log_stdout.push_str(&write_msg);
    genome_genomegenerate_l433_genome_writegenomesequence(genome, &genome.p_ge.g_dir)?;
    result
        .files_written
        .push(format!("{}/Genome", genome.p_ge.g_dir));

    result.log_main.push_str(&format!(
        "SA size in bytes: {}\n",
        genome.sa_packed.length_byte
    ));
    let sa_file = format!("{}/SA", genome.p_ge.g_dir);
    let mut sa_out = streamfuns_l91_ofstropen(&sa_file, "ERROR_OUT")?;
    streamfuns_l51_fstreamwritebig(
        &mut sa_out,
        &genome.sa_packed.char_array,
        genome.sa_packed.length_byte,
    )
    .map_err(|e| e.to_string())?;
    result.files_written.push(sa_file);

    let sa_index_file = format!("{}/SAindex", genome.p_ge.g_dir);
    let mut sai_out = streamfuns_l91_ofstropen(&sa_index_file, "ERROR_OUT")?;
    sai_out
        .write_all(&(genome.p_ge.g_saindex_nbases as u64).to_ne_bytes())
        .map_err(|e| e.to_string())?;
    for value in genome
        .genome_sa_index_start
        .iter()
        .take(genome.p_ge.g_saindex_nbases as usize + 1)
    {
        sai_out
            .write_all(&(*value as u64).to_ne_bytes())
            .map_err(|e| e.to_string())?;
    }
    streamfuns_l51_fstreamwritebig(
        &mut sai_out,
        &genome.sai_packed.char_array,
        genome.sai_packed.length_byte,
    )
    .map_err(|e| e.to_string())?;
    result.files_written.push(sa_index_file);

    let raw_time_finish = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    let finish_msg = format!(
        "{} ..... finished successfully\n",
        timefunctions_l14_timemonthdaytime(raw_time_finish)
    );
    result.log_main.push_str(&finish_msg);
    result.log_stdout.push_str(&finish_msg);
    Ok(result)
}

#[doc = "Original `Genome::writeChrInfo` at STAR/source/Genome_genomeGenerate.cpp:417. Args: dirOut: string"]
pub fn genome_genomegenerate_l417_genome_writechrinfo(
    genome: &crate::genome::Genome,
    dir_out: &str,
) -> Result<(), String> {
    let mut chr_n = String::new();
    let mut chr_s = String::new();
    let mut chr_l = String::new();
    let mut chr_nl = String::new();

    for ii in 0..genome.n_chr_real as usize {
        chr_n.push_str(&genome.chr_name[ii]);
        chr_n.push('\n');
        chr_s.push_str(&format!("{}\n", genome.chr_start[ii]));
        chr_l.push_str(&format!("{}\n", genome.chr_length[ii]));
        chr_nl.push_str(&format!(
            "{}\t{}\n",
            genome.chr_name[ii], genome.chr_length[ii]
        ));
    }
    chr_s.push_str(&format!(
        "{}\n",
        genome.chr_start[genome.n_chr_real as usize]
    ));

    std::fs::write(format!("{}/chrName.txt", dir_out), chr_n).map_err(|err| err.to_string())?;
    std::fs::write(format!("{}/chrStart.txt", dir_out), chr_s).map_err(|err| err.to_string())?;
    std::fs::write(format!("{}/chrLength.txt", dir_out), chr_l).map_err(|err| err.to_string())?;
    std::fs::write(format!("{}/chrNameLength.txt", dir_out), chr_nl)
        .map_err(|err| err.to_string())?;
    Ok(())
}

#[doc = "Original `Genome::writeGenomeSequence` at STAR/source/Genome_genomeGenerate.cpp:433. Args: dirOut: string"]
pub fn genome_genomegenerate_l433_genome_writegenomesequence(
    genome: &crate::genome::Genome,
    dir_out: &str,
) -> Result<(), String> {
    std::fs::write(
        format!("{}/Genome", dir_out),
        &genome.g[..genome.n_genome as usize],
    )
    .map_err(|err| err.to_string())
}
