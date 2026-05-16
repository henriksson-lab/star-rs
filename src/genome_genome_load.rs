#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Genome::genomeLoad` at STAR/source/Genome_genomeLoad.cpp:18. Args: "]
pub fn genome_genomeload_l18_genome_genomeload(
    genome: &mut crate::genome::Genome,
    p: &mut crate::parameters_chimeric::Parameters,
    genome_parameters_contents: Option<&str>,
    chr_name_contents: &str,
    chr_length_contents: &str,
    chr_start_contents: &str,
    genome_contents: Vec<u8>,
    sa_contents: Vec<u8>,
    mut sa_index_contents: Vec<u8>,
    sjdb_info_contents: Option<&str>,
    raw_time: libc::time_t,
) -> Result<crate::parameters_chimeric::GenomeLoadResult, String> {
    let mut result = crate::parameters_chimeric::GenomeLoadResult::default();
    result.log_stdout.push_str(&format!(
        "{} ..... loading genome\n",
        timefunctions_l14_timemonthdaytime(raw_time)
    ));

    if genome.p_ge.g_load == "Remove" {
        result.removed_shared_memory = true;
        result
            .log_main
            .push_str("DONE: removed the genome from shared memory\n");
        return Ok(result);
    }

    genome.gstrand_bit = 0;
    let genome_parameters_contents = genome_parameters_contents.ok_or_else(|| {
        format!(
            "EXITING because of FATAL ERROR: could not open genome file {}/genomeParameters.txt\nSOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permsissions\n",
            genome.p_ge.g_dir
        )
    })?;

    result
        .log_main
        .push_str("Reading genome generation parameters:\n");
    let mut p1 = crate::parameters_chimeric::Parameters::default();
    for line in genome_parameters_contents.lines() {
        let mut words = line.split_whitespace();
        let Some(word1) = words.next() else {
            continue;
        };
        if word1 == "###" {
            let Some(word2) = words.next() else {
                continue;
            };
            if word2 == "GstrandBit" {
                genome.gstrand_bit = words
                    .next()
                    .and_then(|w| w.parse::<u32>().ok())
                    .unwrap_or(0);
                result
                    .log_main
                    .push_str(&format!("### GstrandBit={}\n", genome.gstrand_bit));
            } else {
                result.log_main.push_str("### ");
                result.log_main.push_str(word2);
                let rest = words.collect::<Vec<_>>().join(" ");
                if !rest.is_empty() {
                    result.log_main.push(' ');
                    result.log_main.push_str(&rest);
                }
                result.log_main.push('\n');
            }
            continue;
        }

        match word1 {
            "versionGenome" => {
                p1.version_genome = words.next().unwrap_or("").to_string();
            }
            "genomeType" => {
                p1.p_ge.g_type_string = words.next().unwrap_or("").to_string();
            }
            "genomeSAindexNbases" => {
                p1.p_ge.g_saindex_nbases = words
                    .next()
                    .and_then(|w| w.parse::<u32>().ok())
                    .unwrap_or(0);
            }
            "genomeChrBinNbits" => {
                p1.p_ge.g_chr_bin_nbits = words
                    .next()
                    .and_then(|w| w.parse::<u32>().ok())
                    .unwrap_or(0);
            }
            "genomeSAsparseD" => {
                p1.p_ge.g_sasparse_d = words
                    .next()
                    .and_then(|w| w.parse::<u32>().ok())
                    .unwrap_or(0);
            }
            "sjdbOverhang" => {
                p1.p_ge.sjdb_overhang = words
                    .next()
                    .and_then(|w| w.parse::<u32>().ok())
                    .unwrap_or(0);
            }
            "sjdbInsertSave" => {
                p1.p_ge.sjdb_insert_save = words.next().unwrap_or("").to_string();
            }
            "genomeTransformType" => {
                p1.p_ge.transform.type_string = words.next().unwrap_or("").to_string();
            }
            "genomeFileSizes" => {
                p1.p_ge.g_file_sizes = words.filter_map(|w| w.parse::<u64>().ok()).collect();
            }
            _ => {}
        }
    }

    if p1.version_genome.is_empty() {
        return Err("EXITING because of FATAL ERROR: read no value for the versionGenome parameter from genomeParameters.txt file\nSOLUTION: please re-generate genome from scratch with the latest version of STAR\n".to_string());
    } else if p1.version_genome == p.version_genome {
        result
            .log_main
            .push_str("Genome version is compatible with current STAR\n");
    } else {
        return Err(format!(
            "EXITING because of FATAL ERROR: Genome version: {} is INCOMPATIBLE with running STAR version: {}\nSOLUTION: please re-generate genome from scratch with running version of STAR, or with version: {}\n",
            p1.version_genome, p.version_genome, p.version_genome
        ));
    }

    genome.chr_name.clear();
    for line in chr_name_contents.lines() {
        if line.is_empty() {
            break;
        }
        genome.chr_name.push(line.to_string());
    }
    genome.n_chr_real = genome.chr_name.len() as u32;
    genome.chr_length = chr_length_contents
        .split_whitespace()
        .take(genome.n_chr_real as usize)
        .map(|w| w.parse::<u64>().unwrap_or(0))
        .collect();
    genome.chr_start = chr_start_contents
        .split_whitespace()
        .take(genome.n_chr_real as usize + 1)
        .map(|w| w.parse::<u64>().unwrap_or(0))
        .collect();
    result.log_main.push_str(&format!(
        "Number of real (reference) chromosomes= {}\n",
        genome.n_chr_real
    ));
    genome.chr_name_index.clear();
    for ii in 0..genome.n_chr_real as usize {
        result.log_main.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            ii + 1,
            genome.chr_name[ii],
            genome.chr_length[ii],
            genome.chr_start[ii]
        ));
        genome
            .chr_name_index
            .insert(genome.chr_name[ii].clone(), ii as u64);
    }

    if p.sjdb_insert_yes && sjdb_info_contents.is_some() && p1.p_ge.sjdb_insert_save.is_empty() {
        return Err("EXITING because of FATAL ERROR: old Genome is INCOMPATIBLE with on the fly junction insertion\nSOLUTION: please re-generate genome from scratch with the latest version of STAR\n".to_string());
    }

    genome.p_ge.g_saindex_nbases = p1.p_ge.g_saindex_nbases;
    genome.p_ge.g_chr_bin_nbits = p1.p_ge.g_chr_bin_nbits;
    genome.genome_chr_bin_nbases = 1_u32.checked_shl(genome.p_ge.g_chr_bin_nbits).unwrap_or(0);
    genome.p_ge.g_sasparse_d = p1.p_ge.g_sasparse_d;
    p.p_ge.g_saindex_nbases = p1.p_ge.g_saindex_nbases;
    p.p_ge.g_chr_bin_nbits = p1.p_ge.g_chr_bin_nbits;
    p.p_ge.g_sasparse_d = p1.p_ge.g_sasparse_d;
    if !p1.p_ge.g_file_sizes.is_empty() {
        genome.p_ge.g_file_sizes = p1.p_ge.g_file_sizes.clone();
        p.p_ge.g_file_sizes = p1.p_ge.g_file_sizes.clone();
    }
    let sjdb_overhang_is_default = sjdb_info_contents.is_some() && genome.p_ge.sjdb_overhang == 100;
    if (genome.p_ge.sjdb_overhang == 0 || sjdb_overhang_is_default) && p1.p_ge.sjdb_overhang > 0 {
        genome.p_ge.sjdb_overhang = p1.p_ge.sjdb_overhang;
        p.p_ge.sjdb_overhang = p1.p_ge.sjdb_overhang;
        result.log_main.push_str(&format!(
            "--sjdbOverhang = {} taken from the generated genome\n",
            genome.p_ge.sjdb_overhang
        ));
    } else if sjdb_info_contents.is_some()
        && genome.p_ge.sjdb_overhang > 0
        && genome.p_ge.sjdb_overhang != p1.p_ge.sjdb_overhang
    {
        return Err(format!(
            "EXITING because of fatal PARAMETERS error: present --sjdbOverhang={} is not equal to the value at the genome generation step ={}\nSOLUTION: \n",
            genome.p_ge.sjdb_overhang, p1.p_ge.sjdb_overhang
        ));
    }
    genome.sjdb_overhang = genome.p_ge.sjdb_overhang;
    genome.sjdb_length = if genome.p_ge.sjdb_overhang == 0 {
        0
    } else {
        genome.p_ge.sjdb_overhang * 2 + 1
    };
    genome.p_ge.g_type_string = if p1.p_ge.g_type_string.is_empty() {
        "Full".to_string()
    } else {
        p1.p_ge.g_type_string.clone()
    };
    p.p_ge.g_type_string = genome.p_ge.g_type_string.clone();

    result.log_main.push_str("Started loading the genome: \n");
    genome.n_genome = genome_contents.len() as u64;
    genome.n_sa_byte = sa_contents.len() as u32;

    if sa_index_contents.len() < 8 {
        return Err("EXITING because of FATAL ERROR: failed reading from genome file: SAindex\nSOLUTION: re-generate the genome index\n".to_string());
    }
    let mut first8 = [0u8; 8];
    first8.copy_from_slice(&sa_index_contents[..8]);
    genome.p_ge.g_saindex_nbases = u64::from_ne_bytes(first8) as u32;
    p.p_ge.g_saindex_nbases = genome.p_ge.g_saindex_nbases;
    let index_n = genome.p_ge.g_saindex_nbases as usize + 1;
    let legacy_u32_saindex = if sa_index_contents.len() >= 8 + index_n * 8 {
        let mut start0 = [0u8; 8];
        start0.copy_from_slice(&sa_index_contents[8..16]);
        u64::from_ne_bytes(start0) != 0
    } else {
        true
    };
    let (sa_index_header_bytes, sa_index_word_bytes) = if legacy_u32_saindex {
        (4usize, 4usize)
    } else {
        (8usize, 8usize)
    };
    if sa_index_contents.len() < sa_index_header_bytes + index_n * sa_index_word_bytes {
        return Err("EXITING because of FATAL ERROR: failed reading from genome file: SAindex\nSOLUTION: re-generate the genome index\n".to_string());
    }
    genome.genome_sa_index_start = vec![0; index_n];
    for ii in 0..index_n {
        let start = sa_index_header_bytes + ii * sa_index_word_bytes;
        genome.genome_sa_index_start[ii] = if legacy_u32_saindex {
            let mut b = [0u8; 4];
            b.copy_from_slice(&sa_index_contents[start..start + 4]);
            u32::from_ne_bytes(b)
        } else {
            let mut b = [0u8; 8];
            b.copy_from_slice(&sa_index_contents[start..start + 8]);
            u64::from_ne_bytes(b) as u32
        };
    }
    let n_sai = *genome.genome_sa_index_start.last().unwrap_or(&0);
    result.log_main.push_str(&format!(
        "Read from SAindex: pGe.gSAindexNbases={}  nSAi={}\n",
        genome.p_ge.g_saindex_nbases, n_sai
    ));

    if genome.gstrand_bit == 0 {
        genome.gstrand_bit = if genome.n_genome == 0 {
            32
        } else {
            ((genome.n_genome as f64).log2().floor() as u32 + 1).max(32)
        };
    }
    genome.gstrand_mask = !(1_u64 << genome.gstrand_bit);
    genome.n_sa = if genome.gstrand_bit + 1 == 0 {
        0
    } else {
        (genome.n_sa_byte as u64 * 8 / (genome.gstrand_bit as u64 + 1)) as u32
    };
    packedarray_l8_packedarray_definebits(
        &mut genome.sa_packed,
        genome.gstrand_bit as u64 + 1,
        genome.n_sa as u64,
    );
    genome.sa_packed.char_array = sa_contents;
    if genome.sa_packed.char_array.len() < genome.sa_packed.length_byte as usize {
        genome
            .sa_packed
            .char_array
            .resize(genome.sa_packed.length_byte as usize, 0);
    }
    genome.sa_packed.array_allocated = true;
    genome.sa.clear();

    let sai_bytes_start = sa_index_header_bytes + index_n * sa_index_word_bytes;
    packedarray_l8_packedarray_definebits(
        &mut genome.sai_packed,
        genome.gstrand_bit as u64 + 3,
        n_sai as u64,
    );
    genome.sai_packed.char_array = sa_index_contents.split_off(sai_bytes_start);
    if genome.sai_packed.char_array.len() < genome.sai_packed.length_byte as usize {
        genome
            .sai_packed
            .char_array
            .resize(genome.sai_packed.length_byte as usize, 0);
    }
    genome.sai_packed.array_allocated = true;
    genome.sai = (0..n_sai as u64)
        .map(|ii| packedarray_h18_packedarray_index(&genome.sai_packed, ii))
        .collect();
    genome.sai_mark_nmask_c = 1_u64 << (genome.gstrand_bit + 1);
    genome.sai_mark_nmask = !genome.sai_mark_nmask_c;
    genome.sai_mark_absent_mask_c = 1_u64 << (genome.gstrand_bit + 2);

    result.log_main.push_str(&format!(
        "nGenome={};  nSAbyte={}\nGstrandBit={}   SA number of indices={}\n",
        genome.n_genome, genome.n_sa_byte, genome.gstrand_bit, genome.n_sa
    ));

    genome.g = genome_contents;
    result.log_main.push_str(&format!(
        "Genome file size: {} bytes; state: good=1 eof=0 fail=0 bad=0\nLoading Genome ... done! state: good=1 eof=0 fail=0 bad=0; loaded {} bytes\n",
        genome.n_genome, genome.n_genome
    ));
    result.log_main.push_str(&format!(
        "SA file size: {} bytes; state: good=1 eof=0 fail=0 bad=0\nLoading SA ... done! state: good=1 eof=0 fail=0 bad=0; loaded {} bytes\nLoading SAindex ... done: {} bytes\n",
        genome.sa_packed.length_byte,
        genome.n_sa_byte,
        sa_index_header_bytes + index_n * sa_index_word_bytes + genome.sai_packed.char_array.len()
    ));

    result.log_main.push_str("Finished loading the genome: \n");
    if genome.p_ge.g_load == "LoadAndExit" {
        result.load_and_exit = true;
        result.log_main.push_str(
            "pGe.gLoad=LoadAndExit: completed, the genome is loaded and kept in RAM, EXITING now.\n",
        );
        return Ok(result);
    }

    if genome.p_ge.g_fasta_files.first().is_some_and(|f| f != "-") {
        genome_insertsequences_l9_genome_insertsequences(genome, p)?;
    }
    genome_l209_genome_chrbinfill(genome);
    result
        .log_main
        .push_str(&genome_genomeload_l471_genome_loadsjdb(
            genome,
            sjdb_info_contents,
        )?);

    if p.align_intron_max == 0 && p.align_mates_gap_max == 0 {
        result.log_main.push_str(&format!(
            "alignIntronMax=alignMatesGapMax=0, the max intron size will be approximately determined by (2^winBinNbits)*winAnchorDistNbins={}\n",
            (1_u64 << p.win_bin_nbits) * p.win_anchor_dist_nbins as u64
        ));
    } else {
        let mates_gap = if p.align_mates_gap_max == 0 {
            1000
        } else {
            p.align_mates_gap_max
        };
        let max_gap = std::cmp::max(std::cmp::max(4, p.align_intron_max), mates_gap);
        p.win_bin_nbits = ((max_gap as f64 / 4.0).log2().floor() + 0.5) as u32;
        p.win_bin_nbits = p
            .win_bin_nbits
            .max(((genome.n_genome / 40000 + 1) as f64).log2().floor() as u32);
        result.log_main.push_str(&format!(
            "To accommodate alignIntronMax={} redefined winBinNbits={}\n",
            p.align_intron_max, p.win_bin_nbits
        ));
    }
    if p.win_bin_nbits > genome.p_ge.g_chr_bin_nbits {
        result.log_main.push_str(&format!(
            "winBinNbits={} > pGe.gChrBinNbits={}   redefining:\n",
            p.win_bin_nbits, genome.p_ge.g_chr_bin_nbits
        ));
        p.win_bin_nbits = genome.p_ge.g_chr_bin_nbits;
        result
            .log_main
            .push_str(&format!("winBinNbits={}\n", p.win_bin_nbits));
    }
    if p.align_intron_max != 0 || p.align_mates_gap_max != 0 {
        p.win_flank_nbins = std::cmp::max(p.align_intron_max, p.align_mates_gap_max)
            / (1_u32 << p.win_bin_nbits)
            + 1;
        p.win_anchor_dist_nbins = 2 * p.win_flank_nbins;
        result.log_main.push_str(&format!(
            "To accommodate alignIntronMax={} and alignMatesGapMax={}, redefined winFlankNbins={} and winAnchorDistNbins={}\n",
            p.align_intron_max, p.align_mates_gap_max, p.win_flank_nbins, p.win_anchor_dist_nbins
        ));
    }
    p.win_bin_chr_nbits = genome.p_ge.g_chr_bin_nbits - p.win_bin_nbits;
    p.win_bin_n = (genome.n_genome / (1_u64 << p.win_bin_nbits)) as u32 + 1;

    genome.p_ge.transform.type_ = match p1.p_ge.transform.type_string.as_str() {
        "Haploid" => 1,
        "Diploid" => 2,
        _ => 0,
    };
    p.p_ge.transform.type_ = genome.p_ge.transform.type_;
    if genome.p_ge.transform.out_yes && genome.p_ge.transform.type_ == 0 {
        return Err(format!(
            "EXITING because of FATAL INPUT ERROR: outTransformOutput is set, but the genome {} was generated without transformation\nSOLUTION: use the default--outTransformOutput None, or re-generate the genome with transformation options.\n",
            genome.p_ge.g_dir
        ));
    }

    Ok(result)
}

#[doc = "Original `Genome::loadSJDB` at STAR/source/Genome_genomeLoad.cpp:471. Args: genDir: string"]
pub fn genome_genomeload_l471_genome_loadsjdb(
    genome: &mut crate::genome::Genome,
    sjdb_info_contents: Option<&str>,
) -> Result<String, String> {
    if genome.n_genome == genome.chr_start[genome.n_chr_real as usize] {
        genome.sjdb_n = 0;
        genome.sj_gstart = genome.chr_start[genome.n_chr_real as usize] + 1;
        return Ok(String::new());
    }

    let sjdb_info_contents = sjdb_info_contents.ok_or_else(|| {
        format!(
            "EXITING because of FATAL error, could not open file {}/sjdbInfo.txt\nSOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permsissions\n",
            genome.p_ge.g_dir
        )
    })?;
    let mut words = sjdb_info_contents.split_whitespace();
    genome.sjdb_n = words
        .next()
        .and_then(|w| w.parse::<u32>().ok())
        .unwrap_or(0);
    genome.p_ge.sjdb_overhang = words
        .next()
        .and_then(|w| w.parse::<u32>().ok())
        .unwrap_or(0);
    let log_main = format!(
        "Processing splice junctions database sjdbN={},   pGe.sjdbOverhang={} \n",
        genome.sjdb_n, genome.p_ge.sjdb_overhang
    );

    genome.sj_gstart = genome.chr_start[genome.n_chr_real as usize];
    genome.sj_dstart = vec![0; genome.sjdb_n as usize];
    genome.sj_astart = vec![0; genome.sjdb_n as usize];
    genome.sjdb_start = vec![0; genome.sjdb_n as usize];
    genome.sjdb_end = vec![0; genome.sjdb_n as usize];
    genome.sjdb_motif = vec![0; genome.sjdb_n as usize];
    genome.sjdb_shift_left = vec![0; genome.sjdb_n as usize];
    genome.sjdb_shift_right = vec![0; genome.sjdb_n as usize];
    genome.sjdb_strand = vec![0; genome.sjdb_n as usize];

    for ii in 0..genome.sjdb_n as usize {
        genome.sjdb_start[ii] = words
            .next()
            .and_then(|w| w.parse::<u32>().ok())
            .unwrap_or(0);
        genome.sjdb_end[ii] = words
            .next()
            .and_then(|w| w.parse::<u32>().ok())
            .unwrap_or(0);
        genome.sjdb_motif[ii] = words.next().and_then(|w| w.parse::<u8>().ok()).unwrap_or(0);
        genome.sjdb_shift_left[ii] = words.next().and_then(|w| w.parse::<u8>().ok()).unwrap_or(0);
        genome.sjdb_shift_right[ii] = words.next().and_then(|w| w.parse::<u8>().ok()).unwrap_or(0);
        genome.sjdb_strand[ii] = words.next().and_then(|w| w.parse::<u8>().ok()).unwrap_or(0);

        genome.sj_dstart[ii] = genome.sjdb_start[ii] as u64 - genome.p_ge.sjdb_overhang as u64;
        genome.sj_astart[ii] = genome.sjdb_end[ii] as u64 + 1;
        if genome.sjdb_motif[ii] == 0 {
            genome.sj_dstart[ii] += genome.sjdb_shift_left[ii] as u64;
            genome.sj_astart[ii] += genome.sjdb_shift_left[ii] as u64;
        }
    }

    Ok(log_main)
}
