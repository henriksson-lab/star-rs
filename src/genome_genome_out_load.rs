#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Genome::genomeOutLoad` at STAR/source/Genome_genomeOutLoad.cpp:8. Args: "]
pub fn genome_genomeoutload_l8_genome_genomeoutload(
    genome: &mut crate::genome::Genome,
    genome_parameters_contents: Option<&str>,
    chr_name_contents: &str,
    chr_length_contents: &str,
    chr_start_contents: &str,
    genome_contents: &[u8],
    sjdb_info_contents: Option<&str>,
    transform_blocks_contents: &str,
) -> Result<String, String> {
    let genome_parameters_contents = genome_parameters_contents.ok_or_else(|| {
        format!(
            "EXITING because of FATAL ERROR: could not open genome file {}/genomeParameters.txt\nSOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permsissions\n",
            genome.p_ge.g_dir
        )
    })?;

    let mut log_main = String::from("Reading output genome generation parameters:\n");
    let mut words = genome_parameters_contents.split_whitespace();
    while let Some(key) = words.next() {
        let Some(value) = words.next() else {
            break;
        };
        match key {
            "genomeSAindexNbases" => {
                genome.p_ge.g_saindex_nbases = value.parse::<u32>().unwrap_or(0);
            }
            "genomeChrBinNbits" => {
                genome.p_ge.g_chr_bin_nbits = value.parse::<u32>().unwrap_or(0);
            }
            "genomeSAsparseD" => {
                genome.p_ge.g_sasparse_d = value.parse::<u32>().unwrap_or(0);
            }
            _ => {}
        }
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
        .map(|v| v.parse::<u64>().unwrap_or(0))
        .collect();
    genome.chr_start = chr_start_contents
        .split_whitespace()
        .take(genome.n_chr_real as usize + 1)
        .map(|v| v.parse::<u64>().unwrap_or(0))
        .collect();

    log_main.push_str(&format!(
        "Number of real (reference) chromosomes= {}\n",
        genome.n_chr_real
    ));
    genome.chr_name_index.clear();
    for ii in 0..genome.n_chr_real as usize {
        log_main.push_str(&format!(
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
    genome.p_ge.chr_set_mito.clear();
    for cm in &genome.p_ge.chr_set_mito_strings {
        let ind1 = genome
            .chr_name
            .iter()
            .position(|chr_name| chr_name == cm)
            .unwrap_or(genome.chr_name.len()) as u64;
        genome.p_ge.chr_set_mito.insert(ind1);
    }

    genome.n_genome = genome_contents.len() as u64;
    genome.g = genome_contents.to_vec();
    log_main.push_str(&genome_genomeload_l471_genome_loadsjdb(
        genome,
        sjdb_info_contents,
    )?);

    genome.genome_chr_bin_nbases = 1_u32.checked_shl(genome.p_ge.g_chr_bin_nbits).unwrap_or(0);
    genome_l209_genome_chrbinfill(genome);

    let mut conv_words = transform_blocks_contents.split_whitespace();
    let nconv = conv_words
        .next()
        .and_then(|w| w.parse::<u32>().ok())
        .unwrap_or(0);
    genome.genome_out.n_minus_strand_offset = conv_words
        .next()
        .and_then(|w| w.parse::<i128>().ok())
        .map(|w| w as u64)
        .unwrap_or(0);
    genome.genome_out.conv_blocks = vec![[0; 3]; nconv as usize + 1];
    for ii in 0..nconv as usize {
        genome.genome_out.conv_blocks[ii][0] = conv_words
            .next()
            .and_then(|w| w.parse::<u64>().ok())
            .unwrap_or(0);
        genome.genome_out.conv_blocks[ii][1] = conv_words
            .next()
            .and_then(|w| w.parse::<u64>().ok())
            .unwrap_or(0);
        genome.genome_out.conv_blocks[ii][2] = conv_words
            .next()
            .and_then(|w| w.parse::<u64>().ok())
            .unwrap_or(0);
    }
    if nconv > 0 {
        genome.genome_out.conv_blocks[nconv as usize - 1][1] += 1;
    }
    genome.genome_out.conv_blocks[nconv as usize][0] = u64::MAX;

    Ok(log_main)
}
