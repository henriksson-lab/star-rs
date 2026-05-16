#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `genomeScanFastaFiles` at STAR/source/genomeScanFastaFiles.cpp:5. Args: P: Parameters, G: char, flagRun: bool, mapGen: Genome"]
pub fn genomescanfastafiles_l5_genomescanfastafiles(
    _p: &crate::parameters_chimeric::Parameters,
    g: &mut [u8],
    flag_run: bool,
    map_gen: &mut crate::genome::Genome,
) -> Result<u32, String> {
    let mut n: u32 = 0;
    if !flag_run && !map_gen.chr_length.is_empty() {
        map_gen.chr_start.pop();
        n = (map_gen.chr_start.last().copied().unwrap_or(0)
            + map_gen.chr_length.last().copied().unwrap_or(0)) as u32;
        map_gen.chr_length.pop();
    }

    for fasta_file in map_gen.p_ge.g_fasta_files.clone() {
        let contents = std::fs::read_to_string(&fasta_file).map_err(|_| {
            format!(
                "EXITING because of INPUT ERROR: could not open genomeFastaFile: {}\n",
                fasta_file
            )
        })?;
        let first = contents.as_bytes().first().ok_or_else(|| {
            format!(
                "EXITING because of INPUT ERROR: could not read from genomeFastaFile: {}\n",
                fasta_file
            )
        })?;
        if *first != b'>' {
            return Err(format!(
                "EXITING because of INPUT ERROR: the file format of the genomeFastaFile: {} is not fasta: the first character is '{}' ({}), not '>'.\n Solution: check formatting of the fasta file. Make sure the file is uncompressed (unzipped).\n",
                fasta_file, *first as char, *first
            ));
        }

        for line in contents.lines() {
            if line.as_bytes().first() == Some(&b'>') {
                if !flag_run {
                    let chr_name = line[1..]
                        .split_whitespace()
                        .next()
                        .unwrap_or("")
                        .to_string();
                    map_gen.chr_name.push(chr_name);
                }

                if !flag_run && !map_gen.chr_start.is_empty() {
                    let last_start = *map_gen.chr_start.last().unwrap();
                    map_gen.chr_length.push(n as u64 - last_start);
                }

                if n > 0 {
                    n = ((n + 1) / map_gen.genome_chr_bin_nbases + 1)
                        * map_gen.genome_chr_bin_nbases;
                }

                if !flag_run {
                    map_gen.chr_start.push(n as u64);
                }
            } else if flag_run {
                let wrote = sequencefuns_l170_convertnucleotidestonumbersremovecontrols(
                    line.as_bytes(),
                    &mut g[n as usize..],
                    line.len() as u32,
                );
                n += wrote;
            } else {
                for &b in line.as_bytes() {
                    if b >= 32 {
                        n += 1;
                    }
                }
            }
        }
    }

    if !flag_run {
        let last_start = *map_gen.chr_start.last().unwrap_or(&0);
        map_gen.chr_length.push(n as u64 - last_start);
    }

    n = ((n + 1) / map_gen.genome_chr_bin_nbases + 1) * map_gen.genome_chr_bin_nbases;

    if !flag_run {
        map_gen.n_chr_real = map_gen.chr_start.len() as u32;
        map_gen.chr_start.push(n as u64);
        for ii in 0..map_gen.n_chr_real as usize {
            map_gen
                .chr_name_index
                .insert(map_gen.chr_name[ii].clone(), ii as u64);
        }
    }

    Ok(n)
}
