#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Genome::insertSequences` at STAR/source/Genome_insertSequences.cpp:9. Args: "]
pub fn genome_insertsequences_l9_genome_insertsequences(
    genome: &mut crate::genome::Genome,
    p: &crate::parameters_chimeric::Parameters,
) -> Result<u64, String> {
    if genome.p_ge.g_fasta_files.first().map(String::as_str) == Some("-") {
        return Ok(0);
    }

    let chr_start_back = *genome.chr_start.last().ok_or_else(|| {
        "EXITING because of FATAL ERROR: cannot insert sequences without chrStart entries\n"
            .to_string()
    })?;
    let source_start = chr_start_back
        .checked_sub(genome.genome_insert_l)
        .ok_or_else(|| {
            "EXITING because of FATAL ERROR: genomeInsertL exceeds inserted chromosome start\n"
                .to_string()
        })?;
    let sjdb_len = genome.n_genome.checked_sub(source_start).ok_or_else(|| {
        "EXITING because of FATAL ERROR: inconsistent genome size for sequence insertion\n"
            .to_string()
    })?;

    let dest_start_usize = chr_start_back as usize;
    let source_start_usize = source_start as usize;
    let sjdb_len_usize = sjdb_len as usize;
    if genome.g.len() < dest_start_usize + sjdb_len_usize {
        genome
            .g
            .resize(dest_start_usize + sjdb_len_usize, GENOME_SPACING_CHAR);
    }
    genome.g.copy_within(
        source_start_usize..source_start_usize + sjdb_len_usize,
        dest_start_usize,
    );
    for gg in
        &mut genome.g[source_start_usize..source_start_usize + genome.genome_insert_l as usize]
    {
        *gg = GENOME_SPACING_CHAR;
    }

    let mut g = std::mem::take(&mut genome.g);
    genomescanfastafiles_l5_genomescanfastafiles(p, &mut g[source_start_usize..], true, genome)?;
    genome.g = g;

    let n_genome_old = genome.n_genome;
    genome.n_genome = chr_start_back + sjdb_len;
    let g1_start = source_start_usize;
    let g1_end = g1_start + genome.genome_insert_l as usize;
    let g1 = genome.g[g1_start..g1_end].to_vec();
    let map_gen_snapshot = genome.clone();

    insertseqsa_l18_insertseqsa(
        &mut genome.sa_packed,
        &mut genome.sa_insert,
        &mut genome.sai_packed,
        &genome.g,
        &g1,
        n_genome_old - sjdb_len,
        genome.genome_insert_l,
        sjdb_len,
        p,
        &map_gen_snapshot,
    )
}
