#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `genomeParametersWrite` at STAR/source/genomeParametersWrite.cpp:4. Args: fileName: string, P: Parameters, errorOut: string, mapGen: Genome"]
pub fn genomeparameterswrite_l4_genomeparameterswrite(
    file_name: &str,
    p: &crate::parameters_chimeric::Parameters,
    error_out: &str,
    map_gen: &crate::genome::Genome,
) -> Result<(), String> {
    use std::io::Write;

    let mut genome_par = streamfuns_l91_ofstropen(file_name, error_out)?;

    writeln!(genome_par, "### {}", p.command_line_full).map_err(|e| e.to_string())?;
    writeln!(genome_par, "### GstrandBit {}", map_gen.gstrand_bit).map_err(|e| e.to_string())?;

    writeln!(genome_par, "versionGenome\t{}", p.version_genome).map_err(|e| e.to_string())?;
    writeln!(genome_par, "genomeType\t{}", map_gen.p_ge.g_type_string)
        .map_err(|e| e.to_string())?;

    write!(genome_par, "genomeFastaFiles\t").map_err(|e| e.to_string())?;
    for fasta_file in &map_gen.p_ge.g_fasta_files {
        write!(genome_par, "{} ", fasta_file).map_err(|e| e.to_string())?;
    }
    writeln!(genome_par).map_err(|e| e.to_string())?;
    writeln!(
        genome_par,
        "genomeSAindexNbases\t{}",
        map_gen.p_ge.g_saindex_nbases
    )
    .map_err(|e| e.to_string())?;
    writeln!(
        genome_par,
        "genomeChrBinNbits\t{}",
        map_gen.p_ge.g_chr_bin_nbits
    )
    .map_err(|e| e.to_string())?;
    writeln!(genome_par, "genomeSAsparseD\t{}", map_gen.p_ge.g_sasparse_d)
        .map_err(|e| e.to_string())?;

    writeln!(
        genome_par,
        "genomeTransformType\t{}",
        map_gen.p_ge.transform.type_string
    )
    .map_err(|e| e.to_string())?;
    writeln!(
        genome_par,
        "genomeTransformVCF\t{}",
        map_gen.p_ge.transform.vcf_file
    )
    .map_err(|e| e.to_string())?;

    writeln!(genome_par, "sjdbOverhang\t{}", map_gen.sjdb_overhang).map_err(|e| e.to_string())?;
    write!(genome_par, "sjdbFileChrStartEnd\t").map_err(|e| e.to_string())?;
    for sjdb_file in &map_gen.p_ge.sjdb_file_chr_start_end {
        write!(genome_par, "{} ", sjdb_file).map_err(|e| e.to_string())?;
    }
    writeln!(genome_par).map_err(|e| e.to_string())?;

    writeln!(genome_par, "sjdbGTFfile\t{}", map_gen.p_ge.sjdb_gtf_file)
        .map_err(|e| e.to_string())?;
    writeln!(
        genome_par,
        "sjdbGTFchrPrefix\t{}",
        map_gen.p_ge.sjdb_gtf_chr_prefix
    )
    .map_err(|e| e.to_string())?;
    writeln!(
        genome_par,
        "sjdbGTFfeatureExon\t{}",
        map_gen.p_ge.sjdb_gtf_feature_exon
    )
    .map_err(|e| e.to_string())?;
    writeln!(
        genome_par,
        "sjdbGTFtagExonParentTranscript\t{}",
        map_gen.p_ge.sjdb_gtf_tag_exon_parent_transcript
    )
    .map_err(|e| e.to_string())?;
    writeln!(
        genome_par,
        "sjdbGTFtagExonParentGene\t{}",
        map_gen.p_ge.sjdb_gtf_tag_exon_parent_gene
    )
    .map_err(|e| e.to_string())?;

    writeln!(
        genome_par,
        "sjdbInsertSave\t{}",
        map_gen.p_ge.sjdb_insert_save
    )
    .map_err(|e| e.to_string())?;

    write!(
        genome_par,
        "genomeFileSizes\t{}",
        map_gen.p_ge.g_file_sizes[0]
    )
    .map_err(|e| e.to_string())?;
    for file_size in map_gen.p_ge.g_file_sizes.iter().skip(1) {
        write!(genome_par, " {}", file_size).map_err(|e| e.to_string())?;
    }
    writeln!(genome_par).map_err(|e| e.to_string())?;

    Ok(())
}
