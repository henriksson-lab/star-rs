#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ParametersGenome` at STAR/source/ParametersGenome.h:9."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParametersGenome {
    pub g_dir: String,
    pub transform: ParametersGenomeTransform,
    pub g_type_string: String,
    pub g_fasta_files: Vec<String>,
    pub g_chain_files: Vec<String>,
    pub g_load: String,
    pub g_chr_bin_nbits: u32,
    pub g_saindex_nbases: u32,
    pub g_sasparse_d: u32,
    pub g_suffix_length_max: u32,
    pub sjdb_overhang: u32,
    pub sjdb_file_chr_start_end: Vec<String>,
    pub sjdb_gtf_file: String,
    pub sjdb_gtf_chr_prefix: String,
    pub sjdb_gtf_feature_exon: String,
    pub sjdb_gtf_tag_exon_parent_transcript: String,
    pub sjdb_gtf_tag_exon_parent_gene: String,
    pub sjdb_gtf_tag_exon_parent_gene_name: Vec<String>,
    pub sjdb_gtf_tag_exon_parent_gene_type: Vec<String>,
    pub sjdb_insert_save: String,
    pub g_file_sizes: Vec<u64>,
    pub sjdb_score: i32,
    pub chr_set_mito_strings: Vec<String>,
    pub chr_set_mito: std::collections::BTreeSet<u64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParametersGenomeTransform {
    pub type_string: String,
    pub type_: i32,
    pub vcf_file: String,
    pub output: Vec<String>,
    pub out_yes: bool,
    pub out_sam: bool,
    pub out_sj: bool,
    pub out_quant: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SamAttrPresent {
    pub nh: bool,
    pub hi: bool,
    pub as_: bool,
    pub nm: bool,
    pub md: bool,
    pub n_m: bool,
    pub j_m: bool,
    pub j_i: bool,
    pub rg: bool,
    pub mc: bool,
    pub xs: bool,
    pub ch: bool,
    pub v_a: bool,
    pub v_g: bool,
    pub v_w: bool,
    pub r_b: bool,
    pub ha: bool,
    pub cr: bool,
    pub cy: bool,
    pub ur: bool,
    pub uy: bool,
    pub cb: bool,
    pub ub: bool,
    pub gx: bool,
    pub gn: bool,
    pub gx_lower: bool,
    pub gn_lower: bool,
    pub s_m: bool,
    pub s_s: bool,
    pub s_q: bool,
    pub s_f: bool,
    pub c_n: bool,
}

#[doc = "Original `ParametersGenome::initialize` at STAR/source/ParametersGenome.cpp:5. Args: pPin: Parameters"]
pub fn parametersgenome_l5_parametersgenome_initialize(
    pg: &mut crate::parameters_genome::ParametersGenome,
) -> Result<(), String> {
    if !pg.g_dir.ends_with('/') {
        pg.g_dir.push('/');
    }

    pg.transform.type_ = match pg.transform.type_string.as_str() {
        "None" => 0,
        "Haploid" => 1,
        "Diploid" => 2,
        other => {
            return Err(format!(
                "EXITING because of FATAL PARAMETER ERROR: unrecognized option in --outTransformType = {}\nSOLUTION: use one of the allowed values for --outTransformType : 'None' or 'Haploid' or 'Diploid' \n",
                other
            ));
        }
    };

    pg.transform.out_yes = false;
    pg.transform.out_sam = false;
    pg.transform.out_sj = false;
    pg.transform.out_quant = false;
    if pg.transform.output.first().map(|s| s.as_str()) != Some("None") {
        for ot in pg.transform.output.iter() {
            match ot.as_str() {
                "SAM" => {
                    pg.transform.out_yes = true;
                    pg.transform.out_sam = true;
                }
                "SJ" => {
                    pg.transform.out_yes = true;
                    pg.transform.out_sj = true;
                }
                "Quant" => {
                    pg.transform.out_yes = true;
                    pg.transform.out_quant = true;
                }
                other => {
                    return Err(format!(
                        "EXITING because of FATAL PARAMETER ERROR: unrecognized option in --outTransformOutput = {}\nSOLUTION: use allowed values for --outTransformOutput: None or SAM and/or SJ\n",
                        other
                    ));
                }
            }
        }
    }

    if pg.g_type_string != "Full"
        && pg.g_type_string != "Transcriptome"
        && pg.g_type_string != "SuperTranscriptome"
    {
        return Err(format!(
            "EXITING because of FATAL parameter error: --genomeType={}\nSOLUTION: use one of the allowed values for --genomeLoad : Full OR Transcriptome OR SuperTranscriptome\n",
            pg.g_type_string
        ));
    }

    if pg.g_load != "LoadAndKeep"
        && pg.g_load != "LoadAndRemove"
        && pg.g_load != "Remove"
        && pg.g_load != "LoadAndExit"
        && pg.g_load != "NoSharedMemory"
    {
        return Err(format!(
            "EXITING because of FATAL INPUT ERROR: --genomeLoad={}\nSOLUTION: use one of the allowed values for --genomeLoad : NoSharedMemory,LoadAndKeep,LoadAndRemove,LoadAndExit,Remove.\n",
            pg.g_load
        ));
    }
    Ok(())
}
