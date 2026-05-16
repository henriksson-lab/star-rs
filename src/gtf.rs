#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `GTF` at STAR/source/GTF.h:10."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GTF {
    pub gtf_yes: bool,
    pub exon_n: u64,
    pub exon_loci: Vec<[u64; crate::include_define::GTF_EX_L]>,
    pub transcript_strand: Vec<u32>,
    pub transcript_id: Vec<String>,
    pub gene_id: Vec<String>,
    pub gene_attr: Vec<[String; 2]>,
    pub transcript_seq: Vec<Vec<u8>>,
    pub transcript_start_end: Vec<[u64; 2]>,
    pub super_trome: SuperTranscriptome,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GtfTranscriptGeneSjOutput {
    pub n_junctions_added: u64,
    pub exon_ge_tr_info_tab: String,
    pub gene_info_tab: String,
    pub transcript_info_tab: String,
    pub exon_info_tab: String,
    pub sjdb_list_from_gtf_out_tab: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GtfSuperTranscriptOutput {
    pub transcript_sequences_fasta: String,
    pub super_transcript_sequences_fasta: String,
    pub super_transcript_sj_tsv: String,
    pub conversion_to_full_genome_tsv: String,
    pub full_genome_chr_name_txt: String,
    pub full_genome_chr_start_txt: String,
    pub full_genome_chr_length_txt: String,
    pub full_genome_chr_name_length_txt: String,
    pub full_genome_sequence: Vec<u8>,
    pub log_main: String,
}

#[doc = "Original `GTF::GTF` at STAR/source/GTF.cpp:7. Args: genome: Genome, P: Parameters, dirOut: string, sjdbLoci: SjdbClass"]
pub fn gtf_l7_gtf_gtf(
    genome: &mut crate::genome::Genome,
    _p: &crate::parameters_chimeric::Parameters,
    _dir_out: &str,
    gtf_contents: Option<&str>,
) -> Result<(crate::gtf::GTF, String), String> {
    let mut gtf = crate::gtf::GTF::default();
    let mut log_main = String::new();

    if genome.sjdb_overhang == 0 || genome.p_ge.sjdb_gtf_file == "-" {
        gtf.gtf_yes = false;
        return Ok((gtf, log_main));
    }
    gtf.gtf_yes = true;
    log_main.push_str(" ..... processing annotations GTF\n");

    let contents = gtf_contents.ok_or_else(|| {
        format!(
            "FATAL error, could not open file pGe.sjdbGTFfile={}\n",
            genome.p_ge.sjdb_gtf_file
        )
    })?;

    if genome.chr_name_index.is_empty() {
        for ii in 0..genome.n_chr_real as usize {
            genome
                .chr_name_index
                .insert(genome.chr_name[ii].clone(), ii as u64);
        }
    }

    let feature_exon = &genome.p_ge.sjdb_gtf_feature_exon;
    let exon_estimate = contents
        .lines()
        .filter(|line| {
            let mut fields = line.split_whitespace();
            let chr1 = fields.next().unwrap_or("");
            let _source = fields.next();
            let feature = fields.next().unwrap_or("");
            !chr1.starts_with('#') && feature == feature_exon
        })
        .count();

    if exon_estimate == 0 {
        return Err(format!(
            "Fatal INPUT FILE error, no exon lines in the GTF file: {}\nSolution: check the formatting of the GTF file, it must contain some lines with exon in the 3rd column.\n          Make sure the GTF file is unzipped.\n          If exons are marked with a different word, use --sjdbGTFfeatureExon .\n",
            genome.p_ge.sjdb_gtf_file
        ));
    }

    let mut transcript_id_number = std::collections::BTreeMap::<String, u64>::new();
    let mut gene_id_number = std::collections::BTreeMap::<String, u64>::new();

    for one_line in contents.lines() {
        let mut fields = one_line.split_whitespace();
        let Some(mut chr1) = fields.next().map(|s| s.to_string()) else {
            continue;
        };
        let _source = fields.next();
        let feature_type = fields.next().unwrap_or("");
        if chr1.starts_with('#') || feature_type != feature_exon {
            continue;
        }

        if genome.p_ge.sjdb_gtf_chr_prefix != "-" {
            chr1 = format!("{}{}", genome.p_ge.sjdb_gtf_chr_prefix, chr1);
        }

        let Some(&chr_index_u64) = genome.chr_name_index.get(&chr1) else {
            log_main.push_str(&format!(
                "WARNING: while processing sjdbGTFfile={}: chromosome '{}' not found in Genome fasta files for line:\n{}\n",
                genome.p_ge.sjdb_gtf_file, chr1, one_line
            ));
            continue;
        };
        let chr_index = chr_index_u64 as usize;

        let ex1 = fields
            .next()
            .and_then(|s| s.parse::<u64>().ok())
            .unwrap_or(0);
        let ex2 = fields
            .next()
            .and_then(|s| s.parse::<u64>().ok())
            .unwrap_or(0);
        let _score = fields.next();
        let str1 = fields.next().and_then(|s| s.chars().next()).unwrap_or('.');
        let _frame = fields.next();
        if ex2 > genome.chr_length[chr_index] {
            log_main.push_str(&format!(
                "WARNING: while processing sjdbGTFfile={}, line:\n{}\n exon end = {} is larger than the chromosome {} length = {} , will skip this exon\n",
                genome.p_ge.sjdb_gtf_file,
                one_line,
                ex2,
                chr1,
                genome.chr_length[chr_index]
            ));
            continue;
        }

        let attr_start = one_line
            .splitn(9, char::is_whitespace)
            .collect::<Vec<_>>()
            .get(8)
            .copied()
            .unwrap_or("");
        let mut attrs = String::with_capacity(attr_start.len());
        for ch in attr_start.chars() {
            attrs.push(match ch {
                ';' | '=' | '\t' | '"' => ' ',
                _ => ch,
            });
        }
        let attrs = format!(" {} ", attrs);

        let attr_names = [
            vec![genome.p_ge.sjdb_gtf_tag_exon_parent_transcript.clone()],
            vec![genome.p_ge.sjdb_gtf_tag_exon_parent_gene.clone()],
            genome.p_ge.sjdb_gtf_tag_exon_parent_gene_name.clone(),
            genome.p_ge.sjdb_gtf_tag_exon_parent_gene_type.clone(),
        ];
        let mut ex_attr = vec![String::new(); attr_names.len()];
        for ii in 0..attr_names.len() {
            for attr1 in &attr_names[ii] {
                let needle = format!(" {} ", attr1);
                if let Some(pos1) = attrs.find(&needle) {
                    let value_start = pos1 + attr1.len() + 2;
                    if let Some(rel) = attrs[value_start..].find(|c: char| c != ' ') {
                        let pos2 = value_start + rel;
                        let pos3 = attrs[pos2..]
                            .find(' ')
                            .map(|p| pos2 + p)
                            .unwrap_or(attrs.len());
                        ex_attr[ii] = attrs[pos2..pos3].to_string();
                    }
                }
            }
        }

        if ex_attr[0].is_empty() {
            log_main.push_str(&format!(
                "WARNING: while processing pGe.sjdbGTFfile={}: no transcript_id for line:\n{}\n",
                genome.p_ge.sjdb_gtf_file, one_line
            ));
            ex_attr[0] = format!("tr_{}_{}_{}_{}", chr1, ex1, ex2, gtf.exon_n);
        }

        if ex_attr[1].is_empty() {
            log_main.push_str(&format!(
                "WARNING: while processing pGe.sjdbGTFfile={}: no gene_id for line:\n{}\n",
                genome.p_ge.sjdb_gtf_file, one_line
            ));
            ex_attr[1] = "MissingGeneID".to_string();
        }

        if ex_attr[2].is_empty() {
            ex_attr[2] = ex_attr[1].clone();
        }
        if ex_attr[3].is_empty() {
            ex_attr[3] = "MissingGeneType".to_string();
        }

        if !transcript_id_number.contains_key(&ex_attr[0]) {
            let next = transcript_id_number.len() as u64;
            transcript_id_number.insert(ex_attr[0].clone(), next);
            gtf.transcript_id.push(ex_attr[0].clone());
            gtf.transcript_strand.push(match str1 {
                '+' => 1,
                '-' => 2,
                _ => 0,
            });
        }

        if !gene_id_number.contains_key(&ex_attr[1]) {
            let next = gene_id_number.len() as u64;
            gene_id_number.insert(ex_attr[1].clone(), next);
            gtf.gene_id.push(ex_attr[1].clone());
            gtf.gene_attr.push([ex_attr[2].clone(), ex_attr[3].clone()]);
        }

        gtf.exon_loci.push([
            transcript_id_number[&ex_attr[0]],
            ex1 + genome.chr_start[chr_index] - 1,
            ex2 + genome.chr_start[chr_index] - 1,
            gene_id_number[&ex_attr[1]],
        ]);
        gtf.exon_n += 1;
    }

    if gtf.exon_n == 0 {
        return Err(format!(
            "Fatal INPUT FILE error, no valid exon lines in the GTF file: {}\nSolution: check the formatting of the GTF file. One likely cause is the difference in chromosome naming between GTF and FASTA file.\n",
            genome.p_ge.sjdb_gtf_file
        ));
    }

    Ok((gtf, log_main))
}
