#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `Transcriptome` at STAR/source/Transcriptome.h:13."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Transcriptome {
    pub n_tr: u32,
    pub n_ge: u32,
    pub ge_id: Vec<String>,
    pub ge_name: Vec<String>,
    pub ge_biotype: Vec<String>,
    pub tr_id: Vec<String>,
    pub tr_s: Vec<u32>,
    pub tr_e: Vec<u32>,
    pub tr_e_max: Vec<u32>,
    pub tr_ex_n: Vec<u16>,
    pub tr_ex_i: Vec<u32>,
    pub tr_str: Vec<u8>,
    pub tr_gene: Vec<u32>,
    pub tr_len: Vec<u32>,
    pub ex_se: Vec<u32>,
    pub ex_len_cum: Vec<u32>,
    pub ex_g: TranscriptomeExG,
    pub gene_full: TranscriptomeGeneFull,
    pub quants: Quantifications,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TranscriptomeExG {
    pub n_ex: u64,
    pub s: Vec<u64>,
    pub e: Vec<u64>,
    pub e_max: Vec<u64>,
    pub str_: Vec<u8>,
    pub g: Vec<u32>,
    pub t: Vec<u32>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TranscriptomeGeneFull {
    pub s: Vec<u64>,
    pub e: Vec<u64>,
    pub e_max: Vec<u64>,
    pub str_: Vec<u8>,
    pub g: Vec<u32>,
}

#[doc = "Original `Transcriptome::Transcriptome` at STAR/source/Transcriptome.cpp:7. Args: Pin: Parameters"]
pub fn transcriptome_l7_transcriptome_transcriptome(
    quant_yes: bool,
    transform_out_quant: bool,
    sjdb_gtf_file: &str,
    genome_dir: &str,
    sjdb_insert_out_dir: &str,
    transformed_genome_dir: &str,
    load_transcript_exon_info: bool,
    load_exon_gene_info: bool,
    load_gene_full_info: bool,
    gene_info_tab: &str,
    transcript_info_tab: Option<&str>,
    exon_info_tab: Option<&str>,
    exon_ge_tr_info_tab: Option<&str>,
) -> Result<(crate::transcriptome::Transcriptome, String, String), String> {
    let mut transcriptome = crate::transcriptome::Transcriptome::default();
    let mut log_main = String::new();
    if !quant_yes {
        return Ok((transcriptome, String::new(), log_main));
    }

    let tr_info_dir = if !transform_out_quant {
        if sjdb_gtf_file == "-" {
            genome_dir.to_string()
        } else {
            sjdb_insert_out_dir.to_string()
        }
    } else {
        transformed_genome_dir.to_string()
    };

    let mut gene_tokens = gene_info_tab.split_whitespace();
    transcriptome.n_ge = gene_tokens
        .next()
        .ok_or_else(|| "geneInfo.tab is empty".to_string())?
        .parse::<u32>()
        .map_err(|_| "geneInfo.tab has invalid gene count".to_string())?;
    transcriptome
        .ge_id
        .resize(transcriptome.n_ge as usize, String::new());
    transcriptome
        .ge_name
        .resize(transcriptome.n_ge as usize, String::new());
    transcriptome
        .ge_biotype
        .resize(transcriptome.n_ge as usize, String::new());
    for ii in 0..transcriptome.n_ge as usize {
        transcriptome.ge_id[ii] = gene_tokens
            .next()
            .ok_or_else(|| "geneInfo.tab ended before geID".to_string())?
            .to_string();
        transcriptome.ge_name[ii] = gene_tokens
            .next()
            .ok_or_else(|| "geneInfo.tab ended before geName".to_string())?
            .to_string();
        transcriptome.ge_biotype[ii] = gene_tokens
            .next()
            .ok_or_else(|| "geneInfo.tab ended before geBiotype".to_string())?
            .to_string();
    }

    if load_transcript_exon_info {
        let transcript_info_tab =
            transcript_info_tab.ok_or_else(|| "transcriptInfo.tab is required".to_string())?;
        let mut tr_tokens = transcript_info_tab.split_whitespace();
        transcriptome.n_tr = tr_tokens
            .next()
            .ok_or_else(|| "transcriptInfo.tab is empty".to_string())?
            .parse::<u32>()
            .map_err(|_| "transcriptInfo.tab has invalid transcript count".to_string())?;
        let n_tr = transcriptome.n_tr as usize;
        transcriptome.tr_id.resize(n_tr, String::new());
        transcriptome.tr_s.resize(n_tr, 0);
        transcriptome.tr_e.resize(n_tr, 0);
        transcriptome.tr_e_max.resize(n_tr, 0);
        transcriptome.tr_ex_i.resize(n_tr, 0);
        transcriptome.tr_ex_n.resize(n_tr, 0);
        transcriptome.tr_str.resize(n_tr, 0);
        transcriptome.tr_gene.resize(n_tr, 0);
        transcriptome.tr_len.resize(n_tr, 0);
        for itr in 0..n_tr {
            transcriptome.tr_id[itr] = tr_tokens
                .next()
                .ok_or_else(|| "transcriptInfo.tab ended before trID".to_string())?
                .to_string();
            transcriptome.tr_s[itr] = tr_tokens
                .next()
                .ok_or_else(|| "transcriptInfo.tab ended before trS".to_string())?
                .parse::<u32>()
                .map_err(|_| "transcriptInfo.tab has invalid trS".to_string())?;
            transcriptome.tr_e[itr] = tr_tokens
                .next()
                .ok_or_else(|| "transcriptInfo.tab ended before trE".to_string())?
                .parse::<u32>()
                .map_err(|_| "transcriptInfo.tab has invalid trE".to_string())?;
            transcriptome.tr_e_max[itr] = tr_tokens
                .next()
                .ok_or_else(|| "transcriptInfo.tab ended before trEmax".to_string())?
                .parse::<u32>()
                .map_err(|_| "transcriptInfo.tab has invalid trEmax".to_string())?;
            transcriptome.tr_str[itr] = tr_tokens
                .next()
                .ok_or_else(|| "transcriptInfo.tab ended before trStr".to_string())?
                .parse::<u16>()
                .map_err(|_| "transcriptInfo.tab has invalid trStr".to_string())?
                as u8;
            transcriptome.tr_ex_n[itr] = tr_tokens
                .next()
                .ok_or_else(|| "transcriptInfo.tab ended before trExN".to_string())?
                .parse::<u16>()
                .map_err(|_| "transcriptInfo.tab has invalid trExN".to_string())?;
            transcriptome.tr_ex_i[itr] = tr_tokens
                .next()
                .ok_or_else(|| "transcriptInfo.tab ended before trExI".to_string())?
                .parse::<u32>()
                .map_err(|_| "transcriptInfo.tab has invalid trExI".to_string())?;
            transcriptome.tr_gene[itr] = tr_tokens
                .next()
                .ok_or_else(|| "transcriptInfo.tab ended before trGene".to_string())?
                .parse::<u32>()
                .map_err(|_| "transcriptInfo.tab has invalid trGene".to_string())?;
        }
        log_main.push_str(&format!(
            "Loaded transcript database, nTr={}\n",
            transcriptome.n_tr
        ));

        let exon_info_tab = exon_info_tab.ok_or_else(|| "exonInfo.tab is required".to_string())?;
        let mut ex_tokens = exon_info_tab.split_whitespace();
        let n_ex = ex_tokens
            .next()
            .ok_or_else(|| "exonInfo.tab is empty".to_string())?
            .parse::<u32>()
            .map_err(|_| "exonInfo.tab has invalid exon count".to_string())?;
        transcriptome.ex_se.resize(2 * n_ex as usize, 0);
        transcriptome.ex_len_cum.resize(n_ex as usize, 0);
        for iex in 0..n_ex as usize {
            transcriptome.ex_se[2 * iex] = ex_tokens
                .next()
                .ok_or_else(|| "exonInfo.tab ended before exon start".to_string())?
                .parse::<u32>()
                .map_err(|_| "exonInfo.tab has invalid exon start".to_string())?;
            transcriptome.ex_se[2 * iex + 1] = ex_tokens
                .next()
                .ok_or_else(|| "exonInfo.tab ended before exon end".to_string())?
                .parse::<u32>()
                .map_err(|_| "exonInfo.tab has invalid exon end".to_string())?;
            transcriptome.ex_len_cum[iex] = ex_tokens
                .next()
                .ok_or_else(|| "exonInfo.tab ended before cumulative exon length".to_string())?
                .parse::<u32>()
                .map_err(|_| "exonInfo.tab has invalid cumulative exon length".to_string())?;
        }
        for ii in 0..n_tr {
            let iex1 = transcriptome.tr_ex_i[ii] as usize + transcriptome.tr_ex_n[ii] as usize - 1;
            transcriptome.tr_len[ii] = transcriptome.ex_len_cum[iex1]
                + transcriptome.ex_se[2 * iex1 + 1]
                - transcriptome.ex_se[2 * iex1]
                + 1;
        }
        log_main.push_str(&format!("Loaded exon database, nEx={}\n", n_ex));
    }

    if load_exon_gene_info {
        let exon_ge_tr_info_tab =
            exon_ge_tr_info_tab.ok_or_else(|| "exonGeTrInfo.tab is required".to_string())?;
        let mut exg_tokens = exon_ge_tr_info_tab.split_whitespace();
        transcriptome.ex_g.n_ex = exg_tokens
            .next()
            .ok_or_else(|| "exonGeTrInfo.tab is empty".to_string())?
            .parse::<u64>()
            .map_err(|_| "exonGeTrInfo.tab has invalid exon count".to_string())?;
        let n_ex = transcriptome.ex_g.n_ex as usize;
        transcriptome.ex_g.s.resize(n_ex, 0);
        transcriptome.ex_g.e.resize(n_ex, 0);
        transcriptome.ex_g.e_max.resize(n_ex, 0);
        transcriptome.ex_g.str_.resize(n_ex, 0);
        transcriptome.ex_g.g.resize(n_ex, 0);
        transcriptome.ex_g.t.resize(n_ex, 0);
        for ii in 0..n_ex {
            transcriptome.ex_g.s[ii] = exg_tokens
                .next()
                .ok_or_else(|| "exonGeTrInfo.tab ended before exon start".to_string())?
                .parse::<u64>()
                .map_err(|_| "exonGeTrInfo.tab has invalid exon start".to_string())?;
            transcriptome.ex_g.e[ii] = exg_tokens
                .next()
                .ok_or_else(|| "exonGeTrInfo.tab ended before exon end".to_string())?
                .parse::<u64>()
                .map_err(|_| "exonGeTrInfo.tab has invalid exon end".to_string())?;
            transcriptome.ex_g.str_[ii] = exg_tokens
                .next()
                .ok_or_else(|| "exonGeTrInfo.tab ended before strand".to_string())?
                .parse::<i32>()
                .map_err(|_| "exonGeTrInfo.tab has invalid strand".to_string())?
                as u8;
            transcriptome.ex_g.g[ii] = exg_tokens
                .next()
                .ok_or_else(|| "exonGeTrInfo.tab ended before gene".to_string())?
                .parse::<u32>()
                .map_err(|_| "exonGeTrInfo.tab has invalid gene".to_string())?;
            transcriptome.ex_g.t[ii] = exg_tokens
                .next()
                .ok_or_else(|| "exonGeTrInfo.tab ended before transcript".to_string())?
                .parse::<u32>()
                .map_err(|_| "exonGeTrInfo.tab has invalid transcript".to_string())?;
        }
        if n_ex > 0 {
            transcriptome.ex_g.e_max[0] = transcriptome.ex_g.e[0];
            for iex in 1..n_ex {
                transcriptome.ex_g.e_max[iex] =
                    transcriptome.ex_g.e_max[iex - 1].max(transcriptome.ex_g.e[iex]);
            }
        }
    }

    if load_gene_full_info {
        let exon_ge_tr_info_tab =
            exon_ge_tr_info_tab.ok_or_else(|| "exonGeTrInfo.tab is required".to_string())?;
        let mut exg_tokens = exon_ge_tr_info_tab.split_whitespace();
        let n_ex = exg_tokens
            .next()
            .ok_or_else(|| "exonGeTrInfo.tab is empty".to_string())?
            .parse::<u64>()
            .map_err(|_| "exonGeTrInfo.tab has invalid exon count".to_string())?
            as usize;
        let n_ge = transcriptome.n_ge as usize;
        let mut gene_full_rows = vec![[u64::MAX, 0, 0, 0]; n_ge];
        for (ig, row) in gene_full_rows.iter_mut().enumerate() {
            row[3] = ig as u64;
        }
        for _ in 0..n_ex {
            let s1 = exg_tokens
                .next()
                .ok_or_else(|| "exonGeTrInfo.tab ended before exon start".to_string())?
                .parse::<u64>()
                .map_err(|_| "exonGeTrInfo.tab has invalid exon start".to_string())?;
            let e1 = exg_tokens
                .next()
                .ok_or_else(|| "exonGeTrInfo.tab ended before exon end".to_string())?
                .parse::<u64>()
                .map_err(|_| "exonGeTrInfo.tab has invalid exon end".to_string())?;
            let str1 = exg_tokens
                .next()
                .ok_or_else(|| "exonGeTrInfo.tab ended before strand".to_string())?
                .parse::<u64>()
                .map_err(|_| "exonGeTrInfo.tab has invalid strand".to_string())?;
            let g1 = exg_tokens
                .next()
                .ok_or_else(|| "exonGeTrInfo.tab ended before gene".to_string())?
                .parse::<u64>()
                .map_err(|_| "exonGeTrInfo.tab has invalid gene".to_string())?
                as usize;
            let _t1 = exg_tokens
                .next()
                .ok_or_else(|| "exonGeTrInfo.tab ended before transcript".to_string())?;
            gene_full_rows[g1][0] = gene_full_rows[g1][0].min(s1);
            gene_full_rows[g1][1] = gene_full_rows[g1][1].max(e1);
            gene_full_rows[g1][2] = str1;
        }
        gene_full_rows.sort_by(|a, b| a[0..2].cmp(&b[0..2]));
        transcriptome.gene_full.s = gene_full_rows.iter().map(|row| row[0]).collect();
        transcriptome.gene_full.e = gene_full_rows.iter().map(|row| row[1]).collect();
        transcriptome.gene_full.str_ = gene_full_rows.iter().map(|row| row[2] as u8).collect();
        transcriptome.gene_full.g = gene_full_rows.iter().map(|row| row[3] as u32).collect();
        transcriptome.gene_full.e_max.resize(n_ge, 0);
        if n_ge > 0 {
            transcriptome.gene_full.e_max[0] = transcriptome.gene_full.e[0];
            for ig in 1..n_ge {
                transcriptome.gene_full.e_max[ig] =
                    transcriptome.gene_full.e_max[ig - 1].max(transcriptome.gene_full.e[ig]);
            }
        }
    }

    Ok((transcriptome, tr_info_dir, log_main))
}

#[doc = "Original `Transcriptome::quantsAllocate` at STAR/source/Transcriptome.cpp:150. Args: "]
pub fn transcriptome_l150_transcriptome_quantsallocate(
    transcriptome: &mut crate::transcriptome::Transcriptome,
    ge_count_yes: bool,
) {
    if ge_count_yes {
        transcriptome.quants =
            quantifications_l3_quantifications_quantifications(transcriptome.n_ge);
    }
}

#[doc = "Original `Transcriptome::quantsOutput` at STAR/source/Transcriptome.cpp:156. Args: "]
pub fn transcriptome_l156_transcriptome_quantsoutput(
    transcriptome: &crate::transcriptome::Transcriptome,
    out_file: &str,
    stats_all: &crate::stats::Stats,
) -> Result<(), String> {
    use std::io::Write;

    let mut q_out = std::fs::OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out_file)
        .map_err(|e| e.to_string())?;

    let gc = &transcriptome.quants.gene_counts;
    let unmapped = stats_all.unmapped_mismatch
        + stats_all.unmapped_short
        + stats_all.unmapped_other
        + stats_all.unmapped_multi;

    write!(q_out, "N_unmapped").map_err(|e| e.to_string())?;
    for _ in 0..gc.n_type {
        write!(q_out, "\t{}", unmapped).map_err(|e| e.to_string())?;
    }
    writeln!(q_out).map_err(|e| e.to_string())?;

    write!(q_out, "N_multimapping").map_err(|e| e.to_string())?;
    for _ in 0..gc.n_type {
        write!(q_out, "\t{}", gc.c_multi).map_err(|e| e.to_string())?;
    }
    writeln!(q_out).map_err(|e| e.to_string())?;

    write!(q_out, "N_noFeature").map_err(|e| e.to_string())?;
    for itype in 0..gc.n_type as usize {
        write!(q_out, "\t{}", gc.c_none[itype]).map_err(|e| e.to_string())?;
    }
    writeln!(q_out).map_err(|e| e.to_string())?;

    write!(q_out, "N_ambiguous").map_err(|e| e.to_string())?;
    for itype in 0..gc.n_type as usize {
        write!(q_out, "\t{}", gc.c_ambig[itype]).map_err(|e| e.to_string())?;
    }
    writeln!(q_out).map_err(|e| e.to_string())?;

    for ig in 0..transcriptome.n_ge as usize {
        write!(q_out, "{}", transcriptome.ge_id[ig]).map_err(|e| e.to_string())?;
        for itype in 0..gc.n_type as usize {
            write!(q_out, "\t{}", gc.g_count[itype][ig]).map_err(|e| e.to_string())?;
        }
        writeln!(q_out).map_err(|e| e.to_string())?;
    }

    Ok(())
}
