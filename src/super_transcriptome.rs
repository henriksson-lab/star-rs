#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SuperTranscript` at STAR/source/SuperTranscriptome.h:14."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SuperTranscript {
    pub seq: Vec<u8>,
    pub length: u32,
    pub sj_c: Vec<[u32; 3]>,
    pub sj_donor: Vec<u32>,
}

#[doc = "Original class `SuperTranscriptome` at STAR/source/SuperTranscriptome.h:23."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SuperTranscriptome {
    pub sj_nmax: u32,
    pub sj_donor_nmax: u32,
    pub n: u32,
    pub seq_concat: Vec<u8>,
    pub seq: Vec<Vec<u8>>,
    pub tr_index: Vec<u64>,
    pub tr_start_end: Vec<[u64; 2]>,
    pub sj: Vec<sjInfo>,
    pub super_trs: Vec<SuperTranscript>,
}

#[doc = "Original struct `sjInfo` at STAR/source/SuperTranscriptome.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct sjInfo {
    pub start: u32,
    pub end: u32,
    pub tr: u32,
    pub super_: u32,
}

#[doc = "Original `SuperTranscriptome::sjCollapse` at STAR/source/SuperTranscriptome.cpp:4. Args: "]
pub fn supertranscriptome_l4_supertranscriptome_sjcollapse(
    super_transcriptome: &mut crate::super_transcriptome::SuperTranscriptome,
) -> (String, String) {
    super_transcriptome.sj.sort_by(|sj1, sj2| {
        (sj1.super_, sj1.start, sj1.end).cmp(&(sj2.super_, sj2.start, sj2.end))
    });

    let mut sj_collapsed = vec![Vec::<[u32; 2]>::new(); super_transcriptome.seq.len()];
    for i in 0..super_transcriptome.sj.len() {
        if i == 0
            || super_transcriptome.sj[i].start != super_transcriptome.sj[i - 1].start
            || super_transcriptome.sj[i].end != super_transcriptome.sj[i - 1].end
            || super_transcriptome.sj[i].super_ != super_transcriptome.sj[i - 1].super_
        {
            sj_collapsed[super_transcriptome.sj[i].super_ as usize].push([
                super_transcriptome.sj[i].start,
                super_transcriptome.sj[i].end,
            ]);
        }
    }

    let mut super_tr_sj_stream = String::new();
    for (i, sjs) in sj_collapsed.iter().enumerate() {
        for sj1 in sjs {
            super_tr_sj_stream.push_str(&format!("{}\t{}\t{}\n", i, sj1[0], sj1[1]));
        }
    }

    let log_main = format!(
        "Number of splice junctions in superTranscripts = {}\nNumber of collapsed splice junctions in superTranscripts = {}\n",
        super_transcriptome.sj.len(),
        sj_collapsed.len()
    );

    (super_tr_sj_stream, log_main)
}

#[doc = "Original `SuperTranscriptome::load` at STAR/source/SuperTranscriptome.cpp:32. Args: G: char, chrStart: vector<uint64>, chrLength: vector<uint64>"]
pub fn supertranscriptome_l32_supertranscriptome_load(
    super_transcriptome: &mut crate::super_transcriptome::SuperTranscriptome,
    g: &[u8],
    chr_start: &[u64],
    chr_length: &[u64],
    super_tr_sj_collapsed_tsv: &str,
) -> Result<String, String> {
    if chr_start.len() < chr_length.len() {
        return Err(format!(
            "SuperTranscriptome::load requires chrStart for {} chromosomes, got {}",
            chr_length.len(),
            chr_start.len()
        ));
    }
    super_transcriptome.n = chr_length.len() as u32;
    super_transcriptome
        .super_trs
        .resize(super_transcriptome.n as usize, Default::default());
    for ii in 0..super_transcriptome.n as usize {
        super_transcriptome.super_trs[ii].length = chr_length[ii] as u32;
        let start = chr_start[ii] as usize;
        let end = start
            .checked_add(chr_length[ii] as usize)
            .ok_or_else(|| "SuperTranscriptome::load chromosome interval overflow".to_string())?;
        if end > g.len() {
            return Err(format!(
                "SuperTranscriptome::load chromosome interval {}..{} exceeds genome length {}",
                start,
                end,
                g.len()
            ));
        }
        super_transcriptome.super_trs[ii].seq = g[start..end].to_vec();
    }

    let mut records = Vec::new();
    for chunk in super_tr_sj_collapsed_tsv
        .split_whitespace()
        .collect::<Vec<_>>()
        .chunks(3)
    {
        if chunk.len() != 3 {
            continue;
        }
        records.push([
            chunk[0].parse::<u32>().map_err(|_| {
                format!(
                    "malformed superTranscriptome junction record: invalid super transcript id '{}'",
                    chunk[0]
                )
            })?,
            chunk[1].parse::<u32>().map_err(|_| {
                format!(
                    "malformed superTranscriptome junction record: invalid donor '{}'",
                    chunk[1]
                )
            })?,
            chunk[2].parse::<u32>().map_err(|_| {
                format!(
                    "malformed superTranscriptome junction record: invalid acceptor '{}'",
                    chunk[2]
                )
            })?,
        ]);
    }

    let mut sutr1 = 0u32;
    let mut sj_c1: Vec<[u32; 3]> = Vec::new();
    super_transcriptome.sj_nmax = 0;
    super_transcriptome.sj_donor_nmax = 0;
    let mut sj_donor1: Vec<u32> = Vec::new();

    for rec in records.iter().chain(std::iter::once(&[u32::MAX, 0, 0])) {
        let sutr = rec[0];
        let in_good = sutr != u32::MAX;
        if sutr != sutr1 || !in_good {
            sj_c1.sort_by(|sj1, sj2| (sj1[1], sj1[0]).cmp(&(sj2[1], sj2[0])));
            if (sutr1 as usize) < super_transcriptome.super_trs.len() {
                super_transcriptome.super_trs[sutr1 as usize].sj_c = sj_c1.clone();
                super_transcriptome.super_trs[sutr1 as usize].sj_donor = sj_donor1.clone();
            }
            if super_transcriptome.sj_nmax < sj_c1.len() as u32 {
                super_transcriptome.sj_nmax = sj_c1.len() as u32;
            }
            if super_transcriptome.sj_donor_nmax < sj_donor1.len() as u32 {
                super_transcriptome.sj_donor_nmax = sj_donor1.len() as u32;
            }
            sj_c1.clear();
            sj_donor1.clear();
            sutr1 = sutr;
            if !in_good {
                break;
            }
        }

        let sjd = rec[1];
        let sja = rec[2];
        if sj_donor1.is_empty() || *sj_donor1.last().unwrap() < sjd {
            sj_donor1.push(sjd);
        }
        sj_c1.push([sjd, sja, sj_donor1.len() as u32 - 1]);
    }

    Ok(format!(
        "Max number of splice junctions in a superTranscript = {}\nMax number of donor sites in a superTranscript = {}\n",
        super_transcriptome.sj_nmax, super_transcriptome.sj_donor_nmax
    ))
}
