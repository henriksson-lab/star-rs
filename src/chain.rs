#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `Chain` at STAR/source/Chain.h:16."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Chain {
    pub chain_file_name: String,
    pub chr_chains: std::collections::BTreeMap<String, OneChain>,
}

#[doc = "Original class `OneChain` at STAR/source/Chain.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OneChain {
    pub b_n: u32,
    pub chr1: String,
    pub chr2: String,
    pub b_start1: Vec<u32>,
    pub b_start2: Vec<u32>,
    pub b_len: Vec<u32>,
}

#[doc = "Original `Chain::Chain` at STAR/source/Chain.cpp:5. Args: Pin: Parameters, chainFileNameIn: string"]
pub fn chain_l5_chain_chain(
    chain_file_name: &str,
) -> Result<crate::chain::Chain, String> {
    let mut chain = crate::chain::Chain {
        chain_file_name: chain_file_name.to_string(),
        chr_chains: std::collections::BTreeMap::new(),
    };
    chain_l10_chain_chainload(&mut chain)?;
    Ok(chain)
}

#[doc = "Original `Chain::chainLoad` at STAR/source/Chain.cpp:10. Args: "]
pub fn chain_l10_chain_chainload(
    chain: &mut crate::chain::Chain,
) -> Result<(), String> {
    let contents = std::fs::read_to_string(&chain.chain_file_name).map_err(|err| {
        format!(
            "SOLUTION: check path and permission for the chain file{}: {}",
            chain.chain_file_name, err
        )
    })?;

    let mut chr1 = String::new();
    for (line_no0, line1) in contents.lines().enumerate() {
        let line_no = line_no0 + 1;
        let fields: Vec<&str> = line1.split_whitespace().collect();
        if fields.is_empty() {
        } else if fields.len() == 1 {
            if chr1.is_empty() {
                return Err(format!(
                    "EXITING because of FATAL ERROR in chain file {} line {}: block line appears before chain header\n",
                    chain.chain_file_name, line_no
                ));
            }
            let ch = chain.chr_chains.entry(chr1.clone()).or_default();
            ch.b_len.push(fields[0].parse::<u32>().map_err(|_| {
                format!(
                    "EXITING because of FATAL ERROR in chain file {} line {}: invalid block length {}\n",
                    chain.chain_file_name, line_no, fields[0]
                )
            })?);
            ch.b_n = ch.b_len.len() as u32;
        } else if fields.len() == 3 {
            if chr1.is_empty() {
                return Err(format!(
                    "EXITING because of FATAL ERROR in chain file {} line {}: block line appears before chain header\n",
                    chain.chain_file_name, line_no
                ));
            }
            let ch = chain.chr_chains.entry(chr1.clone()).or_default();
            let block_len = fields[0].parse::<u32>().map_err(|_| {
                format!(
                    "EXITING because of FATAL ERROR in chain file {} line {}: invalid block length {}\n",
                    chain.chain_file_name, line_no, fields[0]
                )
            })?;
            let gap1 = fields[1].parse::<u32>().map_err(|_| {
                format!(
                    "EXITING because of FATAL ERROR in chain file {} line {}: invalid source gap {}\n",
                    chain.chain_file_name, line_no, fields[1]
                )
            })?;
            let gap2 = fields[2].parse::<u32>().map_err(|_| {
                format!(
                    "EXITING because of FATAL ERROR in chain file {} line {}: invalid target gap {}\n",
                    chain.chain_file_name, line_no, fields[2]
                )
            })?;
            ch.b_len.push(block_len);
            let s1 = *ch.b_start1.last().unwrap() + *ch.b_len.last().unwrap() + gap1;
            ch.b_start1.push(s1);
            let s2 = *ch.b_start2.last().unwrap() + *ch.b_len.last().unwrap() + gap2;
            ch.b_start2.push(s2);
        } else {
            if fields.len() < 11 {
                return Err(format!(
                    "EXITING because of FATAL ERROR in chain file {} line {}: malformed chain header\n",
                    chain.chain_file_name, line_no
                ));
            }
            chr1 = fields[2].to_string();
            let ch = chain.chr_chains.entry(chr1.clone()).or_default();
            ch.chr1 = chr1.clone();
            ch.chr2 = fields[7].to_string();
            ch.b_start1.push(fields[5].parse::<u32>().map_err(|_| {
                format!(
                    "EXITING because of FATAL ERROR in chain file {} line {}: invalid source start {}\n",
                    chain.chain_file_name, line_no, fields[5]
                )
            })?);
            ch.b_start2.push(fields[10].parse::<u32>().map_err(|_| {
                format!(
                    "EXITING because of FATAL ERROR in chain file {} line {}: invalid target start {}\n",
                    chain.chain_file_name, line_no, fields[10]
                )
            })?);
        }
    }
    Ok(())
}

#[doc = "Original `Chain::liftOverGTF` at STAR/source/Chain.cpp:58. Args: gtfFileName: string, outFileName: string"]
pub fn chain_l58_chain_liftovergtf(
    chain: &crate::chain::Chain,
    gtf_file_name: &str,
    out_file_name: &str,
) -> Result<(), String> {
    let contents = crate::io_utils::read_to_string_auto_gzip(gtf_file_name).map_err(|err| {
        format!(
            "SOLUTION: check path and permission for the GTF file{}: {}",
            gtf_file_name, err
        )
    })?;
    let mut stream_out = String::new();
    let mut stream_out_unlifted = String::new();

    for (line_no0, line1) in contents.lines().enumerate() {
        let line_no = line_no0 + 1;
        let mut ranges: Vec<(usize, usize)> = Vec::new();
        let mut start = None;
        for (idx, ch) in line1.char_indices() {
            if ch.is_whitespace() {
                if let Some(st) = start.take() {
                    ranges.push((st, idx));
                }
            } else if start.is_none() {
                start = Some(idx);
            }
        }
        if let Some(st) = start {
            ranges.push((st, line1.len()));
        }

        if ranges.is_empty() {
            continue;
        }
        let chr1 = &line1[ranges[0].0..ranges[0].1];
        if chr1.starts_with('#') {
            continue;
        }
        if ranges.len() < 5 {
            return Err(format!(
                "GTF line {} in {} has too few fields for liftOver",
                line_no, gtf_file_name
            ));
        }

        let ch1 = chain.chr_chains.get(chr1).ok_or_else(|| {
            format!(
                "GTF contains chromosome {} not present in the chain file {}",
                chr1, chain.chain_file_name
            )
        })?;

        let str1 = &line1[ranges[1].0..ranges[1].1];
        let str2 = &line1[ranges[2].0..ranges[2].1];
        let mut c2 = [u32::MAX; 2];

        for ii in 0..2 {
            let c1: u32 = line1[ranges[3 + ii].0..ranges[3 + ii].1]
                .parse()
                .map_err(|_| {
                    format!(
                        "GTF line {} in {} has invalid coordinate {}",
                        line_no,
                        gtf_file_name,
                        &line1[ranges[3 + ii].0..ranges[3 + ii].1]
                    )
                })?;
            let i1 = servicefuns_l239_binarysearch1a(c1, &ch1.b_start1, ch1.b_n as i32);

            if i1 >= 0 && c1 < ch1.b_start1[i1 as usize] + ch1.b_len[i1 as usize] {
                c2[ii] = ch1.b_start2[i1 as usize] + c1 - ch1.b_start1[i1 as usize];
            } else if ii == 0 && i1 < ch1.b_n as i32 - 1 {
                c2[ii] = ch1.b_start2[(i1 + 1) as usize];
            } else if ii == 1 && i1 >= 0 {
                c2[ii] = ch1.b_start2[i1 as usize] + ch1.b_len[i1 as usize] - 1;
            }
        }

        if c2[0] != u32::MAX && c2[1] != u32::MAX && c2[1] >= c2[0] {
            let rest = &line1[ranges[4].1..];
            stream_out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}{}\n",
                ch1.chr2, str1, str2, c2[0], c2[1], rest
            ));
        } else {
            stream_out_unlifted.push_str(line1);
            stream_out_unlifted.push('\n');
        }
    }

    std::fs::write(out_file_name, stream_out).map_err(|err| err.to_string())?;
    std::fs::write(format!("{}.unlifted", out_file_name), stream_out_unlifted)
        .map_err(|err| err.to_string())?;
    Ok(())
}
