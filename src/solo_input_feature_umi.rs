#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `soloInputFeatureUMI` at STAR/source/soloInputFeatureUMI.cpp:5. Args: strIn: fstream, featureType: int32, readInfoYes: bool, sjAll: array<vector<uint64>,2>, iread: uint64, cbmatch: int32, feature: uint32, umi: uint64, featVecU32: vector<uint32>, readFlagCounts: SoloReadFlagClass"]
pub fn soloinputfeatureumi_l5_soloinputfeatureumi<'a>(
    tokens: &mut std::str::SplitWhitespace<'a>,
    feature_type: i32,
    read_info_yes: bool,
    sj_all: &[Vec<u64>; 2],
    iread: &mut u64,
    cbmatch: &mut i32,
    feature: &mut u32,
    umi: &mut u64,
    feat_vec_u32: &mut Vec<u32>,
    read_flag_counts: &mut crate::solo_common::SoloReadFlagClass,
) -> Result<bool, String> {
    let Some(umi_token) = tokens.next() else {
        return Ok(false);
    };
    *umi = umi_token
        .parse()
        .map_err(|_| format!("Malformed STARsolo feature record: invalid UMI {umi_token}"))?;

    if read_info_yes {
        let iread_token = tokens
            .next()
            .ok_or_else(|| "Malformed STARsolo feature record: missing read index".to_string())?;
        *iread = iread_token.parse().map_err(|_| {
            format!("Malformed STARsolo feature record: invalid read index {iread_token}")
        })?;
        let flag_token = tokens
            .next()
            .ok_or_else(|| "Malformed STARsolo feature record: missing read flag".to_string())?;
        read_flag_counts.flag = flag_token.parse().map_err(|_| {
            format!("Malformed STARsolo feature record: invalid read flag {flag_token}")
        })?;
    }

    match feature_type {
        SOLO_FEATURE_GENE
        | SOLO_FEATURE_GENE_FULL
        | SOLO_FEATURE_GENE_FULL_EX50P_AS
        | SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON => {
            let feature_token = tokens
                .next()
                .ok_or_else(|| "Malformed STARsolo feature record: missing feature".to_string())?;
            *feature = feature_token.parse().map_err(|_| {
                format!("Malformed STARsolo feature record: invalid feature {feature_token}")
            })?;
        }
        SOLO_FEATURE_SJ => {
            let sj0_token = tokens.next().ok_or_else(|| {
                "Malformed STARsolo feature record: missing splice junction start".to_string()
            })?;
            let sj0: u32 = sj0_token.parse().map_err(|_| {
                format!(
                    "Malformed STARsolo feature record: invalid splice junction start {sj0_token}"
                )
            })?;
            let sj1_token = tokens.next().ok_or_else(|| {
                "Malformed STARsolo feature record: missing splice junction end".to_string()
            })?;
            let sj1: u32 = sj1_token.parse().map_err(|_| {
                format!(
                    "Malformed STARsolo feature record: invalid splice junction end {sj1_token}"
                )
            })?;
            *feature = binarysearch2_l3_binarysearch2(
                sj0,
                sj1,
                &sj_all[0].iter().map(|v| *v as u32).collect::<Vec<_>>(),
                &sj_all[1].iter().map(|v| *v as u32).collect::<Vec<_>>(),
                sj_all[0].len() as i32,
            ) as u32;
        }
        SOLO_FEATURE_TRANSCRIPT3P => {
            *feature = 0;
            let ntr_token = tokens.next().ok_or_else(|| {
                "Malformed STARsolo feature record: missing transcript count".to_string()
            })?;
            let ntr: u32 = ntr_token.parse().map_err(|_| {
                format!("Malformed STARsolo feature record: invalid transcript count {ntr_token}")
            })?;
            feat_vec_u32.resize(2 * ntr as usize, 0);
            for ii in 0..2 * ntr as usize {
                let value_token = tokens.next().ok_or_else(|| {
                    "Malformed STARsolo feature record: missing transcript value".to_string()
                })?;
                feat_vec_u32[ii] = value_token.parse().map_err(|_| {
                    format!(
                        "Malformed STARsolo feature record: invalid transcript value {value_token}"
                    )
                })?;
            }
        }
        _ => {}
    }

    let cbmatch_token = tokens
        .next()
        .ok_or_else(|| "Malformed STARsolo feature record: missing CB match".to_string())?;
    *cbmatch = cbmatch_token.parse().map_err(|_| {
        format!("Malformed STARsolo feature record: invalid CB match {cbmatch_token}")
    })?;
    Ok(true)
}
