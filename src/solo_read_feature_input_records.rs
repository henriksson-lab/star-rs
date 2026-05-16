#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloReadFeature::inputRecords` at STAR/source/SoloReadFeature_inputRecords.cpp:8. Args: cbP: uint32, cbPstride: uint32, cbReadCountTotal: vector<uint32>, readInfo: vector<readInfoStruct>, readFlagCounts: SoloReadFlagClass, nReadPerCBunique1: vector<uint32>, nReadPerCBmulti1: vector<uint32>"]
pub fn soloreadfeature_inputrecords_l8_soloreadfeature_inputrecords(
    rf: &mut crate::solo_read_feature::SoloReadFeature,
    p_solo: &crate::parameters_solo::ParametersSolo,
    sj_all: &[Vec<u64>; 2],
    cb_p: &mut [Vec<u32>],
    cb_p_stride: u32,
    cb_read_count_total: &[u32],
    read_info: &mut [crate::solo_feature::SoloFeatureReadInfo],
    read_flag_counts: &mut crate::solo_common::SoloReadFlagClass,
    n_read_per_cb_unique1: &mut [u32],
    n_read_per_cb_multi1: &mut [u32],
) -> Result<(), String> {
    let mut prev_iread = u64::MAX;
    let mut tr_id_dist = Vec::<u32>::new();

    for line in rf.stream_reads.lines() {
        let mut tokens = line.split_whitespace();
        let mut feature = 0_u32;
        let mut umi = 0_u64;
        let mut iread = 0_u64;
        let mut cbmatch = 0_i32;

        if !soloinputfeatureumi_l5_soloinputfeatureumi(
            &mut tokens,
            rf.feature_type,
            rf.read_index_yes,
            sj_all,
            &mut iread,
            &mut cbmatch,
            &mut feature,
            &mut umi,
            &mut tr_id_dist,
            read_flag_counts,
        )? {
            break;
        }

        if feature == u32::MAX && !rf.read_index_yes {
            continue;
        }

        let mut read_is_counted = false;
        let feat_good = feature != u32::MAX;
        let mut no_mm_to_wl_without_exact = false;
        let mut no_too_many_wl_matches = false;
        let mut cb = -1_i64;

        if cbmatch <= 1 {
            let cb_token = tokens.next().ok_or_else(|| {
                "Malformed STARsolo feature record: missing cell barcode index".to_string()
            })?;
            cb = cb_token.parse::<i64>().map_err(|_| {
                format!("Malformed STARsolo feature record: invalid cell barcode index {cb_token}")
            })?;
            if cb < 0 || cb as usize >= cb_read_count_total.len() {
                return Err(format!(
                    "Malformed STARsolo feature record: cell barcode index {cb} is outside whitelist"
                ));
            }
            if p_solo.cb_match_wl.one_exact && cbmatch == 1 && cb_read_count_total[cb as usize] == 0
            {
                no_mm_to_wl_without_exact = true;
            } else {
                if !p_solo.cb_wl_yes {
                    cb = servicefuns_l294_binarysearchexact(
                        cb as u64,
                        &p_solo.cb_wl,
                        p_solo.cb_wl_size as u64,
                    );
                    if cb + 1 == 0 {
                        continue;
                    }
                }

                if feat_good {
                    read_is_counted = true;
                    let cb_index = cb as usize;
                    let old_len = cb_p[cb_index].len();
                    cb_p[cb_index].resize(old_len + cb_p_stride as usize, 0);
                    cb_p[cb_index][old_len] = feature;
                    cb_p[cb_index][old_len + 1] = umi as u32;
                    if rf.read_index_yes {
                        cb_p[cb_index][old_len + 2] = iread as u32;
                    }
                } else if rf.read_info_yes {
                    read_info[iread as usize].cb = cb;
                    read_info[iread as usize].umi = umi;
                }
            }
        } else {
            let mut ptot = 0.0_f64;
            let mut pmax = 0.0_f64;
            for _ in 0..cbmatch as u32 {
                let cbin_token = tokens.next().ok_or_else(|| {
                    "Malformed STARsolo feature record: missing candidate cell barcode index"
                        .to_string()
                })?;
                let cbin = cbin_token.parse::<u32>().map_err(|_| {
                    format!(
                        "Malformed STARsolo feature record: invalid candidate cell barcode index {cbin_token}"
                    )
                })?;
                let qin_token = tokens.next().ok_or_else(|| {
                    "Malformed STARsolo feature record: missing candidate cell barcode quality"
                        .to_string()
                })?;
                if cbin as usize >= cb_read_count_total.len() {
                    return Err(format!(
                        "Malformed STARsolo feature record: candidate cell barcode index {cbin} is outside whitelist"
                    ));
                }
                if qin_token.is_empty() {
                    return Err(
                        "Malformed STARsolo feature record: empty candidate cell barcode quality"
                            .to_string(),
                    );
                }
                let mut qin = qin_token.as_bytes()[0] as i32 - p_solo.qs_base as i32;
                qin = std::cmp::min(qin, p_solo.qs_max as i32);
                if cb_read_count_total[cbin as usize] > 0 {
                    let pin = cb_read_count_total[cbin as usize] as f64
                        * 10.0_f64.powf(-qin as f64 / 10.0);
                    ptot += pin;
                    if pin > pmax {
                        cb = cbin as i64;
                        pmax = pin;
                    }
                }
            }

            if ptot > 0.0 && pmax >= p_solo.cb_min_p * ptot {
                if feat_good {
                    read_is_counted = true;
                    let cb_index = cb as usize;
                    let old_len = cb_p[cb_index].len();
                    cb_p[cb_index].resize(old_len + cb_p_stride as usize, 0);
                    cb_p[cb_index][old_len] = feature;
                    cb_p[cb_index][old_len + 1] = umi as u32;
                    if rf.read_index_yes {
                        cb_p[cb_index][old_len + 2] = iread as u32;
                    }
                } else if rf.read_info_yes {
                    read_info[iread as usize].cb = cb;
                    read_info[iread as usize].umi = umi;
                }
            } else {
                no_too_many_wl_matches = true;
            }
        }

        if !rf.read_index_yes || iread != prev_iread {
            prev_iread = iread;
            if feat_good {
                if cbmatch == 0 {
                    rf.stats.v[SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_EXACT] += 1;
                } else if no_mm_to_wl_without_exact {
                    rf.stats.v[SOLO_READ_FEATURE_STAT_NO_MM_TO_WL_WITHOUT_EXACT] += 1;
                } else if no_too_many_wl_matches {
                    rf.stats.v[SOLO_READ_FEATURE_STAT_NO_TOO_MANY_WL_MATCHES] += 1;
                }
            }

            if read_is_counted {
                if feature < (1_u32 << 31) {
                    n_read_per_cb_unique1[cb as usize] += 1;
                } else {
                    n_read_per_cb_multi1[cb as usize] += 1;
                }
            }

            let read_stats_yes = p_solo
                .read_stats_yes
                .get(rf.feature_type as usize)
                .copied()
                .unwrap_or(false);
            if read_stats_yes {
                if read_is_counted {
                    if (read_flag_counts.flag >> SOLO_READ_FLAG_FEATURE_U) & 1 != 0 {
                        read_flag_counts.flag |= 1_u32 << SOLO_READ_FLAG_COUNTED_U;
                    }
                    if (read_flag_counts.flag >> SOLO_READ_FLAG_FEATURE_M) & 1 != 0 {
                        read_flag_counts.flag |= 1_u32 << SOLO_READ_FLAG_COUNTED_M;
                    }
                }

                read_flag_counts.flag |= 1_u32 << SOLO_READ_FLAG_CB_MATCH;
                if cbmatch == 0 {
                    read_flag_counts.flag |= 1_u32 << SOLO_READ_FLAG_CB_PERFECT;
                    let counts = read_flag_counts
                        .flag_counts
                        .entry(cb as u64)
                        .or_insert([0; SOLO_READ_FLAG_N_BITS]);
                    for ibit in 0..SOLO_READ_FLAG_N_BITS {
                        counts[ibit] += ((read_flag_counts.flag >> ibit) & 1_u32) as u64;
                    }
                } else if cbmatch == 1 && !no_mm_to_wl_without_exact {
                    read_flag_counts.flag |= 1_u32 << SOLO_READ_FLAG_CB_MM_UNIQUE;
                    let counts = read_flag_counts
                        .flag_counts
                        .entry(cb as u64)
                        .or_insert([0; SOLO_READ_FLAG_N_BITS]);
                    for ibit in 0..SOLO_READ_FLAG_N_BITS {
                        counts[ibit] += ((read_flag_counts.flag >> ibit) & 1_u32) as u64;
                    }
                } else if cbmatch > 1 && !no_too_many_wl_matches {
                    read_flag_counts.flag |= 1_u32 << SOLO_READ_FLAG_CB_MM_MULTIPLE;
                    let counts = read_flag_counts
                        .flag_counts
                        .entry(cb as u64)
                        .or_insert([0; SOLO_READ_FLAG_N_BITS]);
                    for ibit in 0..SOLO_READ_FLAG_N_BITS {
                        counts[ibit] += ((read_flag_counts.flag >> ibit) & 1_u32) as u64;
                    }
                } else {
                    for ibit in 0..SOLO_READ_FLAG_N_BITS {
                        read_flag_counts.flag_counts_no_cb[ibit] +=
                            ((read_flag_counts.flag >> ibit) & 1_u32) as u64;
                    }
                }
            }
        }
    }
    Ok(())
}
