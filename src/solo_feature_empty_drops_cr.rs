#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::emptyDrops_CR` at STAR/source/SoloFeature_emptyDrops_CR.cpp:10. Args: "]
pub fn solofeature_emptydrops_cr_l10_solofeature_emptydrops_cr(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    p_solo: &crate::parameters_solo::ParametersSolo,
) -> crate::solo_filtered_cells::SoloEmptyDropsCrResult {
    use std::collections::{BTreeMap, BTreeSet};

    let mut result = crate::solo_filtered_cells::SoloEmptyDropsCrResult::default();
    let ed = &p_solo.cell_filter.ed_cr;
    if solo_feature.n_cb <= ed.ind_min {
        let msg = format!(
            "emptyDrops_CR filtering: total number of cells: nCB={} is smaller than emptyCellMinIndex={}, which is the starting index for the *true empty* cells. The additional non-empty cells will not be detected.\n",
            solo_feature.n_cb, ed.ind_min
        );
        result.log_main.push_str(&msg);
        result.early_return = Some(msg);
        return result;
    }

    result.log_main.push_str(&format!(
        "{} ... starting emptyDrops_CR filtering\n",
        timefunctions_l4_timemonthdaytime()
    ));

    let main_shift = p_solo.umi_dedup.count_ind_main as usize;
    let stride = solo_feature.count_mat_stride as usize;

    let mut feat_det = BTreeSet::<u32>::new();
    for icb in 0..solo_feature.n_cb as usize {
        for ig in 0..solo_feature.n_gene_per_cb[icb] as usize {
            let irec = solo_feature.count_cell_gene_umi_index[icb] as usize + ig * stride;
            if solo_feature.count_cell_gene_umi[irec + main_shift] > 0 {
                feat_det.insert(solo_feature.count_cell_gene_umi[irec]);
            }
        }
    }
    result.feat_det_n = feat_det.len() as u32;

    let mut ind_count: Vec<(u32, u32)> = (0..solo_feature.n_cb)
        .map(|ii| (ii, solo_feature.n_umi_per_cb[ii as usize]))
        .collect();
    ind_count.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let mut amb_count = vec![0_u32; solo_feature.features_number.max(0) as usize];
    let amb_end = solo_feature.n_cb.min(ed.ind_max);
    for icb in ed.ind_min..amb_end {
        let icb1 = ind_count[icb as usize].0 as usize;
        for ig in 0..solo_feature.n_gene_per_cb[icb1] as usize {
            let irec = solo_feature.count_cell_gene_umi_index[icb1] as usize + ig * stride;
            let gene = solo_feature.count_cell_gene_umi[irec] as usize;
            amb_count[gene] += solo_feature.count_cell_gene_umi[irec + main_shift];
        }
    }
    result.log_main.push_str(&format!(
        "{} ... finished ambient cells counting\n",
        timefunctions_l4_timemonthdaytime()
    ));

    let mut amb_count_freq = BTreeMap::<u32, u32>::new();
    for &ac in &amb_count {
        *amb_count_freq.entry(ac).or_insert(0) += 1;
    }
    if amb_count_freq.len() <= 1 {
        let msg = "emptyDrops_CR filtering: empty cells contain no genes\n".to_string();
        result.log_main.push_str(&msg);
        result.early_return = Some(msg);
        return result;
    }
    if let Some(freq0) = amb_count_freq.get_mut(&0) {
        *freq0 =
            freq0.saturating_sub(solo_feature.features_number.max(0) as u32 - result.feat_det_n);
    }
    let max_freq = *amb_count_freq.keys().next_back().unwrap_or(&0);

    let mut amb_count_freq_sgt = vec![0.0_f64; max_freq as usize + 1];
    {
        let mut data: BTreeMap<u32, (u32, f64)> = BTreeMap::new();
        let mut total_obs = 0_u32;
        for (&obs, &freq) in &amb_count_freq {
            if obs != 0 {
                data.insert(obs, (freq, 0.0));
                total_obs = total_obs.saturating_add(obs.saturating_mul(freq));
            }
        }
        let mut p_zero = 0.0_f64;
        if data.len() >= 5 {
            let big_n: u32 = data.iter().map(|(&obs, &(freq, _))| obs * freq).sum();
            if let Some(&(freq1, _)) = data.get(&1) {
                p_zero = freq1 as f64 / big_n as f64;
            }
            let rows = data.len();
            let mut log_obs = Vec::with_capacity(rows);
            let mut log_z = Vec::with_capacity(rows);
            let mut prev_obs = 0_u32;
            let keys: Vec<u32> = data.keys().copied().collect();
            for (i, &obs) in keys.iter().enumerate() {
                let freq = data[&obs].0;
                let k = if i + 1 == keys.len() {
                    2 * obs - prev_obs
                } else {
                    keys[i + 1]
                };
                let z = 2.0 * freq as f64 / (k - prev_obs) as f64;
                log_obs.push((obs as f64).ln());
                log_z.push(z.ln());
                prev_obs = obs;
            }
            let mean_x = log_obs.iter().sum::<f64>() / rows as f64;
            let mean_y = log_z.iter().sum::<f64>() / rows as f64;
            let mut xys = 0.0;
            let mut xsquares = 0.0;
            for i in 0..rows {
                xys += (log_obs[i] - mean_x) * (log_z[i] - mean_y);
                xsquares += (log_obs[i] - mean_x) * (log_obs[i] - mean_x);
            }
            let slope = xys / xsquares;
            let intercept = mean_y - slope * mean_x;
            let smoothed = |i: u32| (intercept + slope * (i as f64).ln()).exp();
            let mut r_star = vec![0.0; rows];
            let mut indiff_vals_seen = false;
            for (r, &obs) in keys.iter().enumerate() {
                let obs1 = obs + 1;
                let y = obs1 as f64 * smoothed(obs1) / smoothed(obs);
                if let Some(&(next_n, _)) = data.get(&obs1) {
                    if !indiff_vals_seen {
                        let freq = data[&obs].0;
                        let x = obs1 as f64 * next_n as f64 / freq as f64;
                        let sd = ((obs1 as f64 * obs1 as f64) * next_n as f64
                            / (freq as f64 * freq as f64)
                            * (1.0 + next_n as f64 / freq as f64))
                            .sqrt();
                        if (x - y).abs() <= 1.96 * sd {
                            indiff_vals_seen = true;
                        } else {
                            r_star[r] = x;
                        }
                    }
                } else {
                    indiff_vals_seen = true;
                }
                if indiff_vals_seen {
                    r_star[r] = y;
                }
            }
            let big_n_prime: f64 = keys
                .iter()
                .enumerate()
                .map(|(r, obs)| data[obs].0 as f64 * r_star[r])
                .sum();
            for (r, obs) in keys.iter().enumerate() {
                if let Some(value) = data.get_mut(obs) {
                    value.1 = (1.0 - p_zero) * r_star[r] / big_n_prime;
                }
            }
        } else if total_obs > 0 {
            for value in data.values_mut() {
                value.1 = value.0 as f64 / total_obs as f64;
            }
            p_zero = data
                .get(&1)
                .map_or(0.0, |&(freq, _)| freq as f64 / total_obs as f64);
        }

        for freq in 0..=max_freq {
            amb_count_freq_sgt[freq as usize] = if freq == 0 {
                p_zero
            } else {
                data.get(&freq).map_or(0.0, |&(_, estimate)| estimate)
            };
        }
        if let Some(freq0) = amb_count_freq.get(&0) {
            if *freq0 > 0 {
                amb_count_freq_sgt[0] /= *freq0 as f64;
            }
        }
    }
    result.log_main.push_str(&format!(
        "{} ... finished SGT\n",
        timefunctions_l4_timemonthdaytime()
    ));

    let mut amb_profile_log_p = vec![0.0_f64; solo_feature.features_number.max(0) as usize];
    let mut amb_profile_p_non0 = Vec::new();
    let mut amb_profile_log_p_non0 = Vec::new();
    for ig in 0..amb_profile_log_p.len() {
        if feat_det.contains(&(ig as u32)) {
            amb_profile_log_p[ig] = amb_count_freq_sgt[amb_count[ig] as usize];
        }
    }
    let norm1: f64 = amb_profile_log_p.iter().sum();
    if norm1 <= 0.0 {
        let msg =
            "emptyDrops_CR filtering: ambient profile has zero total probability\n".to_string();
        result.log_main.push_str(&msg);
        result.early_return = Some(msg);
        return result;
    }
    for cf in &mut amb_profile_log_p {
        if *cf > 0.0 {
            *cf /= norm1;
            amb_profile_p_non0.push(*cf);
            *cf = cf.ln();
            amb_profile_log_p_non0.push(*cf);
        }
    }
    result.log_main.push_str(&format!(
        "{} ... finished ambient profile\n",
        timefunctions_l4_timemonthdaytime()
    ));

    let i_cand_first = solo_feature.filtered_cells.n_cells_simple as u32;
    result.candidate_first = i_cand_first;
    if i_cand_first as usize >= ind_count.len() {
        return result;
    }
    let median_index = (solo_feature.filtered_cells.n_cells_simple / 2) as usize;
    let median_count = solo_feature
        .n_umi_per_cb_sorted
        .get(median_index)
        .copied()
        .unwrap_or(0);
    let min_umi = ed
        .umi_min
        .max((ed.umi_min_frac_median * median_count as f64) as u32);
    result.min_umi = min_umi;
    let mut i_cand_last = i_cand_first;
    let cand_limit = i_cand_first.saturating_add(ed.cand_max_n);
    while i_cand_last < cand_limit && (i_cand_last as usize) < ind_count.len() {
        if ind_count[i_cand_last as usize].1 < min_umi {
            break;
        }
        i_cand_last += 1;
    }
    if i_cand_last == i_cand_first {
        result.log_main.push_str(&format!(
            "{} ... candidate cells: minUMI={}; number of candidate cells=0\n",
            timefunctions_l4_timemonthdaytime(),
            min_umi
        ));
        return result;
    }
    i_cand_last -= 1;
    result.candidate_last = i_cand_last;
    let n_candidates = i_cand_last - i_cand_first + 1;
    result.log_main.push_str(&format!(
        "{} ... candidate cells: minUMI={}; number of candidate cells={}\n",
        timefunctions_l4_timemonthdaytime(),
        min_umi,
        n_candidates
    ));

    let max_count = ind_count[i_cand_first as usize].1 as usize;
    let mut log_factorial = vec![0.0_f64; max_count + 1];
    for cc in 2..log_factorial.len() {
        log_factorial[cc] = log_factorial[cc - 1] + (cc as f64).ln();
    }

    let mut obs_log_prob = vec![0.0_f64; n_candidates as usize];
    for icand in 0..n_candidates as usize {
        let icell = ind_count[icand + i_cand_first as usize].0 as usize;
        obs_log_prob[icand] = solofeature_emptydrops_cr_l219_logmultinomialpdfsparse(
            &amb_profile_log_p,
            &solo_feature.count_cell_gene_umi,
            solo_feature.count_mat_stride,
            p_solo.umi_dedup.count_ind_main,
            solo_feature.count_cell_gene_umi_index[icell] as i64,
            solo_feature.n_gene_per_cb[icell],
            &log_factorial,
        );
    }
    result.log_main.push_str(&format!(
        "{} ... finished observed logProb\n",
        timefunctions_l4_timemonthdaytime()
    ));

    let mut sim_log_prob = vec![vec![0.0_f64; max_count + 1]; ed.sim_n as usize];
    if !amb_profile_p_non0.is_empty() {
        let mut cumulative = Vec::with_capacity(amb_profile_p_non0.len());
        let mut acc = 0.0_f64;
        for p1 in &amb_profile_p_non0 {
            acc += *p1;
            cumulative.push(acc);
        }
        for isim in 0..ed.sim_n as usize {
            let mut mt = [0_u32; 624];
            mt[0] = 19760110_u32.wrapping_mul(isim as u32 + 1);
            for i in 1..624 {
                mt[i] = 1812433253_u32
                    .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                    .wrapping_add(i as u32);
            }
            let mut mt_index = 624usize;
            let mut next_u32 = || {
                if mt_index >= 624 {
                    for i in 0..624 {
                        let y = (mt[i] & 0x8000_0000) | (mt[(i + 1) % 624] & 0x7fff_ffff);
                        mt[i] = mt[(i + 397) % 624] ^ (y >> 1);
                        if y & 1 != 0 {
                            mt[i] ^= 0x9908_b0df;
                        }
                    }
                    mt_index = 0;
                }
                let mut y = mt[mt_index];
                mt_index += 1;
                y ^= y >> 11;
                y ^= (y << 7) & 0x9d2c_5680;
                y ^= (y << 15) & 0xefc6_0000;
                y ^= y >> 18;
                y
            };
            let mut curr_counts = vec![0_u32; amb_profile_p_non0.len()];
            for ic in 1..=max_count {
                let u = (next_u32() as f64) / (u32::MAX as f64 + 1.0) * acc;
                let ig1 = cumulative
                    .binary_search_by(|probe| probe.partial_cmp(&u).unwrap())
                    .unwrap_or_else(|idx| idx)
                    .min(cumulative.len() - 1);
                curr_counts[ig1] += 1;
                sim_log_prob[isim][ic] =
                    sim_log_prob[isim][ic - 1] + amb_profile_log_p_non0[ig1] + (ic as f64).ln()
                        - (curr_counts[ig1] as f64).ln();
            }
        }
    }
    result.log_main.push_str(&format!(
        "{} ... finished simulations\n",
        timefunctions_l4_timemonthdaytime()
    ));

    let mut p_values: Vec<(u32, f64, f64)> = Vec::with_capacity(n_candidates as usize);
    for icand in 0..n_candidates as usize {
        let index = ind_count[icand + i_cand_first as usize].0;
        let count1 = ind_count[icand + i_cand_first as usize].1 as usize;
        let mut n_lower_p = 0_u32;
        for sp in &sim_log_prob {
            if sp[count1] < obs_log_prob[icand] {
                n_lower_p += 1;
            }
        }
        let pval = (1 + n_lower_p) as f64 / (1 + sim_log_prob.len()) as f64;
        p_values.push((index, pval, 0.0));
    }
    p_values.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let p_values_len = p_values.len();
    for (rank0, item) in p_values.iter_mut().enumerate() {
        let rank = rank0 + 1;
        item.2 = item.1 * p_values_len as f64 / rank as f64;
    }
    if p_values.len() > 1 {
        for ii in (0..p_values.len() - 1).rev() {
            p_values[ii].2 = p_values[ii].2.min(p_values[ii + 1].2);
        }
    }

    for &(index, _, padj) in &p_values {
        if padj <= ed.fdr {
            result.extra_cells += 1;
            if let Some(value) = solo_feature
                .filtered_cells
                .filt_vec_bool
                .get_mut(index as usize)
            {
                *value = true;
            }
        }
    }
    result.p_values = p_values;
    result.log_main.push_str(&format!(
        "{} ... finished emptyDrops_CR filtering: number of additional non-ambient cells={}\n",
        timefunctions_l4_timemonthdaytime(),
        result.extra_cells
    ));
    result
}

#[doc = "Original `logMultinomialPDFsparse` at STAR/source/SoloFeature_emptyDrops_CR.cpp:219. Args: ambProfileLogP: vector<double>, countCellGeneUMI: vector<uint32>, stride: uint32, shift: uint32, start: int64, nGenes: uint32, logFactorial: vector<double>"]
pub fn solofeature_emptydrops_cr_l219_logmultinomialpdfsparse(
    amb_profile_log_p: &[f64],
    count_cell_gene_umi: &[u32],
    stride: u32,
    shift: u32,
    start: i64,
    n_genes: u32,
    log_factorial: &[f64],
) -> f64 {
    let mut sum_count = 0u32;
    let mut sum_log_fac = 0.0;
    let mut sum_count_log_p = 0.0;
    for ig in 0..n_genes {
        let count1 = count_cell_gene_umi[(start + (ig * stride + shift) as i64) as usize];
        sum_count += count1;
        sum_log_fac += log_factorial[count1 as usize];
        sum_count_log_p += amb_profile_log_p
            [count_cell_gene_umi[(start + (ig * stride) as i64) as usize] as usize]
            * count1 as f64;
    }

    log_factorial[sum_count as usize] - sum_log_fac + sum_count_log_p
}
