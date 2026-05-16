#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::quantTranscript` at STAR/source/SoloFeature_quantTranscript.cpp:12. Args: "]
pub fn solofeature_quanttranscript_l12_solofeature_quanttranscript(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    p_solo: &crate::parameters_solo::ParametersSolo,
    trans: &crate::transcriptome::Transcriptome,
    run_thread_n: i32,
    cluster_contents: &str,
    time_finished_input: &str,
    time_finished_cluster_prefix: &str,
    time_finished_quantification: &str,
) -> crate::solo_filtered_cells::SoloFeatureQuantTranscriptResult {
    let mut result = crate::solo_filtered_cells::SoloFeatureQuantTranscriptResult::default();
    if p_solo.cluster_cb_file == "-" {
        result.returned_no_cluster_file = true;
        return result;
    }

    let mut cluster_cb_ind = std::collections::BTreeMap::<u32, u32>::new();
    let mut cluster_ind = std::collections::BTreeSet::<u32>::new();
    let mut words = cluster_contents.split_whitespace();
    while let Some(seq1) = words.next() {
        let Some(icl1s) = words.next() else {
            break;
        };
        let icl1 = icl1s.parse::<u32>().unwrap_or(0);
        let mut cb1 = 0_u64;
        if sequencefuns_l249_convertnuclstrtoint64(seq1, &mut cb1) != 0 {
            let cb1ind = p_solo.cb_wl.partition_point(|x| *x < cb1) as u32;
            if cb1ind < p_solo.cb_wl.len() as u32 {
                cluster_cb_ind.insert(cb1ind, icl1);
                cluster_ind.insert(icl1);
            } else {
                result.log_main.push_str(&format!(
                    "WARNING: cluster CB sequence not present in whitelist and is ignored: {}\n",
                    seq1
                ));
            }
        } else {
            result.log_main.push_str(&format!(
                "WARNING: cluster CB sequence contains non-ACGT base and is ignored: {}\n",
                seq1
            ));
        }
    }

    let tr_dist_count = solo_feature
        .read_feat_sum
        .as_ref()
        .map(|rf| rf.transcript_dist_count.clone())
        .unwrap_or_default();
    let mut tr_dist_fun = vec![0.0_f64; tr_dist_count.len()];
    let mut tr_dist_fun_tr_factor = vec![0.0_f64; trans.n_tr as usize];

    let run_aver_n = 50_i32;
    let run_aver_start = 0_i32;
    for ii in 0..run_aver_start as usize {
        tr_dist_fun[ii] = tr_dist_count[ii] as f64;
    }
    let upper = tr_dist_count.len() as i32 - run_aver_n - 1;
    for ii in run_aver_start..upper {
        let start = std::cmp::max(run_aver_start, ii - run_aver_n) as usize;
        let end = (ii + run_aver_n + 1) as usize;
        let sum: u32 = tr_dist_count[start..end].iter().sum();
        let denom = std::cmp::min(2 * run_aver_n + 1, ii - run_aver_start + run_aver_n);
        tr_dist_fun[ii as usize] = sum as f64 / denom as f64;
    }

    let mut imax = 1000_usize;
    while imax + 1 < tr_dist_fun.len() && tr_dist_fun[imax + 1] > tr_dist_fun[imax] {
        imax += 1;
    }
    result.log_main.push_str(&format!(
        "SoloQuant: distance distribution past maximum = {}\n",
        imax
    ));
    while imax + 1 < tr_dist_fun.len() && tr_dist_fun[imax + 1] < tr_dist_fun[imax] {
        imax += 1;
    }
    result.log_main.push_str(&format!(
        "SoloQuant: distance distribution cutoff = {}\n",
        imax
    ));
    tr_dist_fun.truncate(imax);

    let norm1: f64 = tr_dist_fun.iter().sum();
    let mut dist_file = String::new();
    if norm1 != 0.0 {
        for ff in &mut tr_dist_fun {
            *ff /= norm1;
            dist_file.push_str(&format!("{}\n", *ff));
        }
    }
    result.files.insert(
        format!(
            "{}transcriptEndDistanceDistribution.txt",
            solo_feature.output_prefix
        ),
        dist_file,
    );

    let mut tr_dist_fun_cum = tr_dist_fun.clone();
    for ii in 1..tr_dist_fun_cum.len() {
        tr_dist_fun_cum[ii] += tr_dist_fun_cum[ii - 1];
    }
    for ii in 0..trans.n_tr as usize {
        if trans.tr_len[ii] as usize <= tr_dist_fun_cum.len() && trans.tr_len[ii] > 0 {
            tr_dist_fun_tr_factor[ii] = -tr_dist_fun_cum[trans.tr_len[ii] as usize - 1].ln();
        }
    }
    for ff in &mut tr_dist_fun {
        *ff = ff.ln();
    }

    let mut map_tr_dist =
        std::collections::BTreeMap::<u32, std::collections::BTreeMap<u64, Vec<(u32, f64)>>>::new();
    for i_thread in 0..run_thread_n as usize {
        let Some(read_feat) = solo_feature.read_feat_all.get(i_thread) else {
            continue;
        };
        for line1 in read_feat.stream_reads.lines() {
            let mut fields = line1.split_whitespace();
            let Some(cb_s) = fields.next() else {
                continue;
            };
            let Some(umi_s) = fields.next() else {
                continue;
            };
            let Some(ntr_s) = fields.next() else {
                continue;
            };
            let cb = cb_s.parse::<u32>().unwrap_or(0);
            let mut umi = umi_s.parse::<u64>().unwrap_or(0);
            let n_tr = ntr_s.parse::<u32>().unwrap_or(0);
            let Some(cb_cl) = cluster_cb_ind.get(&cb).copied() else {
                for _ in 0..n_tr {
                    let _ = fields.next();
                    let _ = fields.next();
                }
                continue;
            };
            umi += (cb as u64) << 32;

            let mut td = Vec::<(u32, f64)>::with_capacity(n_tr as usize);
            for _ in 0..n_tr {
                let Some(tr1_s) = fields.next() else {
                    break;
                };
                let Some(d1_s) = fields.next() else {
                    break;
                };
                let tr1 = tr1_s.parse::<u32>().unwrap_or(0);
                let d1 = d1_s.parse::<usize>().unwrap_or(0);
                if d1 >= tr_dist_fun.len() || tr1 as usize >= tr_dist_fun_tr_factor.len() {
                    continue;
                }
                td.push((tr1, tr_dist_fun[d1] + tr_dist_fun_tr_factor[tr1 as usize]));
            }
            if td.is_empty() {
                continue;
            }
            td.sort_by_key(|x| x.0);

            let cl_map = map_tr_dist.entry(cb_cl).or_default();
            if let Some(old_td) = cl_map.get(&umi).cloned() {
                let mut inew = 0usize;
                let mut td1 = Vec::<(u32, f64)>::with_capacity(old_td.len());
                for old in old_td {
                    while inew < td.len() && old.0 > td[inew].0 {
                        inew += 1;
                    }
                    if inew == td.len() {
                        break;
                    }
                    if old.0 == td[inew].0 {
                        td1.push((td[inew].0, old.1 + td[inew].1));
                    }
                }
                cl_map.insert(umi, td1);
            } else {
                cl_map.insert(umi, td);
            }
        }
    }
    result.log_main.push_str(&format!(
        "{} ... Transcript3p counting: finished input\n",
        time_finished_input
    ));

    let mut cluster_expression = std::collections::BTreeMap::<u32, Vec<f64>>::new();
    for (cluster, mut cl_tr_dist) in map_tr_dist {
        let mut tr_unique = vec![0.0_f64; trans.n_tr as usize];
        let mut tr_initial = vec![0.0_f64; trans.n_tr as usize];
        let mut n_umi_tot = 0_u64;
        let mut n_umi0 = 0_u64;
        let mut n_umi1 = 0_u64;

        let keys: Vec<u64> = cl_tr_dist.keys().copied().collect();
        for key in keys {
            let Some(tr_dist) = cl_tr_dist.get_mut(&key) else {
                continue;
            };
            if tr_dist.is_empty() {
                cl_tr_dist.remove(&key);
                n_umi0 += 1;
            } else if tr_dist.len() == 1 {
                let tr = tr_dist[0].0 as usize;
                tr_unique[tr] += 1.0;
                tr_initial[tr] += 1.0;
                cl_tr_dist.remove(&key);
                n_umi1 += 1;
                n_umi_tot += 1;
            } else {
                let max1 = tr_dist
                    .iter()
                    .map(|td| td.1)
                    .fold(f64::NEG_INFINITY, f64::max);
                let tr_dist_len = tr_dist.len() as f64;
                for tt in tr_dist.iter_mut() {
                    tr_initial[tt.0 as usize] += 1.0 / tr_dist_len;
                    tt.1 = (tt.1 - max1).exp();
                }
                n_umi_tot += 1;
            }
        }
        result.log_main.push_str(&format!(
            "{} ... Transcript3p counting: cluster {} nUMItot={} nUMI0={} nUMI1={}\n",
            time_finished_input, cluster, n_umi_tot, n_umi0, n_umi1
        ));

        let mut th_old_new = [tr_initial, vec![0.0_f64; trans.n_tr as usize]];
        let mut old_i = 0usize;
        let mut new_i = 1usize;
        let mut tr_converged = vec![false; trans.n_tr as usize];
        for iteration in 0..10000_u32 {
            let (old_slice, new_slice): (Vec<f64>, &mut Vec<f64>) = {
                let old = th_old_new[old_i].clone();
                let new = &mut th_old_new[new_i];
                (old, new)
            };
            new_slice.copy_from_slice(&tr_unique);
            for tr_dist in cl_tr_dist.values() {
                let mut denom1 = 0.0_f64;
                for td in tr_dist {
                    denom1 += td.1 * old_slice[td.0 as usize];
                }
                for td in tr_dist {
                    if !tr_converged[td.0 as usize] && denom1 != 0.0 {
                        new_slice[td.0 as usize] += td.1 * old_slice[td.0 as usize] / denom1;
                    }
                }
            }

            let diff_threshold_max = 1e-5_f64;
            let diff_threshold_one = diff_threshold_max * 0.1;
            let expr_threshold = 1e-8_f64 * n_umi_tot as f64;
            let mut diff_max = 0.0_f64;
            let mut diff_sum = 0.0_f64;
            let mut above_thr_n = 0_u64;
            let mut above_thr_expr_sum = 0.0_f64;
            let mut above_thr_one_n = 0_u64;
            for itr in 0..new_slice.len() {
                if tr_converged[itr] || old_slice[itr] == 0.0 {
                    continue;
                }
                let diff1 = (new_slice[itr] - old_slice[itr]).abs() / old_slice[itr];
                diff_sum += diff1;
                diff_max = diff_max.max(diff1);
                if diff1 > diff_threshold_max {
                    above_thr_n += 1;
                    above_thr_expr_sum += new_slice[itr];
                }
                if new_slice[itr] < expr_threshold {
                    tr_converged[itr] = true;
                    tr_unique[itr] = 0.0;
                }
                if diff1 < diff_threshold_one {
                    tr_converged[itr] = true;
                    tr_unique[itr] = new_slice[itr];
                } else {
                    above_thr_one_n += 1;
                }
            }
            result.stdout.push_str(&format!(
                "{} {} {} {} {} {}\n",
                iteration, diff_max, diff_sum, above_thr_n, above_thr_expr_sum, above_thr_one_n
            ));
            if diff_max < diff_threshold_max {
                break;
            }
            std::mem::swap(&mut old_i, &mut new_i);
        }

        let mut th_out = th_old_new[new_i].clone();
        let mut norm = 0.0_f64;
        for itr in 0..th_out.len() {
            th_out[itr] *= tr_dist_fun_tr_factor[itr].exp();
            norm += th_out[itr];
        }
        if norm != 0.0 {
            norm = n_umi_tot as f64 / norm;
            for tt in &mut th_out {
                *tt *= norm;
            }
        }
        cluster_expression.insert(cluster, th_out);
        result.log_main.push_str(&format!(
            "{} ... Transcript3p counting: finished cluster{}\n",
            time_finished_cluster_prefix, cluster
        ));
    }

    let mut matrix = String::new();
    matrix.push_str("%%MatrixMarket matrix coordinate real general\n%\n");
    let mut n_cell_gene_entries = 0_u32;
    for ctpm in cluster_expression.values() {
        for tt in ctpm {
            if *tt > 0.0 {
                n_cell_gene_entries += 1;
            }
        }
    }
    let max_cluster = cluster_ind.iter().next_back().copied().unwrap_or(0);
    matrix.push_str(&format!(
        "{} {} {}\n",
        trans.n_tr, max_cluster, n_cell_gene_entries
    ));
    for (cluster, ctpm) in &cluster_expression {
        for (itr, value) in ctpm.iter().enumerate() {
            if *value > 0.0 {
                matrix.push_str(&format!("{} {} {}\n", itr + 1, cluster, value));
            }
        }
    }
    result.files.insert(
        format!(
            "{}{}",
            solo_feature.output_prefix,
            p_solo
                .out_file_names
                .get(3)
                .cloned()
                .unwrap_or_else(|| "matrix.mtx".to_string())
        ),
        matrix,
    );

    let mut features = String::new();
    for ii in 0..trans.n_tr as usize {
        let gene = trans.tr_gene.get(ii).copied().unwrap_or(0) as usize;
        features.push_str(&format!(
            "{}\t{}\t{}\n",
            trans.tr_id.get(ii).cloned().unwrap_or_default(),
            trans.tr_len.get(ii).copied().unwrap_or(0),
            trans.ge_name.get(gene).cloned().unwrap_or_default()
        ));
    }
    result.files.insert(
        format!("{}/features.tsv", solo_feature.output_prefix),
        features,
    );
    result.log_main.push_str(&format!(
        "{} ... Transcript3p counting: finished transcript quantification\n",
        time_finished_quantification
    ));
    result.cluster_expression = cluster_expression;
    result
}
