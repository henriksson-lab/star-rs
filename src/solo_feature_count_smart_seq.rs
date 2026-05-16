#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::countSmartSeq` at STAR/source/SoloFeature_countSmartSeq.cpp:9. Args: "]
pub fn solofeature_countsmartseq_l9_solofeature_countsmartseq(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    p: &crate::parameters_chimeric::Parameters,
    time_finished_redistribution: &str,
    time_finished_collapsing: &str,
    time_finished_counting: &str,
) -> Result<String, String> {
    for ii in 0..p.run_thread_n as usize {
        let rf_in = solo_feature.read_feat_all[ii].clone();
        let read_feat_sum = solo_feature
            .read_feat_sum
            .as_mut()
            .expect("SoloFeature::countSmartSeq requires readFeatSum");
        soloreadfeature_l47_soloreadfeature_addstats(read_feat_sum, &rf_in);
    }

    let mut log_main = solofeature_redistributereadsbycb_l8_solofeature_redistributereadsbycb(
        solo_feature,
        &p.p_solo,
        p.run_thread_n,
    );
    log_main.push_str(time_finished_redistribution);
    log_main.push_str(" ... Finished redistribution of reads from Solo read files\n");

    solo_feature
        .n_read_per_cb
        .resize(solo_feature.n_cb as usize, 0);
    solo_feature
        .cb_feature_umi_map
        .resize(solo_feature.n_cb as usize, ());

    let mut v_cell_feature_count: Vec<Vec<[u32; 3]>> = vec![Vec::new(); solo_feature.n_cb as usize];

    for ired in 0..solo_feature.redistr_files_cb_first.len().saturating_sub(1) {
        let i_cb1 = solo_feature.redistr_files_cb_first[ired] as usize;
        let i_cb2 = solo_feature.redistr_files_cb_first[ired + 1] as usize;
        let mut v_feature_umi: Vec<Vec<(u32, u64)>> = vec![Vec::new(); i_cb2 - i_cb1];

        let mut tokens = solo_feature.redistr_files_streams[ired].split_whitespace();
        loop {
            let mut feature = 0_u32;
            let mut umi = 0_u64;
            let mut iread = 0_u64;
            let mut cbmatch = 0_i32;
            let mut tr_id_dist = Vec::new();
            let mut read_flag_counts = crate::solo_common::SoloReadFlagClass::default();
            if !soloinputfeatureumi_l5_soloinputfeatureumi(
                &mut tokens,
                solo_feature.feature_type,
                false,
                &p.sj_all,
                &mut iread,
                &mut cbmatch,
                &mut feature,
                &mut umi,
                &mut tr_id_dist,
                &mut read_flag_counts,
            )? {
                break;
            }

            let cb_token = tokens.next().ok_or_else(|| {
                "Malformed STARsolo SmartSeq record: missing cell barcode index".to_string()
            })?;
            let cb: i64 = cb_token.parse().map_err(|_| {
                format!("Malformed STARsolo SmartSeq record: invalid cell barcode index {cb_token}")
            })?;
            if cb < 0 || cb as usize >= solo_feature.ind_cb_wl.len() {
                return Err(format!(
                    "Malformed STARsolo SmartSeq record: cell barcode index {cb} is outside whitelist"
                ));
            }
            if feature == u32::MAX {
                continue;
            }
            let icb = solo_feature.ind_cb_wl[cb as usize] as usize;
            if icb < i_cb1 || icb >= i_cb2 {
                return Err(format!(
                    "Malformed STARsolo SmartSeq record: redistributed cell barcode index {icb} is outside chunk"
                ));
            }
            v_feature_umi[icb - i_cb1].push((feature, umi));
            solo_feature.n_read_per_cb[icb] += 1;
        }

        for icb in i_cb1..i_cb2 {
            if solo_feature.n_read_per_cb[icb] == 0 {
                continue;
            }

            let records = &mut v_feature_umi[icb - i_cb1];
            records.sort_by(|a, b| {
                if a.0 == b.0 {
                    a.1.cmp(&b.1)
                } else {
                    a.0.cmp(&b.0)
                }
            });

            v_cell_feature_count[icb].reserve(8192);
            v_cell_feature_count[icb].push([records[0].0, 1, 1]);
            for fu in 1..records.len() {
                if records[fu].0 != records[fu - 1].0 {
                    v_cell_feature_count[icb].push([records[fu].0, 1, 1]);
                } else {
                    let last = v_cell_feature_count[icb].len() - 1;
                    v_cell_feature_count[icb][last][1] += 1;
                    if records[fu].1 != records[fu - 1].1 {
                        v_cell_feature_count[icb][last][2] += 1;
                    }
                }
            }
        }
    }

    log_main.push_str(time_finished_collapsing);
    log_main.push_str(" ... Finished reading / collapsing\n");

    solo_feature
        .n_umi_per_cb
        .resize(solo_feature.n_cb as usize, 0);
    solo_feature
        .n_gene_per_cb
        .resize(solo_feature.n_cb as usize, 0);

    solo_feature.count_mat_stride = p.p_solo.umi_dedup.yes_n + 1;
    let ccg_n: usize = v_cell_feature_count.iter().map(|cbf| cbf.len()).sum();
    solo_feature
        .count_cell_gene_umi
        .resize(ccg_n * solo_feature.count_mat_stride as usize, 0);
    solo_feature
        .count_cell_gene_umi_index
        .resize(solo_feature.n_cb as usize + 1, 0);
    solo_feature.count_cell_gene_umi_index[0] = 0;

    for (icb, cbf) in v_cell_feature_count
        .iter()
        .enumerate()
        .take(solo_feature.n_cb as usize)
    {
        solo_feature.n_gene_per_cb[icb] = cbf.len() as u32;
        solo_feature.count_cell_gene_umi_index[icb + 1] = solo_feature.count_cell_gene_umi_index
            [icb]
            + solo_feature.n_gene_per_cb[icb] * solo_feature.count_mat_stride;

        let mut ig = 0usize;
        let mut ic = solo_feature.count_cell_gene_umi_index[icb] as usize;
        while ic < solo_feature.count_cell_gene_umi_index[icb + 1] as usize {
            solo_feature.count_cell_gene_umi[ic] = cbf[ig][0];
            if p.p_solo.umi_dedup.yes_b[0] {
                solo_feature.count_cell_gene_umi[ic + p.p_solo.umi_dedup.count_ind_i[0] as usize] =
                    cbf[ig][1];
            }
            if p.p_solo.umi_dedup.yes_b[1] {
                solo_feature.count_cell_gene_umi[ic + p.p_solo.umi_dedup.count_ind_i[1] as usize] =
                    cbf[ig][2];
            }
            solo_feature.n_umi_per_cb[icb] +=
                solo_feature.count_cell_gene_umi[ic + p.p_solo.umi_dedup.count_ind_main as usize];
            ic += solo_feature.count_mat_stride as usize;
            ig += 1;
        }
    }

    solo_feature.n_read_per_cb_total = solo_feature.n_read_per_cb.clone();
    solo_feature.n_read_per_cb_unique = solo_feature.n_read_per_cb.clone();

    if let Some(read_feat_sum) = solo_feature.read_feat_sum.as_mut() {
        for icb in 0..solo_feature.n_cb as usize {
            read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_WL_MATCH] +=
                solo_feature.n_read_per_cb_total[icb] as u64;
            read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE] +=
                solo_feature.n_read_per_cb_unique[icb] as u64;
            read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_UMIS] +=
                solo_feature.n_umi_per_cb[icb] as u64;
            if solo_feature.n_gene_per_cb[icb] > 0 {
                read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_CELL_BARCODES] += 1;
            }
        }
        read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_EXACT] =
            read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_WL_MATCH];
    }

    log_main.push_str(time_finished_counting);
    log_main.push_str(" ... Finished SmartSeq counting\n");
    Ok(log_main)
}
