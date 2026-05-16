#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::countCBgeneUMI` at STAR/source/SoloFeature_countCBgeneUMI.cpp:7. Args: "]
pub fn solofeature_countcbgeneumi_l7_solofeature_countcbgeneumi(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    p: &crate::parameters_chimeric::Parameters,
    p_solo: &crate::parameters_solo::ParametersSolo,
    raw_time: libc::time_t,
) -> Result<String, String> {
    use std::fmt::Write;

    let mut log_main = String::new();
    let feature_type = solo_feature.feature_type as usize;
    solo_feature.rgu_stride = 2;
    if p_solo
        .read_index_yes
        .get(feature_type)
        .copied()
        .unwrap_or(false)
    {
        solo_feature.rgu_stride = 3;
    }

    if p_solo
        .read_info_yes
        .get(feature_type)
        .copied()
        .unwrap_or(false)
    {
        solo_feature.read_info = vec![
            crate::solo_feature::SoloFeatureReadInfo {
                cb: -1,
                umi: u32::MAX as u64,
            };
            solo_feature.n_reads_input as usize
        ];
        writeln!(
            log_main,
            "{} ... Allocated and initialized readInfo array, nReadsInput = {}",
            timefunctions_l14_timemonthdaytime(raw_time),
            solo_feature.n_reads_input
        )
        .unwrap();
    }

    let read_feat_sum_counts = solo_feature
        .read_feat_sum
        .as_ref()
        .ok_or_else(|| "SoloFeature::countCBgeneUMI requires readFeatSum".to_string())?
        .cb_read_count
        .clone();

    let mut r_cbpa = vec![Vec::<u32>::new(); p_solo.cb_wl_size as usize + 1];
    solo_feature.r_cbp.clear();
    solo_feature.r_cbp.push(Vec::new());
    solo_feature.n_cb = 0;
    for ii in 0..p_solo.cb_wl_size as usize {
        if read_feat_sum_counts.get(ii).copied().unwrap_or(0) > 0 {
            solo_feature.r_cbp.push(Vec::new());
            solo_feature.n_cb += 1;
        }
        r_cbpa[ii + 1] = solo_feature.r_cbp[solo_feature.n_cb as usize].clone();
    }

    writeln!(
        log_main,
        "{} ... Finished allocating arrays for Solo {} GiB",
        timefunctions_l14_timemonthdaytime(raw_time),
        solo_feature.n_reads_mapped as f64 * solo_feature.rgu_stride as f64 * 4.0
            / 1024.0
            / 1024.0
            / 1024.0
    )
    .unwrap();

    solo_feature.read_flag_counts.flag_counts.clear();
    solo_feature.read_flag_counts.flag_counts_no_cb = [0; SOLO_READ_FLAG_N_BITS];
    let mut n_read_per_cb_unique1 = vec![0_u32; p_solo.cb_wl_size as usize];
    let mut n_read_per_cb_multi1 = vec![0_u32; p_solo.cb_wl_size as usize];
    let cb_read_count_exact = solo_feature
        .read_bar_sum
        .as_ref()
        .ok_or_else(|| "SoloFeature::countCBgeneUMI requires readBarSum".to_string())?
        .cb_read_count_exact
        .clone();

    let n_threads = p.run_thread_n.max(0) as usize;
    for ii in 0..n_threads {
        if ii >= solo_feature.read_feat_all.len() {
            break;
        }
        soloreadfeature_inputrecords_l8_soloreadfeature_inputrecords(
            &mut solo_feature.read_feat_all[ii],
            p_solo,
            &solo_feature.sj_all,
            &mut r_cbpa[..p_solo.cb_wl_size as usize],
            solo_feature.rgu_stride,
            &cb_read_count_exact,
            &mut solo_feature.read_info,
            &mut solo_feature.read_flag_counts,
            &mut n_read_per_cb_unique1,
            &mut n_read_per_cb_multi1,
        )?;
        let rf_in = solo_feature.read_feat_all[ii].clone();
        let read_feat_sum = solo_feature
            .read_feat_sum
            .as_mut()
            .ok_or_else(|| "SoloFeature::countCBgeneUMI requires readFeatSum".to_string())?;
        soloreadfeature_l47_soloreadfeature_addstats(read_feat_sum, &rf_in);
    }

    {
        let read_feat_sum = solo_feature
            .read_feat_sum
            .as_ref()
            .ok_or_else(|| "SoloFeature::countCBgeneUMI requires readFeatSum".to_string())?;
        for ii in 0..SOLO_READ_FLAG_N_BITS {
            solo_feature.read_flag_counts.flag_counts_no_cb[ii] +=
                read_feat_sum.read_flag.flag_counts_no_cb[ii];
        }
    }

    solo_feature
        .n_read_per_cb_total
        .resize(solo_feature.n_cb as usize, 0);
    solo_feature
        .n_read_per_cb_unique
        .resize(solo_feature.n_cb as usize, 0);
    for icb in 0..solo_feature.n_cb as usize {
        let wl = solo_feature.ind_cb[icb] as usize;
        solo_feature.n_read_per_cb_unique[icb] = n_read_per_cb_unique1[wl];
        solo_feature.n_read_per_cb_total[icb] =
            solo_feature.n_read_per_cb_unique[icb] + n_read_per_cb_multi1[wl];
    }

    solo_feature.r_cbp.clear();
    solo_feature
        .n_read_per_cb
        .resize(solo_feature.n_cb as usize, 0);
    let mut n_read_per_cb_max = 0_u32;
    for icb in 0..solo_feature.n_cb as usize {
        let wl = solo_feature.ind_cb[icb] as usize;
        let records = r_cbpa[wl].clone();
        let n_read = records.len() as u32 / solo_feature.rgu_stride;
        solo_feature.n_read_per_cb[icb] = n_read;
        n_read_per_cb_max = n_read_per_cb_max.max(n_read);
        solo_feature.r_cbp.push(records);
    }

    write!(
        log_main,
        "{} ... Finished reading reads from Solo files nCB={}, nReadPerCBmax={}",
        timefunctions_l14_timemonthdaytime(raw_time),
        solo_feature.n_cb,
        n_read_per_cb_max
    )
    .unwrap();
    if let Some(read_feat_sum) = &solo_feature.read_feat_sum {
        let yes_wl = read_feat_sum
            .stats
            .v
            .get(SOLO_READ_FEATURE_STAT_YES_WL_MATCH)
            .copied()
            .unwrap_or(0);
        writeln!(log_main, ", yesWLmatch={}", yes_wl).unwrap();
    } else {
        writeln!(log_main, ", yesWLmatch=0").unwrap();
    }

    solo_feature
        .n_umi_per_cb
        .resize(solo_feature.n_cb as usize, 0);
    solo_feature
        .n_gene_per_cb
        .resize(solo_feature.n_cb as usize, 0);
    solo_feature.count_mat_stride = p_solo.umi_dedup.yes_n + 1;
    solo_feature.count_cell_gene_umi.resize(
        solo_feature.n_reads_mapped as usize * solo_feature.count_mat_stride as usize / 5 + 16,
        0,
    );
    solo_feature
        .count_cell_gene_umi_index
        .resize(solo_feature.n_cb as usize + 1, 0);

    if p_solo.multi_map.yes_multi {
        solo_feature.count_mat_mult_s = 1 + p_solo.multi_map.yes_n * p_solo.umi_dedup.yes_n;
        solo_feature.count_mat_mult_m.resize(
            solo_feature.n_reads_mapped as usize * solo_feature.count_mat_mult_s as usize / 5 + 16,
            0.0,
        );
        solo_feature
            .count_mat_mult_i
            .resize(solo_feature.n_cb as usize + 1, 0);
    }

    solofeature_collapseumiall_l11_solofeature_collapseumiall(solo_feature, p_solo)?;
    log_main.push_str(&format!(
        "RAM for solo feature {}\n{}",
        solo_feature.feature_type,
        systemfunctions_l6_linuxprocmemory()
    ));
    writeln!(
        log_main,
        "{} ... Finished collapsing UMIs",
        timefunctions_l14_timemonthdaytime(raw_time)
    )
    .unwrap();

    solo_feature.r_cbp.clear();
    Ok(log_main)
}
