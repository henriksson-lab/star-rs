#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::sumThreads` at STAR/source/SoloFeature_sumThreads.cpp:8. Args: "]
pub fn solofeature_sumthreads_l8_solofeature_sumthreads(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    p: &crate::parameters_chimeric::Parameters,
    p_solo: &mut crate::parameters_solo::ParametersSolo,
    read_bar_sum: &mut crate::solo_read_barcode::SoloReadBarcode,
    read_feat_all: &[crate::solo_read_feature::SoloReadFeature],
    stats_all_read_n: u64,
) -> Result<(), String> {
    solo_feature.n_reads_input = stats_all_read_n + 1;
    solo_feature.read_feat_all_len = read_feat_all.len();

    let read_feat_sum = solo_feature
        .read_feat_sum
        .as_mut()
        .expect("SoloFeature::sumThreads requires readFeatSum");

    for read_feat in read_feat_all {
        soloreadfeature_l29_soloreadfeature_addcounts(read_feat_sum, read_feat);
    }

    if !p_solo.cb_wl_yes {
        p_solo.cb_wl_size = read_feat_sum.cb_read_count_map.len() as u32;
        p_solo.cb_wl.resize(p_solo.cb_wl_size as usize, 0);
        p_solo
            .cb_wl_str
            .resize(p_solo.cb_wl_size as usize, String::new());
        read_feat_sum
            .cb_read_count
            .resize(p_solo.cb_wl_size as usize, 0);
        read_bar_sum
            .cb_read_count_exact
            .resize(p_solo.cb_wl_size as usize, 0);

        if p_solo.cb_type_type == 1 {
            let mut icb = 0usize;
            for (cb, count) in read_feat_sum.cb_read_count_map.iter() {
                p_solo.cb_wl[icb] = *cb;
                p_solo.cb_wl_str[icb] =
                    sequencefuns_l267_convertnuclint64tostring(p_solo.cb_wl[icb], p_solo.cb_l);
                read_feat_sum.cb_read_count[icb] = *count;
                read_bar_sum.cb_read_count_exact[icb] = *count;
                icb += 1;
            }
        } else if p_solo.cb_type_type == 2 {
            let mut cbiter = vec![String::new(); p_solo.cb_type_str_map.len()];
            for (cb_string, cb_index) in p_solo.cb_type_str_map.iter() {
                cbiter[*cb_index as usize] = cb_string.clone();
            }

            let mut icb = 0usize;
            for (cb, count) in read_feat_sum.cb_read_count_map.iter() {
                p_solo.cb_wl[icb] = *cb;
                p_solo.cb_wl_str[icb] = cbiter[*cb as usize].clone();
                read_feat_sum.cb_read_count[icb] = *count;
                read_bar_sum.cb_read_count_exact[icb] = *count;
                icb += 1;
            }
        }

        if p_solo.cb_match_wl.mm1_multi_pc {
            for ii in 0..p_solo.cb_wl_size as usize {
                read_bar_sum.cb_read_count_exact[ii] += 1;
            }
        }
    }

    if p.run_restart_type == 1 {
        for (read_feat_i, read_feat) in read_feat_all.iter().enumerate() {
            for (line_i, line1) in read_feat.stream_reads.lines().enumerate() {
                let mut line1stream = line1.split_whitespace();
                for field_i in 0..3 {
                    let token = line1stream.next().ok_or_else(|| {
                        format!(
                            "Malformed STARsolo restart record: missing field {field_i} in read feature {read_feat_i} line {}",
                            line_i + 1
                        )
                    })?;
                    token.parse::<usize>().map_err(|_| {
                        format!(
                            "Malformed STARsolo restart record: invalid field {field_i} value {token} in read feature {read_feat_i} line {}",
                            line_i + 1
                        )
                    })?;
                }
                if solo_feature.feature_type == SOLO_FEATURE_SJ {
                    let token = line1stream.next().ok_or_else(|| {
                        format!(
                            "Malformed STARsolo restart record: missing SJ field in read feature {read_feat_i} line {}",
                            line_i + 1
                        )
                    })?;
                    token.parse::<usize>().map_err(|_| {
                        format!(
                            "Malformed STARsolo restart record: invalid SJ field value {token} in read feature {read_feat_i} line {}",
                            line_i + 1
                        )
                    })?;
                }
                let cb_token = line1stream.next().ok_or_else(|| {
                    format!(
                        "Malformed STARsolo restart record: missing cell barcode index in read feature {read_feat_i} line {}",
                        line_i + 1
                    )
                })?;
                let cb1 = cb_token.parse::<usize>().map_err(|_| {
                    format!(
                        "Malformed STARsolo restart record: invalid cell barcode index {cb_token} in read feature {read_feat_i} line {}",
                        line_i + 1
                    )
                })?;
                if cb1 >= read_feat_sum.cb_read_count.len() {
                    return Err(format!(
                        "Malformed STARsolo restart record: cell barcode index {cb1} is outside whitelist"
                    ));
                }
                read_feat_sum.cb_read_count[cb1] += 1;
            }
        }
    }

    solo_feature.n_cb = 0;
    solo_feature.n_reads_mapped = 0;
    for ii in 0..p_solo.cb_wl_size as usize {
        if read_feat_sum.cb_read_count[ii] > 0 {
            solo_feature.n_cb += 1;
            solo_feature.n_reads_mapped += read_feat_sum.cb_read_count[ii] as u64;
        }
    }

    solo_feature
        .ind_cb_wl
        .resize(p_solo.cb_wl_size as usize, u32::MAX);
    solo_feature.ind_cb.resize(solo_feature.n_cb as usize, 0);
    solo_feature.n_cb = 0;
    for ii in 0..p_solo.cb_wl_size as usize {
        if read_feat_sum.cb_read_count[ii] > 0 {
            solo_feature.ind_cb[solo_feature.n_cb as usize] = ii as u32;
            solo_feature.ind_cb_wl[ii] = solo_feature.n_cb;
            solo_feature.n_cb += 1;
        }
    }
    Ok(())
}
