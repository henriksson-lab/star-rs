#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::redistributeReadsByCB` at STAR/source/SoloFeature_redistributeReadsByCB.cpp:8. Args: "]
pub fn solofeature_redistributereadsbycb_l8_solofeature_redistributereadsbycb(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    p_solo: &crate::parameters_solo::ParametersSolo,
    run_thread_n: i32,
) -> String {
    let read_feat_sum = solo_feature
        .read_feat_sum
        .as_ref()
        .expect("SoloFeature::redistributeReadsByCB requires readFeatSum");
    let n_read_rec = read_feat_sum
        .cb_read_count
        .iter()
        .fold(0_u64, |acc, x| acc + *x as u64);
    let n_read_rec_bin = n_read_rec / p_solo.redistr_reads_nfiles as u64;
    let log_main = format!(
        "     Redistributing reads into {}files; nReadRec={};   nReadRecBin={}\n",
        p_solo.redistr_reads_nfiles, n_read_rec, n_read_rec_bin
    );

    solo_feature.redistr_files_cb_first.clear();
    solo_feature.redistr_files_cb_first.push(0);
    solo_feature
        .redistr_files_cb_index
        .resize(solo_feature.n_cb as usize, 0);
    solo_feature.redistr_files_nreads.clear();

    let mut nreads = 0_u64;
    let mut ind = 0_u32;
    for icb in 0..solo_feature.n_cb as usize {
        solo_feature.redistr_files_cb_index[icb] = ind;
        nreads += read_feat_sum.cb_read_count[solo_feature.ind_cb[icb] as usize] as u64;
        if nreads >= n_read_rec_bin {
            ind += 1;
            solo_feature.redistr_files_cb_first.push(icb as u32 + 1);
            solo_feature.redistr_files_nreads.push(nreads);
            nreads = 0;
        }
    }
    if nreads > 0 {
        solo_feature.redistr_files_cb_first.push(solo_feature.n_cb);
        solo_feature.redistr_files_nreads.push(nreads);
    }

    solo_feature
        .redistr_files_streams
        .resize(solo_feature.redistr_files_nreads.len(), String::new());

    for ii in 0..run_thread_n as usize {
        for line1 in solo_feature.read_feat_all[ii].stream_reads.lines() {
            if line1.is_empty() {
                break;
            }

            let mut line1stream = line1.split_whitespace();
            let _umi = line1stream.next();
            let mut cb1 = line1stream
                .next()
                .and_then(|x| x.parse::<u64>().ok())
                .unwrap();
            cb1 = line1stream
                .next()
                .and_then(|x| x.parse::<u64>().ok())
                .unwrap_or(cb1);
            if solo_feature.feature_type == SOLO_FEATURE_SJ {
                cb1 = line1stream
                    .next()
                    .and_then(|x| x.parse::<u64>().ok())
                    .unwrap_or(cb1);
            }
            cb1 = line1stream
                .next()
                .and_then(|x| x.parse::<u64>().ok())
                .unwrap_or(cb1);

            let cb_compact = solo_feature.ind_cb_wl[cb1 as usize] as usize;
            let file_index = solo_feature.redistr_files_cb_index[cb_compact] as usize;
            solo_feature.redistr_files_streams[file_index].push_str(line1);
            solo_feature.redistr_files_streams[file_index].push('\n');
        }
    }

    log_main
}
