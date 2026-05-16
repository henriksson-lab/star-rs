#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `Solo` at STAR/source/Solo.h:11."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Solo {
    pub p_solo: ParametersSolo,
    pub solo_feat: Vec<SoloFeature>,
    pub read_bar_sum: Option<SoloReadBarcode>,
}

#[doc = "Original `Solo::Solo` at STAR/source/Solo.cpp:5. Args: RAchunkIn: ReadAlignChunk, Pin: Parameters, inTrans: Transcriptome"]
pub fn solo_l5_solo_solo(
    p: &crate::parameters_chimeric::Parameters,
    trans: &crate::transcriptome::Transcriptome,
) -> crate::solo::Solo {
    let mut solo = crate::solo::Solo {
        p_solo: p.p_solo.clone(),
        ..Default::default()
    };

    if solo.p_solo.solo_type == SOLO_TYPE_NONE {
        return solo;
    }

    solo.read_bar_sum = Some(soloreadbarcode_l4_soloreadbarcode_soloreadbarcode(
        p.p_solo.solo_type,
        p.p_solo.cb_wl_yes,
        p.p_solo.cb_wl_size,
        p.p_solo.umi_l,
    ));

    if solo.p_solo.solo_type == SOLO_TYPE_CB_SAM_TAG_OUT {
        return solo;
    }

    solo.solo_feat = Vec::with_capacity(solo.p_solo.n_features as usize);
    for ii in 0..solo.p_solo.n_features as usize {
        solo.solo_feat.push(solofeature_l4_solofeature_solofeature(
            p,
            solo.p_solo.features[ii] as i32,
            trans.n_ge,
        ));
    }
    solo
}

#[doc = "Original `Solo::Solo` at STAR/source/Solo.cpp:23. Args: Pin: Parameters, inTrans: Transcriptome"]
pub fn solo_l23_solo_solo(
    p: &crate::parameters_chimeric::Parameters,
    trans: &crate::transcriptome::Transcriptome,
    matrix_contents: &str,
    barcodes_contents: &str,
    features_contents: &str,
    current_dir: &str,
) -> Result<crate::solo_filtered_cells::SoloConstructorCellFilteringResult, String> {
    let mut result = crate::solo_filtered_cells::SoloConstructorCellFilteringResult {
        solo: crate::solo::Solo {
            p_solo: p.p_solo.clone(),
            ..Default::default()
        },
        ..Default::default()
    };

    if p.run_mode_in.first().map(|s| s.as_str()) != Some("soloCellFiltering") {
        return Ok(result);
    }

    result.log_stdout.push_str(&format!(
        "{} ..... starting SoloCellFiltering\n",
        timefunctions_l4_timemonthdaytime()
    ));

    result.solo.solo_feat = vec![solofeature_l4_solofeature_solofeature(p, -1, trans.n_ge)];
    let mut p_solo = result.solo.p_solo.clone();
    solofeature_loadrawmatrix_l7_solofeature_loadrawmatrix(
        &mut result.solo.solo_feat[0],
        p,
        &mut p_solo,
        matrix_contents,
        barcodes_contents,
        features_contents,
    )?;
    result.solo.p_solo = p_solo;

    let cell_filtering = solofeature_cellfiltering_l5_solofeature_cellfiltering(
        &mut result.solo.solo_feat[0],
        &result.solo.p_solo,
        None,
        p,
        trans,
        current_dir,
    )?;
    result.log_main.push_str(&cell_filtering.log_main);
    result.cell_filtering = Some(cell_filtering);

    result.log_stdout.push_str(&format!(
        "{} ..... finished successfully\n",
        timefunctions_l4_timemonthdaytime()
    ));
    result.log_main.push_str("ALL DONE!\n");
    result.exited = true;
    Ok(result)
}

#[doc = "Original `Solo::processAndOutput` at STAR/source/Solo.cpp:48. Args: "]
pub fn solo_l48_solo_processandoutput<FCount, FQuant>(
    solo: &mut crate::solo::Solo,
    p: &mut crate::parameters_chimeric::Parameters,
    trans: &crate::transcriptome::Transcriptome,
    ra_chunks: &mut [crate::read_align_chunk::ReadAlignChunk],
    g_stats_all: &crate::stats::Stats,
    current_dir: &str,
    sj_start_gap_rows: &[(u64, u64)],
    time_counting_start: &str,
    time_counting_finish: &str,
    time_process_start: &str,
    time_writing_raw_matrix: &str,
    time_cell_filtering: &str,
    time_finished_redistribution: &str,
    time_finished_collapsing: &str,
    time_finished_counting: &str,
    linux_proc_memory: &str,
    mut count_cb_gene_umi: FCount,
    mut quant_transcript: FQuant,
) -> Result<crate::solo_filtered_cells::SoloProcessAndOutputResult, String>
where
    FCount: FnMut(usize, &mut crate::solo_feature::SoloFeature) -> String,
    FQuant: FnMut(usize, &mut crate::solo_feature::SoloFeature) -> String,
{
    let mut result = crate::solo_filtered_cells::SoloProcessAndOutputResult::default();
    if solo.p_solo.solo_type == SOLO_TYPE_NONE {
        return Ok(result);
    }

    {
        let read_bar_sum = solo
            .read_bar_sum
            .as_mut()
            .expect("Solo::processAndOutput requires readBarSum");
        if solo.p_solo.cb_wl_yes {
            for ii in 0..p.run_thread_n as usize {
                if let Some(read_bar) = ra_chunks
                    .get_mut(ii)
                    .and_then(|chunk| chunk.ra.solo_read.read_bar.take())
                {
                    soloreadbarcode_l26_soloreadbarcode_addcounts(read_bar_sum, &read_bar);
                    soloreadbarcode_l38_soloreadbarcode_addstats(read_bar_sum, &read_bar);
                }
            }
        }

        result.files.insert(
            format!(
                "{}{}Barcodes.stats",
                p.out_file_name_prefix,
                solo.p_solo
                    .out_file_names
                    .first()
                    .cloned()
                    .unwrap_or_default()
            ),
            soloreadbarcode_l44_soloreadbarcode_statsout(read_bar_sum),
        );

        if solo.p_solo.cb_match_wl.mm1_multi_pc {
            for ii in 0..solo.p_solo.cb_wl_size as usize {
                read_bar_sum.cb_read_count_exact[ii] += 1;
            }
        }
    }

    if solo.p_solo.solo_type == SOLO_TYPE_CB_SAM_TAG_OUT {
        result.returned_after_barcode_output = true;
        return Ok(result);
    }

    result.log_stdout.push_str(&format!(
        "{} ..... started Solo counting\n",
        time_counting_start
    ));
    result.log_main.push_str(&format!(
        "{} ..... started Solo counting\n",
        time_counting_start
    ));

    let ra_reads: Vec<crate::read_align::ReadAlign> =
        ra_chunks.iter().map(|chunk| chunk.ra.clone()).collect();
    let mut p_solo = solo.p_solo.clone();
    let mut read_bar_sum = solo
        .read_bar_sum
        .clone()
        .expect("Solo::processAndOutput requires readBarSum");

    for ii in 0..p_solo.n_features as usize {
        let feature_type = solo.solo_feat[ii].feature_type;
        let read_feat_all: Vec<crate::solo_read_feature::SoloReadFeature> = ra_chunks
            .iter()
            .filter_map(|chunk| {
                chunk
                    .ra
                    .solo_read
                    .read_feat
                    .iter()
                    .find(|rf| rf.feature_type == feature_type)
                    .cloned()
            })
            .collect();
        let gene_solo_feature = solo
            .solo_feat
            .iter()
            .find(|sf| sf.feature_type == SOLO_FEATURE_GENE)
            .cloned();

        let feature_result = solofeature_processrecords_l8_solofeature_processrecords(
            &mut solo.solo_feat[ii],
            p,
            &mut p_solo,
            &mut read_bar_sum,
            &read_feat_all,
            gene_solo_feature.as_ref(),
            trans,
            g_stats_all,
            &ra_reads,
            current_dir,
            sj_start_gap_rows,
            time_process_start,
            time_writing_raw_matrix,
            time_cell_filtering,
            time_finished_redistribution,
            time_finished_collapsing,
            time_finished_counting,
            linux_proc_memory,
            |sf| count_cb_gene_umi(ii, sf),
            |sf| quant_transcript(ii, sf),
        )?;
        result
            .created_directories
            .extend(feature_result.created_directories.clone());
        result.files.extend(feature_result.files.clone());
        result.symlinks.extend(feature_result.symlinks.clone());
        result.log_main.push_str(&feature_result.log_main);
        result.feature_results.push(feature_result);
    }

    solo.p_solo = p_solo;
    solo.read_bar_sum = Some(read_bar_sum);

    result.log_stdout.push_str(&format!(
        "{} ..... finished Solo counting\n",
        time_counting_finish
    ));
    result.log_main.push_str(&format!(
        "{} ..... finished Solo counting\n",
        time_counting_finish
    ));
    Ok(result)
}
