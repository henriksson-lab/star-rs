#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::processRecords` at STAR/source/SoloFeature_processRecords.cpp:8. Args: "]
pub fn solofeature_processrecords_l8_solofeature_processrecords<FCount, FQuant>(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    p: &mut crate::parameters_chimeric::Parameters,
    p_solo: &mut crate::parameters_solo::ParametersSolo,
    read_bar_sum: &mut crate::solo_read_barcode::SoloReadBarcode,
    read_feat_all: &[crate::solo_read_feature::SoloReadFeature],
    gene_solo_feature: Option<&crate::solo_feature::SoloFeature>,
    trans: &crate::transcriptome::Transcriptome,
    g_stats_all: &crate::stats::Stats,
    ra_chunks: &[crate::read_align::ReadAlign],
    current_dir: &str,
    sj_start_gap_rows: &[(u64, u64)],
    time_start: &str,
    time_writing_raw_matrix: &str,
    time_cell_filtering: &str,
    time_finished_redistribution: &str,
    time_finished_collapsing: &str,
    time_finished_counting: &str,
    linux_proc_memory: &str,
    mut count_cb_gene_umi: FCount,
    mut quant_transcript: FQuant,
) -> Result<crate::solo_filtered_cells::SoloFeatureProcessRecordsResult, String>
where
    FCount: FnMut(&mut crate::solo_feature::SoloFeature) -> String,
    FQuant: FnMut(&mut crate::solo_feature::SoloFeature) -> String,
{
    let mut result = crate::solo_filtered_cells::SoloFeatureProcessRecordsResult::default();
    if p_solo.solo_type == SOLO_TYPE_NONE {
        return Ok(result);
    }

    let feature_names = [
        "SJ",
        "Transcript3p",
        "GeneFull",
        "GeneFull_ExonOverIntron",
        "GeneFull_Ex50pAS",
        "Gene",
        "VelocytoSimple",
        "Velocyto",
    ];
    let feature_name = if solo_feature.feature_type >= 0 {
        feature_names
            .get(solo_feature.feature_type as usize)
            .copied()
            .unwrap_or("")
    } else {
        ""
    };

    result.log_main.push_str(&format!(
        "{} ... Starting Solo post-map for {}\n",
        time_start, feature_name
    ));

    solo_feature.output_prefix = format!(
        "{}{}{}{}/",
        p.out_file_name_prefix,
        p_solo.out_file_names.first().cloned().unwrap_or_default(),
        feature_name,
        ""
    );
    solo_feature.output_prefix_filtered = format!("{}filtered/", solo_feature.output_prefix);
    result
        .created_directories
        .push(solo_feature.output_prefix.clone());

    if solo_feature.feature_type == SOLO_FEATURE_SJ && p.sj_all[0].is_empty() {
        p.sj_all[0].reserve(10_000_000);
        p.sj_all[1].reserve(10_000_000);
        for (start, gap) in sj_start_gap_rows {
            p.sj_all[0].push(*start);
            p.sj_all[1].push(*gap);
        }
        result.log_main.push_str(&format!(
            "Read splice junctions for Solo SJ feature: {}\n",
            p.sj_all[0].len()
        ));
    }

    solofeature_sumthreads_l8_solofeature_sumthreads(
        solo_feature,
        p,
        p_solo,
        read_bar_sum,
        read_feat_all,
        g_stats_all.read_n as u64,
    )?;

    if solo_feature.feature_type == SOLO_FEATURE_VELOCYTO {
        let gene_read_info = gene_solo_feature
            .map(|sf| sf.read_info.as_slice())
            .unwrap_or(&[]);
        result
            .log_main
            .push_str(&solofeature_countvelocyto_l12_solofeature_countvelocyto(
                solo_feature,
                gene_read_info,
                trans,
                p.run_thread_n,
                time_start,
                time_finished_counting,
                time_finished_collapsing,
                linux_proc_memory,
            )?);
    } else if solo_feature.feature_type == SOLO_FEATURE_TRANSCRIPT3P {
        result.quant_transcript_called = true;
        result.log_main.push_str(&quant_transcript(solo_feature));
        result.returned_after_quant_transcript = true;
        return Ok(result);
    } else if p_solo.solo_type == SOLO_TYPE_SMART_SEQ {
        result
            .log_main
            .push_str(&solofeature_countsmartseq_l9_solofeature_countsmartseq(
                solo_feature,
                p,
                time_finished_redistribution,
                time_finished_collapsing,
                time_finished_counting,
            )?);
    } else {
        result.count_cb_gene_umi_called = true;
        result.log_main.push_str(&count_cb_gene_umi(solo_feature));
    }

    let read_feat_sum = solo_feature
        .read_feat_sum
        .as_ref()
        .expect("SoloFeature::processRecords requires readFeatSum");
    result.files.insert(
        format!("{}Features.stats", solo_feature.output_prefix),
        soloreadfeature_l56_soloreadfeature_statsout(read_feat_sum),
    );

    result.log_main.push_str(&format!(
        "{} ... Solo: writing raw matrix\n",
        time_writing_raw_matrix
    ));
    let raw = solofeature_outputresults_l12_solofeature_outputresults(
        solo_feature,
        false,
        &format!("{}/raw/", solo_feature.output_prefix),
        p,
        p_solo,
        trans,
        current_dir,
    )?;
    result.created_directories.push(raw.created_directory);
    result.files.extend(raw.files);
    result.symlinks.extend(raw.symlinks);

    result.log_main.push_str(&format!(
        "{} ... Solo: cell filtering\n",
        time_cell_filtering
    ));
    let cell_filtering = solofeature_cellfiltering_l5_solofeature_cellfiltering(
        solo_feature,
        p_solo,
        gene_solo_feature,
        p,
        trans,
        current_dir,
    )?;
    if let Some(out) = &cell_filtering.output_results {
        result
            .created_directories
            .push(out.created_directory.clone());
        result.files.extend(out.files.clone());
        result.symlinks.extend(out.symlinks.clone());
    }
    result.log_main.push_str(&cell_filtering.log_main);
    result.cell_filtering = Some(cell_filtering);

    let stats_output = solofeature_statsoutput_l6_solofeature_statsoutput(
        solo_feature,
        p,
        p_solo,
        g_stats_all,
        ra_chunks,
    );
    result.files.extend(stats_output.files.clone());
    result.stats_output = Some(stats_output);

    solofeature_l29_solofeature_clearlarge(solo_feature);
    result.log_main.push_str("RAM after completing solo:\n");
    result.log_main.push_str(linux_proc_memory);
    Ok(result)
}
