#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::cellFiltering` at STAR/source/SoloFeature_cellFiltering.cpp:5. Args: "]
pub fn solofeature_cellfiltering_l5_solofeature_cellfiltering(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    p_solo: &crate::parameters_solo::ParametersSolo,
    gene_solo_feature: Option<&crate::solo_feature::SoloFeature>,
    p: &crate::parameters_chimeric::Parameters,
    trans: &crate::transcriptome::Transcriptome,
    current_dir: &str,
) -> Result<crate::solo_filtered_cells::SoloCellFilteringResult, String> {
    let mut result = crate::solo_filtered_cells::SoloCellFilteringResult::default();

    if p_solo
        .cell_filter
        .type_
        .first()
        .map(|s| s.as_str())
        .unwrap_or("None")
        == "None"
        || solo_feature.n_cb < 1
    {
        return Ok(result);
    }

    match solo_feature.feature_type {
        SOLO_FEATURE_VELOCYTO => {
            solo_feature.filtered_cells = crate::solo_filtered_cells::SoloFilteredCells {
                filt_vec_bool: vec![false; solo_feature.n_cb as usize],
                ..Default::default()
            };

            if let Some(so_fe_ge) = gene_solo_feature {
                for ic in 0..so_fe_ge.n_cb as usize {
                    let gene_cb = so_fe_ge.ind_cb[ic] as usize;
                    if so_fe_ge.filtered_cells.filt_vec_bool[ic]
                        && solo_feature
                            .ind_cb_wl
                            .get(gene_cb)
                            .copied()
                            .unwrap_or(u32::MAX)
                            != u32::MAX
                    {
                        let velo_cb = solo_feature.ind_cb_wl[gene_cb] as usize;
                        solo_feature.filtered_cells.filt_vec_bool[velo_cb] = true;
                    }
                }
            }

            solo_feature.n_umi_per_cb_sorted = solo_feature.n_umi_per_cb.clone();
            solo_feature
                .n_umi_per_cb_sorted
                .sort_by(|u1, u2| u2.cmp(u1));
        }
        SOLO_FEATURE_GENE
        | SOLO_FEATURE_GENE_FULL
        | SOLO_FEATURE_GENE_FULL_EX50P_AS
        | SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON
        | -1 => {
            solo_feature.n_umi_per_cb_sorted = solo_feature.n_umi_per_cb.clone();
            solo_feature
                .n_umi_per_cb_sorted
                .sort_by(|u1, u2| u2.cmp(u1));

            let mut n_umi_max = 0u32;
            let n_umi_min = if p_solo.cell_filter.type_[0] == "TopCells" {
                let index = std::cmp::min(
                    solo_feature.n_cb.saturating_sub(1),
                    p_solo.cell_filter.top_cells,
                ) as usize;
                solo_feature.n_umi_per_cb_sorted[index]
            } else {
                let maxind = (p_solo.cell_filter.knee.n_expected_cells
                    * (1.0 - p_solo.cell_filter.knee.max_percentile))
                    .round() as u32;
                let index = std::cmp::min(solo_feature.n_cb.saturating_sub(1), maxind) as usize;
                n_umi_max = solo_feature.n_umi_per_cb_sorted[index];
                (n_umi_max as f64 / p_solo.cell_filter.knee.max_min_ratio).round() as u32
            }
            .max(1);

            solo_feature.filtered_cells = crate::solo_filtered_cells::SoloFilteredCells {
                filt_vec_bool: vec![false; solo_feature.n_cb as usize],
                ..Default::default()
            };

            for icb in 0..solo_feature.n_cb as usize {
                if solo_feature.n_umi_per_cb[icb] >= n_umi_min {
                    solo_feature.filtered_cells.filt_vec_bool[icb] = true;
                    solo_feature.filtered_cells.n_cells_simple += 1;
                }
            }

            result.log_main.push_str(&format!(
                "cellFiltering: simple: nUMImax={}; nUMImin={}; nCellsSimple={}\n",
                n_umi_max, n_umi_min, solo_feature.filtered_cells.n_cells_simple
            ));

            if p_solo.cell_filter.type_[0] == "EmptyDrops_CR" {
                result.empty_drops_requested = true;
                let empty_drops =
                    solofeature_emptydrops_cr_l10_solofeature_emptydrops_cr(solo_feature, p_solo);
                result.log_main.push_str(&empty_drops.log_main);
                result.empty_drops = Some(empty_drops);
            }
        }
        _ => return Ok(result),
    }

    let mut gene_detected = vec![0u32; solo_feature.features_number.max(0) as usize];
    for icb in 0..solo_feature.n_cb as usize {
        if solo_feature.filtered_cells.filt_vec_bool[icb] {
            solo_feature.filtered_cells.n_cells += 1;
            solo_feature.filtered_cells.n_umi_in_cells += solo_feature.n_umi_per_cb[icb] as u64;

            if !solo_feature.n_read_per_cb_unique.is_empty() {
                solo_feature.filtered_cells.n_read_in_cells +=
                    solo_feature.n_read_per_cb_total[icb] as u64;
                solo_feature.filtered_cells.n_read_in_cells_unique +=
                    solo_feature.n_read_per_cb_unique[icb] as u64;
                solo_feature
                    .filtered_cells
                    .n_read_per_cell_unique
                    .push(solo_feature.n_read_per_cb_unique[icb]);
            }

            let mut ng1 = 0u32;
            for ig in 0..solo_feature.n_gene_per_cb[icb] as usize {
                let ind_g1 = solo_feature.count_cell_gene_umi_index[icb] as usize
                    + ig * solo_feature.count_mat_stride as usize;
                if solo_feature.count_cell_gene_umi
                    [ind_g1 + p_solo.umi_dedup.count_ind_main as usize]
                    > 0
                {
                    let gene = solo_feature.count_cell_gene_umi[ind_g1] as usize;
                    gene_detected[gene] = 1;
                    ng1 += 1;
                }
            }
            solo_feature.filtered_cells.n_gene_in_cells += ng1 as u64;
            solo_feature.filtered_cells.n_gene_per_cell.push(ng1);
        }
    }

    if solo_feature.filtered_cells.n_cells == 0 {
        return Ok(result);
    }

    solo_feature.filtered_cells.n_gene_detected =
        gene_detected.iter().filter(|&&ii| ii > 0).count() as u64;
    solo_feature.filtered_cells.mean_umi_per_cell =
        solo_feature.filtered_cells.n_umi_in_cells / solo_feature.filtered_cells.n_cells as u64;
    solo_feature.filtered_cells.mean_read_per_cell_unique =
        solo_feature.filtered_cells.n_read_in_cells_unique
            / solo_feature.filtered_cells.n_cells as u64;
    solo_feature.filtered_cells.mean_gene_per_cell =
        solo_feature.filtered_cells.n_gene_in_cells / solo_feature.filtered_cells.n_cells as u64;

    solo_feature.filtered_cells.n_read_per_cell_unique.sort();
    solo_feature.filtered_cells.n_gene_per_cell.sort();

    let median_index = solo_feature.filtered_cells.n_cells as usize / 2;
    solo_feature.filtered_cells.median_umi_per_cell =
        solo_feature.n_umi_per_cb_sorted[median_index] as u64;
    solo_feature.filtered_cells.median_gene_per_cell =
        solo_feature.filtered_cells.n_gene_per_cell[median_index] as u64;
    solo_feature.filtered_cells.median_read_per_cell_unique = solo_feature
        .filtered_cells
        .n_read_per_cell_unique
        .get(median_index)
        .copied()
        .unwrap_or(0) as u64;

    result.output_results = Some(solofeature_outputresults_l12_solofeature_outputresults(
        solo_feature,
        true,
        &solo_feature.output_prefix_filtered.clone(),
        p,
        p_solo,
        trans,
        current_dir,
    )?);
    Ok(result)
}
