#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::outputResults` at STAR/source/SoloFeature_outputResults.cpp:12. Args: cellFilterYes: bool, outputPrefixMat: string"]
pub fn solofeature_outputresults_l12_solofeature_outputresults(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    cell_filter_yes: bool,
    output_prefix_mat: &str,
    p: &crate::parameters_chimeric::Parameters,
    p_solo: &crate::parameters_solo::ParametersSolo,
    trans: &crate::transcriptome::Transcriptome,
    current_dir: &str,
) -> Result<crate::solo_filtered_cells::SoloOutputResults, String> {
    let mut out = crate::solo_filtered_cells::SoloOutputResults {
        created_directory: output_prefix_mat.to_string(),
        ..Default::default()
    };

    match solo_feature.feature_type {
        SOLO_FEATURE_GENE
        | SOLO_FEATURE_GENE_FULL
        | SOLO_FEATURE_GENE_FULL_EX50P_AS
        | SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON
        | SOLO_FEATURE_VELOCYTO
        | SOLO_FEATURE_VELOCYTO_SIMPLE => {
            let mut gene_str = String::new();
            for ii in 0..trans.n_ge as usize {
                gene_str.push_str(&trans.ge_id[ii]);
                gene_str.push('\t');
                if trans.ge_name[ii].is_empty() {
                    gene_str.push_str(&trans.ge_id[ii]);
                } else {
                    gene_str.push_str(&trans.ge_name[ii]);
                }
                if p_solo.out_format_features_gene_field3 != "-" {
                    gene_str.push('\t');
                    gene_str.push_str(&p_solo.out_format_features_gene_field3);
                }
                gene_str.push('\n');
            }
            out.files.insert(
                format!("{}{}", output_prefix_mat, p_solo.out_file_names[1]),
                gene_str,
            );
        }
        SOLO_FEATURE_SJ => {
            let sjout = if p.out_file_name_prefix.starts_with('/') {
                format!("{}SJ.out.tab", p.out_file_name_prefix)
            } else {
                format!("{}/{}SJ.out.tab", current_dir, p.out_file_name_prefix)
            };
            out.symlinks.push((
                sjout,
                format!("{}{}", output_prefix_mat, p_solo.out_file_names[1]),
            ));
        }
        _ => {}
    }

    let mut cb_str = String::new();
    let mut n_cell_gene_entries = 0u64;
    if cell_filter_yes {
        for icb in 0..solo_feature.n_cb as usize {
            if solo_feature.filtered_cells.filt_vec_bool[icb] {
                cb_str.push_str(&p_solo.cb_wl_str[solo_feature.ind_cb[icb] as usize]);
                cb_str.push('\n');
                n_cell_gene_entries += solo_feature.n_gene_per_cb[icb] as u64;
            }
        }
    } else {
        for ii in 0..p_solo.cb_wl_size as usize {
            cb_str.push_str(&p_solo.cb_wl_str[ii]);
            cb_str.push('\n');
        }
        for icb in 0..solo_feature.n_cb as usize {
            n_cell_gene_entries += solo_feature.n_gene_per_cb[icb] as u64;
        }
    }
    out.files.insert(
        format!("{}{}", output_prefix_mat, p_solo.out_file_names[2]),
        cb_str,
    );

    let velo_names = ["spliced.mtx", "unspliced.mtx", "ambiguous.mtx"];
    let umi_type_names = [
        "NoDedup",
        "Exact",
        "1MM_All",
        "1MM_Directional",
        "1MM_CR",
        "1MM_Directional_UMItools",
    ];
    for i_col in 1..solo_feature.count_mat_stride as usize {
        let matrix_file_name = if solo_feature.feature_type == SOLO_FEATURE_VELOCYTO {
            format!("{}{}", output_prefix_mat, velo_names[i_col - 1])
        } else if i_col > 1 && cell_filter_yes {
            break;
        } else if p_solo.umi_dedup.types.len() > 1 {
            format!(
                "{}umiDedup-{}.mtx",
                output_prefix_mat,
                umi_type_names[p_solo.umi_dedup.types[i_col - 1] as usize]
            )
        } else {
            format!("{}{}", output_prefix_mat, p_solo.out_file_names[3])
        };

        let mut count_matrix_stream = String::new();
        count_matrix_stream.push_str("%%MatrixMarket matrix coordinate integer general\n");
        count_matrix_stream.push_str("%\n");
        count_matrix_stream.push_str(&format!(
            "{} {} {}\n",
            solo_feature.features_number,
            if cell_filter_yes {
                solo_feature.filtered_cells.n_cells
            } else {
                p_solo.cb_wl_size
            },
            n_cell_gene_entries
        ));

        let mut cb_ind1 = 0u32;
        for icb in 0..solo_feature.n_cb as usize {
            if cell_filter_yes {
                if solo_feature.filtered_cells.filt_vec_bool[icb] {
                    cb_ind1 += 1;
                } else {
                    continue;
                }
            } else {
                cb_ind1 = solo_feature.ind_cb[icb] + 1;
            }
            for ig in 0..solo_feature.n_gene_per_cb[icb] as usize {
                let ind_g1 = solo_feature.count_cell_gene_umi_index[icb] as usize
                    + ig * solo_feature.count_mat_stride as usize;
                count_matrix_stream.push_str(&format!(
                    "{} {} {}\n",
                    solo_feature.count_cell_gene_umi[ind_g1] + 1,
                    cb_ind1,
                    solo_feature.count_cell_gene_umi[ind_g1 + i_col]
                ));
            }
        }
        out.files.insert(matrix_file_name, count_matrix_stream);
    }

    if p_solo.multi_map.yes_multi
        && !cell_filter_yes
        && (solo_feature.feature_type == SOLO_FEATURE_GENE
            || solo_feature.feature_type == SOLO_FEATURE_GENE_FULL
            || solo_feature.feature_type == SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON
            || solo_feature.feature_type == SOLO_FEATURE_GENE_FULL_EX50P_AS)
    {
        let multi_type_names = ["Unique", "Uniform", "Rescue", "PropUnique", "EM"];
        solo_feature
            .n_umi_per_cb_multi
            .resize(solo_feature.n_cb as usize, 0);
        solo_feature
            .n_gene_per_cb_multi
            .resize(solo_feature.n_cb as usize, 0);
        let mut n_gene_umi_per_cb_multi_fill = true;

        for &i_mult in p_solo.multi_map.types.iter() {
            for i_ded in 0..p_solo.umi_dedup.yes_n as usize {
                let mut matrix_file_name = format!(
                    "{}UniqueAndMult-{}",
                    output_prefix_mat, multi_type_names[i_mult as usize]
                );
                if p_solo.umi_dedup.types.len() > 1 {
                    matrix_file_name.push_str(&format!(
                        "_umiDedup-{}",
                        umi_type_names[p_solo.umi_dedup.types[i_ded] as usize]
                    ));
                }
                matrix_file_name.push_str(".mtx");

                let m_index = p_solo.multi_map.count_ind_i[i_mult as usize] as usize + i_ded;
                let mut mat_out_string_stream = String::new();
                n_cell_gene_entries = 0;

                for icb in 0..solo_feature.n_cb as usize {
                    let cb_ind1 = solo_feature.ind_cb[icb] + 1;
                    let mut igm1 = solo_feature.count_cell_gene_umi_index[icb] as usize;
                    let mut igm2 = solo_feature.count_mat_mult_i[icb] as usize;
                    while igm1 < solo_feature.count_cell_gene_umi_index[icb + 1] as usize
                        || igm2 < solo_feature.count_mat_mult_i[icb + 1] as usize
                    {
                        let (g1, c1) =
                            if igm1 < solo_feature.count_cell_gene_umi_index[icb + 1] as usize {
                                (
                                    solo_feature.count_cell_gene_umi[igm1],
                                    solo_feature.count_cell_gene_umi[igm1 + 1 + i_ded],
                                )
                            } else {
                                (u32::MAX, 0)
                            };
                        let (g2, c2) = if igm2 < solo_feature.count_mat_mult_i[icb + 1] as usize {
                            (
                                solo_feature.count_mat_mult_m[igm2] as u32,
                                solo_feature.count_mat_mult_m[igm2 + m_index],
                            )
                        } else {
                            (u32::MAX, 0.0)
                        };

                        if g1 < g2 {
                            mat_out_string_stream.push_str(&format!(
                                "{} {} {}\n",
                                g1 + 1,
                                cb_ind1,
                                c1
                            ));
                            igm1 += solo_feature.count_mat_stride as usize;
                        } else if g1 > g2 {
                            mat_out_string_stream.push_str(&format!(
                                "{} {} {}\n",
                                g2 + 1,
                                cb_ind1,
                                c2
                            ));
                            igm2 += solo_feature.count_mat_mult_s as usize;
                            if n_gene_umi_per_cb_multi_fill {
                                solo_feature.n_umi_per_cb_multi[icb] += c2 as u32;
                                solo_feature.n_gene_per_cb_multi[icb] += 1;
                            }
                        } else {
                            mat_out_string_stream.push_str(&format!(
                                "{} {} {}\n",
                                g1 + 1,
                                cb_ind1,
                                c1 as f64 + c2
                            ));
                            igm1 += solo_feature.count_mat_stride as usize;
                            igm2 += solo_feature.count_mat_mult_s as usize;
                            if n_gene_umi_per_cb_multi_fill {
                                solo_feature.n_umi_per_cb_multi[icb] += c2 as u32;
                            }
                        }
                        n_cell_gene_entries += 1;
                    }
                }
                n_gene_umi_per_cb_multi_fill = false;

                let mut count_matrix_stream = String::new();
                count_matrix_stream.push_str("%%MatrixMarket matrix coordinate real general\n");
                count_matrix_stream.push_str("%\n");
                count_matrix_stream.push_str(&format!(
                    "{} {} {}\n",
                    solo_feature.features_number, p_solo.cb_wl_size, n_cell_gene_entries
                ));
                count_matrix_stream.push_str(&mat_out_string_stream);
                out.files.insert(matrix_file_name, count_matrix_stream);
            }
        }
    }

    Ok(out)
}
