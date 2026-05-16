#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::loadRawMatrix` at STAR/source/SoloFeature_loadRawMatrix.cpp:7. Args: "]
pub fn solofeature_loadrawmatrix_l7_solofeature_loadrawmatrix(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    p: &crate::parameters_chimeric::Parameters,
    p_solo: &mut crate::parameters_solo::ParametersSolo,
    matrix_contents: &str,
    barcodes_contents: &str,
    features_contents: &str,
) -> Result<(), String> {
    if p.run_mode_in.len() < 3 {
        return Err("Exiting because of fatal PARAMETER error: --runMode soloCellFiltering should contain paths to count matrix input directorry and output prefix.\nSOLUTION: re-run with --runMode soloCellFiltering </path/to/raw/count/dir/> </path/to/output/prefix>\n".to_string());
    }

    let input_prefix = format!("{}/", p.run_mode_in[1]);
    solo_feature.output_prefix = p.run_mode_in[2].clone();
    solo_feature.output_prefix_filtered = solo_feature.output_prefix.clone();

    let matrix_file_name = format!("{}{}", input_prefix, p_solo.out_file_names[3]);
    let mut mat_lines = matrix_contents.lines();
    let mut data_tokens = Vec::new();
    for line in mat_lines.by_ref() {
        if line.starts_with('%') {
            continue;
        }
        data_tokens.extend(line.split_whitespace().map(|s| s.to_string()));
        break;
    }
    for line in mat_lines {
        data_tokens.extend(line.split_whitespace().map(|s| s.to_string()));
    }

    let mut it = data_tokens.iter();
    solo_feature.features_number = it
        .next()
        .ok_or_else(|| format!("missing feature count in {}", matrix_file_name))?
        .parse::<i64>()
        .map_err(|e| e.to_string())?;
    let n_cb1 = it
        .next()
        .ok_or_else(|| format!("missing barcode count in {}", matrix_file_name))?
        .parse::<u32>()
        .map_err(|e| e.to_string())?;
    let n_tot = it
        .next()
        .ok_or_else(|| format!("missing entry count in {}", matrix_file_name))?
        .parse::<u64>()
        .map_err(|e| e.to_string())?;

    if n_tot == 0 {
        return Err(format!(
            "Exiting because of fatal INPUT FILE error: no counts detected in {}\nSOLUTION: check the formatting of the matrix file.\n",
            matrix_file_name
        ));
    }

    solo_feature.count_mat_stride = 3;
    let stride = solo_feature.count_mat_stride as usize;
    solo_feature
        .count_cell_gene_umi
        .resize(n_tot as usize * stride, 0);

    for ii in 0..n_tot as usize {
        let gene = it
            .next()
            .ok_or_else(|| format!("missing gene entry {} in {}", ii, matrix_file_name))?
            .parse::<u32>()
            .map_err(|e| e.to_string())?;
        let cell = it
            .next()
            .ok_or_else(|| format!("missing cell entry {} in {}", ii, matrix_file_name))?
            .parse::<u32>()
            .map_err(|e| e.to_string())?;
        let count1 = it
            .next()
            .ok_or_else(|| format!("missing count entry {} in {}", ii, matrix_file_name))?
            .parse::<f64>()
            .map_err(|e| e.to_string())?;

        solo_feature.count_cell_gene_umi[ii * stride] = gene - 1;
        solo_feature.count_cell_gene_umi[ii * stride + 1] = cell - 1;
        solo_feature.count_cell_gene_umi[ii * stride + 2] = count1.round() as u32;
    }

    let mut entries = Vec::with_capacity(n_tot as usize);
    for ii in 0..n_tot as usize {
        entries.push([
            solo_feature.count_cell_gene_umi[ii * stride],
            solo_feature.count_cell_gene_umi[ii * stride + 1],
            solo_feature.count_cell_gene_umi[ii * stride + 2],
        ]);
    }
    entries.sort_by(|a, b| a[1].cmp(&b[1]).then(a[0].cmp(&b[0])));
    for (ii, entry) in entries.iter().enumerate() {
        solo_feature.count_cell_gene_umi[ii * stride..ii * stride + stride].copy_from_slice(entry);
    }

    solo_feature.n_cb = 0;
    let mut ci_prev = u32::MAX;
    for ii in 0..n_tot as usize {
        let ci1 = solo_feature.count_cell_gene_umi[ii * stride + 1];
        if ci1 != ci_prev {
            ci_prev = ci1;
            solo_feature.n_cb += 1;
        }
    }

    solo_feature.ind_cb.resize(solo_feature.n_cb as usize, 0);
    solo_feature
        .count_cell_gene_umi_index
        .resize(solo_feature.n_cb as usize, 0);
    solo_feature
        .n_umi_per_cb
        .resize(solo_feature.n_cb as usize, 0);
    solo_feature
        .n_gene_per_cb
        .resize(solo_feature.n_cb as usize, 0);
    solo_feature
        .n_read_per_cb
        .resize(solo_feature.n_cb as usize, 0);

    solo_feature.n_cb = u32::MAX;
    ci_prev = u32::MAX;
    for ii in 0..n_tot as usize {
        let ci1 = solo_feature.count_cell_gene_umi[ii * stride + 1];
        if ci1 != ci_prev {
            ci_prev = ci1;
            solo_feature.n_cb = solo_feature.n_cb.wrapping_add(1);
            solo_feature.ind_cb[solo_feature.n_cb as usize] = ci1;
            solo_feature.count_cell_gene_umi_index[solo_feature.n_cb as usize] =
                (ii * stride) as u32;
        }
        let cb = solo_feature.n_cb as usize;
        solo_feature.n_gene_per_cb[cb] += 1;
        solo_feature.n_umi_per_cb[cb] += solo_feature.count_cell_gene_umi[ii * stride + 2];
        solo_feature.count_cell_gene_umi[ii * stride + 1] =
            solo_feature.count_cell_gene_umi[ii * stride + 2];
    }

    p_solo.cb_wl_str = barcodes_contents
        .lines()
        .take(n_cb1 as usize)
        .map(|s| s.to_string())
        .collect();
    p_solo.cb_wl_str.resize(n_cb1 as usize, String::new());
    solo_feature.copied_features_tsv = features_contents.to_string();

    Ok(())
}
