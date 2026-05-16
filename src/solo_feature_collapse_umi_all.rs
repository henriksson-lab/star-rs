#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::collapseUMIall` at STAR/source/SoloFeature_collapseUMIall.cpp:11. Args: "]
pub fn solofeature_collapseumiall_l11_solofeature_collapseumiall(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    p_solo: &crate::parameters_solo::ParametersSolo,
) -> Result<(), String> {
    let n_read_per_cb_max = solo_feature
        .n_read_per_cb
        .iter()
        .copied()
        .max()
        .unwrap_or(0);
    let mut umi_array = vec![0_u32; n_read_per_cb_max as usize * 3];
    let gene_slots = std::cmp::min(
        2_u32.saturating_mul(solo_feature.features_number.max(0) as u32),
        n_read_per_cb_max,
    )
    .saturating_add(1) as usize;
    let mut g_id = vec![0_u32; gene_slots];
    let mut g_read_s = vec![0_u32; gene_slots];

    for icb in 0..solo_feature.n_cb {
        solofeature_collapseumiall_l30_solofeature_collapseumipercb(
            solo_feature,
            p_solo,
            icb,
            &mut umi_array,
            &mut g_id,
            &mut g_read_s,
        )?;

        let read_feat_sum = solo_feature
            .read_feat_sum
            .as_mut()
            .ok_or_else(|| "SoloFeature::collapseUMIall requires readFeatSum".to_string())?;
        if read_feat_sum.stats.v.len() > SOLO_READ_FEATURE_STAT_YES_UMIS as usize {
            read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_UMIS as usize] +=
                solo_feature.n_umi_per_cb[icb as usize] as u64;
        }
        if solo_feature.n_gene_per_cb[icb as usize] > 0
            && read_feat_sum.stats.v.len() > SOLO_READ_FEATURE_STAT_YES_CELL_BARCODES as usize
        {
            read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_CELL_BARCODES as usize] += 1;
        }
        if read_feat_sum.stats.v.len() > SOLO_READ_FEATURE_STAT_YES_WL_MATCH as usize {
            read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_WL_MATCH as usize] +=
                solo_feature.n_read_per_cb_total[icb as usize] as u64;
        }
        if read_feat_sum.stats.v.len()
            > SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE as usize
        {
            read_feat_sum.stats.v
                [SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE as usize] +=
                solo_feature.n_read_per_cb_unique[icb as usize] as u64;
        }
    }
    Ok(())
}

#[doc = "Original `SoloFeature::collapseUMIperCB` at STAR/source/SoloFeature_collapseUMIall.cpp:30. Args: iCB: uint32, umiArray: vector<uint32>, gID: vector<uint32>, gReadS: vector<uint32>"]
pub fn solofeature_collapseumiall_l30_solofeature_collapseumipercb(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    p_solo: &crate::parameters_solo::ParametersSolo,
    i_cb: u32,
    umi_array: &mut [u32],
    g_id: &mut [u32],
    g_read_s: &mut [u32],
) -> Result<(), String> {
    const UMI_ARRAY_STRIDE: usize = 3;
    const RGU_G: usize = 0;
    const RGU_U: usize = 1;
    const RGU_R: usize = 2;
    const GENE_MULT_MARK: u32 = 1_u32 << 31;
    const UMI_DEDUP_NODEDUP: usize = 0;
    const UMI_DEDUP_EXACT: usize = 1;
    const UMI_DEDUP_ALL: usize = 2;
    const UMI_DEDUP_DIRECTIONAL: usize = 3;
    const UMI_DEDUP_CR: usize = 4;
    const UMI_DEDUP_DIRECTIONAL_UMITOOLS: usize = 5;
    const MULTI_UNIFORM: usize = 1;
    const MULTI_RESCUE: usize = 2;
    const MULTI_PROP_UNIQUE: usize = 3;
    const MULTI_EM: usize = 4;

    let stride = solo_feature.rgu_stride as usize;
    if stride < 2 {
        return Err("SoloFeature::collapseUMIperCB requires rguStride >= 2".to_string());
    }
    if stride < 3 && !solo_feature.read_info.is_empty() {
        return Err(
            "SoloFeature::collapseUMIperCB requires read index column for readInfo".to_string(),
        );
    }

    let icb = i_cb as usize;
    let mut r_gu = solo_feature.r_cbp.get(icb).cloned().unwrap_or_default();
    let r_n = solo_feature.n_read_per_cb.get(icb).copied().unwrap_or(0);

    let mut records = Vec::<Vec<u32>>::with_capacity(r_n as usize);
    for rec in r_gu.chunks(stride).take(r_n as usize) {
        records.push(rec.to_vec());
    }
    records.sort_by(|a, b| {
        a[RGU_G]
            .cmp(&b[RGU_G])
            .then_with(|| a[RGU_U].cmp(&b[RGU_U]))
    });
    r_gu.clear();
    for rec in records {
        r_gu.extend_from_slice(&rec);
    }

    let mut gid1 = u32::MAX;
    let mut n_genes = 0_usize;
    let mut n_genes_mult = 0_usize;
    for i_r in (0..r_n as usize * stride).step_by(stride) {
        if r_gu[i_r + RGU_G] != gid1 {
            if n_genes >= g_read_s.len() || n_genes >= g_id.len() {
                return Err(
                    "SoloFeature::collapseUMIperCB gene scratch array is too small".to_string(),
                );
            }
            g_read_s[n_genes] = i_r as u32;
            gid1 = r_gu[i_r + RGU_G];
            g_id[n_genes] = gid1;
            n_genes += 1;
            if p_solo.multi_map.yes_multi && (gid1 & GENE_MULT_MARK) != 0 {
                n_genes_mult += 1;
            }
        }
    }
    if n_genes >= g_read_s.len() {
        return Err("SoloFeature::collapseUMIperCB gene scratch array is too small".to_string());
    }
    g_read_s[n_genes] = (stride * r_n as usize) as u32;
    n_genes -= n_genes_mult;

    let mut umi_gene_map_count =
        std::collections::BTreeMap::<u32, std::collections::BTreeMap<u32, u32>>::new();
    let mut umi_gene_map_count0 =
        std::collections::BTreeMap::<u32, std::collections::BTreeMap<u32, u32>>::new();

    if p_solo.umi_filtering.multi_gene_umi {
        for i_r in (0..g_read_s[n_genes] as usize).step_by(stride) {
            *umi_gene_map_count
                .entry(r_gu[i_r + RGU_U])
                .or_default()
                .entry(r_gu[i_r + RGU_G])
                .or_insert(0) += 1;
        }
        for counts in umi_gene_map_count.values_mut() {
            if counts.len() == 1 {
                continue;
            }
            let mut maxu = counts.values().copied().max().unwrap_or(0);
            if maxu == 1 {
                maxu = 2;
            }
            for count in counts.values_mut() {
                if maxu > *count {
                    *count = 0;
                }
            }
        }
    }

    if p_solo.umi_filtering.multi_gene_umi_all {
        for i_r in (0..g_read_s[n_genes] as usize).step_by(stride) {
            *umi_gene_map_count
                .entry(r_gu[i_r + RGU_U])
                .or_default()
                .entry(r_gu[i_r + RGU_G])
                .or_insert(0) += 1;
        }
        for counts in umi_gene_map_count.values_mut() {
            if counts.len() > 1 {
                for count in counts.values_mut() {
                    *count = 0;
                }
            }
        }
    }

    let mut umi_corrected = vec![std::collections::BTreeMap::<u32, u32>::new(); n_genes];

    let min_needed = solo_feature.count_cell_gene_umi_index[icb] as usize
        + n_genes * solo_feature.count_mat_stride as usize;
    if solo_feature.count_cell_gene_umi.len() < min_needed {
        solo_feature.count_cell_gene_umi.resize(
            (solo_feature.count_cell_gene_umi.len()
                + n_genes * solo_feature.count_mat_stride as usize)
                * 2,
            0,
        );
    }

    solo_feature.n_gene_per_cb[icb] = 0;
    solo_feature.n_umi_per_cb[icb] = 0;
    solo_feature.count_cell_gene_umi_index[icb + 1] = solo_feature.count_cell_gene_umi_index[icb];

    for i_g in 0..n_genes {
        let start = g_read_s[i_g] as usize;
        let end = g_read_s[i_g + 1] as usize;
        let n_r0 = (end - start) / stride;
        if n_r0 == 0 {
            continue;
        }

        let mut gene_records = Vec::<Vec<u32>>::with_capacity(n_r0);
        for rec in r_gu[start..end].chunks(stride) {
            gene_records.push(rec.to_vec());
        }
        gene_records.sort_by(|a, b| a[RGU_U].cmp(&b[RGU_U]));
        for (irec, rec) in gene_records.iter().enumerate() {
            let pos = start + irec * stride;
            r_gu[pos..pos + stride].copy_from_slice(rec);
        }

        let mut i_r1 = 0_usize;
        let mut u1 = u32::MAX;
        let mut n_u0 = 0_usize;
        for rel in (RGU_U..end - start).step_by(stride) {
            let pos = start + rel;
            if p_solo.umi_filtering.multi_gene_umi
                && umi_gene_map_count
                    .get(&r_gu[pos])
                    .and_then(|m| m.get(&g_id[i_g]))
                    .copied()
                    .unwrap_or(0)
                    == 0
            {
                if p_solo.umi_dedup.type_main != UMI_DEDUP_NODEDUP as i32 {
                    r_gu[pos] = u32::MAX;
                }
                continue;
            }

            if r_gu[pos] != u1 {
                i_r1 = n_u0 * UMI_ARRAY_STRIDE;
                u1 = r_gu[pos];
                umi_array[i_r1] = u1;
                umi_array[i_r1 + 1] = 0;
                umi_array[i_r1 + 2] = 0;
                n_u0 += 1;
            }
            umi_array[i_r1 + 1] += 1;
        }

        if p_solo.umi_filtering.multi_gene_umi_cr {
            if n_u0 == 0 {
                continue;
            }
            for iu in (0..n_u0 * UMI_ARRAY_STRIDE).step_by(UMI_ARRAY_STRIDE) {
                *umi_gene_map_count0
                    .entry(umi_array[iu])
                    .or_default()
                    .entry(i_g as u32)
                    .or_insert(0) += umi_array[iu + 1];
            }
            solofeature_collapseumiall_l580_solofeature_umiarraycorrect_cr(
                n_u0 as u32,
                umi_array,
                UMI_ARRAY_STRIDE as u32,
                !solo_feature.read_info.is_empty(),
                false,
                &mut umi_corrected[i_g],
            );
            for iu in (0..n_u0 * UMI_ARRAY_STRIDE).step_by(UMI_ARRAY_STRIDE) {
                *umi_gene_map_count
                    .entry(umi_array[iu + 2])
                    .or_default()
                    .entry(i_g as u32)
                    .or_insert(0) += umi_array[iu + 1];
            }
            continue;
        }

        let rec_index = solo_feature.count_cell_gene_umi_index[icb + 1] as usize;
        let stride_count = solo_feature.count_mat_stride as usize;
        if solo_feature.count_cell_gene_umi.len() < rec_index + stride_count {
            solo_feature
                .count_cell_gene_umi
                .resize((rec_index + stride_count) * 2, 0);
        }
        for value in &mut solo_feature.count_cell_gene_umi[rec_index..rec_index + stride_count] {
            *value = 0;
        }

        if p_solo.umi_dedup.yes_b[UMI_DEDUP_NODEDUP] {
            solo_feature.count_cell_gene_umi
                [rec_index + p_solo.umi_dedup.count_ind_i[UMI_DEDUP_NODEDUP] as usize] =
                n_r0 as u32;
        }

        if n_u0 > 0 {
            if p_solo.umi_dedup.yes_b[UMI_DEDUP_EXACT] {
                solo_feature.count_cell_gene_umi
                    [rec_index + p_solo.umi_dedup.count_ind_i[UMI_DEDUP_EXACT] as usize] =
                    n_u0 as u32;
            }
            if p_solo.umi_dedup.yes_b[UMI_DEDUP_CR] {
                solo_feature.count_cell_gene_umi
                    [rec_index + p_solo.umi_dedup.count_ind_i[UMI_DEDUP_CR] as usize] =
                    solofeature_collapseumiall_l580_solofeature_umiarraycorrect_cr(
                        n_u0 as u32,
                        umi_array,
                        UMI_ARRAY_STRIDE as u32,
                        !solo_feature.read_info.is_empty()
                            && p_solo.umi_dedup.type_main == UMI_DEDUP_CR as i32,
                        true,
                        &mut umi_corrected[i_g],
                    );
            }
            if p_solo.umi_dedup.yes_b[UMI_DEDUP_DIRECTIONAL] {
                solo_feature.count_cell_gene_umi
                    [rec_index + p_solo.umi_dedup.count_ind_i[UMI_DEDUP_DIRECTIONAL] as usize] =
                    solofeature_collapseumiall_l617_solofeature_umiarraycorrect_directional(
                        n_u0 as u32,
                        umi_array,
                        UMI_ARRAY_STRIDE as u32,
                        !solo_feature.read_info.is_empty()
                            && p_solo.umi_dedup.type_main == UMI_DEDUP_DIRECTIONAL as i32,
                        true,
                        &mut umi_corrected[i_g],
                        0,
                    );
            }
            if p_solo.umi_dedup.yes_b[UMI_DEDUP_DIRECTIONAL_UMITOOLS] {
                solo_feature.count_cell_gene_umi[rec_index
                    + p_solo.umi_dedup.count_ind_i[UMI_DEDUP_DIRECTIONAL_UMITOOLS] as usize] =
                    solofeature_collapseumiall_l617_solofeature_umiarraycorrect_directional(
                        n_u0 as u32,
                        umi_array,
                        UMI_ARRAY_STRIDE as u32,
                        !solo_feature.read_info.is_empty()
                            && p_solo.umi_dedup.type_main == UMI_DEDUP_DIRECTIONAL_UMITOOLS as i32,
                        true,
                        &mut umi_corrected[i_g],
                        -1,
                    );
            }
            if p_solo.umi_dedup.yes_b[UMI_DEDUP_ALL] {
                solo_feature.count_cell_gene_umi
                    [rec_index + p_solo.umi_dedup.count_ind_i[UMI_DEDUP_ALL] as usize] =
                    solofeature_collapseumi_graph_l16_solofeature_umiarraycorrect_graph(
                        p_solo,
                        n_u0 as u32,
                        umi_array,
                        UMI_ARRAY_STRIDE as u32,
                        !solo_feature.read_info.is_empty()
                            && p_solo.umi_dedup.type_main == UMI_DEDUP_ALL as i32,
                        true,
                        &mut umi_corrected[i_g],
                    );
            }
        }

        let mut totcount = 0_u32;
        for ii in rec_index + 1..rec_index + stride_count {
            totcount += solo_feature.count_cell_gene_umi[ii];
        }
        if totcount > 0 {
            solo_feature.count_cell_gene_umi[rec_index] = g_id[i_g];
            solo_feature.n_gene_per_cb[icb] += 1;
            solo_feature.n_umi_per_cb[icb] += solo_feature.count_cell_gene_umi
                [rec_index + p_solo.umi_dedup.count_ind_main as usize];
            solo_feature.count_cell_gene_umi_index[icb + 1] += solo_feature.count_mat_stride;
        }

        if !solo_feature.read_info.is_empty() {
            for i_r in (0..end - start).step_by(stride) {
                let pos = start + i_r;
                let iread1 = r_gu[pos + RGU_R] as usize;
                if let Some(info) = solo_feature.read_info.get_mut(iread1) {
                    info.cb = solo_feature.ind_cb[icb] as i64;
                    let mut umi = r_gu[pos + RGU_U];
                    if let Some(corrected) = umi_corrected[i_g].get(&umi) {
                        umi = *corrected;
                    }
                    info.umi = umi as u64;
                }
            }
        }
    }

    if p_solo.umi_filtering.multi_gene_umi_cr {
        let mut gene_counts = vec![0_u32; n_genes];
        let mut gene_umi_hash = if solo_feature.read_info.is_empty() {
            Vec::new()
        } else {
            vec![std::collections::BTreeSet::<u32>::new(); n_genes]
        };

        for (umi, counts) in umi_gene_map_count.iter() {
            let mut maxu = 0_u32;
            let mut maxg = u32::MAX;
            for (gene, count) in counts.iter() {
                if *count > maxu {
                    maxu = *count;
                    maxg = *gene;
                } else if *count == maxu {
                    maxg = u32::MAX;
                }
            }
            if maxg == u32::MAX {
                continue;
            }
            if let Some(counts0) = umi_gene_map_count0.get(umi) {
                for count in counts0.values() {
                    if *count > *counts0.get(&maxg).unwrap_or(&0) {
                        maxg = u32::MAX;
                        break;
                    }
                }
            }
            if maxg != u32::MAX {
                gene_counts[maxg as usize] += 1;
                if !solo_feature.read_info.is_empty() {
                    gene_umi_hash[maxg as usize].insert(*umi);
                }
            }
        }

        for ig in 0..n_genes {
            if gene_counts[ig] == 0 {
                continue;
            }
            let rec_index = solo_feature.count_cell_gene_umi_index[icb + 1] as usize;
            let stride_count = solo_feature.count_mat_stride as usize;
            if solo_feature.count_cell_gene_umi.len() < rec_index + stride_count {
                solo_feature
                    .count_cell_gene_umi
                    .resize((rec_index + stride_count) * 2, 0);
            }
            for value in &mut solo_feature.count_cell_gene_umi[rec_index..rec_index + stride_count]
            {
                *value = 0;
            }
            solo_feature.n_gene_per_cb[icb] += 1;
            solo_feature.n_umi_per_cb[icb] += gene_counts[ig];
            solo_feature.count_cell_gene_umi[rec_index] = g_id[ig];
            solo_feature.count_cell_gene_umi
                [rec_index + p_solo.umi_dedup.count_ind_i[UMI_DEDUP_CR] as usize] = gene_counts[ig];
            solo_feature.count_cell_gene_umi_index[icb + 1] += solo_feature.count_mat_stride;
        }

        if !solo_feature.read_info.is_empty() {
            for i_g in 0..n_genes {
                let start = g_read_s[i_g] as usize;
                let end = g_read_s[i_g + 1] as usize;
                for i_r in (0..end - start).step_by(stride) {
                    let pos = start + i_r;
                    let iread1 = r_gu[pos + RGU_R] as usize;
                    if let Some(info) = solo_feature.read_info.get_mut(iread1) {
                        info.cb = solo_feature.ind_cb[icb] as i64;
                        let mut umi = r_gu[pos + RGU_U];
                        if let Some(corrected) = umi_corrected[i_g].get(&umi) {
                            umi = *corrected;
                        }
                        info.umi = if gene_umi_hash[i_g].contains(&umi) {
                            umi as u64
                        } else {
                            u32::MAX as u64
                        };
                    }
                }
            }
        }
    }

    if p_solo.multi_map.yes_multi {
        solo_feature.count_mat_mult_i[icb + 1] = solo_feature.count_mat_mult_i[icb];
    }

    if n_genes_mult > 0 {
        if !solo_feature.read_info.is_empty() {
            for i_r in (g_read_s[n_genes] as usize..g_read_s[n_genes + n_genes_mult] as usize)
                .step_by(stride)
            {
                let iread1 = r_gu[i_r + RGU_R] as usize;
                if let Some(info) = solo_feature.read_info.get_mut(iread1) {
                    info.cb = solo_feature.ind_cb[icb] as i64;
                    info.umi = r_gu[i_r + RGU_U] as u64;
                }
            }
        }

        let mut umi_genes = Vec::<Vec<u32>>::new();
        let m_start = g_read_s[n_genes] as usize;
        let m_end = g_read_s[n_genes + n_genes_mult] as usize;
        let n_rm = (m_end - m_start) / stride;
        let mut multi_records = Vec::<Vec<u32>>::with_capacity(n_rm);
        for rec in r_gu[m_start..m_end].chunks(stride) {
            multi_records.push(rec.to_vec());
        }
        multi_records.sort_by(|a, b| {
            a[RGU_U]
                .cmp(&b[RGU_U])
                .then_with(|| {
                    a.get(RGU_R)
                        .copied()
                        .unwrap_or(0)
                        .cmp(&b.get(RGU_R).copied().unwrap_or(0))
                })
                .then_with(|| a[RGU_G].cmp(&b[RGU_G]))
        });
        for (irec, rec) in multi_records.iter().enumerate() {
            let pos = m_start + irec * stride;
            r_gu[pos..pos + stride].copy_from_slice(rec);
        }

        let mut gene_read_count = std::collections::BTreeMap::<u32, u32>::new();
        let mut n_rumi = 0_u32;
        let mut skip_umi = false;
        let mut umi_prev = u32::MAX;
        let mut read_prev = u32::MAX;
        for i_r in (0..n_rm * stride).step_by(stride) {
            let pos = m_start + i_r;
            let umi1 = r_gu[pos + RGU_U];
            if umi1 != umi_prev {
                umi_prev = umi1;
                if umi_gene_map_count.contains_key(&umi1) {
                    skip_umi = true;
                } else {
                    skip_umi = false;
                    gene_read_count.clear();
                    n_rumi = 0;
                    read_prev = u32::MAX;
                }
            }
            if skip_umi {
                continue;
            }
            let read1 = r_gu.get(pos + RGU_R).copied().unwrap_or(0);
            if read1 != read_prev {
                n_rumi += 1;
                read_prev = read1;
            }
            let g1 = r_gu[pos + RGU_G] ^ GENE_MULT_MARK;
            *gene_read_count.entry(g1).or_insert(0) += 1;
            let next_is_new = i_r == n_rm * stride - stride || umi1 != r_gu[pos + RGU_U + stride];
            if next_is_new {
                let mut vg = Vec::new();
                for (gene, count) in gene_read_count.iter() {
                    if *count == n_rumi {
                        vg.push(*gene);
                    }
                }
                umi_genes.push(vg);
            }
        }

        let mut genes_m = std::collections::BTreeMap::<u32, u32>::new();
        let mut ng = 0_u32;
        for ug in umi_genes.iter_mut() {
            for gg in ug.iter_mut() {
                let entry = *genes_m.entry(*gg).or_insert_with(|| {
                    let out = ng;
                    ng += 1;
                    out
                });
                *gg = entry;
            }
        }

        let mut g_euniform = vec![0.0_f64; genes_m.len()];
        for ug in umi_genes.iter() {
            if ug.is_empty() {
                continue;
            }
            for gg in ug.iter() {
                g_euniform[*gg as usize] += 1.0 / ug.len() as f64;
            }
        }

        let mut g_erescue = vec![Vec::<f64>::new(); p_solo.umi_dedup.yes_n as usize];
        if p_solo.multi_map.yes_b[MULTI_RESCUE] {
            for ind_dedup in 0..p_solo.umi_dedup.yes_n as usize {
                let mut g_eu = vec![0.0_f64; genes_m.len()];
                for igm in (solo_feature.count_cell_gene_umi_index[icb] as usize
                    ..solo_feature.count_cell_gene_umi_index[icb + 1] as usize)
                    .step_by(solo_feature.count_mat_stride as usize)
                {
                    let g1 = solo_feature.count_cell_gene_umi[igm];
                    if let Some(&idx) = genes_m.get(&g1) {
                        g_eu[idx as usize] =
                            solo_feature.count_cell_gene_umi[igm + 1 + ind_dedup] as f64;
                    }
                }
                g_erescue[ind_dedup] = vec![0.0; genes_m.len()];
                for ug in umi_genes.iter() {
                    let norm: f64 = ug
                        .iter()
                        .map(|gg| g_euniform[*gg as usize] + g_eu[*gg as usize])
                        .sum();
                    if norm == 0.0 {
                        continue;
                    }
                    for gg in ug.iter() {
                        g_erescue[ind_dedup][*gg as usize] +=
                            (g_euniform[*gg as usize] + g_eu[*gg as usize]) / norm;
                    }
                }
            }
        }

        let mut g_eprop_unique = vec![Vec::<f64>::new(); p_solo.umi_dedup.yes_n as usize];
        if p_solo.multi_map.yes_b[MULTI_PROP_UNIQUE] {
            for ind_dedup in 0..p_solo.umi_dedup.yes_n as usize {
                let mut g_eu = vec![0.0_f64; genes_m.len()];
                for igm in (solo_feature.count_cell_gene_umi_index[icb] as usize
                    ..solo_feature.count_cell_gene_umi_index[icb + 1] as usize)
                    .step_by(solo_feature.count_mat_stride as usize)
                {
                    let g1 = solo_feature.count_cell_gene_umi[igm];
                    if let Some(&idx) = genes_m.get(&g1) {
                        g_eu[idx as usize] =
                            solo_feature.count_cell_gene_umi[igm + 1 + ind_dedup] as f64;
                    }
                }
                g_eprop_unique[ind_dedup] = vec![0.0; genes_m.len()];
                for ug in umi_genes.iter() {
                    let norm: f64 = ug.iter().map(|gg| g_eu[*gg as usize]).sum();
                    if norm == 0.0 {
                        if !ug.is_empty() {
                            for gg in ug.iter() {
                                g_eprop_unique[ind_dedup][*gg as usize] += 1.0 / ug.len() as f64;
                            }
                        }
                    } else {
                        for gg in ug.iter() {
                            g_eprop_unique[ind_dedup][*gg as usize] += g_eu[*gg as usize] / norm;
                        }
                    }
                }
            }
        }

        let mut g_eem = vec![Vec::<f64>::new(); p_solo.umi_dedup.yes_n as usize];
        if p_solo.multi_map.yes_b[MULTI_EM] {
            for ind_dedup in 0..p_solo.umi_dedup.yes_n as usize {
                let mut g_eu = vec![0.0_f64; genes_m.len()];
                for igm in (solo_feature.count_cell_gene_umi_index[icb] as usize
                    ..solo_feature.count_cell_gene_umi_index[icb + 1] as usize)
                    .step_by(solo_feature.count_mat_stride as usize)
                {
                    let g1 = solo_feature.count_cell_gene_umi[igm];
                    if let Some(&idx) = genes_m.get(&g1) {
                        g_eu[idx as usize] =
                            solo_feature.count_cell_gene_umi[igm + 1 + ind_dedup] as f64;
                    }
                }
                let mut g_em_old = g_euniform.clone();
                for ii in 0..g_em_old.len() {
                    g_em_old[ii] += g_eu[ii];
                }
                let mut g_em_new = vec![0.0; genes_m.len()];
                let mut iter_i = 0_u32;
                loop {
                    iter_i += 1;
                    g_em_new.copy_from_slice(&g_eu);
                    for value in g_em_old.iter_mut() {
                        if *value < 0.01 {
                            *value = 0.0;
                        }
                    }
                    for ug in umi_genes.iter() {
                        let norm: f64 = ug.iter().map(|gg| g_em_old[*gg as usize]).sum();
                        if norm == 0.0 {
                            continue;
                        }
                        for gg in ug.iter() {
                            g_em_new[*gg as usize] += g_em_old[*gg as usize] / norm;
                        }
                    }
                    let mut max_abs_change = 0.0_f64;
                    for ii in 0..g_em_new.len() {
                        max_abs_change = max_abs_change.max((g_em_new[ii] - g_em_old[ii]).abs());
                    }
                    if max_abs_change < 0.01 || iter_i > 100 {
                        g_eem[ind_dedup] = g_em_new.clone();
                        break;
                    }
                    std::mem::swap(&mut g_em_old, &mut g_em_new);
                }
                for ii in 0..g_eem[ind_dedup].len() {
                    g_eem[ind_dedup][ii] -= g_eu[ii];
                }
            }
        }

        if solo_feature.count_mat_mult_m.len()
            < solo_feature.count_mat_mult_i[icb + 1] as usize
                + genes_m.len()
                    * solo_feature.count_mat_mult_s as usize
                    * p_solo.umi_dedup.yes_n as usize
                + 100
        {
            let new_len = (solo_feature.count_mat_mult_i[icb + 1] as usize
                + genes_m.len()
                    * solo_feature.count_mat_mult_s as usize
                    * p_solo.umi_dedup.yes_n as usize
                + 100)
                * 2;
            solo_feature.count_mat_mult_m.resize(new_len, 0.0);
        }

        for (gene, idx) in genes_m.iter() {
            let base = solo_feature.count_mat_mult_i[icb + 1] as usize;
            solo_feature.count_mat_mult_m[base] = *gene as f64;
            for ind_dedup in 0..p_solo.umi_dedup.yes_n as usize {
                let ind1 = solo_feature.count_mat_mult_i[icb + 1] as usize + ind_dedup;
                if p_solo.multi_map.yes_b[MULTI_UNIFORM] {
                    solo_feature.count_mat_mult_m
                        [ind1 + p_solo.multi_map.count_ind_i[MULTI_UNIFORM] as usize] =
                        g_euniform[*idx as usize];
                }
                if p_solo.multi_map.yes_b[MULTI_RESCUE] {
                    solo_feature.count_mat_mult_m
                        [ind1 + p_solo.multi_map.count_ind_i[MULTI_RESCUE] as usize] =
                        g_erescue[ind_dedup][*idx as usize];
                }
                if p_solo.multi_map.yes_b[MULTI_PROP_UNIQUE] {
                    solo_feature.count_mat_mult_m
                        [ind1 + p_solo.multi_map.count_ind_i[MULTI_PROP_UNIQUE] as usize] =
                        g_eprop_unique[ind_dedup][*idx as usize];
                }
                if p_solo.multi_map.yes_b[MULTI_EM] {
                    solo_feature.count_mat_mult_m
                        [ind1 + p_solo.multi_map.count_ind_i[MULTI_EM] as usize] =
                        g_eem[ind_dedup][*idx as usize];
                }
                solo_feature.count_mat_mult_i[icb + 1] += solo_feature.count_mat_mult_s;
            }
        }
    }

    if solo_feature.r_cbp.len() <= icb {
        solo_feature.r_cbp.resize(icb + 1, Vec::new());
    }
    solo_feature.r_cbp[icb] = r_gu;
    Ok(())
}

#[doc = "Original `funCompareSolo1` at STAR/source/SoloFeature_collapseUMIall.cpp:540. Args: a: void, b: void"]
pub fn solofeature_collapseumiall_l540_funcomparesolo1(a: &[u32], b: &[u32]) -> i32 {
    if a[1] > b[1] {
        1
    } else if a[1] < b[1] {
        -1
    } else if a[0] > b[0] {
        1
    } else if a[0] < b[0] {
        -1
    } else {
        0
    }
}

#[doc = "Original `funCompare_uint32_1_2_0` at STAR/source/SoloFeature_collapseUMIall.cpp:557. Args: a: void, b: void"]
pub fn solofeature_collapseumiall_l557_funcompare_uint32_1_2_0(a: &[u32], b: &[u32]) -> i32 {
    if a[1] > b[1] {
        1
    } else if a[1] < b[1] {
        -1
    } else if a[2] > b[2] {
        1
    } else if a[2] < b[2] {
        -1
    } else if a[0] > b[0] {
        1
    } else if a[0] < b[0] {
        -1
    } else {
        0
    }
}

#[doc = "Original `SoloFeature::umiArrayCorrect_CR` at STAR/source/SoloFeature_collapseUMIall.cpp:580. Args: nU0: uint32, umiArr: uintUMI, readInfoRec: bool, nUMIyes: bool, umiCorr: unordered_map <uintUMI,uintUMI>"]
pub fn solofeature_collapseumiall_l580_solofeature_umiarraycorrect_cr(
    n_u0: u32,
    umi_arr: &mut [u32],
    umi_array_stride: u32,
    read_info_rec: bool,
    n_umi_yes: bool,
    umi_corr: &mut std::collections::BTreeMap<u32, u32>,
) -> u32 {
    let stride = umi_array_stride as usize;
    let n_records = n_u0 as usize;
    let mut records = Vec::<Vec<u32>>::with_capacity(n_records);
    for i in 0..n_records {
        records.push(umi_arr[i * stride..(i + 1) * stride].to_vec());
    }
    records.sort_by(
        |a, b| match solofeature_collapseumiall_l540_funcomparesolo1(a, b) {
            x if x < 0 => std::cmp::Ordering::Less,
            x if x > 0 => std::cmp::Ordering::Greater,
            _ => std::cmp::Ordering::Equal,
        },
    );
    for i in 0..n_records {
        umi_arr[i * stride..(i + 1) * stride].copy_from_slice(&records[i]);
    }

    for iu in (0..n_records * stride).step_by(stride) {
        umi_arr[iu + 2] = umi_arr[iu];
        for iuu in ((iu + stride)..n_records * stride).step_by(stride).rev() {
            let uu_xor = umi_arr[iu] ^ umi_arr[iuu];
            if uu_xor != 0 && (uu_xor >> ((uu_xor.trailing_zeros() / 2) * 2)) <= 3 {
                umi_arr[iu + 2] = umi_arr[iuu];
                break;
            }
        }
    }

    if read_info_rec {
        for iu in (0..n_records * stride).step_by(stride) {
            if umi_arr[iu] != umi_arr[iu + 2] {
                umi_corr.insert(umi_arr[iu], umi_arr[iu + 2]);
            }
        }
    }

    if !n_umi_yes {
        0
    } else {
        let mut umi_c = std::collections::BTreeSet::new();
        for iu in (0..n_records * stride).step_by(stride) {
            umi_c.insert(umi_arr[iu + 2]);
        }
        umi_c.len() as u32
    }
}

#[doc = "Original `SoloFeature::umiArrayCorrect_Directional` at STAR/source/SoloFeature_collapseUMIall.cpp:617. Args: nU0: uint32, umiArr: uintUMI, readInfoRec: bool, nUMIyes: bool, umiCorr: unordered_map <uintUMI,uintUMI>, dirCountAdd: int32"]
pub fn solofeature_collapseumiall_l617_solofeature_umiarraycorrect_directional(
    n_u0: u32,
    umi_arr: &mut [u32],
    umi_array_stride: u32,
    read_info_rec: bool,
    n_umi_yes: bool,
    umi_corr: &mut std::collections::BTreeMap<u32, u32>,
    dir_count_add: i32,
) -> u32 {
    let stride = umi_array_stride as usize;
    let n_records = n_u0 as usize;
    let mut records = Vec::<Vec<u32>>::with_capacity(n_records);
    for i in 0..n_records {
        records.push(umi_arr[i * stride..(i + 1) * stride].to_vec());
    }
    records.sort_by(
        |a, b| match servicefuns_l39_funcomparenumbersreverseshift::<u32, 1>(a, b) {
            x if x < 0 => std::cmp::Ordering::Less,
            x if x > 0 => std::cmp::Ordering::Greater,
            _ => std::cmp::Ordering::Equal,
        },
    );
    for i in 0..n_records {
        umi_arr[i * stride..(i + 1) * stride].copy_from_slice(&records[i]);
    }

    for iu in (0..n_records * stride).step_by(stride) {
        umi_arr[iu + 2] = umi_arr[iu];
    }

    let mut n_u1 = n_u0;
    for iu in (stride..n_records * stride).step_by(stride) {
        for iuu in (0..iu).step_by(stride) {
            let uu_xor = umi_arr[iu] ^ umi_arr[iuu];
            if uu_xor != 0
                && (uu_xor >> ((uu_xor.trailing_zeros() / 2) * 2)) <= 3
                && (umi_arr[iuu + 1] as i64)
                    >= 2_i64 * umi_arr[iu + 1] as i64 + dir_count_add as i64
            {
                umi_arr[iu + 2] = umi_arr[iuu + 2];
                n_u1 -= 1;
                break;
            }
        }
    }

    if read_info_rec {
        for iu in (0..n_records * stride).step_by(stride) {
            if umi_arr[iu] != umi_arr[iu + 2] {
                umi_corr.insert(umi_arr[iu], umi_arr[iu + 2]);
            }
        }
    }

    if !n_umi_yes {
        0
    } else {
        let mut umi_c = std::collections::BTreeSet::new();
        for iu in (0..n_records * stride).step_by(stride) {
            umi_c.insert(umi_arr[iu + 2]);
        }
        let _ = n_u1;
        umi_c.len() as u32
    }
}
