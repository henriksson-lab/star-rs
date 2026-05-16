#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::countVelocyto` at STAR/source/SoloFeature_countVelocyto.cpp:12. Args: "]
pub fn solofeature_countvelocyto_l12_solofeature_countvelocyto(
    solo_feature: &mut crate::solo_feature::SoloFeature,
    gene_feature_read_info: &[crate::solo_feature::SoloFeatureReadInfo],
    transcriptome: &crate::transcriptome::Transcriptome,
    run_thread_n: i32,
    time_allocated: &str,
    time_finished_input: &str,
    time_finished_collapsing: &str,
    linux_proc_memory: &str,
) -> Result<String, String> {
    use std::collections::{BTreeMap, HashMap};
    use std::fmt::Write;

    solo_feature
        .n_read_per_cb
        .resize(solo_feature.n_cb as usize, 0);

    let read_feat_sum = solo_feature
        .read_feat_sum
        .as_mut()
        .expect("SoloFeature::countVelocyto requires readFeatSum");

    let mut cu_tr_types: Vec<HashMap<u64, Vec<crate::solo_read_feature_record::TrTypeStruct>>> =
        Vec::with_capacity(solo_feature.n_cb as usize);
    for ii in 0..solo_feature.n_cb as usize {
        let cap = if read_feat_sum.cb_read_count[ii] > 100 {
            read_feat_sum.cb_read_count[ii]
        } else {
            read_feat_sum.cb_read_count[ii] / 5
        };
        cu_tr_types.push(HashMap::with_capacity(cap as usize));
    }

    let mut log_main = String::new();
    writeln!(
        log_main,
        "{} ... Velocyto counting: allocated arrays",
        time_allocated
    )
    .unwrap();

    for i_thread in 0..run_thread_n as usize {
        let Some(read_feat) = solo_feature.read_feat_all.get(i_thread) else {
            continue;
        };
        let mut tokens = read_feat.stream_reads.split_whitespace();
        while let Some(iread_token) = tokens.next() {
            let iread = iread_token.parse::<usize>().map_err(|_| {
                format!("Malformed STARsolo Velocyto record: invalid read index {iread_token}")
            })?;
            if iread >= gene_feature_read_info.len() {
                return Err(format!(
                    "Malformed STARsolo Velocyto record: read index {iread} is outside read info"
                ));
            }
            let cb = gene_feature_read_info[iread].cb;
            let umi = gene_feature_read_info[iread].umi;
            if cb == -1 || umi == u64::MAX {
                let n_tr = tokens
                    .next()
                    .and_then(|x| x.parse::<usize>().ok())
                    .unwrap_or(0);
                for _ in 0..2 * n_tr {
                    let _ = tokens.next();
                }
                continue;
            }

            if cb < 0 || cb as usize >= solo_feature.ind_cb_wl.len() {
                return Err(format!(
                    "Malformed STARsolo Velocyto record: cell barcode index {cb} is outside whitelist"
                ));
            }
            let i_cb = solo_feature.ind_cb_wl[cb as usize] as usize;
            if i_cb >= cu_tr_types.len() {
                return Err(format!(
                    "Malformed STARsolo Velocyto record: internal cell barcode index {i_cb} is outside allocated cells"
                ));
            }
            solo_feature.n_read_per_cb[i_cb] += 1;

            if cu_tr_types[i_cb].get(&umi).is_some_and(|v| v.is_empty()) {
                let n_tr = tokens
                    .next()
                    .and_then(|x| x.parse::<usize>().ok())
                    .unwrap_or(0);
                for _ in 0..2 * n_tr {
                    let _ = tokens.next();
                }
                continue;
            }

            let n_tr_token = tokens.next().ok_or_else(|| {
                "Malformed STARsolo Velocyto record: missing transcript count".to_string()
            })?;
            let n_tr = n_tr_token.parse::<usize>().map_err(|_| {
                format!("Malformed STARsolo Velocyto record: invalid transcript count {n_tr_token}")
            })?;
            let mut tr_t = Vec::with_capacity(n_tr);
            for _ in 0..n_tr {
                let tr_token = tokens.next().ok_or_else(|| {
                    "Malformed STARsolo Velocyto record: missing transcript id".to_string()
                })?;
                let tr = tr_token.parse::<u32>().map_err(|_| {
                    format!("Malformed STARsolo Velocyto record: invalid transcript id {tr_token}")
                })?;
                let type_token = tokens.next().ok_or_else(|| {
                    "Malformed STARsolo Velocyto record: missing transcript type".to_string()
                })?;
                let type_ = type_token.parse::<u32>().map_err(|_| {
                    format!(
                        "Malformed STARsolo Velocyto record: invalid transcript type {type_token}"
                    )
                })? as u8;
                tr_t.push(crate::solo_read_feature_record::TrTypeStruct { tr, type_ });
            }

            if !cu_tr_types[i_cb].contains_key(&umi) {
                cu_tr_types[i_cb].insert(umi, tr_t);
                continue;
            }

            let old = cu_tr_types[i_cb].get(&umi).cloned().unwrap();
            let mut inew = 0usize;
            let mut tr_t1 = Vec::with_capacity(old.len());
            for old_tt in old {
                while inew < tr_t.len() && old_tt.tr > tr_t[inew].tr {
                    inew += 1;
                }
                if inew == tr_t.len() {
                    break;
                }
                if old_tt.tr == tr_t[inew].tr {
                    tr_t1.push(crate::solo_read_feature_record::TrTypeStruct {
                        tr: tr_t[inew].tr,
                        type_: old_tt.type_ | tr_t[inew].type_,
                    });
                }
            }
            cu_tr_types[i_cb].insert(umi, tr_t1);
        }
    }

    writeln!(
        log_main,
        "{} ... Velocyto counting: finished input",
        time_finished_input
    )
    .unwrap();

    solo_feature
        .n_umi_per_cb
        .resize(solo_feature.n_cb as usize, 0);
    solo_feature
        .n_gene_per_cb
        .resize(solo_feature.n_cb as usize, 0);
    solo_feature.count_mat_stride = 4;
    solo_feature.count_cell_gene_umi.resize(
        solo_feature.n_reads_mapped as usize * solo_feature.count_mat_stride as usize / 5 + 16,
        0,
    );
    solo_feature
        .count_cell_gene_umi_index
        .resize(solo_feature.n_cb as usize + 1, 0);
    solo_feature.count_cell_gene_umi_index[0] = 0;

    for (i_cb, umi_map) in cu_tr_types
        .iter()
        .enumerate()
        .take(solo_feature.n_cb as usize)
    {
        let mut gene_c: BTreeMap<u32, [u32; 3]> = BTreeMap::new();
        for tr_types in umi_map.values() {
            if tr_types.is_empty() {
                continue;
            }
            if tr_types[0].tr as usize >= transcriptome.tr_gene.len() {
                return Err(format!(
                    "Malformed STARsolo Velocyto record: transcript id {} is outside transcriptome",
                    tr_types[0].tr
                ));
            }
            let mut gene_i = transcriptome.tr_gene[tr_types[0].tr as usize];
            let mut exon_model = false;
            let mut intron_model = false;
            let mut span_model = true;
            let mut mixed_model = false;

            for tt in tr_types {
                if tt.tr as usize >= transcriptome.tr_gene.len() {
                    return Err(format!(
                        "Malformed STARsolo Velocyto record: transcript id {} is outside transcriptome",
                        tt.tr
                    ));
                }
                if transcriptome.tr_gene[tt.tr as usize] != gene_i {
                    gene_i = u32::MAX;
                    break;
                }
                let has_intron = (tt.type_ & (1_u8 << ALIGN_VS_TRANSCRIPT_INTRON as u32)) != 0;
                let has_exon_intron =
                    (tt.type_ & (1_u8 << ALIGN_VS_TRANSCRIPT_EXON_INTRON as u32)) != 0;
                let has_span =
                    (tt.type_ & (1_u8 << ALIGN_VS_TRANSCRIPT_EXON_INTRON_SPAN as u32)) != 0;
                let has_concordant =
                    (tt.type_ & (1_u8 << ALIGN_VS_TRANSCRIPT_CONCORDANT as u32)) != 0;

                mixed_model |= ((has_intron && has_concordant) || has_exon_intron) && !has_span;
                span_model &= has_span;
                exon_model |= has_concordant && !has_intron && !has_exon_intron;
                intron_model |= has_intron && !has_exon_intron && !has_concordant;
            }

            if gene_i == u32::MAX {
                continue;
            }

            let counts = gene_c.entry(gene_i).or_insert([0; 3]);
            if exon_model && !intron_model && !mixed_model {
                counts[0] += 1;
            } else if span_model || ((intron_model || mixed_model) && !exon_model) {
                counts[1] += 1;
            } else {
                counts[2] += 1;
            }
            solo_feature.n_umi_per_cb[i_cb] += 1;
        }

        solo_feature.count_cell_gene_umi_index[i_cb + 1] =
            solo_feature.count_cell_gene_umi_index[i_cb];
        if solo_feature.n_umi_per_cb[i_cb] == 0 {
            continue;
        }

        solo_feature.n_gene_per_cb[i_cb] += gene_c.len() as u32;
        read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_UMIS] +=
            solo_feature.n_umi_per_cb[i_cb] as u64;
        read_feat_sum.stats.v[SOLO_READ_FEATURE_STAT_YES_CELL_BARCODES] += 1;

        let needed = solo_feature.count_cell_gene_umi_index[i_cb + 1] as usize
            + gene_c.len() * solo_feature.count_mat_stride as usize;
        if solo_feature.count_cell_gene_umi.len() < needed {
            let new_len = (solo_feature.count_cell_gene_umi.len() * 2).max(needed);
            solo_feature.count_cell_gene_umi.resize(new_len, 0);
        }

        for (gene, counts) in gene_c {
            let idx = solo_feature.count_cell_gene_umi_index[i_cb + 1] as usize;
            solo_feature.count_cell_gene_umi[idx] = gene;
            for ii in 0..3 {
                solo_feature.count_cell_gene_umi[idx + 1 + ii] = counts[ii];
            }
            solo_feature.count_cell_gene_umi_index[i_cb + 1] += solo_feature.count_mat_stride;
        }
    }

    solo_feature.n_read_per_cb_total = solo_feature.n_read_per_cb.clone();
    solo_feature.n_read_per_cb_unique = solo_feature.n_read_per_cb.clone();

    writeln!(
        log_main,
        "{} ... Velocyto counting: finished collapsing UMIs",
        time_finished_collapsing
    )
    .unwrap();
    write!(
        log_main,
        "RAM for solo feature Velocyto\n{}",
        linux_proc_memory
    )
    .unwrap();
    Ok(log_main)
}
