#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::statsOutput` at STAR/source/SoloFeature_statsOutput.cpp:6. Args: "]
pub fn solofeature_statsoutput_l6_solofeature_statsoutput(
    solo_feature: &crate::solo_feature::SoloFeature,
    p: &crate::parameters_chimeric::Parameters,
    p_solo: &crate::parameters_solo::ParametersSolo,
    g_stats_all: &crate::stats::Stats,
    ra_chunks: &[crate::read_align::ReadAlign],
) -> crate::solo_filtered_cells::SoloStatsOutput {
    use std::collections::BTreeMap;
    use std::fmt::Write;

    let read_bar_sum = solo_feature
        .read_bar_sum
        .as_ref()
        .expect("SoloFeature::statsOutput requires readBarSum");
    let read_feat_sum = solo_feature
        .read_feat_sum
        .as_ref()
        .expect("SoloFeature::statsOutput requires readFeatSum");

    let mut summary = String::new();
    writeln!(summary, "Number of Reads,{}", g_stats_all.read_n).unwrap();
    writeln!(
        summary,
        "Reads With Valid Barcodes,{}",
        1.0_f64
            - (read_bar_sum.stats.v.iter().take(9).sum::<u64>()
                + read_feat_sum
                    .stats
                    .v
                    .get(SOLO_READ_FEATURE_STAT_NO_TOO_MANY_WL_MATCHES)
                    .copied()
                    .unwrap_or(0)
                + read_feat_sum
                    .stats
                    .v
                    .get(SOLO_READ_FEATURE_STAT_NO_MM_TO_WL_WITHOUT_EXACT)
                    .copied()
                    .unwrap_or(0)) as f64
                / g_stats_all.read_n as f64
    )
    .unwrap();
    writeln!(
        summary,
        "Sequencing Saturation,{}",
        1.0_f64
            - read_feat_sum
                .stats
                .v
                .get(SOLO_READ_FEATURE_STAT_YES_UMIS)
                .copied()
                .unwrap_or(0) as f64
                / read_feat_sum
                    .stats
                    .v
                    .get(SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE)
                    .copied()
                    .unwrap_or(0) as f64
    )
    .unwrap();

    if p_solo.solo_type != SOLO_TYPE_SMART_SEQ {
        let mut q30 = 0_u64;
        let mut ntot = 0_u64;
        for ix in 0..256_usize {
            let n = read_bar_sum.qual_hist.get(ix).copied().unwrap_or(0);
            ntot += n;
            if ix >= p.read_quality_score_base as usize + 30 {
                q30 += n;
            }
        }
        writeln!(summary, "Q30 Bases in CB+UMI,{}", q30 as f64 / ntot as f64).unwrap();
    }

    {
        let mut q30 = 0_u64;
        let mut ntot = 0_u64;
        for chunk in ra_chunks.iter().take(p.run_thread_n as usize) {
            for imate in 0..p.read_nmates as usize {
                for ix in 0..256_usize {
                    let n = chunk
                        .qual_hist
                        .get(imate)
                        .and_then(|mate_hist| mate_hist.get(ix))
                        .copied()
                        .unwrap_or(0);
                    ntot += n;
                    if ix >= p.read_quality_score_base as usize + 30 {
                        q30 += n;
                    }
                }
            }
        }
        writeln!(
            summary,
            "Q30 Bases in RNA read,{}",
            q30 as f64 / ntot as f64
        )
        .unwrap();
    }

    writeln!(
        summary,
        "Reads Mapped to Genome: Unique+Multiple,{}",
        (g_stats_all.mapped_reads_u + g_stats_all.mapped_reads_m) as f64
            / g_stats_all.read_n as f64
    )
    .unwrap();
    writeln!(
        summary,
        "Reads Mapped to Genome: Unique,{}",
        g_stats_all.mapped_reads_u as f64 / g_stats_all.read_n as f64
    )
    .unwrap();

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
    let mapfeat = feature_names[solo_feature.feature_type as usize];
    write!(
        summary,
        "Reads Mapped to {}: Unique+Multiple {},",
        mapfeat, mapfeat
    )
    .unwrap();
    if p_solo.multi_map.yes_multi {
        writeln!(
            summary,
            "{}",
            read_feat_sum
                .stats
                .v
                .get(SOLO_READ_FEATURE_STAT_YES_WL_MATCH)
                .copied()
                .unwrap_or(0) as f64
                / g_stats_all.read_n as f64
        )
        .unwrap();
    } else {
        summary.push_str("NoMulti\n");
    }

    writeln!(
        summary,
        "Reads Mapped to {}: Unique {},{}",
        mapfeat,
        mapfeat,
        read_feat_sum
            .stats
            .v
            .get(SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE)
            .copied()
            .unwrap_or(0) as f64
            / g_stats_all.read_n as f64
    )
    .unwrap();

    let mut umi_per_cell_sorted = String::new();
    if p_solo
        .cell_filter
        .type_
        .first()
        .map(|x| x.as_str())
        .unwrap_or("None")
        != "None"
        && (solo_feature.feature_type == SOLO_FEATURE_GENE
            || solo_feature.feature_type == SOLO_FEATURE_GENE_FULL
            || solo_feature.feature_type == SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON
            || solo_feature.feature_type == SOLO_FEATURE_GENE_FULL_EX50P_AS)
    {
        writeln!(
            summary,
            "Estimated Number of Cells,{}",
            solo_feature.filtered_cells.n_cells
        )
        .unwrap();
        writeln!(
            summary,
            "Unique Reads in Cells Mapped to {},{}",
            mapfeat, solo_feature.filtered_cells.n_read_in_cells_unique
        )
        .unwrap();
        writeln!(
            summary,
            "Fraction of Unique Reads in Cells,{}",
            solo_feature.filtered_cells.n_read_in_cells_unique as f64
                / read_feat_sum
                    .stats
                    .v
                    .get(SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE)
                    .copied()
                    .unwrap_or(0) as f64
        )
        .unwrap();
        writeln!(
            summary,
            "Mean Reads per Cell,{}",
            solo_feature.filtered_cells.mean_read_per_cell_unique
        )
        .unwrap();
        writeln!(
            summary,
            "Median Reads per Cell,{}",
            solo_feature.filtered_cells.median_read_per_cell_unique
        )
        .unwrap();
        writeln!(
            summary,
            "UMIs in Cells,{}",
            solo_feature.filtered_cells.n_umi_in_cells
        )
        .unwrap();
        writeln!(
            summary,
            "Mean UMI per Cell,{}",
            solo_feature.filtered_cells.mean_umi_per_cell
        )
        .unwrap();
        writeln!(
            summary,
            "Median UMI per Cell,{}",
            solo_feature.filtered_cells.median_umi_per_cell
        )
        .unwrap();
        writeln!(
            summary,
            "Mean {} per Cell,{}",
            mapfeat, solo_feature.filtered_cells.mean_gene_per_cell
        )
        .unwrap();
        writeln!(
            summary,
            "Median {} per Cell,{}",
            mapfeat, solo_feature.filtered_cells.median_gene_per_cell
        )
        .unwrap();
        writeln!(
            summary,
            "Total {} Detected,{}",
            mapfeat, solo_feature.filtered_cells.n_gene_detected
        )
        .unwrap();

        for n in &solo_feature.n_umi_per_cb_sorted {
            if *n == 0 {
                break;
            }
            writeln!(umi_per_cell_sorted, "{n}").unwrap();
        }
    }

    let mut cell_reads_stats = String::new();
    if p_solo
        .read_stats_yes
        .get(solo_feature.feature_type as usize)
        .copied()
        .unwrap_or(false)
    {
        let stat_names = [
            "cbMatch",
            "cbPerfect",
            "cbMMunique",
            "cbMMmultiple",
            "genomeU",
            "genomeM",
            "featureU",
            "featureM",
            "exonic",
            "intronic",
            "exonicAS",
            "intronicAS",
            "mito",
            "countedU",
            "countedM",
        ];
        cell_reads_stats.push_str("CB");
        for sn in stat_names {
            write!(cell_reads_stats, "\t{sn}").unwrap();
        }
        cell_reads_stats.push_str("\tnUMIunique\tnGenesUnique\tnUMImulti\tnGenesMulti\n");

        cell_reads_stats.push_str("CBnotInPasslist");
        for cc in solo_feature.read_flag_counts.flag_counts_no_cb {
            write!(cell_reads_stats, "\t{cc}").unwrap();
        }
        cell_reads_stats.push_str("\t0\t0\t0\t0\n");

        for (cb, counts) in &solo_feature.read_flag_counts.flag_counts {
            cell_reads_stats.push_str(&p_solo.cb_wl_str[*cb as usize]);
            for cc in counts.iter().take(SOLO_READ_FLAG_N_BITS) {
                write!(cell_reads_stats, "\t{cc}").unwrap();
            }

            let ind = solo_feature
                .ind_cb_wl
                .get(*cb as usize)
                .copied()
                .unwrap_or(u32::MAX);
            if ind == u32::MAX {
                cell_reads_stats.push_str("\t0\t0\t0\t0");
            } else {
                let ind = ind as usize;
                write!(
                    cell_reads_stats,
                    "\t{}\t{}",
                    solo_feature.n_umi_per_cb[ind], solo_feature.n_gene_per_cb[ind]
                )
                .unwrap();
                if solo_feature.n_umi_per_cb_multi.is_empty() {
                    cell_reads_stats.push_str("\t0\t0");
                } else {
                    write!(
                        cell_reads_stats,
                        "\t{}\t{}",
                        solo_feature.n_umi_per_cb_multi[ind], solo_feature.n_gene_per_cb_multi[ind]
                    )
                    .unwrap();
                }
            }
            cell_reads_stats.push('\n');
        }
    }

    let mut files = BTreeMap::new();
    files.insert(
        format!("{}Summary.csv", solo_feature.output_prefix),
        summary,
    );
    if !umi_per_cell_sorted.is_empty() {
        files.insert(
            format!("{}UMIperCellSorted.txt", solo_feature.output_prefix),
            umi_per_cell_sorted,
        );
    }
    if !cell_reads_stats.is_empty() {
        files.insert(
            format!("{}CellReads.stats", solo_feature.output_prefix),
            cell_reads_stats,
        );
    }
    crate::solo_filtered_cells::SoloStatsOutput { files }
}
