#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ReadSoloFeatures` at STAR/source/SoloReadFeature_record.cpp:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadSoloFeatures {
    pub gene: u32,
    pub gene_mult: Vec<u32>,
    pub sj: Vec<[u64; 2]>,
    pub sj_annot: bool,
    pub ind_annot_tr: usize,
    pub align_out: Vec<Transcript>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TrTypeStruct {
    pub tr: u32,
    pub type_: u8,
}

#[doc = "Original `SoloReadFeature::record` at STAR/source/SoloReadFeature_record.cpp:20. Args: soloBar: SoloReadBarcode, nTr: uint, alignOut: Transcript, iRead: uint64, readAnnot: ReadAnnotations"]
pub fn soloreadfeature_record_l20_soloreadfeature_record(
    rf: &mut crate::solo_read_feature::SoloReadFeature,
    p: &crate::parameters_chimeric::Parameters,
    solo_bar: &mut crate::solo_read_barcode::SoloReadBarcode,
    n_tr: u32,
    align_out: &[crate::transcript::Transcript],
    i_read: u64,
    read_annot: &crate::read_annotations::ReadAnnotations,
) {
    if p.p_solo.solo_type == 0 {
        return;
    }

    if rf.stats.v.len() < SOLO_READ_FEATURE_N_STATS {
        rf.stats.v.resize(SOLO_READ_FEATURE_N_STATS, 0);
    }

    let feature_index = rf.feature_type as usize;
    let ann_feature_default = crate::read_annotations::ReadAnnotFeature::default();
    let ann_feature = read_annot
        .annot_features
        .get(feature_index)
        .unwrap_or(&ann_feature_default);
    let read_stats_yes = p
        .p_solo
        .read_stats_yes
        .get(feature_index)
        .copied()
        .unwrap_or(false);

    if read_stats_yes {
        if n_tr == 1 {
            rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_GENOME_U;
        } else if n_tr > 1 {
            rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_GENOME_M;
        }

        for tr in align_out.iter().take(n_tr as usize) {
            if p.p_ge.chr_set_mito.contains(&(tr.chr as u64)) {
                rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_MITO;
            }
        }

        match ann_feature.ov_type {
            1 | 3 => rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_EXONIC,
            5 => rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_INTRONIC,
            2 | 4 => rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_EXONIC_AS,
            6 => rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_INTRONIC_AS,
            _ => {}
        }

        if solo_bar.cb_match < 0 && p.p_solo.cb_wl_yes {
            if ann_feature.f_set.len() == 1 {
                rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_FEATURE_U;
            } else if ann_feature.f_set.len() > 1 {
                rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_FEATURE_M;
            }
            rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_CB_MATCH;
            for ibit in 0..SOLO_READ_FLAG_N_BITS {
                rf.read_flag.flag_counts_no_cb[ibit] +=
                    ((rf.read_flag.flag >> ibit) & 1_u32) as u64;
            }
        }
    }

    if solo_bar.cb_match < 0 {
        return;
    }

    let mut re_fe = crate::solo_read_feature_record::ReadSoloFeatures {
        align_out: align_out.to_vec(),
        ..Default::default()
    };

    let mut n_feat = 0_u32;
    if n_tr == 0 {
        rf.stats.v[SOLO_READ_FEATURE_STAT_NO_UNMAPPED] += 1;
    } else {
        match rf.feature_type {
            SOLO_FEATURE_GENE
            | SOLO_FEATURE_GENE_FULL
            | SOLO_FEATURE_GENE_FULL_EX50P_AS
            | SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON => {
                if solo_bar.solo_type == 4 {
                    for itr in (0..n_tr as usize).rev() {
                        if ann_feature
                            .f_align
                            .get(itr)
                            .map(|genes| !genes.is_empty())
                            .unwrap_or(false)
                        {
                            re_fe.ind_annot_tr = itr;
                            break;
                        }
                    }
                }

                if ann_feature.f_set.is_empty() {
                    rf.stats.v[SOLO_READ_FEATURE_STAT_NO_NO_FEATURE] += 1;
                } else if ann_feature.f_set.len() > 1 {
                    rf.stats.v[SOLO_READ_FEATURE_STAT_MULTI_FEATURE] += 1;
                    rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_FEATURE_M;
                    if n_tr > 1 {
                        rf.stats.v[SOLO_READ_FEATURE_STAT_SUB_MULTI_FEATURE_MULTI_GENOMIC] += 1;
                    }
                    if p.p_solo.multi_map.yes_multi {
                        re_fe.gene_mult.reserve(ann_feature.f_set.len());
                        for g in &ann_feature.f_set {
                            re_fe.gene_mult.push(*g | (1_u32 << 31));
                        }
                        n_feat = soloreadfeature_record_l206_outputreadcb(
                            &mut rf.stream_reads,
                            i_read,
                            rf.feature_type,
                            solo_bar,
                            &re_fe,
                            read_annot,
                            &rf.read_flag,
                        );
                    }
                } else {
                    re_fe.gene = *ann_feature.f_set.iter().next().unwrap();
                    rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_FEATURE_U;
                    let out_i_read = if rf.read_index_yes { i_read } else { u64::MAX };
                    n_feat = soloreadfeature_record_l206_outputreadcb(
                        &mut rf.stream_reads,
                        out_i_read,
                        rf.feature_type,
                        solo_bar,
                        &re_fe,
                        read_annot,
                        &rf.read_flag,
                    );
                }
            }
            SOLO_FEATURE_SJ => {
                if n_tr > 1 {
                    rf.stats.v[SOLO_READ_FEATURE_STAT_SUB_MULTI_FEATURE_MULTI_GENOMIC] += 1;
                    rf.stats.v[SOLO_READ_FEATURE_STAT_MULTI_FEATURE] += 1;
                } else if let Some(tr) = align_out.first() {
                    transcript_l38_transcript_extractsplicejunctions(
                        tr,
                        &mut re_fe.sj,
                        &mut re_fe.sj_annot,
                    );
                    if re_fe.sj.is_empty() {
                        rf.stats.v[SOLO_READ_FEATURE_STAT_NO_NO_FEATURE] += 1;
                    } else {
                        rf.read_flag.flag |= 1_u32 << SOLO_READ_FLAG_FEATURE_U;
                        let out_i_read = if rf.read_index_yes { i_read } else { u64::MAX };
                        n_feat = soloreadfeature_record_l206_outputreadcb(
                            &mut rf.stream_reads,
                            out_i_read,
                            rf.feature_type,
                            solo_bar,
                            &re_fe,
                            read_annot,
                            &rf.read_flag,
                        );
                    }
                }
            }
            SOLO_FEATURE_TRANSCRIPT3P => {
                if read_annot.transcript_concordant.is_empty() || solo_bar.cb_match > 1 {
                    rf.stats.v[SOLO_READ_FEATURE_STAT_NO_NO_FEATURE] += 1;
                } else {
                    n_feat = soloreadfeature_record_l206_outputreadcb(
                        &mut rf.stream_reads,
                        i_read,
                        rf.feature_type,
                        solo_bar,
                        &re_fe,
                        read_annot,
                        &rf.read_flag,
                    );
                }
                if read_annot.transcript_concordant.len() == 1 {
                    let dist = read_annot.transcript_concordant[0][1] as usize;
                    if dist < rf.transcript_dist_count.len() {
                        rf.transcript_dist_count[dist] += 1;
                    }
                }
            }
            SOLO_FEATURE_VELOCYTO => {
                if read_annot.tr_velocyto_type.is_empty() {
                    rf.stats.v[SOLO_READ_FEATURE_STAT_NO_NO_FEATURE] += 1;
                } else {
                    let mut tr_velocyto_type = read_annot.tr_velocyto_type.clone();
                    tr_velocyto_type.sort_by_key(|tt| tt.tr);
                    rf.stream_reads
                        .push_str(&format!("{} {}", i_read, tr_velocyto_type.len()));
                    for tt in &tr_velocyto_type {
                        rf.stream_reads
                            .push_str(&format!(" {} {}", tt.tr, tt.type_ as u32));
                    }
                    rf.stream_reads.push('\n');
                    n_feat = 1;
                }
            }
            _ => {}
        }
    }

    if n_feat == 0 && (rf.read_info_yes || read_stats_yes) {
        soloreadfeature_record_l206_outputreadcb(
            &mut rf.stream_reads,
            i_read,
            -1,
            solo_bar,
            &re_fe,
            read_annot,
            &rf.read_flag,
        );
    }

    if n_feat == 0 {
        return;
    }

    if p.p_solo.cb_wl_yes {
        for cbi in &solo_bar.cb_match_ind {
            let cbi = *cbi as usize;
            if cbi >= rf.cb_read_count.len() {
                rf.cb_read_count.resize(cbi + 1, 0);
            }
            rf.cb_read_count[cbi] += n_feat;
        }
    } else if let Some(cbi) = solo_bar.cb_match_ind.first() {
        *rf.cb_read_count_map.entry(*cbi).or_insert(0) += n_feat;
    }
}

#[doc = "Original `outputReadCB` at STAR/source/SoloReadFeature_record.cpp:206. Args: streamOut: fstream, iRead: uint64, featureType: int32, soloBar: SoloReadBarcode, reFe: ReadSoloFeatures, readAnnot: ReadAnnotations, readFlag: SoloReadFlagClass"]
pub fn soloreadfeature_record_l206_outputreadcb(
    stream_out: &mut String,
    i_read: u64,
    feature_type: i32,
    solo_bar: &mut crate::solo_read_barcode::SoloReadBarcode,
    re_fe: &crate::solo_read_feature_record::ReadSoloFeatures,
    read_annot: &crate::read_annotations::ReadAnnotations,
    read_flag: &crate::solo_common::SoloReadFlagClass,
) -> u32 {
    if solo_bar.solo_type == 4 && feature_type != -1 {
        solo_bar.umi_b =
            transcript_l53_transcript_chrstartlengthextended(&re_fe.align_out[re_fe.ind_annot_tr]);
    }

    let mut nout = 1_u32;

    match feature_type {
        -1 => {
            stream_out.push_str(&format!(
                "{} {} {} -1 {} {}\n",
                solo_bar.umi_b, i_read, read_flag.flag, solo_bar.cb_match, solo_bar.cb_match_string
            ));
        }
        SOLO_FEATURE_GENE
        | SOLO_FEATURE_GENE_FULL
        | SOLO_FEATURE_GENE_FULL_EX50P_AS
        | SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON => {
            if re_fe.gene_mult.is_empty() {
                stream_out.push_str(&format!("{} ", solo_bar.umi_b));
                if i_read != u64::MAX {
                    stream_out.push_str(&format!("{} {} ", i_read, read_flag.flag));
                }
                stream_out.push_str(&format!(
                    "{} {} {}\n",
                    re_fe.gene, solo_bar.cb_match, solo_bar.cb_match_string
                ));
            } else {
                for g in &re_fe.gene_mult {
                    stream_out.push_str(&format!(
                        "{} {} {} {} {} {}\n",
                        solo_bar.umi_b,
                        i_read,
                        read_flag.flag,
                        g,
                        solo_bar.cb_match,
                        solo_bar.cb_match_string
                    ));
                }
                nout = re_fe.gene_mult.len() as u32;
            }
        }
        SOLO_FEATURE_SJ => {
            for sj in &re_fe.sj {
                stream_out.push_str(&format!("{} ", solo_bar.umi_b));
                if i_read != u64::MAX {
                    stream_out.push_str(&format!("{} {} ", i_read, read_flag.flag));
                }
                stream_out.push_str(&format!(
                    "{} {} {} {}\n",
                    sj[0], sj[1], solo_bar.cb_match, solo_bar.cb_match_string
                ));
            }
            nout = re_fe.sj.len() as u32;
        }
        SOLO_FEATURE_TRANSCRIPT3P => {
            stream_out.push_str(&format!(
                "{} {} {}",
                solo_bar.cb_match_string,
                solo_bar.umi_b,
                read_annot.transcript_concordant.len()
            ));
            for tt in &read_annot.transcript_concordant {
                stream_out.push_str(&format!(" {} {}", tt[0], tt[1]));
            }
            if i_read != u64::MAX {
                stream_out.push_str(&format!(" {}", i_read));
            }
            stream_out.push('\n');
            nout = 1;
        }
        _ => {}
    }

    nout
}
