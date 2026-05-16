#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::outputAlignments` at STAR/source/ReadAlign_outputAlignments.cpp:5. Args: "]
pub fn readalign_outputalignments_l5_readalign_outputalignments(
    read_align: &mut crate::read_align::ReadAlign,
    p: &mut crate::parameters_chimeric::Parameters,
    map_gen: &crate::genome::Genome,
    transcriptome: &mut crate::transcriptome::Transcriptome,
    solo_read: &mut crate::solo_read::SoloRead,
    tr_mult: &[crate::transcript::Transcript],
    aligns_gen_out: &crate::quantifications::ReadAlignGenomeTransformResult,
    chunk_out_sj: &mut crate::out_sj::OutSJ,
    chunk_out_sj1: &mut crate::out_sj::OutSJ,
    chunk_out_filter_by_sjout_files: &mut [String],
    chunk_out_unmapped_reads_stream: &mut [String],
    read_name_mates: &[String],
    read_name_extra: &[String],
    read0: &[String],
    qual0: &[String],
    read_file_type: i32,
    out_sam_stream: &mut String,
    primary_pick_fraction: f64,
) -> Result<crate::quantifications::OutputAlignmentsResult, String> {
    let mut result = crate::quantifications::OutputAlignmentsResult {
        unmap_type: read_align.unmap_type,
        out_filter_by_sjout_pass: true,
        ..Default::default()
    };
    let mut read_annot = crate::read_annotations::ReadAnnotations::default();

    if map_gen.p_ge.g_type_string == "SuperTranscriptome" {
        result.out_bam_bytes = readalign_outputalignments_l277_readalign_splicegraphwritesam(
            tr_mult,
            read_align.n_tr,
            out_sam_stream,
            read_align.read_filter as u8,
            read_align.unmap_type,
            &read_align.read_name,
            read0,
            read_file_type,
            qual0,
            p,
            read_align.read_files_index,
            read_name_extra,
            read_align.l_read,
            map_gen,
        )?;
        result.sam = out_sam_stream.clone();
        return Ok(result);
    }

    readalign_outputalignments_l90_readalign_outfilterbysjout(
        &mut read_align.unmap_type,
        &mut result.out_filter_by_sjout_pass,
        p.out_filter_by_sjout_stage,
        read_align.n_tr,
        tr_mult,
        &mut read_align.stats_ra,
        &read_align.read_length,
        p.read_nends,
        chunk_out_filter_by_sjout_files,
        read_name_mates,
        read_align.i_read_all,
        read_align.read_filter,
        read_align.read_files_index,
        read_name_extra,
        read0,
        read_file_type,
        qual0,
        p.out_sj,
        &p.out_sj_filter_reads,
        chunk_out_sj1,
    );

    if result.out_filter_by_sjout_pass {
        if read_align.unmap_type < 0 {
            let (n_tr1, tr_out1) = if p.p_ge.transform.out_yes {
                (aligns_gen_out.al_n as u64, aligns_gen_out.al_mult.first())
            } else {
                (read_align.n_tr, tr_mult.first())
            };

            if n_tr1 > 1 {
                read_align.stats_ra.mapped_reads_m += 1;
                read_align.unmap_type = -2;
            } else if n_tr1 == 1 {
                let tr_out1 = tr_out1.ok_or_else(|| {
                    "ReadAlign::outputAlignments expected one transcript for unique mapper"
                        .to_string()
                })?;
                read_align.stats_ra.mapped_reads_u += 1;
                stats_l35_stats_transcriptstats(
                    &mut read_align.stats_ra,
                    tr_out1,
                    read_align.l_read,
                );
            }

            if p.p_ge.transform.out_sam && (!p.two_pass_yes || p.two_pass_pass2) {
                readalign_outputalignments_l76_readalign_recordsj(
                    aligns_gen_out.al_n as u64,
                    &aligns_gen_out.al_mult,
                    chunk_out_sj,
                    p.out_sj,
                    &p.out_sj_filter_reads,
                );
            } else {
                readalign_outputalignments_l76_readalign_recordsj(
                    read_align.n_tr,
                    tr_mult,
                    chunk_out_sj,
                    p.out_sj,
                    &p.out_sj_filter_reads,
                );
            }

            readalign_outputalignments_l298_readalign_alignedannotation(
                transcriptome,
                p,
                read_align.n_tr,
                tr_mult,
                aligns_gen_out,
                &mut read_annot,
            );
        }

        if let Some(read_bar) = solo_read.read_bar.as_mut() {
            let read_len: Vec<u64> = read_align
                .read_length_original
                .iter()
                .map(|len| *len as u64)
                .collect();
            soloreadbarcode_getcbandumi_l147_soloreadbarcode_getcbandumi(
                read_bar,
                p,
                read0,
                qual0,
                &read_len,
                read_name_extra.first().map(String::as_str).unwrap_or(""),
                read_align.read_files_index,
                &read_align.read_name,
            )?;
        }

        if p.quant_tr_sam_yes && read_align.unmap_type < 0 {
            let quant_aligns = if p.p_ge.transform.out_quant {
                &aligns_gen_out.al_mult[..aligns_gen_out.al_n as usize]
            } else {
                &tr_mult[..read_align.n_tr as usize]
            };
            result.quant_transcriptome = Some(
                readalign_quanttranscriptome_l7_readalign_quanttranscriptome(
                    read_align,
                    p,
                    map_gen,
                    transcriptome,
                    quant_aligns,
                    primary_pick_fraction,
                ),
            );
        }

        soloread_record_l3_soloread_record(
            solo_read,
            p,
            if read_align.unmap_type < 0 {
                read_align.n_tr
            } else {
                0
            },
            tr_mult,
            read_align.i_read_all,
            &read_annot,
        );

        let (n_tr_out_sam, tr_out_sam, tr_best_sam) = if p.p_ge.transform.out_sam {
            (
                aligns_gen_out.al_n as u64,
                aligns_gen_out.al_mult.as_slice(),
                &aligns_gen_out.al_best,
            )
        } else {
            (read_align.n_tr, tr_mult, &read_align.tr_best)
        };

        result.write_sam = readalign_outputalignments_l132_readalign_writesam(
            n_tr_out_sam,
            tr_out_sam,
            tr_best_sam,
            read_align,
            p,
            map_gen,
            read_align.unmap_type,
            result.out_filter_by_sjout_pass,
            read0,
            read_file_type,
            qual0,
            read_name_extra,
        )?;
        read_align.unmap_type = result.write_sam.unmap_type;
        result.out_bam_bytes = result.write_sam.out_bam_bytes;
        result.mate_mapped = result.write_sam.mate_mapped;
        out_sam_stream.push_str(&result.write_sam.sam);
        result.sam = result.write_sam.sam.clone();
    }

    if read_align.unmap_type >= 0 {
        read_align.stats_ra.unmapped_all += 1;
        readalign_outputalignments_l259_readalign_outreadsunmapped(
            &p.out_reads_unmapped,
            p.read_nends,
            p.read_nmates,
            read_name_mates,
            read_align.read_filter,
            read_name_extra,
            &result.mate_mapped,
            read0,
            read_file_type,
            qual0,
            chunk_out_unmapped_reads_stream,
        );
    }

    result.unmap_type = read_align.unmap_type;
    Ok(result)
}

#[doc = "Original `ReadAlign::recordSJ` at STAR/source/ReadAlign_outputAlignments.cpp:76. Args: nTrO: uint64, trO: Transcript, cSJ: OutSJ"]
pub fn readalign_outputalignments_l76_readalign_recordsj(
    n_tr_o: u64,
    tr_o: &[crate::transcript::Transcript],
    c_sj: &mut crate::out_sj::OutSJ,
    out_sj_yes: bool,
    out_sj_filter_reads: &str,
) {
    if !out_sj_yes {
        return;
    }

    if out_sj_filter_reads == "All" || n_tr_o == 1 {
        let sj_read_start_n = c_sj.n;
        for i_tr in 0..n_tr_o as usize {
            readalign_outputtranscriptsj_l4_readalign_outputtranscriptsj(
                &tr_o[i_tr],
                n_tr_o as u32,
                c_sj,
                sj_read_start_n,
            );
        }
    }
}

#[doc = "Original `ReadAlign::outFilterBySJout` at STAR/source/ReadAlign_outputAlignments.cpp:90. Args: "]
pub fn readalign_outputalignments_l90_readalign_outfilterbysjout(
    unmap_type: &mut i32,
    out_filter_by_sjout_pass: &mut bool,
    out_filter_by_sjout_stage: i32,
    n_tr: u64,
    tr_mult: &[crate::transcript::Transcript],
    stats_ra: &mut crate::stats::Stats,
    read_length: &[u32],
    read_nends: u32,
    chunk_out_filter_by_sjout_files: &mut [String],
    read_name_mates: &[String],
    i_read_all: u64,
    read_filter: i32,
    read_files_index: u32,
    read_name_extra: &[String],
    read0: &[String],
    read_file_type: i32,
    qual0: &[String],
    out_sj_yes: bool,
    out_sj_filter_reads: &str,
    chunk_out_sj1: &mut crate::out_sj::OutSJ,
) {
    *out_filter_by_sjout_pass = true;

    if *unmap_type > 0 || out_filter_by_sjout_stage != 1 {
        return;
    }

    for i_tr in 0..n_tr as usize {
        for iex in 0..tr_mult[i_tr].n_exons.saturating_sub(1) as usize {
            if tr_mult[i_tr].canon_sj[iex] >= 0 && tr_mult[i_tr].sj_annot[iex] == 0 {
                *out_filter_by_sjout_pass = false;
                break;
            }
        }
        if !*out_filter_by_sjout_pass {
            break;
        }
    }

    if !*out_filter_by_sjout_pass {
        *unmap_type = -3;
        stats_ra.read_n = stats_ra.read_n.wrapping_sub(1);
        stats_ra.read_bases = stats_ra
            .read_bases
            .wrapping_sub(read_length[0].wrapping_add(read_length[1]));

        for im in 0..read_nends as usize {
            chunk_out_filter_by_sjout_files[im].push_str(&format!(
                "{} {} {} {}",
                read_name_mates[im], i_read_all, read_filter, read_files_index
            ));
            if !read_name_extra[im].is_empty() {
                chunk_out_filter_by_sjout_files[im].push(' ');
                chunk_out_filter_by_sjout_files[im].push_str(&read_name_extra[im]);
            }
            chunk_out_filter_by_sjout_files[im].push('\n');
            chunk_out_filter_by_sjout_files[im].push_str(&read0[im]);
            chunk_out_filter_by_sjout_files[im].push('\n');
            if read_file_type == 2 {
                chunk_out_filter_by_sjout_files[im].push_str("+\n");
                chunk_out_filter_by_sjout_files[im].push_str(&qual0[im]);
                chunk_out_filter_by_sjout_files[im].push('\n');
            }
        }
    }

    readalign_outputalignments_l76_readalign_recordsj(
        n_tr,
        tr_mult,
        chunk_out_sj1,
        out_sj_yes,
        out_sj_filter_reads,
    );
}

#[doc = "Original `ReadAlign::writeSAM` at STAR/source/ReadAlign_outputAlignments.cpp:132. Args: nTrOutSAM: uint64, trOutSAM: Transcript, trBestSAM: Transcript"]
pub fn readalign_outputalignments_l132_readalign_writesam(
    mut n_tr_out_sam: u64,
    tr_out_sam: &[crate::transcript::Transcript],
    tr_best_sam: &crate::transcript::Transcript,
    read_align: &crate::read_align::ReadAlign,
    p: &crate::parameters_chimeric::Parameters,
    map_gen: &crate::genome::Genome,
    unmap_type_in: i32,
    out_filter_by_sjout_pass: bool,
    read0: &[String],
    read_file_type: i32,
    qual0: &[String],
    read_name_extra: &[String],
) -> Result<crate::quantifications::WriteSamResult, String> {
    let mut result = crate::quantifications::WriteSamResult {
        unmap_type: unmap_type_in,
        ..Default::default()
    };

    if result.unmap_type < 0 && out_filter_by_sjout_pass {
        let mut tr_out = tr_out_sam.to_vec();

        if p.out_sam_filter_yes {
            if p.out_sam_filter_keep_only_added_references {
                for tr in tr_out.iter().take(n_tr_out_sam as usize) {
                    if tr.chr < p.genome_insert_chr_ind_first {
                        return Ok(result);
                    }
                }
            } else if p.out_sam_filter_keep_all_added_references {
                let mut filtered = Vec::new();
                for mut tr in tr_out.into_iter().take(n_tr_out_sam as usize) {
                    if tr.chr >= p.genome_insert_chr_ind_first {
                        tr.primary_flag = false;
                        filtered.push(tr);
                    }
                }
                if filtered.is_empty() {
                    return Ok(result);
                }
                filtered[0].primary_flag = true;
                n_tr_out_sam = filtered.len() as u64;
                tr_out = filtered;
            }
        }

        let n_tr_out_write = if p.out_sam_mult_nmax == 0 {
            0
        } else {
            n_tr_out_sam.min(p.out_sam_mult_nmax)
        };

        for i_tr in 0..n_tr_out_write as usize {
            let tr = &tr_out[i_tr];
            let mut mate_mapped1 = [false; 2];
            mate_mapped1[tr.exons[0][EX_IFRAG] as usize] = true;
            mate_mapped1[tr.exons[tr.n_exons as usize - 1][EX_IFRAG] as usize] = true;

            if p.out_sam_bool {
                result.out_bam_bytes +=
                    readalign_outputtranscriptsam_l5_readalign_outputtranscriptsam(
                        tr,
                        n_tr_out_sam as u32,
                        i_tr as u32,
                        u32::MAX,
                        u32::MAX,
                        0,
                        -1,
                        None,
                        &mut result.sam,
                        read_align.read_filter as u8,
                        &read_align.read_name,
                        read0,
                        read_file_type,
                        qual0,
                        p,
                        read_align.read_files_index,
                        read_name_extra,
                        read_align.l_read,
                        &read_align.read_length,
                        &read_align.read_length_original,
                        &read_align.clip_mates,
                        &read_align.read1,
                        map_gen,
                    )?;
                if p.out_sam_unmapped_keep_pairs
                    && p.read_nmates > 1
                    && (!mate_mapped1[0] || !mate_mapped1[1])
                {
                    result.out_bam_bytes +=
                        readalign_outputtranscriptsam_l5_readalign_outputtranscriptsam(
                            tr,
                            0,
                            0,
                            u32::MAX,
                            u32::MAX,
                            0,
                            4,
                            Some(&mate_mapped1),
                            &mut result.sam,
                            read_align.read_filter as u8,
                            &read_align.read_name,
                            read0,
                            read_file_type,
                            qual0,
                            p,
                            read_align.read_files_index,
                            read_name_extra,
                            read_align.l_read,
                            &read_align.read_length,
                            &read_align.read_length_original,
                            &read_align.clip_mates,
                            &read_align.read1,
                            map_gen,
                        )?;
                }
            }

            if p.out_bam_unsorted || p.out_bam_coord {
                result
                    .bam_requests
                    .push(crate::quantifications::AlignBamRequest {
                        transcript: tr.clone(),
                        n_tr_out: n_tr_out_sam as u32,
                        i_tr_out: i_tr as u32,
                        tr_chr_start: map_gen.chr_start[tr.chr as usize],
                        mate_chr: u32::MAX,
                        mate_start: u32::MAX,
                        mate_strand: 0,
                        align_type: -1,
                        mate_map: None,
                    });
                if p.out_sam_unmapped_keep_pairs
                    && p.read_nmates > 1
                    && (!mate_mapped1[0] || !mate_mapped1[1])
                {
                    result
                        .bam_requests
                        .push(crate::quantifications::AlignBamRequest {
                            transcript: tr.clone(),
                            n_tr_out: 0,
                            i_tr_out: 0,
                            tr_chr_start: map_gen.chr_start[tr.chr as usize],
                            mate_chr: u32::MAX,
                            mate_start: u32::MAX,
                            mate_strand: 0,
                            align_type: 4,
                            mate_map: Some(mate_mapped1),
                        });
                }
            }
        }

        result.mate_mapped[tr_best_sam.exons[0][EX_IFRAG] as usize] = true;
        result.mate_mapped
            [tr_best_sam.exons[tr_best_sam.n_exons as usize - 1][EX_IFRAG] as usize] = true;

        if p.read_nmates > 1 && !(result.mate_mapped[0] && result.mate_mapped[1]) {
            result.unmap_type = 4;
        }

        if result.unmap_type == 4 && p.out_sam_unmapped_within {
            if p.out_sam_bool && !p.out_sam_unmapped_keep_pairs {
                result.out_bam_bytes +=
                    readalign_outputtranscriptsam_l5_readalign_outputtranscriptsam(
                        tr_best_sam,
                        0,
                        0,
                        u32::MAX,
                        u32::MAX,
                        0,
                        result.unmap_type,
                        Some(&result.mate_mapped),
                        &mut result.sam,
                        read_align.read_filter as u8,
                        &read_align.read_name,
                        read0,
                        read_file_type,
                        qual0,
                        p,
                        read_align.read_files_index,
                        read_name_extra,
                        read_align.l_read,
                        &read_align.read_length,
                        &read_align.read_length_original,
                        &read_align.clip_mates,
                        &read_align.read1,
                        map_gen,
                    )?;
            }

            if p.out_bam_coord || (p.out_bam_unsorted && !p.out_sam_unmapped_keep_pairs) {
                result
                    .bam_requests
                    .push(crate::quantifications::AlignBamRequest {
                        transcript: tr_best_sam.clone(),
                        n_tr_out: 0,
                        i_tr_out: 0,
                        tr_chr_start: map_gen.chr_start[tr_best_sam.chr as usize],
                        mate_chr: u32::MAX,
                        mate_start: u32::MAX,
                        mate_strand: 0,
                        align_type: result.unmap_type,
                        mate_map: Some(result.mate_mapped),
                    });
            }
        }
    } else if result.unmap_type >= 0 && p.out_sam_unmapped_within {
        if p.out_bam_coord || p.out_bam_unsorted || p.quant_tr_sam_bam_yes {
            result
                .bam_requests
                .push(crate::quantifications::AlignBamRequest {
                    transcript: tr_best_sam.clone(),
                    n_tr_out: 0,
                    i_tr_out: 0,
                    tr_chr_start: map_gen.chr_start[tr_best_sam.chr as usize],
                    mate_chr: u32::MAX,
                    mate_start: u32::MAX,
                    mate_strand: 0,
                    align_type: result.unmap_type,
                    mate_map: Some(result.mate_mapped),
                });
        }

        if p.out_sam_bool {
            result.out_bam_bytes += readalign_outputtranscriptsam_l5_readalign_outputtranscriptsam(
                tr_best_sam,
                0,
                0,
                u32::MAX,
                u32::MAX,
                0,
                result.unmap_type,
                Some(&result.mate_mapped),
                &mut result.sam,
                read_align.read_filter as u8,
                &read_align.read_name,
                read0,
                read_file_type,
                qual0,
                p,
                read_align.read_files_index,
                read_name_extra,
                read_align.l_read,
                &read_align.read_length,
                &read_align.read_length_original,
                &read_align.clip_mates,
                &read_align.read1,
                map_gen,
            )?;
        }
    }

    Ok(result)
}

#[doc = "Original `ReadAlign::outReadsUnmapped` at STAR/source/ReadAlign_outputAlignments.cpp:259. Args: "]
pub fn readalign_outputalignments_l259_readalign_outreadsunmapped(
    out_reads_unmapped: &str,
    read_nends: u32,
    read_nmates: u32,
    read_name_mates: &[String],
    read_filter: i32,
    read_name_extra: &[String],
    mate_mapped: &[bool; 2],
    read0: &[String],
    read_file_type: i32,
    qual0: &[String],
    chunk_out_unmapped_reads_stream: &mut [String],
) {
    if out_reads_unmapped == "Fastx" {
        for im in 0..read_nends as usize {
            chunk_out_unmapped_reads_stream[im].push_str(&read_name_mates[im]);
            chunk_out_unmapped_reads_stream[im]
                .push_str(&format!(" {}:{}: {}", im, read_filter, read_name_extra[im]));
            if read_nmates > 1 {
                chunk_out_unmapped_reads_stream[im].push_str(&format!(
                    " {}{}",
                    i32::from(mate_mapped[0]),
                    i32::from(mate_mapped[1])
                ));
            }
            chunk_out_unmapped_reads_stream[im].push('\n');
            chunk_out_unmapped_reads_stream[im].push_str(&read0[im]);
            chunk_out_unmapped_reads_stream[im].push('\n');
            if read_file_type == 2 {
                chunk_out_unmapped_reads_stream[im].push_str("+\n");
                chunk_out_unmapped_reads_stream[im].push_str(&qual0[im]);
                chunk_out_unmapped_reads_stream[im].push('\n');
            }
        }
    }
}

#[doc = "Original `ReadAlign::spliceGraphWriteSAM` at STAR/source/ReadAlign_outputAlignments.cpp:277. Args: "]
pub fn readalign_outputalignments_l277_readalign_splicegraphwritesam(
    tr_mult: &[crate::transcript::Transcript],
    n_tr: u64,
    out_sam_stream: &mut String,
    read_filter: u8,
    unmap_type: i32,
    read_name: &str,
    read0: &[String],
    read_file_type: i32,
    qual0: &[String],
    p: &crate::parameters_chimeric::Parameters,
    read_files_index: u32,
    read_name_extra: &[String],
    l_read: u32,
    map_gen: &crate::genome::Genome,
) -> Result<u64, String> {
    let mut out_bam_bytes = 0_u64;
    let mut n_tr_out_sam = n_tr;
    if map_gen.genome_out.conv_yes {
        n_tr_out_sam = 0;
        for _i_tr in 0..n_tr_out_sam as usize {}
    }

    for i_tr in 0..n_tr_out_sam as usize {
        let (sam, bytes) = readalign_outputsplicegraphsam_l5_readalign_outputsplicegraphsam(
            &tr_mult[i_tr],
            n_tr_out_sam as u32,
            i_tr as u32,
            read_filter,
            unmap_type,
            read_name,
            read0,
            read_file_type,
            qual0,
            p,
            read_files_index,
            read_name_extra,
            l_read,
            map_gen,
        )?;
        out_sam_stream.push_str(&sam);
        out_bam_bytes += bytes;
    }

    Ok(out_bam_bytes)
}

#[doc = "Original `ReadAlign::alignedAnnotation` at STAR/source/ReadAlign_outputAlignments.cpp:298. Args: "]
pub fn readalign_outputalignments_l298_readalign_alignedannotation(
    transcriptome: &mut crate::transcriptome::Transcriptome,
    p: &crate::parameters_chimeric::Parameters,
    n_tr: u64,
    tr_mult: &[crate::transcript::Transcript],
    aligns_gen_out: &crate::quantifications::ReadAlignGenomeTransformResult,
    read_annot: &mut crate::read_annotations::ReadAnnotations,
) {
    let min_features = (SOLO_FEATURE_VELOCYTO as usize) + 1;
    if read_annot.annot_features.len() < min_features {
        read_annot
            .annot_features
            .resize(min_features, Default::default());
    }

    if p.quant_ge_count_yes {
        if p.p_ge.transform.out_quant {
            transcriptome_genecountsaddalign_l4_transcriptome_genecountsaddalign(
                transcriptome,
                aligns_gen_out.al_n,
                &aligns_gen_out.al_mult,
                &mut read_annot.gene_exon_overlap,
            );
        } else {
            transcriptome_genecountsaddalign_l4_transcriptome_genecountsaddalign(
                transcriptome,
                n_tr as u32,
                tr_mult,
                &mut read_annot.gene_exon_overlap,
            );
        }
    }

    if p.quant_gene_full_yes {
        transcriptome_genefullalignoverlap_l5_transcriptome_genefullalignoverlap(
            transcriptome,
            n_tr as u32,
            tr_mult,
            p.p_solo.strand,
            &mut read_annot.annot_features[SOLO_FEATURE_GENE_FULL as usize],
        );
    }

    if p.quant_gene_yes {
        transcriptome_classifyalign_l177_transcriptome_classifyalign(
            transcriptome,
            &p.p_solo,
            tr_mult,
            read_annot,
        );
    }

    if p.quant_gene_full_exon_over_intron_yes {
        let gene_concordant = read_annot.annot_features[SOLO_FEATURE_GENE as usize].clone();
        transcriptome_genefullalignoverlap_exonoverintron_l5_transcriptome_genefullalignoverlap_exonoverintron(
            transcriptome,
            n_tr as u32,
            tr_mult,
            p.p_solo.strand,
            &mut read_annot.annot_features[SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON as usize],
            &gene_concordant,
        );
    }

    if p.quant_gene_full_ex50p_as_yes {
        transcriptome_alignexonoverlap_l10_transcriptome_alignexonoverlap(
            transcriptome,
            n_tr as u32,
            tr_mult,
            p.p_solo.strand,
            &mut read_annot.annot_features[SOLO_FEATURE_GENE_FULL_EX50P_AS as usize],
        );
    }
}
