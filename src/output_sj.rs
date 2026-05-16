#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `compareUint` at STAR/source/outputSJ.cpp:7. Args: i1: void, i2: void"]
pub fn outputsj_l7_compareuint(i1: u32, i2: u32) -> i32 {
    if i1 > i2 {
        1
    } else if i1 < i2 {
        -1
    } else {
        0
    }
}

#[doc = "Original `outputSJ` at STAR/source/outputSJ.cpp:20. Args: RAchunk: ReadAlignChunk, P: Parameters"]
pub fn outputsj_l20_outputsj(
    ra_chunk: &[crate::read_align_chunk::ReadAlignChunk],
    p: &mut crate::parameters_chimeric::Parameters,
    gen_out: &crate::genome::Genome,
) -> Result<crate::parameters_chimeric::OutputSjResult, String> {
    use std::fmt::Write;

    let mut result = crate::parameters_chimeric::OutputSjResult::default();
    let mut one_sj = crate::out_sj::Junction {
        gen_out: genome_output_sj_snapshot(gen_out),
        ..Default::default()
    };
    let mut all_sj = outsj_l4_outsj_outsj(p.limit_out_sj_collapsed * 2);
    let run_thread_n = p.run_thread_n as usize;
    let mut sj_chunks =
        Vec::<Vec<crate::out_sj::JunctionRecord>>::with_capacity(run_thread_n);
    let mut sj_index = vec![0usize; run_thread_n];

    for chunk in ra_chunk.iter().take(run_thread_n) {
        let out_sj = if p.out_filter_by_sjout_stage != 1 {
            &chunk.chunk_out_sj
        } else {
            &chunk.chunk_out_sj1
        };
        sj_chunks.push(
            out_sj
                .junctions
                .iter()
                .take(out_sj.n as usize)
                .cloned()
                .collect(),
        );
    }

    loop {
        let mut ic_out: Option<usize> = None;
        for ic in 0..run_thread_n {
            if sj_index[ic] < sj_chunks[ic].len()
                && ic_out.map_or(true, |ic_prev| {
                    sj_chunks[ic][sj_index[ic]]
                        .start
                        .cmp(&sj_chunks[ic_prev][sj_index[ic_prev]].start)
                        .then_with(|| {
                            sj_chunks[ic][sj_index[ic]]
                                .gap
                                .cmp(&sj_chunks[ic_prev][sj_index[ic_prev]].gap)
                        })
                        .is_lt()
                })
            {
                ic_out = Some(ic);
            }
        }

        let Some(ic_out) = ic_out else {
            break;
        };

        let mut collapsed = sj_chunks[ic_out][sj_index[ic_out]].clone();
        for ic in 0..run_thread_n {
            if ic != ic_out
                && sj_index[ic] < sj_chunks[ic].len()
                && sj_chunks[ic][sj_index[ic]].start == collapsed.start
                && sj_chunks[ic][sj_index[ic]].gap == collapsed.gap
            {
                outsj_l92_junction_collapseonesj(&mut collapsed, &sj_chunks[ic][sj_index[ic]])?;
                sj_index[ic] += 1;
            }
        }

        let motif_index = ((collapsed.motif + 1) / 2) as usize;
        let count_total = collapsed.count_multiple + collapsed.count_unique;
        let intron_ok = count_total as usize > p.out_sjfilter_intron_max_vs_read_n.len()
            || collapsed.gap
                <= p.out_sjfilter_intron_max_vs_read_n[count_total.saturating_sub(1) as usize];
        let sj_filter = collapsed.annot > 0
            || ((collapsed.count_unique as i32 >= p.out_sjfilter_count_unique_min[motif_index]
                || count_total as i32 >= p.out_sjfilter_count_total_min[motif_index])
                && collapsed.overhang_left as i32 >= p.out_sjfilter_overhang_min[motif_index]
                && collapsed.overhang_right as i32 >= p.out_sjfilter_overhang_min[motif_index]
                && intron_ok);

        if sj_filter {
            all_sj.junctions.push(collapsed);
            all_sj.n += 1;
            if all_sj.n == all_sj.n_store - 1 {
                outsj_l62_outsj_datasizeincrease(&mut all_sj);
                writeln!(
                    result.log_main,
                    "Increased the size of chunkOutSJ to {}",
                    all_sj.n_store
                )
                .unwrap();
            }
        }

        sj_index[ic_out] += 1;
    }

    let mut sj_filter = vec![true; all_sj.n as usize];
    if p.out_filter_by_sjout_stage != 2 {
        let mut sj_a = Vec::<(u32, usize, u32)>::with_capacity(all_sj.n as usize);
        for ii in 0..all_sj.n as usize {
            let record = &all_sj.junctions[ii];
            let x1 = if ii > 0 {
                all_sj.junctions[ii - 1].start
            } else {
                0
            };
            let x2 = if ii + 1 < all_sj.n as usize {
                all_sj.junctions[ii + 1].start
            } else {
                u32::MAX
            };
            let min_dist = record
                .start
                .wrapping_sub(x1)
                .min(x2.wrapping_sub(record.start));
            let motif_index = ((record.motif + 1) / 2) as usize;
            sj_filter[ii] = min_dist >= p.out_sjfilter_dist_to_other_sj_min[motif_index] as u32;
            let motif_for_acceptor = if record.annot == 0 {
                record.motif as u32
            } else {
                SJ_MOTIF_SIZE as u32 + 1
            };
            sj_a.push((
                record.start.wrapping_add(record.gap),
                ii,
                motif_for_acceptor,
            ));
        }
        sj_a.sort_by(|a, b| a.0.cmp(&b.0));
        for ii in 0..sj_a.len() {
            let (_, original_index, motif_for_acceptor) = sj_a[ii];
            if motif_for_acceptor == SJ_MOTIF_SIZE as u32 + 1 {
                sj_filter[original_index] = true;
            } else {
                let x1 = if ii > 0 { sj_a[ii - 1].0 } else { 0 };
                let x2 = if ii + 1 < sj_a.len() {
                    sj_a[ii + 1].0
                } else {
                    u32::MAX
                };
                let min_dist = sj_a[ii].0.wrapping_sub(x1).min(x2.wrapping_sub(sj_a[ii].0));
                let motif_index = ((motif_for_acceptor + 1) / 2) as usize;
                sj_filter[original_index] = sj_filter[original_index]
                    && min_dist >= p.out_sjfilter_dist_to_other_sj_min[motif_index] as u32;
            }
        }
    }

    p.sj_all[0].reserve(all_sj.n as usize);
    p.sj_all[1].reserve(all_sj.n as usize);

    if p.out_filter_by_sjout_stage != 1 {
        for ii in 0..all_sj.n as usize {
            if p.out_filter_by_sjout_stage == 2 || sj_filter[ii] {
                outsj_l72_junction_junctionpointer(&mut one_sj, &all_sj.junctions, ii as u32)?;
                result
                    .sj_out_tab
                    .push_str(&outsj_l85_junction_outputstream(&one_sj)?);
                writeln!(
                    result.sj_start_gap_tsv,
                    "{}\t{}",
                    all_sj.junctions[ii].start, all_sj.junctions[ii].gap
                )
                .unwrap();
                p.sj_all[0].push(all_sj.junctions[ii].start as u64);
                p.sj_all[1].push(all_sj.junctions[ii].gap as u64);
            }
        }
    } else {
        p.sj_novel_n = 0;
        for ii in 0..all_sj.n as usize {
            if sj_filter[ii] && all_sj.junctions[ii].annot == 0 {
                p.sj_novel_n += 1;
            }
        }
        writeln!(
            result.log_main,
            "Detected {} novel junctions that passed filtering, will proceed to filter reads that contained unannotated junctions",
            p.sj_novel_n
        )
        .unwrap();
        p.sj_novel_start.clear();
        p.sj_novel_end.clear();
        p.sj_novel_start.reserve(p.sj_novel_n as usize);
        p.sj_novel_end.reserve(p.sj_novel_n as usize);
        for ii in 0..all_sj.n as usize {
            if sj_filter[ii] && all_sj.junctions[ii].annot == 0 {
                p.sj_novel_start.push(all_sj.junctions[ii].start);
                p.sj_novel_end.push(
                    all_sj.junctions[ii]
                        .start
                        .wrapping_add(all_sj.junctions[ii].gap)
                        .wrapping_sub(1),
                );
            }
        }
    }

    Ok(result)
}

pub fn genome_output_sj_snapshot(
    genome: &crate::genome::Genome,
) -> crate::genome::Genome {
    let mut p_ge = crate::parameters_genome::ParametersGenome::default();
    p_ge.g_chr_bin_nbits = genome.p_ge.g_chr_bin_nbits;
    crate::genome::Genome {
        p_ge,
        chr_bin: genome.chr_bin.clone(),
        chr_name: genome.chr_name.clone(),
        chr_start: genome.chr_start.clone(),
        ..Default::default()
    }
}
