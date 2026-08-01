#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SpliceGraph::findSuperTr` at STAR/source/SpliceGraph_findSuperTr.cpp:5. Args: readSeq: char, readSeqRevCompl: char, readLen: uint32, readName: string, mapGen: Genome"]
pub fn splicegraph_findsupertr_l5_splicegraph_findsupertr(
    splice_graph: &mut crate::splice_graph::SpliceGraph,
    read_seq: &[u8],
    read_seq_rev_compl: &[u8],
    read_len: u64,
    read_name: &str,
    map_gen: &crate::genome::Genome,
    seed_multimap_nmax: u32,
) -> String {
    let seed_coverage_threshold = 0.02_f32;
    let seed_coverage_min_to_max = 0.5_f32;
    let seed_mult_max = seed_multimap_nmax;
    let seed_spacing = 1_u32;
    let seed_len = map_gen.p_ge.g_saindex_nbases;
    let mut out = String::new();

    let mut seed_super_tr = Vec::with_capacity(seed_mult_max as usize);
    splice_graph
        .super_tr_seed_count
        .resize(2 * splice_graph.super_trome.n as usize, 0);
    splice_graph.super_tr_seed_count.fill(0);

    let mut iseed = 0_u64;
    while iseed < read_len {
        let mut ind1 = 0_u64;
        for ii in iseed..iseed + seed_len as u64 {
            let b = read_seq[ii as usize] as u64;
            if b > 3 {
                continue;
            } else {
                ind1 <<= 2;
                ind1 += b;
            }
        }

        let sai_index = map_gen.genome_sa_index_start[seed_len as usize - 1] as u64 + ind1;
        let i_sa1 = map_gen.sai_value(sai_index);
        if (i_sa1 & map_gen.sai_mark_absent_mask_c as u64) != 0 {
            iseed += seed_spacing as u64;
            continue;
        }

        let i_sa2 = if sai_index + 1 < map_gen.genome_sa_index_start[seed_len as usize] as u64 {
            ((map_gen.sai_value(sai_index + 1) & map_gen.sai_mark_nmask as u64)
                & !(map_gen.sai_mark_absent_mask_c as u64))
                - 1
        } else {
            map_gen.n_sa as u64 - 1
        };

        if i_sa2 - i_sa1 >= seed_mult_max as u64 {
            iseed += seed_spacing as u64;
            continue;
        }

        seed_super_tr.clear();
        seed_super_tr.resize((i_sa2 - i_sa1 + 1) as usize, 0);
        for isa in i_sa1..=i_sa2 {
            let mut a1 = genome_sa_index_value(map_gen, isa);
            let a_str = a1 >> map_gen.gstrand_bit;
            a1 &= map_gen.gstrand_mask as u64;
            if a_str == 1 {
                a1 = map_gen.n_genome - (seed_len as u64 + a1);
            }

            if a1 >= map_gen.sj_gstart {
                let mut a1_d = 0_u64;
                let mut a_length_d = 0_u64;
                let mut a1_a = 0_u64;
                let mut a_length_a = 0_u64;
                let mut sj1 = 0_u64;
                if sjalignsplit_l3_sjalignsplit(
                    a1,
                    seed_len as u64,
                    map_gen,
                    &mut a1_d,
                    &mut a_length_d,
                    &mut a1_a,
                    &mut a_length_a,
                    &mut sj1,
                ) {
                    a1 = a1_d;
                } else {
                    continue;
                }
            }

            seed_super_tr[(isa - i_sa1) as usize] = (a_str as u32 * splice_graph.super_trome.n)
                + map_gen.chr_bin[(a1 >> map_gen.p_ge.g_chr_bin_nbits) as usize];
        }

        seed_super_tr.sort_unstable();
        let mut su1prev = u32::MAX;
        for su1 in &seed_super_tr {
            if *su1 != su1prev {
                splice_graph.super_tr_seed_count[*su1 as usize] += 1;
                su1prev = *su1;
            }
        }

        iseed += seed_spacing as u64;
    }

    let mut count_max = 0_u16;
    for ii in 0..2 * splice_graph.super_trome.n as usize {
        count_max = count_max.max(splice_graph.super_tr_seed_count[ii]);
    }

    if (count_max as f32) < read_len as f32 * seed_coverage_threshold / seed_spacing as f32 {
        return out;
    }

    let mut n_super_tr = 0_usize;
    let mut max_max_score = 0_i32;
    for ii in 0..splice_graph.super_trome.n as usize {
        let sutr1 = ii % splice_graph.super_trome.n as usize;
        let str1 = ii / splice_graph.super_trome.n as usize;
        let seed_count = splice_graph.super_tr_seed_count[ii];

        if (seed_count as f32) < read_len as f32 * seed_coverage_threshold / seed_spacing as f32
            || (seed_count as f32) < count_max as f32 * seed_coverage_min_to_max
        {
            continue;
        }

        if read_len >= 100000
            || splice_graph.super_trome.super_trs[sutr1].length as u64 * read_len as u64
                >= 1_000_000_000
        {
            continue;
        }

        let read_for_score = if str1 == 0 {
            read_seq
        } else {
            read_seq_rev_compl
        };
        let super_tr = splice_graph.super_trome.super_trs[sutr1].clone();
        let mut tr_a = splice_graph
            .ra
            .as_ref()
            .map(|ra| ra.tr_best.clone())
            .unwrap_or_else(transcript_l3_transcript_transcript);
        tr_a.read_name = read_name.to_string();
        tr_a.chr = sutr1 as u64;
        tr_a.str_ = str1 as u64;
        let sw_score = splicegraph_swscorespliced_l8_splicegraph_swscorespliced(
            splice_graph,
            read_for_score,
            read_len as u32,
            &super_tr,
            &mut tr_a.cigar,
        );

        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            read_name,
            sutr1,
            str1,
            super_tr.length,
            read_len,
            seed_count as f32 / read_len as f32 * seed_spacing as f32,
            seed_count as f32 / count_max as f32,
            sw_score,
            sw_score as f32 / read_len as f32,
            splice_graph.align_info.a_start[0],
            splice_graph.align_info.a_end[0],
            splice_graph.align_info.a_start[1],
            splice_graph.align_info.a_end[1]
        ));

        tr_a.max_score = sw_score;
        tr_a.n_match = sw_score as u64;
        tr_a.n_exons = 0;
        tr_a.g_start =
            map_gen.chr_start[tr_a.chr as usize] + splice_graph.align_info.a_start[1] as u64;
        tr_a.n_mm = splice_graph.align_info.n_mm as u64;
        tr_a.l_ins = splice_graph.align_info.n_i as u64;
        tr_a.l_del = splice_graph.align_info.n_d as u64;
        tr_a.r_length = splice_graph.align_info.n_map as u64;

        if let Some(ra) = splice_graph.ra.as_mut() {
            if ra.tr_all.len() <= n_super_tr {
                ra.tr_all.resize(n_super_tr + 1, Vec::new());
            }
            ra.tr_all[n_super_tr] = vec![tr_a.clone()];
            if ra.n_win_tr.len() <= n_super_tr {
                ra.n_win_tr.resize(n_super_tr + 1, 0);
            }
            ra.n_win_tr[n_super_tr] = 1;
            if sw_score > max_max_score {
                max_max_score = sw_score;
                ra.tr_best = tr_a;
            }
        }

        n_super_tr += 1;
    }

    if let Some(ra) = splice_graph.ra.as_mut() {
        ra.n_w = n_super_tr as u64;
    }

    out
}

pub fn format_local_time_month_day_time(raw_time: libc::time_t) -> String {
    format_local_time_month_day_time_impl(raw_time)
}
