#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `sjdbBuildIndex` at STAR/source/sjdbBuildIndex.cpp:16. Args: P: Parameters, Gsj: char, G: char, SA: PackedArray, SA2: PackedArray, SAi: PackedArray, mapGen: Genome, mapGen1: Genome"]
pub fn sjdbbuildindex_l16_sjdbbuildindex(
    p: &crate::parameters_chimeric::Parameters,
    gsj: &mut Vec<u8>,
    g: &mut Vec<u8>,
    map_gen: &mut crate::genome::Genome,
    map_gen1: &crate::genome::Genome,
) -> Result<crate::parameters_chimeric::SjdbBuildIndexResult, String> {
    let mut result = crate::parameters_chimeric::SjdbBuildIndexResult::default();
    if map_gen.sjdb_n == 0 {
        return Ok(result);
    }

    result.log_main.push_str(&format!(
        "{} ..... inserting junctions into the genome indices\n",
        timefunctions_l4_timemonthdaytime()
    ));

    let n_gsj = map_gen.sjdb_length * map_gen.sjdb_n;
    result.n_gsj = n_gsj;
    let needed = n_gsj as usize * 2 + 1;
    if gsj.len() < needed {
        gsj.resize(needed, GENOME_SPACING_CHAR);
    }
    for ii in 1..=map_gen.sjdb_n as usize {
        gsj[ii * map_gen.sjdb_length as usize - 1] = GENOME_SPACING_CHAR;
    }
    gsj[n_gsj as usize * 2] = GENOME_SPACING_CHAR;

    for ii in 0..n_gsj as usize {
        let value = gsj[ii];
        gsj[n_gsj as usize * 2 - 1 - ii] = if value < 4 { 3 - value } else { value };
    }

    let mut g1c = vec![0_u8; n_gsj as usize * 2 + 1];
    sequencefuns_l4_complementseqnumbers(gsj, &mut g1c, (n_gsj * 2 + 1) as u64);

    let mut old_sj_ind = vec![0_u32; map_gen1.sjdb_n as usize];
    let n_indices_sj1 = map_gen.sjdb_length;
    let mut ind_array: Vec<[u64; 2]> = Vec::new();
    let mut sj_new = 0_u64;

    // STAR's sjdbBuildIndex parallelizes this loop with
    // `#pragma omp parallel for schedule(dynamic,1000) reduction(+:sjNew)`.
    // Each outer iteration is independent except for the `sjNew` reduction
    // and an `old_sj_ind[sjdb_ind] = isj1` write that only fires when
    // map_gen1.sjdb_n > 0 (i.e. when re-indexing on top of an existing
    // genome). Collect per-iteration results in isj order and merge.
    let n_outer = 2 * map_gen.sjdb_n as u32;
    let sjdb_n = map_gen.sjdb_n;
    let sjdb_length = map_gen.sjdb_length;
    let n_sa_search_high = map_gen.n_sa.saturating_sub(1);
    let n_threads = std::cmp::max(1, p.run_thread_n as usize);

    type IterOut = (u64, Option<(usize, u32)>, Vec<[u64; 2]>);

    let compute_one = |isj: u32| -> IterOut {
        let isj1 = if isj < sjdb_n { isj } else { 2 * sjdb_n - 1 - isj };
        let sjdb_ind = if map_gen1.sjdb_n == 0 {
            -1
        } else {
            binarysearch2_l3_binarysearch2(
                map_gen.sjdb_start[isj1 as usize],
                map_gen.sjdb_end[isj1 as usize],
                &map_gen1.sjdb_start,
                &map_gen1.sjdb_end,
                map_gen1.sjdb_n as i32,
            )
        };
        let mut entries: Vec<[u64; 2]> = Vec::new();
        if sjdb_ind < 0 {
            for istart1 in 0..n_indices_sj1 {
                let istart = istart1;
                let seq_offset = isj * sjdb_length + istart;
                if gsj[seq_offset as usize] > 3 {
                    continue;
                }
                let mut l = 0_u64;
                let sa_index = suffixarrayfuns_l297_suffixarraysearch1(
                    map_gen,
                    [&gsj[..], &g1c[..]],
                    seq_offset as u64,
                    10000,
                    u64::MAX,
                    true,
                    0,
                    n_sa_search_high,
                    &mut l,
                );
                entries.push([
                    sa_index,
                    (isj as u64) * sjdb_length as u64 + istart as u64,
                ]);
            }
            (1, None, entries)
        } else if (sjdb_ind as usize) < old_sj_ind.len() {
            (0, Some((sjdb_ind as usize, isj1)), entries)
        } else {
            (0, None, entries)
        }
    };

    let per_iter: Vec<IterOut> = if n_threads <= 1 || n_outer <= 1 {
        (0..n_outer).map(compute_one).collect()
    } else {
        use rayon::iter::IntoParallelIterator;
        use rayon::iter::ParallelIterator;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build()
            .map_err(|e| e.to_string())?;
        pool.install(|| (0..n_outer).into_par_iter().map(&compute_one).collect())
    };

    let total_entries: usize = per_iter.iter().map(|r| r.2.len()).sum();
    ind_array.reserve(total_entries);
    for (inc, write, entries) in per_iter.into_iter() {
        sj_new += inc;
        if let Some((sjdb_ind, isj1)) = write {
            old_sj_ind[sjdb_ind] = isj1;
        }
        ind_array.extend(entries);
    }
    sj_new /= 2;
    result.sj_new = sj_new;
    result.log_main.push_str(&format!(
        "{}   Finished SA search: number of new junctions={}, old junctions={}\n",
        timefunctions_l4_timemonthdaytime(),
        sj_new,
        map_gen.sjdb_n as u64 - sj_new
    ));

    ind_array.sort_by(|a, b| {
        let cmp = funcompareuintandsuffixes_l6_funcompareuintandsuffixes(a, b, gsj);
        cmp.cmp(&0)
    });
    let n_ind = ind_array.len() as u32;
    result.n_ind = n_ind;
    result.ind_array = ind_array.clone();
    result.log_main.push_str(&format!(
        "{}   Finished sorting SA indicesL nInd={}\n",
        timefunctions_l4_timemonthdaytime(),
        n_ind
    ));

    map_gen.n_genome = map_gen.chr_start[map_gen.n_chr_real as usize] + n_gsj as u64;
    map_gen.n_sa += n_ind as u64;

    let mut gstrand_bit1 = ((map_gen.n_genome as f64).ln() / 2.0_f64.ln()).floor() as u32 + 1;
    if gstrand_bit1 < 32 {
        gstrand_bit1 = 32;
    }
    result.log_main.push_str(&format!(
        "Genome size with junctions={}  {}   {}\n",
        map_gen.n_genome, map_gen.chr_start[map_gen.n_chr_real as usize], n_gsj
    ));
    result.log_main.push_str(&format!(
        "GstrandBit1={}   GstrandBit={}\n",
        gstrand_bit1, map_gen.gstrand_bit
    ));
    if gstrand_bit1 > map_gen.gstrand_bit {
        return Err("EXITING because of FATAL ERROR: cannot insert junctions on the fly because of strand GstrandBit problem\nSOLUTION: please contact STAR author at https://groups.google.com/forum/#!forum/rna-star\n".to_string());
    }

    packedarray_l8_packedarray_definebits(
        &mut map_gen.sa_insert,
        map_gen.gstrand_bit as u64 + 1,
        map_gen.n_sa as u64,
    );
    packedarray_l31_packedarray_allocatearray(&mut map_gen.sa_insert);

    let n_gsj_new = sj_new * map_gen.sjdb_length as u64;
    let n2bit = 1_u64 << map_gen.gstrand_bit;
    let strand_mask = !n2bit;

    let mut isj = 0_usize;
    let mut isa2 = 0_u64;
    for isa in 0..map_gen1.n_sa as u64 {
        while isj < ind_array.len() && isa == ind_array[isj][0] {
            let mut ind1 = ind_array[isj][1];
            if ind1 < n_gsj as u64 {
                ind1 += map_gen.chr_start[map_gen.n_chr_real as usize];
            } else {
                ind1 = (ind1 - n_gsj as u64) | n2bit;
            }
            packedarray_l17_packedarray_writepacked(&mut map_gen.sa_insert, isa2, ind1);
            isa2 += 1;
            isj += 1;
        }

        let mut ind1 = packedarray_h18_packedarray_index(&map_gen.sa_packed, isa);
        if (ind1 & n2bit) > 0 {
            let mut ind1s = map_gen1.n_genome - (ind1 & strand_mask);
            if ind1s >= map_gen.chr_start[map_gen.n_chr_real as usize] {
                let sj1 = ((ind1s - map_gen.chr_start[map_gen.n_chr_real as usize])
                    / map_gen.sjdb_length as u64) as usize;
                let old = old_sj_ind.get(sj1).copied().unwrap_or(sj1 as u32) as i64 - sj1 as i64;
                ind1s = (ind1s as i64 + old * map_gen.sjdb_length as i64) as u64;
                ind1 = (map_gen.n_genome - ind1s) | n2bit;
            } else {
                ind1 += n_gsj_new;
            }
        } else if ind1 >= map_gen.chr_start[map_gen.n_chr_real as usize] {
            let sj1 = ((ind1 - map_gen.chr_start[map_gen.n_chr_real as usize])
                / map_gen.sjdb_length as u64) as usize;
            let old = old_sj_ind.get(sj1).copied().unwrap_or(sj1 as u32) as i64 - sj1 as i64;
            ind1 = (ind1 as i64 + old * map_gen.sjdb_length as i64) as u64;
        }
        packedarray_l17_packedarray_writepacked(&mut map_gen.sa_insert, isa2, ind1);
        isa2 += 1;
    }

    while isj < ind_array.len() {
        let mut ind1 = ind_array[isj][1];
        if ind1 < n_gsj as u64 {
            ind1 += map_gen.chr_start[map_gen.n_chr_real as usize];
        } else {
            ind1 = (ind1 - n_gsj as u64) | n2bit;
        }
        packedarray_l17_packedarray_writepacked(&mut map_gen.sa_insert, isa2, ind1);
        isa2 += 1;
        isj += 1;
    }
    result.log_main.push_str(&format!(
        "{}   Finished inserting junction indices\n",
        timefunctions_l4_timemonthdaytime()
    ));

    if !map_gen.sai_packed.char_array.is_empty() && !map_gen.genome_sa_index_start.is_empty() {
        for i_l in 0..map_gen.p_ge.g_saindex_nbases {
            let mut i_sj = 0_usize;
            let start = map_gen.genome_sa_index_start[i_l as usize];
            let end = map_gen.genome_sa_index_start[i_l as usize + 1];
            let mut ind0 = start.saturating_sub(1);
            for ii in start..end {
                let i_sa1 = packedarray_h18_packedarray_index(&map_gen.sai_packed, ii as u64);
                let i_sa2 = i_sa1 & map_gen.sai_mark_nmask & !map_gen.sai_mark_absent_mask_c;
                if i_sj < ind_array.len() && (i_sa1 & map_gen.sai_mark_absent_mask_c) > 0 {
                    let i_sj1 = i_sj;
                    let mut ind1 = suffixarrayfuns_l397_funcalcsai(
                        &gsj[ind_array[2 * 0.min(ind_array.len())][1] as usize..],
                        0,
                    );
                    if i_sj < ind_array.len() {
                        ind1 = suffixarrayfuns_l397_funcalcsai(
                            &gsj[ind_array[i_sj][1] as usize..],
                            i_l,
                        );
                    }
                    while i_sj < ind_array.len()
                        && ind1 < (ii - start) as i64
                        && ind_array[i_sj][0].saturating_sub(1) < i_sa2
                    {
                        i_sj += 1;
                        if i_sj < ind_array.len() {
                            ind1 = suffixarrayfuns_l397_funcalcsai(
                                &gsj[ind_array[i_sj][1] as usize..],
                                i_l,
                            );
                        }
                    }
                    if i_sj < ind_array.len() && ind1 == (ii - start) as i64 {
                        let value = ind_array[i_sj][0].saturating_sub(1) + i_sj as u64 + 1;
                        packedarray_l17_packedarray_writepacked(
                            &mut map_gen.sai_packed,
                            ii as u64,
                            value,
                        );
                        for ii0 in ind0 + 1..ii {
                            packedarray_l17_packedarray_writepacked(
                                &mut map_gen.sai_packed,
                                ii0 as u64,
                                value | map_gen.sai_mark_absent_mask_c,
                            );
                        }
                        i_sj += 1;
                        ind0 = ii;
                    } else {
                        i_sj = i_sj1;
                    }
                } else {
                    while i_sj < ind_array.len() && ind_array[i_sj][0].saturating_sub(1) + 1 < i_sa2
                    {
                        i_sj += 1;
                    }
                    while i_sj < ind_array.len()
                        && ind_array[i_sj][0].saturating_sub(1) + 1 == i_sa2
                    {
                        if suffixarrayfuns_l397_funcalcsai(&gsj[ind_array[i_sj][1] as usize..], i_l)
                            >= (ii - start) as i64
                        {
                            break;
                        }
                        i_sj += 1;
                    }
                    packedarray_l17_packedarray_writepacked(
                        &mut map_gen.sai_packed,
                        ii as u64,
                        i_sa1 + i_sj as u64,
                    );
                    for ii0 in ind0 + 1..ii {
                        packedarray_l17_packedarray_writepacked(
                            &mut map_gen.sai_packed,
                            ii0 as u64,
                            (i_sa2 + i_sj as u64) | map_gen.sai_mark_absent_mask_c,
                        );
                    }
                    ind0 = ii;
                }
            }
        }

        for pair in &ind_array {
            let mut ind1 = 0_i64;
            for i_l in 0..map_gen.p_ge.g_saindex_nbases {
                let gg = gsj[pair[1] as usize + i_l as usize] as u32;
                ind1 <<= 2;
                if gg > 3 {
                    for i_l1 in i_l..map_gen.p_ge.g_saindex_nbases {
                        ind1 += 3;
                        let mut ind2 = map_gen.genome_sa_index_start[i_l1 as usize] as i64 + ind1;
                        while ind2 >= 0 {
                            if (packedarray_h18_packedarray_index(&map_gen.sai_packed, ind2 as u64)
                                & map_gen.sai_mark_absent_mask_c)
                                == 0
                            {
                                break;
                            }
                            ind2 -= 1;
                        }
                        if ind2 >= 0 {
                            let old =
                                packedarray_h18_packedarray_index(&map_gen.sai_packed, ind2 as u64);
                            packedarray_l17_packedarray_writepacked(
                                &mut map_gen.sai_packed,
                                ind2 as u64,
                                old | map_gen.sai_mark_nmask_c,
                            );
                        }
                        ind1 <<= 2;
                    }
                    break;
                } else {
                    ind1 += gg as i64;
                }
            }
        }
    }
    result.log_main.push_str(&format!(
        "{}   Finished SAi\n",
        timefunctions_l4_timemonthdaytime()
    ));

    packedarray_l8_packedarray_definebits(
        &mut map_gen.sa_packed,
        map_gen.gstrand_bit as u64 + 1,
        map_gen.n_sa as u64,
    );
    map_gen.sa_packed.char_array = std::mem::take(&mut map_gen.sa_insert.char_array);
    map_gen.sa_packed.array_allocated = true;
    map_gen.sa_insert.array_allocated = false;
    map_gen.n_sa_byte = map_gen.sa_packed.length_byte;
    map_gen.sj_gstart = map_gen.chr_start[map_gen.n_chr_real as usize];
    let g_start = map_gen.chr_start[map_gen.n_chr_real as usize] as usize;
    if g.len() < g_start + n_gsj as usize {
        g.resize(g_start + n_gsj as usize, GENOME_SPACING_CHAR);
    }
    g[g_start..g_start + n_gsj as usize].copy_from_slice(&gsj[..n_gsj as usize]);
    map_gen.g = std::mem::take(g);
    Ok(result)
}
