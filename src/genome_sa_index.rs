#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `genomeSAindex` at STAR/source/genomeSAindex.cpp:6. Args: G: char, SA: PackedArray, P: Parameters, SAi: PackedArray, mapGen: Genome"]
pub fn genomesaindex_l6_genomesaindex(
    g: &[u8],
    sa: &[u64],
    map_gen: &mut crate::genome::Genome,
) -> Result<Vec<u64>, String> {
    map_gen.genome_sa_index_start = vec![0; map_gen.p_ge.g_saindex_nbases as usize + 1];
    for ii in 1..=map_gen.p_ge.g_saindex_nbases as usize {
        map_gen.genome_sa_index_start[ii] =
            map_gen.genome_sa_index_start[ii - 1] + (1u32 << (2 * ii as u32));
    }
    let n_sai = map_gen.genome_sa_index_start[map_gen.p_ge.g_saindex_nbases as usize];

    map_gen.sai_mark_nmask_c = 1u64 << (map_gen.gstrand_bit + 1);
    map_gen.sai_mark_nmask = !map_gen.sai_mark_nmask_c;
    map_gen.sai_mark_absent_mask_c = 1u64 << (map_gen.gstrand_bit + 2);
    if map_gen.n_sa == 0 {
        map_gen.n_sa = sa.len() as u32;
    }

    let mut sai = vec![0u64; n_sai as usize];
    if !sa.is_empty() {
        genomesaindex_l117_genomesaindexchunk(g, sa, &mut sai, 0, sa.len() as u32 - 1, map_gen)?;
    }
    map_gen.sai = sai.clone();
    Ok(sai)
}

pub fn genomesaindex_from_packed(
    map_gen: &mut crate::genome::Genome,
) -> Result<Vec<u64>, String> {
    map_gen.genome_sa_index_start = vec![0; map_gen.p_ge.g_saindex_nbases as usize + 1];
    for ii in 1..=map_gen.p_ge.g_saindex_nbases as usize {
        map_gen.genome_sa_index_start[ii] =
            map_gen.genome_sa_index_start[ii - 1] + (1u32 << (2 * ii as u32));
    }
    let n_sai = map_gen.genome_sa_index_start[map_gen.p_ge.g_saindex_nbases as usize];

    map_gen.sai_mark_nmask_c = 1u64 << (map_gen.gstrand_bit + 1);
    map_gen.sai_mark_nmask = !map_gen.sai_mark_nmask_c;
    map_gen.sai_mark_absent_mask_c = 1u64 << (map_gen.gstrand_bit + 2);

    let mut sai = vec![0u64; n_sai as usize];
    if map_gen.sa_packed.length > 0 {
        genomesaindexchunk_from_packed(&mut sai, 0, map_gen.sa_packed.length as u32 - 1, map_gen)?;
    }
    map_gen.sai = sai.clone();
    Ok(sai)
}

#[doc = "Original `genomeSAindexChunk` at STAR/source/genomeSAindex.cpp:117. Args: G: char, SA: PackedArray, P: Parameters, SAi: PackedArray, iSA1: uint, iSA2: uint, mapGen: Genome"]
pub fn genomesaindex_l117_genomesaindexchunk(
    g: &[u8],
    sa: &[u64],
    sai: &mut [u64],
    i_sa1: u32,
    i_sa2: u32,
    map_gen: &crate::genome::Genome,
) -> Result<(), String> {
    let index_bases = map_gen.p_ge.g_saindex_nbases as usize;
    if index_bases == 0 {
        return Ok(());
    }
    if map_gen.n_sa == 0 {
        return Err("BUG: empty suffix array in genomeSAindexChunk\n".to_string());
    }

    let sai_mark_nmask_c = map_gen.sai_mark_nmask_c;
    let sai_mark_absent_mask_c = map_gen.sai_mark_absent_mask_c;

    let mut ind0 = vec![u32::MAX; index_bases];
    let isa_step = map_gen.n_sa as u64 / (1u64 << (2 * map_gen.p_ge.g_saindex_nbases)) + 1;
    let mut isa = i_sa1 as u64;
    let mut i_l4 = 0;
    let mut ind_full = suffixarrayfuns_l353_funcalcsaifromsa(
        g,
        sa,
        map_gen.n_genome,
        map_gen.gstrand_bit as u8,
        map_gen.gstrand_mask as u64,
        isa,
        map_gen.p_ge.g_saindex_nbases as i32,
        &mut i_l4,
    );

    while isa <= i_sa2 as u64 {
        for i_l in 0..index_bases {
            let shift = 2 * (map_gen.p_ge.g_saindex_nbases as usize - 1 - i_l);
            let ind_pref = (ind_full >> shift) as u32;

            if i_l as i32 == i_l4 {
                for i_l1 in i_l..index_bases {
                    let jj = map_gen.genome_sa_index_start[i_l1].wrapping_add(ind0[i_l1]);
                    let jj_usize = jj as usize;
                    if jj_usize < sai.len() {
                        sai[jj_usize] |= sai_mark_nmask_c;
                    }
                }
                break;
            }

            if ind_pref > ind0[i_l] || isa == 0 {
                let jj = map_gen.genome_sa_index_start[i_l].wrapping_add(ind_pref);
                let jj_usize = jj as usize;
                if jj_usize >= sai.len() {
                    return Err("BUG: SA index write is out of bounds\n".to_string());
                }
                sai[jj_usize] = isa;

                let mut ii = ind0[i_l].wrapping_add(1);
                while ii < ind_pref {
                    let jj_absent = map_gen.genome_sa_index_start[i_l].wrapping_add(ii);
                    let jj_absent_usize = jj_absent as usize;
                    if jj_absent_usize >= sai.len() {
                        return Err("BUG: SA absent-index write is out of bounds\n".to_string());
                    }
                    sai[jj_absent_usize] = isa | sai_mark_absent_mask_c;
                    ii = ii.wrapping_add(1);
                }
                ind0[i_l] = ind_pref;
            } else if ind_pref < ind0[i_l] {
                return Err("BUG: next index is smaller than previous, EXITING\n".to_string());
            }
        }

        genomesaindex_l178_funsaifindnextindex(
            g,
            sa,
            isa_step,
            &mut isa,
            &mut ind_full,
            &mut i_l4,
            map_gen,
        );
    }

    for i_l in 0..index_bases {
        let mut ii = map_gen.genome_sa_index_start[i_l]
            .wrapping_add(ind0[i_l])
            .wrapping_add(1);
        while ii < map_gen.genome_sa_index_start[i_l + 1] {
            let ii_usize = ii as usize;
            if ii_usize >= sai.len() {
                return Err("BUG: SA tail-index write is out of bounds\n".to_string());
            }
            sai[ii_usize] = map_gen.n_sa as u64 | sai_mark_absent_mask_c;
            ii = ii.wrapping_add(1);
        }
    }

    Ok(())
}

pub fn genomesaindexchunk_from_packed(
    sai: &mut [u64],
    i_sa1: u32,
    i_sa2: u32,
    map_gen: &crate::genome::Genome,
) -> Result<(), String> {
    let index_bases = map_gen.p_ge.g_saindex_nbases as usize;
    if index_bases == 0 {
        return Ok(());
    }
    if map_gen.n_sa == 0 {
        return Err("BUG: empty suffix array in genomeSAindexChunk\n".to_string());
    }

    let sai_mark_nmask_c = map_gen.sai_mark_nmask_c;
    let sai_mark_absent_mask_c = map_gen.sai_mark_absent_mask_c;

    let mut ind0 = vec![u32::MAX; index_bases];
    let isa_step = map_gen.n_sa as u64 / (1u64 << (2 * map_gen.p_ge.g_saindex_nbases)) + 1;
    let mut isa = i_sa1 as u64;
    let mut i_l4 = 0;
    let mut ind_full = suffixarrayfuns_l353_funcalcsaifrompacked(
        &map_gen.g,
        &map_gen.sa_packed,
        map_gen.n_genome,
        map_gen.gstrand_bit as u8,
        map_gen.gstrand_mask as u64,
        isa,
        map_gen.p_ge.g_saindex_nbases as i32,
        &mut i_l4,
    );

    while isa <= i_sa2 as u64 {
        for i_l in 0..index_bases {
            let shift = 2 * (map_gen.p_ge.g_saindex_nbases as usize - 1 - i_l);
            let ind_pref = (ind_full >> shift) as u32;

            if i_l as i32 == i_l4 {
                for i_l1 in i_l..index_bases {
                    let jj = map_gen.genome_sa_index_start[i_l1].wrapping_add(ind0[i_l1]);
                    let jj_usize = jj as usize;
                    if jj_usize < sai.len() {
                        sai[jj_usize] |= sai_mark_nmask_c;
                    }
                }
                break;
            }

            if ind_pref > ind0[i_l] || isa == 0 {
                let jj = map_gen.genome_sa_index_start[i_l].wrapping_add(ind_pref);
                let jj_usize = jj as usize;
                if jj_usize >= sai.len() {
                    return Err("BUG: SA index write is out of bounds\n".to_string());
                }
                sai[jj_usize] = isa;

                let mut ii = ind0[i_l].wrapping_add(1);
                while ii < ind_pref {
                    let jj_absent = map_gen.genome_sa_index_start[i_l].wrapping_add(ii);
                    let jj_absent_usize = jj_absent as usize;
                    if jj_absent_usize >= sai.len() {
                        return Err("BUG: SA absent-index write is out of bounds\n".to_string());
                    }
                    sai[jj_absent_usize] = isa | sai_mark_absent_mask_c;
                    ii = ii.wrapping_add(1);
                }
                ind0[i_l] = ind_pref;
            } else if ind_pref < ind0[i_l] {
                return Err("BUG: next index is smaller than previous, EXITING\n".to_string());
            }
        }

        funsaifindnextindex_from_packed(
            &map_gen.sa_packed,
            isa_step,
            &mut isa,
            &mut ind_full,
            &mut i_l4,
            map_gen,
        );
    }

    for i_l in 0..index_bases {
        let mut ii = map_gen.genome_sa_index_start[i_l]
            .wrapping_add(ind0[i_l])
            .wrapping_add(1);
        while ii < map_gen.genome_sa_index_start[i_l + 1] {
            let ii_usize = ii as usize;
            if ii_usize >= sai.len() {
                return Err("BUG: SA tail-index write is out of bounds\n".to_string());
            }
            sai[ii_usize] = map_gen.n_sa as u64 | sai_mark_absent_mask_c;
            ii = ii.wrapping_add(1);
        }
    }

    Ok(())
}

#[doc = "Original `funSAiFindNextIndex` at STAR/source/genomeSAindex.cpp:178. Args: G: char, SA: PackedArray, isaStep: uint, isa: uint, indFull: uint, iL4: int, mapGen: Genome"]
pub fn genomesaindex_l178_funsaifindnextindex(
    g: &[u8],
    sa: &[u64],
    isa_step: u64,
    isa: &mut u64,
    ind_full: &mut u64,
    i_l4: &mut i32,
    map_gen: &crate::genome::Genome,
) {
    let ind_full_prev = *ind_full;
    let i_l4_prev = *i_l4;
    *isa += isa_step;
    while *isa < map_gen.n_sa as u64 {
        *ind_full = suffixarrayfuns_l353_funcalcsaifromsa(
            g,
            sa,
            map_gen.n_genome,
            map_gen.gstrand_bit as u8,
            map_gen.gstrand_mask as u64,
            *isa,
            map_gen.p_ge.g_saindex_nbases as i32,
            i_l4,
        );
        if *ind_full != ind_full_prev || *i_l4 != i_l4_prev {
            break;
        }
        *isa += isa_step;
    }

    if *isa >= map_gen.n_sa as u64 {
        *ind_full = suffixarrayfuns_l353_funcalcsaifromsa(
            g,
            sa,
            map_gen.n_genome,
            map_gen.gstrand_bit as u8,
            map_gen.gstrand_mask as u64,
            map_gen.n_sa as u64 - 1,
            map_gen.p_ge.g_saindex_nbases as i32,
            i_l4,
        );
        if *ind_full == ind_full_prev && *i_l4 == i_l4_prev {
            *isa = map_gen.n_sa as u64;
            return;
        }
    }

    let mut i1 = *isa - isa_step;
    let mut i2 = std::cmp::min(*isa, map_gen.n_sa as u64 - 1);
    while i1 + 1 < i2 {
        *isa = i1 / 2 + i2 / 2 + (i1 % 2 + i2 % 2) / 2;
        *ind_full = suffixarrayfuns_l353_funcalcsaifromsa(
            g,
            sa,
            map_gen.n_genome,
            map_gen.gstrand_bit as u8,
            map_gen.gstrand_mask as u64,
            *isa,
            map_gen.p_ge.g_saindex_nbases as i32,
            i_l4,
        );
        if *ind_full == ind_full_prev && *i_l4 == i_l4_prev {
            i1 = *isa;
        } else {
            i2 = *isa;
        }
    }
    if *isa == i1 {
        *isa = i2;
        *ind_full = suffixarrayfuns_l353_funcalcsaifromsa(
            g,
            sa,
            map_gen.n_genome,
            map_gen.gstrand_bit as u8,
            map_gen.gstrand_mask as u64,
            *isa,
            map_gen.p_ge.g_saindex_nbases as i32,
            i_l4,
        );
    }
}

pub fn funsaifindnextindex_from_packed(
    sa: &crate::packed_array::PackedArray,
    isa_step: u64,
    isa: &mut u64,
    ind_full: &mut u64,
    i_l4: &mut i32,
    map_gen: &crate::genome::Genome,
) {
    let ind_full_prev = *ind_full;
    let i_l4_prev = *i_l4;
    *isa += isa_step;
    while *isa < map_gen.n_sa as u64 {
        *ind_full = suffixarrayfuns_l353_funcalcsaifrompacked(
            &map_gen.g,
            sa,
            map_gen.n_genome,
            map_gen.gstrand_bit as u8,
            map_gen.gstrand_mask as u64,
            *isa,
            map_gen.p_ge.g_saindex_nbases as i32,
            i_l4,
        );
        if *ind_full != ind_full_prev || *i_l4 != i_l4_prev {
            break;
        }
        *isa += isa_step;
    }

    if *isa >= map_gen.n_sa as u64 {
        *ind_full = suffixarrayfuns_l353_funcalcsaifrompacked(
            &map_gen.g,
            sa,
            map_gen.n_genome,
            map_gen.gstrand_bit as u8,
            map_gen.gstrand_mask as u64,
            map_gen.n_sa as u64 - 1,
            map_gen.p_ge.g_saindex_nbases as i32,
            i_l4,
        );
        if *ind_full == ind_full_prev && *i_l4 == i_l4_prev {
            *isa = map_gen.n_sa as u64;
            return;
        }
    }

    let mut i1 = *isa - isa_step;
    let mut i2 = std::cmp::min(*isa, map_gen.n_sa as u64 - 1);
    while i1 + 1 < i2 {
        *isa = i1 / 2 + i2 / 2 + (i1 % 2 + i2 % 2) / 2;
        *ind_full = suffixarrayfuns_l353_funcalcsaifrompacked(
            &map_gen.g,
            sa,
            map_gen.n_genome,
            map_gen.gstrand_bit as u8,
            map_gen.gstrand_mask as u64,
            *isa,
            map_gen.p_ge.g_saindex_nbases as i32,
            i_l4,
        );
        if *ind_full == ind_full_prev && *i_l4 == i_l4_prev {
            i1 = *isa;
        } else {
            i2 = *isa;
        }
    }
    if *isa == i1 {
        *isa = i2;
        *ind_full = suffixarrayfuns_l353_funcalcsaifrompacked(
            &map_gen.g,
            sa,
            map_gen.n_genome,
            map_gen.gstrand_bit as u8,
            map_gen.gstrand_mask as u64,
            *isa,
            map_gen.p_ge.g_saindex_nbases as i32,
            i_l4,
        );
    }
}
