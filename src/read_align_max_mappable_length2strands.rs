#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::maxMappableLength2strands` at STAR/source/ReadAlign_maxMappableLength2strands.cpp:5. Args: pieceStartIn: uint, pieceLengthIn: uint, iDir: uint, iSA1: uint, iSA2: uint, maxLbest: uint, iFrag: uint"]
pub fn readalign_maxmappablelength2strands_l5_readalign_maxmappablelength2strands(
    map_gen: &crate::genome::Genome,
    p_ge: &crate::parameters_genome::ParametersGenome,
    read1: [&[u8]; 2],
    piece_start_in: u64,
    piece_length_in: u64,
    i_dir: u64,
    i_sa1: u64,
    _i_sa2: u64,
    max_lbest: &mut u64,
    i_frag: u64,
    stored_aligns: &mut Vec<crate::read_align::StoredAlign>,
) -> Result<u64, String> {
    let mut nrep: u64 = 0;
    let mut ind_start_end = [0u64; 2];
    let sparse_d = (p_ge.g_sasparse_d as u64).min(piece_length_in);
    let mut nrep_all = vec![0u64; sparse_d as usize];
    let mut ind_start_end_all = vec![[0u64; 2]; sparse_d as usize];
    let mut max_l_all = vec![0u64; sparse_d as usize];
    *max_lbest = 0;

    let dir_r = i_dir == 0;

    for i_dist in 0..sparse_d {
        let piece_length = piece_length_in - i_dist;
        let l_max = (p_ge.g_saindex_nbases as u64).min(piece_length);
        let mut ind1 = 0u32;
        let piece_start;
        if dir_r {
            piece_start = piece_start_in + i_dist;
            for ii in 0..l_max {
                ind1 <<= 2;
                ind1 += read1[0][(piece_start + ii) as usize] as u32;
            }
        } else {
            piece_start = piece_start_in - i_dist;
            for ii in 0..l_max {
                ind1 <<= 2;
                ind1 += 3 - read1[0][(piece_start - ii) as usize] as u32;
            }
        }

        let mut i_sa1_raw = i_sa1 as u64;
        let mut i_sa2_raw: u64;
        let mut l_ind = l_max;
        let mut i_sa2good = true;
        if l_ind == 0 {
            i_sa1_raw = 0;
            i_sa2_raw = map_gen.n_sa.saturating_sub(1) as u64;
        } else {
            while l_ind > 0 {
                i_sa1_raw = map_gen.sai
                    [(map_gen.genome_sa_index_start[(l_ind - 1) as usize] + ind1 as u64) as usize];
                if (i_sa1_raw & map_gen.sai_mark_absent_mask_c) == 0 {
                    break;
                }
                l_ind -= 1;
                ind1 >>= 2;
            }

            if l_ind == 0 {
                i_sa1_raw = 0;
                i_sa2_raw = map_gen.n_sa.saturating_sub(1) as u64;
            } else if map_gen.genome_sa_index_start[(l_ind - 1) as usize] + ind1 as u64 + 1
                < map_gen.genome_sa_index_start[l_ind as usize]
            {
                i_sa2_raw = map_gen.sai
                    [(map_gen.genome_sa_index_start[(l_ind - 1) as usize] + ind1 as u64 + 1) as usize];
                if (i_sa2_raw & map_gen.sai_mark_absent_mask_c) == 0 {
                    i_sa2_raw = (i_sa2_raw & map_gen.sai_mark_nmask).saturating_sub(1);
                } else {
                    i_sa2_raw = (map_gen.n_sa - 1) as u64;
                    i_sa2good = false;
                }
            } else {
                i_sa2_raw = (map_gen.n_sa - 1) as u64;
                i_sa2good = false;
            };
        }

        let mut i_sa1_search = i_sa1_raw & map_gen.sai_mark_nmask;
        let mut i_sa2_search = i_sa2_raw & map_gen.sai_mark_nmask;
        if map_gen.sai_mark_absent_mask_c != 0 {
            i_sa1_search &= !map_gen.sai_mark_absent_mask_c;
            i_sa2_search &= !map_gen.sai_mark_absent_mask_c;
        }
        let search_bounds_valid =
            i_sa1_search < map_gen.n_sa && i_sa2_search < map_gen.n_sa;
        if !search_bounds_valid || i_sa1_search > i_sa2_search {
            i_sa1_search = 0;
            i_sa2_search = map_gen.n_sa.saturating_sub(1);
        }

        let i_sa1no_n = (i_sa1_raw & map_gen.sai_mark_nmask_c) == 0;
        let max_l;
        if l_ind < p_ge.g_saindex_nbases as u64 && i_sa1no_n && i_sa2good && search_bounds_valid {
            ind_start_end[0] = i_sa1_search;
            ind_start_end[1] = i_sa2_search;
            nrep = ind_start_end[1] - ind_start_end[0] + 1;
            max_l = l_ind;
        } else if i_sa1_search == i_sa2_search && i_sa1no_n && i_sa2good && search_bounds_valid {
            if (i_sa1_raw & map_gen.sai_mark_nmask_c) != 0 {
                return Err("BUG: in ReadAlign::maxMappableLength2strands".to_string());
            }
            ind_start_end[0] = i_sa1_search;
            ind_start_end[1] = i_sa1_search;
            nrep = 1;
            let mut compar_res = false;
            max_l = suffixarrayfuns_l10_compareseqtogenome(
                map_gen,
                read1,
                piece_start,
                piece_length,
                l_ind,
                i_sa1_search,
                dir_r,
                &mut compar_res,
            );
        } else {
            let mut max_l_in = if i_sa2good && i_sa1no_n && search_bounds_valid {
                l_ind
            } else {
                0
            };
            nrep = suffixarrayfuns_l133_maxmappablelength(
                map_gen,
                read1,
                piece_start,
                piece_length,
                i_sa1_search,
                i_sa2_search,
                dir_r,
                &mut max_l_in,
                &mut ind_start_end,
            );
            max_l = max_l_in;
        }
        if max_l + i_dist > *max_lbest {
            *max_lbest = max_l + i_dist;
        }
        nrep_all[i_dist as usize] = nrep;
        ind_start_end_all[i_dist as usize] = ind_start_end;
        max_l_all[i_dist as usize] = max_l;
    }

    for i_dist in 0..sparse_d {
        if max_l_all[i_dist as usize] + i_dist == *max_lbest {
            stored_aligns.push(crate::read_align::StoredAlign {
                i_dir,
                shift: if dir_r {
                    piece_start_in + i_dist
                } else {
                    piece_start_in - i_dist
                },
                n_rep: nrep_all[i_dist as usize],
                l: max_l_all[i_dist as usize],
                ind_start_end: ind_start_end_all[i_dist as usize],
                i_frag,
            });
        }
    }
    Ok(nrep)
}
