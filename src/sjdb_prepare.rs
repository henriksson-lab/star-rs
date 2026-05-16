#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `sjdbPrepare` at STAR/source/sjdbPrepare.cpp:5. Args: sjdbLoci: SjdbClass, P: Parameters, nGenomeReal: uint, outDir: string, mapGen: Genome, Gsj: char"]
pub fn sjdbprepare_l5_sjdbprepare(
    sjdb_loci: &crate::sjdb_class::SjdbClass,
    p: &crate::parameters_chimeric::Parameters,
    n_genome_real: u64,
    _out_dir: &str,
    map_gen: &mut crate::genome::Genome,
) -> Result<crate::parameters_chimeric::SjdbPrepareResult, String> {
    let nsj_in = sjdb_loci.chr.len();
    let mut result = crate::parameters_chimeric::SjdbPrepareResult::default();
    if nsj_in == 0 {
        map_gen.sjdb_n = 0;
        result
            .sjdb_info_txt
            .push_str(&format!("0\t{}\n", map_gen.sjdb_overhang));
        return Ok(result);
    }

    let g = &map_gen.g;
    let mut sjdb_s = vec![0u64; nsj_in];
    let mut sjdb_e = vec![0u64; nsj_in];
    let mut sjdb_motif = vec![0u8; nsj_in];
    let mut sjdb_shift_left = vec![0u8; nsj_in];
    let mut sjdb_shift_right = vec![0u8; nsj_in];

    let mut chr_old = String::new();
    let mut i_chr = 0usize;
    for ii in 0..nsj_in {
        if chr_old != sjdb_loci.chr[ii] {
            i_chr = 0;
            while i_chr < map_gen.n_chr_real as usize {
                if sjdb_loci.chr[ii] == map_gen.chr_name[i_chr] {
                    break;
                }
                i_chr += 1;
            }
            if i_chr >= map_gen.n_chr_real as usize {
                return Err(format!(
                    "EXITING because of FATAL error, the sjdb chromosome {} is not found among the genomic chromosomes\nSOLUTION: fix your file(s) --sjdbFileChrStartEnd or --sjdbGTFfile, offending junction:{}\t{}\t{}\n",
                    sjdb_loci.chr[ii], sjdb_loci.chr[ii], sjdb_loci.start[ii], sjdb_loci.end[ii]
                ));
            }
            chr_old = sjdb_loci.chr[ii].clone();
        }

        sjdb_s[ii] = sjdb_loci.start[ii] + map_gen.chr_start[i_chr] - 1;
        sjdb_e[ii] = sjdb_loci.end[ii] + map_gen.chr_start[i_chr] - 1;

        let s = sjdb_s[ii] as usize;
        let e = sjdb_e[ii] as usize;
        sjdb_motif[ii] = if g[s] == 2 && g[s + 1] == 3 && g[e - 1] == 0 && g[e] == 2 {
            1
        } else if g[s] == 1 && g[s + 1] == 3 && g[e - 1] == 0 && g[e] == 1 {
            2
        } else if g[s] == 2 && g[s + 1] == 1 && g[e - 1] == 0 && g[e] == 2 {
            3
        } else if g[s] == 1 && g[s + 1] == 3 && g[e - 1] == 2 && g[e] == 1 {
            4
        } else if g[s] == 0 && g[s + 1] == 3 && g[e - 1] == 0 && g[e] == 1 {
            5
        } else if g[s] == 2 && g[s + 1] == 3 && g[e - 1] == 0 && g[e] == 3 {
            6
        } else {
            0
        };

        let mut jj_l = 0u64;
        while jj_l <= sjdb_s[ii].saturating_sub(1)
            && g[(sjdb_s[ii] - 1 - jj_l) as usize] == g[(sjdb_e[ii] - jj_l) as usize]
            && g[(sjdb_s[ii] - 1 - jj_l) as usize] < 4
            && jj_l < 255
        {
            jj_l += 1;
        }
        sjdb_shift_left[ii] = jj_l as u8;

        let mut jj_r = 0u64;
        while sjdb_s[ii] + jj_r < n_genome_real
            && g[(sjdb_s[ii] + jj_r) as usize] == g[(sjdb_e[ii] + 1 + jj_r) as usize]
            && g[(sjdb_s[ii] + jj_r) as usize] < 4
            && jj_r < 255
        {
            jj_r += 1;
        }
        sjdb_shift_right[ii] = jj_r as u8;

        if jj_r == 255 || jj_l == 255 {
            result.log_main.push_str(&format!(
                "WARNING: long repeat for junction # {} : {} {} {}; left shift = {}; right shift = {}\n",
                ii + 1,
                sjdb_loci.chr[ii],
                sjdb_s[ii] - map_gen.chr_start[i_chr] + 1,
                sjdb_e[ii] - map_gen.chr_start[i_chr] + 1,
                sjdb_shift_left[ii],
                sjdb_shift_right[ii]
            ));
        }

        sjdb_s[ii] -= sjdb_shift_left[ii] as u64;
        sjdb_e[ii] -= sjdb_shift_left[ii] as u64;
    }

    let mut sjdb_sort = Vec::<[u64; 3]>::with_capacity(nsj_in);
    for ii in 0..nsj_in {
        let shift1 = match sjdb_loci.str_[ii] {
            '+' => 0u64,
            '-' => n_genome_real,
            _ => 2 * n_genome_real,
        };
        sjdb_sort.push([sjdb_s[ii] + shift1, sjdb_e[ii] + shift1, ii as u64]);
    }
    sjdb_sort.sort_by_key(|row| (row[0], row[1]));

    let mut keep_i = Vec::<usize>::new();
    for row in &sjdb_sort {
        let isj = row[2] as usize;
        if keep_i.is_empty()
            || sjdb_s[isj] != sjdb_s[*keep_i.last().unwrap()]
            || sjdb_e[isj] != sjdb_e[*keep_i.last().unwrap()]
        {
            keep_i.push(isj);
            continue;
        }

        let isj0 = *keep_i.last().unwrap();
        if sjdb_loci.priority[isj] < sjdb_loci.priority[isj0] {
        } else if sjdb_loci.priority[isj] > sjdb_loci.priority[isj0]
            || (sjdb_motif[isj] > 0 && sjdb_motif[isj0] == 0)
            || ((sjdb_motif[isj] > 0) == (sjdb_motif[isj0] > 0)
                && sjdb_shift_left[isj] < sjdb_shift_left[isj0])
        {
            *keep_i.last_mut().unwrap() = isj;
        }
    }

    sjdb_sort.clear();
    for &isj in &keep_i {
        let left_restore = if sjdb_motif[isj] == 0 {
            0u64
        } else {
            sjdb_shift_left[isj] as u64
        };
        sjdb_sort.push([
            sjdb_s[isj] + left_restore,
            sjdb_e[isj] + left_restore,
            isj as u64,
        ]);
    }
    sjdb_sort.sort_by_key(|row| (row[0], row[1]));

    map_gen.sjdb_start.clear();
    map_gen.sjdb_end.clear();
    map_gen.sjdb_motif.clear();
    map_gen.sjdb_shift_left.clear();
    map_gen.sjdb_shift_right.clear();
    map_gen.sjdb_strand.clear();

    for ii in 0..sjdb_sort.len() {
        let isj = sjdb_sort[ii][2] as usize;
        if !map_gen.sjdb_start.is_empty()
            && *map_gen.sjdb_start.last().unwrap() == sjdb_sort[ii][0]
            && *map_gen.sjdb_end.last().unwrap() == sjdb_sort[ii][1]
        {
            let isj0 = sjdb_sort[ii - 1][2] as usize;
            let last = map_gen.sjdb_start.len() - 1;
            if sjdb_loci.priority[isj] < sjdb_loci.priority[isj0] {
                continue;
            } else if sjdb_loci.priority[isj] > sjdb_loci.priority[isj0] {
                map_gen.sjdb_start.pop();
                map_gen.sjdb_end.pop();
                map_gen.sjdb_motif.pop();
                map_gen.sjdb_shift_left.pop();
                map_gen.sjdb_shift_right.pop();
                map_gen.sjdb_strand.pop();
            } else if map_gen.sjdb_strand[last] > 0 && sjdb_loci.str_[isj] == '.' {
                continue;
            } else if map_gen.sjdb_strand[last] == 0 && sjdb_loci.str_[isj] != '.' {
                map_gen.sjdb_start.pop();
                map_gen.sjdb_end.pop();
                map_gen.sjdb_motif.pop();
                map_gen.sjdb_shift_left.pop();
                map_gen.sjdb_shift_right.pop();
                map_gen.sjdb_strand.pop();
            } else if map_gen.sjdb_motif[last] == 0 && sjdb_motif[isj] == 0 {
                map_gen.sjdb_strand[last] = 0;
                continue;
            } else if (map_gen.sjdb_motif[last] > 0 && sjdb_motif[isj] == 0)
                || (map_gen.sjdb_motif[last] % 2 == (2 - map_gen.sjdb_strand[last]))
            {
                continue;
            } else {
                map_gen.sjdb_start.pop();
                map_gen.sjdb_end.pop();
                map_gen.sjdb_motif.pop();
                map_gen.sjdb_shift_left.pop();
                map_gen.sjdb_shift_right.pop();
                map_gen.sjdb_strand.pop();
            }
        }

        map_gen.sjdb_start.push(sjdb_sort[ii][0]);
        map_gen.sjdb_end.push(sjdb_sort[ii][1]);
        map_gen.sjdb_motif.push(sjdb_motif[isj]);
        map_gen.sjdb_shift_left.push(sjdb_shift_left[isj]);
        map_gen.sjdb_shift_right.push(sjdb_shift_right[isj]);
        let strand = if sjdb_loci.str_[isj] == '+' {
            1
        } else if sjdb_loci.str_[isj] == '-' {
            2
        } else if *map_gen.sjdb_motif.last().unwrap() == 0 {
            0
        } else {
            2 - (*map_gen.sjdb_motif.last().unwrap() % 2)
        };
        map_gen.sjdb_strand.push(strand);
    }

    map_gen.sjdb_n = map_gen.sjdb_start.len() as u32;
    map_gen.sj_dstart.clear();
    map_gen.sj_astart.clear();
    result
        .sjdb_info_txt
        .push_str(&format!("{}\t{}\n", map_gen.sjdb_n, map_gen.sjdb_overhang));

    let strand_char = ['.', '+', '-'];
    let mut sj_gstart = 0usize;
    result.gsj.resize(
        map_gen.sjdb_n as usize * map_gen.sjdb_length as usize,
        GENOME_SPACING_CHAR,
    );
    for ii in 0..map_gen.sjdb_n as usize {
        let mut d_start = map_gen.sjdb_start[ii] - map_gen.sjdb_overhang as u64;
        let mut a_start = map_gen.sjdb_end[ii] + 1;
        if map_gen.sjdb_motif[ii] == 0 {
            d_start += map_gen.sjdb_shift_left[ii] as u64;
            a_start += map_gen.sjdb_shift_left[ii] as u64;
        }
        map_gen.sj_dstart.push(d_start);
        map_gen.sj_astart.push(a_start);

        let over = map_gen.sjdb_overhang as usize;
        result.gsj[sj_gstart..sj_gstart + over]
            .copy_from_slice(&g[d_start as usize..d_start as usize + over]);
        result.gsj[sj_gstart + over..sj_gstart + 2 * over]
            .copy_from_slice(&g[a_start as usize..a_start as usize + over]);
        sj_gstart += map_gen.sjdb_length as usize;
        result.gsj[sj_gstart - 1] = GENOME_SPACING_CHAR;

        result.sjdb_info_txt.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            map_gen.sjdb_start[ii],
            map_gen.sjdb_end[ii],
            map_gen.sjdb_motif[ii],
            map_gen.sjdb_shift_left[ii],
            map_gen.sjdb_shift_right[ii],
            map_gen.sjdb_strand[ii]
        ));

        let chr1 =
            map_gen.chr_bin[(map_gen.sjdb_start[ii] >> p.p_ge.g_chr_bin_nbits) as usize] as usize;
        let restore = if map_gen.sjdb_motif[ii] > 0 {
            0
        } else {
            map_gen.sjdb_shift_left[ii] as u64
        };
        result.sjdb_list_out_tab.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            map_gen.chr_name[chr1],
            map_gen.sjdb_start[ii] - map_gen.chr_start[chr1] + 1 + restore,
            map_gen.sjdb_end[ii] - map_gen.chr_start[chr1] + 1 + restore,
            strand_char[map_gen.sjdb_strand[ii] as usize]
        ));
    }

    Ok(result)
}
