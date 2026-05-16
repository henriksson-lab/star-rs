#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::storeAligns` at STAR/source/ReadAlign_storeAligns.cpp:10. Args: iDir: uint, Shift: uint, Nrep: uint, L: uint, indStartEnd: uint, iFrag: uint"]
pub fn readalign_storealigns_l10_readalign_storealigns(
    read_align: &mut crate::read_align::ReadAlign,
    p: &crate::parameters_chimeric::Parameters,
    i_dir: u64,
    shift: u64,
    nrep: u64,
    l: u64,
    ind_start_end: [u64; 2],
    i_frag: u64,
) -> Result<(), String> {
    if nrep > p.seed_multimap_nmax as u64 {
        if nrep < read_align.mult_nmin || read_align.mult_nmin == 0 {
            read_align.mult_nmin = nrep;
            read_align.mult_nmin_l = l;
        }
        return Ok(());
    }

    read_align.n_um[if nrep == 1 { 0 } else { 1 }] += nrep;
    read_align.n_a += nrep;

    let r_start = if i_dir == 0 {
        shift
    } else {
        shift.wrapping_add(1).wrapping_sub(l)
    };

    let mut ip = read_align.n_p as i32 - 1;
    while ip >= 0 {
        let piece = &read_align.pc[ip as usize];
        if piece[PC_R_START] <= r_start {
            if piece[PC_R_START] == r_start && piece[PC_LENGTH] < l {
                ip -= 1;
                continue;
            }
            if piece[PC_R_START] == r_start && piece[PC_LENGTH] == l {
                return Ok(());
            }
            break;
        }
        ip -= 1;
    }

    let ip = (ip + 1) as usize;
    let n_p = read_align.n_p as usize;
    if read_align.n_p >= p.seed_per_read_nmax as u64 {
        return Err(
            "EXITING because of FATAL error: too many pieces pere read\nSOLUTION: increase input parameter --seedPerReadNmax"
                .to_string(),
        );
    }
    if read_align.pc.len() <= n_p {
        read_align.pc.resize(n_p + 1, [0; PC_SIZE]);
    }
    if ip < n_p {
        read_align.pc.copy_within(ip..n_p, ip + 1);
    }
    read_align.n_p += 1;

    read_align.pc[ip][PC_R_START] = r_start;
    read_align.pc[ip][PC_LENGTH] = l;
    read_align.pc[ip][PC_DIR] = i_dir;
    read_align.pc[ip][PC_NREP] = nrep;
    read_align.pc[ip][PC_SASTART] = ind_start_end[0];
    read_align.pc[ip][PC_SAEND] = ind_start_end[1];
    read_align.pc[ip][PC_IFRAG] = i_frag;

    let best_l = if l < read_align.stored_lmin {
        read_align.stored_lmin
    } else {
        l
    };

    if nrep == 1 {
        if best_l > read_align.uniq_lmax {
            read_align.uniq_lmax = best_l;
            read_align.uniq_lmax_ind = read_align.n_p - 1;
        }
    } else {
        if nrep < read_align.mult_nmin || read_align.mult_nmin == 0 {
            read_align.mult_nmin = nrep;
            read_align.mult_nmin_l = best_l;
        }
        if best_l > read_align.mult_lmax {
            read_align.mult_lmax = best_l;
            read_align.mult_lmax_n = nrep;
        }
        if nrep > read_align.mult_nmax {
            read_align.mult_nmax = nrep;
            read_align.mult_nmax_l = best_l;
        }
    }

    Ok(())
}
