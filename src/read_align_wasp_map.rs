#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ReadAlign::waspMap` at STAR/source/ReadAlign_waspMap.cpp:3. Args: "]
pub fn readalign_waspmap_l3_readalign_waspmap(
    read_align: &mut crate::read_align::ReadAlign,
    wasp_ra: &mut crate::read_align::ReadAlign,
    p: &crate::parameters_chimeric::Parameters,
    map_gen: &crate::genome::Genome,
    gen_out: Option<&crate::genome::Genome>,
    aligns_gen_out: Option<&crate::quantifications::ReadAlignGenomeTransformResult>,
    remap_outcomes: &[crate::quantifications::WaspMapOutcome],
) -> crate::quantifications::WaspMapResult {
    let mut result = crate::quantifications::WaspMapResult::default();

    if !p.wasp_yes {
        read_align.wasp_type = -1;
        result.wasp_type = read_align.wasp_type;
        return result;
    }

    let mut n_tr1 = read_align.n_tr;
    let mut align1 = read_align.tr_best.clone();
    let mut var1 = map_gen.var.clone();
    if p.p_ge.transform.out_yes {
        if let Some(aligns_gen_out) = aligns_gen_out {
            n_tr1 = aligns_gen_out.al_n as u64;
            if let Some(al) = aligns_gen_out.al_mult.first() {
                align1 = al.clone();
            }
        }
        if let Some(gen_out) = gen_out {
            var1 = gen_out.var.clone();
            let read_index = if align1.str_ == 0 { 0 } else { 2 };
            transcript_variationadjust_l4_transcript_variationadjust(
                &mut align1,
                gen_out,
                &read_align.read1[read_index],
            );
        }
    }

    if align1.var_allele.is_empty() {
        read_align.wasp_type = -1;
        result.wasp_type = read_align.wasp_type;
        return result;
    } else if n_tr1 > 1 {
        read_align.wasp_type = 2;
        result.wasp_type = read_align.wasp_type;
        return result;
    } else if align1.var_allele.len() > 10 {
        read_align.wasp_type = 7;
        result.wasp_type = read_align.wasp_type;
        return result;
    }

    readalign_waspmap_l115_readalign_copyread(wasp_ra, read_align);
    let v_a = align1.var_allele.clone();
    for &a in &v_a {
        if a > 3 {
            read_align.wasp_type = 3;
            result.wasp_type = read_align.wasp_type;
            return result;
        }
    }

    let mut vv_a: Vec<Vec<u8>> = vec![Vec::new()];
    for _ in &v_a {
        let mut r = Vec::with_capacity(vv_a.len() * 2);
        for x in &vv_a {
            for y in [1_u8, 2_u8] {
                let mut x1 = x.clone();
                x1.push(y);
                r.push(x1);
            }
        }
        vv_a = r;
    }

    for v_a1 in vv_a {
        if v_a1 == v_a {
            continue;
        }

        for iv in 0..v_a1.len() {
            let isnp = align1.var_ind[iv] as usize;
            let allele = v_a1[iv] as usize;
            let mut nt2 = var1.snp.nt[isnp][allele];
            let mut vr = align1.var_read_coord[iv];
            if align1.str_ == 1 {
                nt2 = 3 - nt2;
                vr = read_align.l_read - 1 - vr;
            }
            let vr_usize = vr as usize;
            wasp_ra.read1[0][vr_usize] = nt2;
            wasp_ra.read1[1][vr_usize] = 3 - nt2;
            let rev_index = (read_align.l_read - 1 - vr) as usize;
            wasp_ra.read1[2][rev_index] = 3 - nt2;
        }

        result
            .requests
            .push(crate::quantifications::WaspMapRequest {
                alleles: v_a1.clone(),
                read1: wasp_ra.read1.clone(),
            });

        let Some(outcome) = remap_outcomes.get(result.requests.len() - 1) else {
            result.wasp_type = read_align.wasp_type;
            return result;
        };

        let (n_tr2, align2) = if p.p_ge.transform.out_yes {
            if let Some(transformed) = &outcome.transformed {
                (
                    transformed.al_n as u64,
                    transformed
                        .al_mult
                        .first()
                        .cloned()
                        .unwrap_or_else(|| outcome.align.clone()),
                )
            } else {
                (outcome.n_tr, outcome.align.clone())
            }
        } else {
            (outcome.n_tr, outcome.align.clone())
        };

        if outcome.unmap_type != -1 {
            read_align.wasp_type = 4;
            result.wasp_type = read_align.wasp_type;
            return result;
        } else if n_tr2 > 1 {
            read_align.wasp_type = 5;
            result.wasp_type = read_align.wasp_type;
            return result;
        } else if align2.n_exons != align1.n_exons {
            read_align.wasp_type = 6;
            result.wasp_type = read_align.wasp_type;
            return result;
        } else {
            for ii in 0..align1.n_exons as usize {
                for jj in 0..=2 {
                    if align1.exons[ii][jj] != align2.exons[ii][jj] {
                        read_align.wasp_type = 6;
                        result.wasp_type = read_align.wasp_type;
                        return result;
                    }
                }
            }
        }
    }

    read_align.wasp_type = 1;
    result.wasp_type = read_align.wasp_type;
    result
}

#[doc = "Original `ReadAlign::copyRead` at STAR/source/ReadAlign_waspMap.cpp:115. Args: r: ReadAlign"]
pub fn readalign_waspmap_l115_readalign_copyread(
    read_align: &mut crate::read_align::ReadAlign,
    r: &crate::read_align::ReadAlign,
) {
    read_align.l_read = r.l_read;
    if read_align.read_length.len() < 2 {
        read_align.read_length.resize(2, 0);
    }
    if read_align.read_length_original.len() < 2 {
        read_align.read_length_original.resize(2, 0);
    }
    read_align.read_length[0] = *r.read_length.first().unwrap_or(&0);
    read_align.read_length[1] = *r.read_length.get(1).unwrap_or(&0);
    read_align.read_length_original[0] = *r.read_length_original.first().unwrap_or(&0);
    read_align.read_length_original[1] = *r.read_length_original.get(1).unwrap_or(&0);
    read_align.read_length_pair_original = r.read_length_pair_original;
    read_align.out_filter_mismatch_nmax_total = r.out_filter_mismatch_nmax_total;
    read_align.read_name = r.read_name.clone();
    read_align.i_read_all = r.i_read_all;
    read_align.read_filter = r.read_filter;
    read_align.read_files_index = r.read_files_index;

    let l_read = r.l_read as usize;
    for ii in 0..=2 {
        read_align.read1[ii].clear();
        read_align.read1[ii].extend_from_slice(&r.read1[ii][..l_read]);
    }
}

pub fn streamfuns_create_dir(path: &str, mode: u32) -> std::io::Result<()> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::DirBuilderExt;
        std::fs::DirBuilder::new().mode(mode).create(path)
    }
    #[cfg(not(unix))]
    {
        let _ = mode;
        std::fs::create_dir(path)
    }
}
