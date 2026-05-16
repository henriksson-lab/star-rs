#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SoloBarcode` at STAR/source/SoloBarcode.h:9."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloBarcode {
    pub anchor_type: [i32; 2],
    pub anchor_dist: [i32; 2],
    pub adapter_length: i32,
    pub wl: Vec<Vec<u64>>,
    pub wl_ed: Vec<Vec<u64>>,
    pub wl_ed_ind: Vec<Vec<u32>>,
    pub wl_factor: u64,
    pub wl_add: Vec<u32>,
    pub min_len: u32,
    pub total_size: u32,
    pub i_cb: u32,
    pub i_len: u32,
}

#[doc = "Original struct `type_cbMMind` at STAR/source/SoloBarcode.cpp:49."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct type_cbMMind {
    pub cb: i64,
    pub ind: u32,
    pub mm: u32,
}

#[doc = "Original `SoloBarcode::sortWhiteList` at STAR/source/SoloBarcode.cpp:9. Args: pSolo: ParametersSolo"]
pub fn solobarcode_l9_solobarcode_sortwhitelist(
    barcode: &mut crate::solo_barcode::SoloBarcode,
    edit_dist_2: bool,
) {
    barcode.total_size = 0;
    barcode.min_len = u32::MAX;
    barcode.wl_add.resize(barcode.wl.len(), 0);
    if edit_dist_2 {
        barcode.wl_ed.resize(barcode.wl.len(), Vec::new());
        barcode.wl_ed_ind.resize(barcode.wl.len(), Vec::new());
    }

    for ilen1 in 1..barcode.wl.len() {
        barcode.wl_add[ilen1] = barcode.total_size;
        if !barcode.wl[ilen1].is_empty() {
            if (ilen1 as u32) < barcode.min_len {
                barcode.min_len = ilen1 as u32;
            }
            barcode.wl[ilen1].sort_unstable();
            barcode.wl[ilen1].dedup();
            barcode.total_size += barcode.wl[ilen1].len() as u32;

            if edit_dist_2 {
                solobarcode_l47_wladdmismatches(
                    2,
                    ilen1 as u32,
                    &barcode.wl[ilen1],
                    &mut barcode.wl_ed[ilen1],
                    &mut barcode.wl_ed_ind[ilen1],
                );
            }
        }
    }
}

#[doc = "Original `SoloBarcode::extractPositionsFromString` at STAR/source/SoloBarcode.cpp:37. Args: strIn: string"]
pub fn solobarcode_l37_solobarcode_extractpositionsfromstring(
    barcode: &mut crate::solo_barcode::SoloBarcode,
    str_in: &str,
) -> Result<(), String> {
    let mut parts = Vec::new();
    servicefuns_l167_splitstring(str_in, '_', &mut parts);
    if parts.len() != 4 {
        return Err(format!(
            "malformed solo barcode position '{}': expected 4 underscore-separated fields",
            str_in
        ));
    }

    barcode.anchor_type[0] = parts[0].parse::<i32>().map_err(|_| {
        format!(
            "malformed solo barcode position '{}': invalid anchor type '{}'",
            str_in, parts[0]
        )
    })?;
    barcode.anchor_type[1] = parts[2].parse::<i32>().map_err(|_| {
        format!(
            "malformed solo barcode position '{}': invalid anchor type '{}'",
            str_in, parts[2]
        )
    })?;
    barcode.anchor_dist[0] = parts[1].parse::<i32>().map_err(|_| {
        format!(
            "malformed solo barcode position '{}': invalid anchor distance '{}'",
            str_in, parts[1]
        )
    })?;
    barcode.anchor_dist[1] = parts[3].parse::<i32>().map_err(|_| {
        format!(
            "malformed solo barcode position '{}': invalid anchor distance '{}'",
            str_in, parts[3]
        )
    })?;
    Ok(())
}

#[doc = "Original `wlAddMismatches` at STAR/source/SoloBarcode.cpp:47. Args: nMM: uint32, cbLen: uint32, wl: vector<uintCB>, wlEd1: vector<uintCB>, wlEdInd1: vector<uint32>"]
pub fn solobarcode_l47_wladdmismatches(
    n_mm: u32,
    cb_len: u32,
    wl: &[u64],
    wl_ed1: &mut Vec<u64>,
    wl_ed_ind1: &mut Vec<u32>,
) {
    #[derive(Clone)]
    struct TypeCbMmInd {
        cb: u64,
        ind: u32,
        mm: u32,
    }

    let ntot = (wl.len() as f64 * (((cb_len * 3) as f64).powi((n_mm + 1) as i32) - 1.0)
        / ((cb_len * 3 - 1) as f64)) as usize;
    let mut cb_mm_ind = Vec::with_capacity(ntot);

    for (icb, cb) in wl.iter().enumerate() {
        cb_mm_ind.push(TypeCbMmInd {
            cb: *cb,
            ind: icb as u32,
            mm: 0,
        });
    }

    let mut ind1 = 0usize;
    let mut ind2 = wl.len();
    for mm in 1..=n_mm {
        let mut ind3 = ind2;
        for ii in ind1..ind2 {
            for ll in (0..cb_len * 2).step_by(2) {
                for jj in 1..4u64 {
                    let cbmm = cb_mm_ind[ii].cb ^ (jj << ll);
                    cb_mm_ind.push(TypeCbMmInd {
                        cb: cbmm,
                        ind: cb_mm_ind[ii].ind,
                        mm,
                    });
                    ind3 += 1;
                }
            }
        }

        if mm == 2 {
            for ii in 0..wl.len() {
                let cbmm = cb_mm_ind[ii].cb;
                for ld in (0..cb_len * 2).step_by(2) {
                    let maskd1 = u64::MAX << (ld + 2);
                    let maskd = (!maskd1) >> 2;
                    let cbmmd = (cbmm & maskd) | ((cbmm & maskd1) >> 2);
                    for ll in (0..cb_len * 2).step_by(2) {
                        let cbmm1 = cbmmd << 2;
                        let mask1 = u64::MAX << (ll + 2);
                        let mask = (!mask1) >> 2;
                        let cbmm2 = (cbmmd & mask) | (cbmm1 & mask1);
                        for jj in 0..4u64 {
                            let cbmm3 = cbmm2 | (jj << ll);
                            cb_mm_ind.push(TypeCbMmInd {
                                cb: cbmm3,
                                ind: cb_mm_ind[ii].ind,
                                mm,
                            });
                            ind3 += 1;
                        }
                    }
                }
            }
        }
        ind1 = ind2;
        ind2 = ind3;
    }

    cb_mm_ind.sort_by(|c1, c2| {
        c1.cb
            .cmp(&c2.cb)
            .then(c1.mm.cmp(&c2.mm))
            .then(c1.ind.cmp(&c2.ind))
    });

    let mut n_cb_out = 0usize;
    let mut prev_cb = u64::MAX;
    for ii in 0..cb_mm_ind.len() {
        if ii < cb_mm_ind.len() - 1
            && cb_mm_ind[ii].cb == cb_mm_ind[ii + 1].cb
            && cb_mm_ind[ii].mm == cb_mm_ind[ii + 1].mm
            && cb_mm_ind[ii].ind == cb_mm_ind[ii + 1].ind
        {
            cb_mm_ind[ii].mm = n_mm + 1;
            continue;
        }

        if (ii > 0 && cb_mm_ind[ii].cb == prev_cb)
            || (ii < cb_mm_ind.len() - 1
                && cb_mm_ind[ii].cb == cb_mm_ind[ii + 1].cb
                && cb_mm_ind[ii].mm == cb_mm_ind[ii + 1].mm)
        {
            cb_mm_ind[ii].mm = n_mm + 1;
        } else {
            n_cb_out += 1;
        }
        prev_cb = cb_mm_ind[ii].cb;
    }

    wl_ed1.resize(n_cb_out, 0);
    wl_ed_ind1.resize(n_cb_out, 0);
    let mut icb = 0usize;
    for cb1 in &cb_mm_ind {
        if cb1.mm <= n_mm {
            wl_ed1[icb] = cb1.cb;
            wl_ed_ind1[icb] = cb1.ind;
            icb += 1;
        }
    }
}
