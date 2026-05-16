#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SoloReadBarcode` at STAR/source/SoloReadBarcode.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloReadBarcode {
    pub solo_type: i32,
    pub cb_wl_yes: bool,
    pub cb_wl_size: u32,
    pub umi_l: u32,
    pub homo_polymer: [u64; 4],
    pub umi_seq: String,
    pub cb_seq: String,
    pub cb_qual: String,
    pub umi_qual: String,
    pub b_seq: String,
    pub b_qual: String,
    pub b_strings: Vec<String>,
    pub cb_seq_corrected: String,
    pub umi_b: u64,
    pub umi_check: i32,
    pub cb_match: i32,
    pub cb_match_string: String,
    pub cb_match_ind: Vec<u64>,
    pub cb_read_count_exact: Vec<u32>,
    pub qual_hist: Vec<u64>,
    pub stats: SoloReadBarcodeStats,
}

#[doc = "Original `SoloReadBarcode::SoloReadBarcode` at STAR/source/SoloReadBarcode.cpp:4. Args: P: Parameters"]
pub fn soloreadbarcode_l4_soloreadbarcode_soloreadbarcode(
    solo_type: i32,
    cb_wl_yes: bool,
    cb_wl_size: u32,
    umi_l: u32,
) -> crate::solo_read_barcode::SoloReadBarcode {
    let mut solo_read_barcode = crate::solo_read_barcode::SoloReadBarcode {
        solo_type,
        cb_wl_yes,
        cb_wl_size,
        umi_l,
        stats: crate::solo_read_barcode_stats::SoloReadBarcodeStats {
            names: vec![
                "noNoAdapter".to_string(),
                "noNoUMI".to_string(),
                "noNoCB".to_string(),
                "noNinCB".to_string(),
                "noNinUMI".to_string(),
                "noUMIhomopolymer".to_string(),
                "noNoWLmatch".to_string(),
                "noTooManyMM".to_string(),
                "noTooManyWLmatches".to_string(),
                "yesWLmatchExact".to_string(),
                "yesOneWLmatchWithMM".to_string(),
                "yesMultWLmatchWithMM".to_string(),
            ],
            v: vec![0; SOLO_READ_BARCODE_N_STATS],
        },
        ..Default::default()
    };

    if solo_read_barcode.solo_type == 0 {
        return solo_read_barcode;
    }

    if solo_read_barcode.cb_wl_yes {
        solo_read_barcode.cb_read_count_exact = vec![0; solo_read_barcode.cb_wl_size as usize];
    }

    for jj in 0..4u64 {
        solo_read_barcode.homo_polymer[jj as usize] = 0;
        for _ in 0..solo_read_barcode.umi_l {
            solo_read_barcode.homo_polymer[jj as usize] =
                (solo_read_barcode.homo_polymer[jj as usize] << 2) + jj;
        }
    }

    solo_read_barcode.qual_hist = vec![0; 256];
    solo_read_barcode
}

#[doc = "Original `SoloReadBarcode::addCounts` at STAR/source/SoloReadBarcode.cpp:26. Args: rfIn: SoloReadBarcode"]
pub fn soloreadbarcode_l26_soloreadbarcode_addcounts(
    solo_read_barcode: &mut crate::solo_read_barcode::SoloReadBarcode,
    rf_in: &crate::solo_read_barcode::SoloReadBarcode,
) {
    if solo_read_barcode.cb_wl_yes {
        for ii in 0..solo_read_barcode.cb_wl_size as usize {
            solo_read_barcode.cb_read_count_exact[ii] += rf_in.cb_read_count_exact[ii];
        }
    }

    for ii in 0..solo_read_barcode.qual_hist.len() {
        solo_read_barcode.qual_hist[ii] += rf_in.qual_hist[ii];
    }
}

#[doc = "Original `SoloReadBarcode::addStats` at STAR/source/SoloReadBarcode.cpp:38. Args: rfIn: SoloReadBarcode"]
pub fn soloreadbarcode_l38_soloreadbarcode_addstats(
    solo_read_barcode: &mut crate::solo_read_barcode::SoloReadBarcode,
    rf_in: &crate::solo_read_barcode::SoloReadBarcode,
) {
    for ii in 0..solo_read_barcode.stats.v.len() {
        solo_read_barcode.stats.v[ii] += rf_in.stats.v[ii];
    }
}

#[doc = "Original `SoloReadBarcode::statsOut` at STAR/source/SoloReadBarcode.cpp:44. Args: streamOut: ofstream"]
pub fn soloreadbarcode_l44_soloreadbarcode_statsout(
    solo_read_barcode: &crate::solo_read_barcode::SoloReadBarcode,
) -> String {
    use std::fmt::Write;

    let mut out = String::new();
    for ii in 0..solo_read_barcode.stats.v.len() {
        writeln!(
            out,
            "{:>50}{:>15}",
            solo_read_barcode.stats.names[ii], solo_read_barcode.stats.v[ii]
        )
        .unwrap();
    }
    out
}
