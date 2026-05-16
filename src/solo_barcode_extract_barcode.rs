#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloBarcode::extractBarcode` at STAR/source/SoloBarcode_extractBarcode.cpp:4. Args: seqIn: string, qualIn: string, adapterStart: uint32, bSeq: string, bQual: string"]
pub fn solobarcode_extractbarcode_l4_solobarcode_extractbarcode(
    barcode: &crate::solo_barcode::SoloBarcode,
    seq_in: &str,
    qual_in: &str,
    adapter_start: u32,
) -> Option<(String, String)> {
    let mut pos = [0i32; 2];
    for ii in 0..2 {
        pos[ii] = match barcode.anchor_type[ii] {
            0 => 0,
            1 => seq_in.len() as i32 - 1,
            2 => adapter_start as i32,
            3 => adapter_start as i32 + barcode.adapter_length - 1,
            _ => 0,
        };
        pos[ii] += barcode.anchor_dist[ii];
    }

    if pos[0] < 0 || pos[1] > seq_in.len() as i32 || pos[0] > pos[1] {
        return None;
    }

    let start = pos[0] as usize;
    let end_exclusive = pos[1] as usize + 1;
    Some((
        seq_in[start..end_exclusive].to_string(),
        qual_in[start..end_exclusive].to_string(),
    ))
}
