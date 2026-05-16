#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `Junction` at STAR/source/OutSJ.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Junction {
    pub gen_out: Genome,
    pub record: Option<JunctionRecord>,
}

#[doc = "Original class `OutSJ` at STAR/source/OutSJ.h:36."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OutSJ {
    pub n: u64,
    pub n_store: u64,
    pub junctions: Vec<JunctionRecord>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct JunctionRecord {
    pub start: u64,
    pub gap: u32,
    pub strand: i8,
    pub motif: i32,
    pub annot: u8,
    pub count_unique: u32,
    pub count_multiple: u32,
    pub overhang_left: u16,
    pub overhang_right: u16,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ExitWithErrorResult {
    pub stream_out1: String,
    pub stream_out2: String,
    pub error_int: i32,
    pub thread_mutex_locked: bool,
    pub in_out_deleted: bool,
}

#[doc = "Original `OutSJ::OutSJ` at STAR/source/OutSJ.cpp:4. Args: nSJmax: uint, Pin: Parameters, genomeIn: Genome"]
pub fn outsj_l4_outsj_outsj(n_sjmax: u64) -> crate::out_sj::OutSJ {
    crate::out_sj::OutSJ {
        n: 0,
        n_store: n_sjmax,
        junctions: Vec::new(),
    }
}

#[doc = "Original `compareSJ` at STAR/source/OutSJ.cpp:15. Args: i1: void, i2: void"]
pub fn outsj_l15_comparesj(i1: &[u8], i2: &[u8]) -> i32 {
    let mut s1_bytes = [0u8; std::mem::size_of::<u64>()];
    let mut s2_bytes = [0u8; std::mem::size_of::<u64>()];
    if let Some(bytes) = i1.get(0..std::mem::size_of::<u64>()) {
        s1_bytes.copy_from_slice(bytes);
    }
    if let Some(bytes) = i2.get(0..std::mem::size_of::<u64>()) {
        s2_bytes.copy_from_slice(bytes);
    }
    let s1 = u64::from_ne_bytes(s1_bytes);
    let s2 = u64::from_ne_bytes(s2_bytes);
    if s1 > s2 {
        1
    } else if s1 < s2 {
        -1
    } else {
        let mut g1_bytes = [0u8; std::mem::size_of::<u32>()];
        let mut g2_bytes = [0u8; std::mem::size_of::<u32>()];
        if let Some(bytes) = i1.get(
            std::mem::size_of::<u64>()..std::mem::size_of::<u64>() + std::mem::size_of::<u32>(),
        ) {
            g1_bytes.copy_from_slice(bytes);
        }
        if let Some(bytes) = i2.get(
            std::mem::size_of::<u64>()..std::mem::size_of::<u64>() + std::mem::size_of::<u32>(),
        ) {
            g2_bytes.copy_from_slice(bytes);
        }
        let g1 = u32::from_ne_bytes(g1_bytes);
        let g2 = u32::from_ne_bytes(g2_bytes);
        if g1 > g2 {
            1
        } else if g1 < g2 {
            -1
        } else {
            0
        }
    }
}

#[doc = "Original `OutSJ::collapseSJ` at STAR/source/OutSJ.cpp:36. Args: "]
pub fn outsj_l36_outsj_collapsesj(out_sj: &mut crate::out_sj::OutSJ) -> Result<(), String> {
    if out_sj.n == 0 {
        return Ok(());
    }

    out_sj
        .junctions
        .sort_by(|a, b| a.start.cmp(&b.start).then_with(|| a.gap.cmp(&b.gap)));

    let mut isj1 = 0usize;
    for isj in 1..out_sj.n as usize {
        let sj = out_sj.junctions[isj].clone();
        if out_sj.junctions[isj1].start == sj.start && out_sj.junctions[isj1].gap == sj.gap {
            outsj_l92_junction_collapseonesj(&mut out_sj.junctions[isj1], &sj)?;
        } else {
            isj1 += 1;
            if isj != isj1 {
                out_sj.junctions[isj1] = sj;
            }
        }
    }
    out_sj.n = isj1 as u64 + 1;
    out_sj.junctions.truncate(out_sj.n as usize);
    Ok(())
}

#[doc = "Original `OutSJ::dataSizeIncrease` at STAR/source/OutSJ.cpp:62. Args: "]
pub fn outsj_l62_outsj_datasizeincrease(out_sj: &mut crate::out_sj::OutSJ) {
    out_sj.n_store *= 2;
    let target = out_sj.n_store as usize;
    if out_sj.junctions.capacity() < target {
        out_sj
            .junctions
            .reserve_exact(target - out_sj.junctions.len());
    }
}

#[doc = "Original `Junction::Junction` at STAR/source/OutSJ.cpp:68. Args: genOut: Genome"]
pub fn outsj_l68_junction_junction() -> crate::out_sj::Junction {
    crate::out_sj::Junction::default()
}

#[doc = "Original `Junction::junctionPointer` at STAR/source/OutSJ.cpp:72. Args: sjPoint: char, isj: uint"]
pub fn outsj_l72_junction_junctionpointer(
    junction: &mut crate::out_sj::Junction,
    sj_point: &[crate::out_sj::JunctionRecord],
    isj: u32,
) -> Result<(), String> {
    let record = sj_point
        .get(isj as usize)
        .ok_or_else(|| format!("junctionPointer index {} out of range", isj))?;
    junction.record = Some(record.clone());
    Ok(())
}

#[doc = "Original `Junction::outputStream` at STAR/source/OutSJ.cpp:85. Args: outStream: ostream"]
pub fn outsj_l85_junction_outputstream(
    junction: &crate::out_sj::Junction,
) -> Result<String, String> {
    let record = junction
        .record
        .as_ref()
        .ok_or_else(|| "junction outputStream called before junctionPointer".to_string())?;
    let bin_index = (record.start >> junction.gen_out.p_ge.g_chr_bin_nbits as u64) as usize;
    let sj_chr = *junction
        .gen_out
        .chr_bin
        .get(bin_index)
        .ok_or_else(|| format!("chrBin index {} out of range", bin_index))?
        as usize;
    let chr_name = junction
        .gen_out
        .chr_name
        .get(sj_chr)
        .ok_or_else(|| format!("chrName index {} out of range", sj_chr))?;
    let chr_start = *junction
        .gen_out
        .chr_start
        .get(sj_chr)
        .ok_or_else(|| format!("chrStart index {} out of range", sj_chr))?;

    Ok(format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        chr_name,
        record.start + 1 - chr_start,
        record.start + record.gap as u64 - chr_start,
        record.strand as i32,
        record.motif,
        record.annot,
        record.count_unique,
        record.count_multiple,
        record.overhang_left
    ))
}

#[doc = "Original `Junction::collapseOneSJ` at STAR/source/OutSJ.cpp:92. Args: isj1P: char, isjP: char, P: Parameters"]
pub fn outsj_l92_junction_collapseonesj(
    isj1: &mut crate::out_sj::JunctionRecord,
    isj: &crate::out_sj::JunctionRecord,
) -> Result<(), String> {
    isj1.count_unique += isj.count_unique;
    isj1.count_multiple += isj.count_multiple;

    if isj1.overhang_left < isj.overhang_left {
        isj1.overhang_left = isj.overhang_left;
    }
    if isj1.overhang_right < isj.overhang_right {
        isj1.overhang_right = isj.overhang_right;
    }

    if isj1.motif != isj.motif {
        return Err(format!(
            "EXITING because of BUG: different motifs for the same junction while collapsing junctions\n{} {} {} {} {} {}\n",
            isj1.start, isj1.gap, isj1.motif, isj.motif, isj1.annot, isj.annot
        ));
    }
    if isj1.annot < isj.annot {
        return Err(format!(
            "EXITING because  of BUG: different annotation status for the same junction while collapsing junctions:{} {} {} {}\n",
            isj1.start, isj1.gap, isj1.annot, isj.annot
        ));
    }
    Ok(())
}
