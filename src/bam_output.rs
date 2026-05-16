#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `BAMoutput` at STAR/source/BAMoutput.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BAMoutput {
    pub n_bins: u32,
    pub bin_total_n: Vec<u64>,
    pub bin_total_bytes: Vec<u64>,
    pub bam_array_size: u64,
    pub bam_array: Vec<u8>,
    pub bin_size: u64,
    pub bin_size1: u64,
    pub bin_bytes: Vec<u64>,
    pub bin_bytes1: u64,
    pub bin_buffers: Vec<Vec<u8>>,
    pub bin_streams: Vec<Vec<u8>>,
    pub bin_stream_by_sjout: Vec<bool>,
    pub bgzf_bam: Vec<u8>,
    pub bam_dir: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BamCore {
    pub tid: i32,
    pub pos: i32,
    pub bin: u32,
    pub qual: u32,
    pub l_qname: u32,
    pub flag: u32,
    pub n_cigar: u32,
    pub l_qseq: i32,
    pub mtid: i32,
    pub mpos: i32,
    pub isize: i32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Bam1 {
    pub core: BamCore,
    pub l_data: i32,
    pub m_data: i32,
    pub data_offset: usize,
}

#[doc = "Original `BAMoutput::BAMoutput` at STAR/source/BAMoutput.cpp:9. Args: iChunk: int, tmpDir: string, Pin: Parameters"]
pub fn bamoutput_l9_bamoutput_bamoutput(
    i_chunk: i32,
    tmp_dir: &str,
    p: &crate::parameters_chimeric::Parameters,
) -> crate::bam_output::BAMoutput {
    let mut n_bins = p.out_bam_coord_nbins;
    let bin_size = p.chunk_out_bam_size_bytes / n_bins as u64;
    let bam_array_size = bin_size * n_bins as u64;
    let bam_array = vec![0u8; bam_array_size as usize];
    let bam_dir = format!("{}{}", tmp_dir, i_chunk as u32);
    let bin_buffers = vec![Vec::new(); n_bins as usize];
    let bin_streams = vec![Vec::new(); n_bins as usize];
    let bin_stream_by_sjout = vec![false; n_bins as usize];
    let bin_bytes = vec![0u64; n_bins as usize];
    let bin_total_n = vec![0u64; n_bins as usize];
    let bin_total_bytes = vec![0u64; n_bins as usize];
    let bin_size1 = bin_size * (n_bins as u64 - 1);
    n_bins = 1;

    crate::bam_output::BAMoutput {
        n_bins,
        bin_total_n,
        bin_total_bytes,
        bam_array_size,
        bam_array,
        bin_size,
        bin_size1,
        bin_bytes,
        bin_bytes1: 0,
        bin_buffers,
        bin_streams,
        bin_stream_by_sjout,
        bgzf_bam: Vec::new(),
        bam_dir,
    }
}

#[doc = "Original `BAMoutput::BAMoutput` at STAR/source/BAMoutput.cpp:36. Args: bgzfBAMin: BGZF, Pin: Parameters"]
pub fn bamoutput_l36_bamoutput_bamoutput(
    bgzf_bam_in: Vec<u8>,
    p: &crate::parameters_chimeric::Parameters,
) -> crate::bam_output::BAMoutput {
    crate::bam_output::BAMoutput {
        bam_array_size: p.chunk_out_bam_size_bytes,
        bam_array: vec![0u8; p.chunk_out_bam_size_bytes as usize],
        bin_bytes1: 0,
        bgzf_bam: bgzf_bam_in,
        bin_size: 0,
        bin_streams: Vec::new(),
        bin_buffers: Vec::new(),
        bin_bytes: Vec::new(),
        bin_total_bytes: Vec::new(),
        bin_total_n: Vec::new(),
        n_bins: 0,
        ..Default::default()
    }
}

#[doc = "Original `BAMoutput::unsortedOneAlign` at STAR/source/BAMoutput.cpp:52. Args: bamIn: char, bamSize: uint, bamSize2: uint"]
pub fn bamoutput_l52_bamoutput_unsortedonealign(
    bam_output: &mut crate::bam_output::BAMoutput,
    bam_in: &[u8],
    bam_size: u64,
    bam_size2: u64,
) -> Result<(), String> {
    if bam_size == 0 {
        return Ok(());
    }

    if bam_size as usize > bam_in.len() {
        return Err(format!(
            "malformed BAM record for unsorted output: record size {} exceeds input buffer {}",
            bam_size,
            bam_in.len()
        ));
    }
    if bam_size > bam_output.bam_array_size {
        return Err(format!(
            "malformed BAM record for unsorted output: record size {} exceeds output buffer {}",
            bam_size, bam_output.bam_array_size
        ));
    }

    let needed = bam_output
        .bin_bytes1
        .checked_add(bam_size2)
        .ok_or_else(|| {
            "malformed BAM unsorted output state: bin byte count overflow".to_string()
        })?;
    if needed > bam_output.bam_array_size {
        bam_output
            .bgzf_bam
            .extend_from_slice(&bam_output.bam_array[..bam_output.bin_bytes1 as usize]);
        bam_output.bin_bytes1 = 0;
    }

    let start = bam_output.bin_bytes1 as usize;
    let end = start.checked_add(bam_size as usize).ok_or_else(|| {
        "malformed BAM unsorted output state: bin byte count overflow".to_string()
    })?;
    if end > bam_output.bam_array.len() {
        return Err(format!(
            "malformed BAM record for unsorted output: record size {} exceeds output buffer {}",
            bam_size,
            bam_output.bam_array.len()
        ));
    }
    bam_output.bam_array[start..end].copy_from_slice(&bam_in[..bam_size as usize]);
    bam_output.bin_bytes1 = bam_output.bin_bytes1.checked_add(bam_size).ok_or_else(|| {
        "malformed BAM unsorted output state: bin byte count overflow".to_string()
    })?;
    Ok(())
}

#[doc = "Original `BAMoutput::unsortedFlush` at STAR/source/BAMoutput.cpp:70. Args: "]
pub fn bamoutput_l70_bamoutput_unsortedflush(
    bam_output: &mut crate::bam_output::BAMoutput,
) {
    bam_output
        .bgzf_bam
        .extend_from_slice(&bam_output.bam_array[..bam_output.bin_bytes1 as usize]);
    bam_output.bin_bytes1 = 0;
}

#[doc = "Original `BAMoutput::coordOneAlign` at STAR/source/BAMoutput.cpp:77. Args: bamIn: char, bamSize: uint, iRead: uint"]
pub fn bamoutput_l77_bamoutput_coordonealign(
    bam_output: &mut crate::bam_output::BAMoutput,
    p: &mut crate::parameters_chimeric::Parameters,
    bam_in: &[u8],
    bam_size: u64,
    i_read: u64,
) -> Result<(), String> {
    let mut i_bin = 0usize;
    if bam_size == 0 {
        return Ok(());
    } else {
        if bam_size as usize > bam_in.len() {
            return Err(format!(
                "malformed BAM record for coordinate output: record size {} exceeds input buffer {}",
                bam_size,
                bam_in.len()
            ));
        }
        if bam_size < 12 {
            return Err(format!(
                "malformed BAM record for coordinate output: record size {} is too small",
                bam_size
            ));
        }
        let bam_in32_1 = u32::from_le_bytes(bam_in[4..8].try_into().unwrap());
        let bam_in32_2 = u32::from_le_bytes(bam_in[8..12].try_into().unwrap());
        let align_g = ((bam_in32_1 as u64) << 32) | bam_in32_2 as u64;
        if bam_in32_1 == u32::MAX {
            i_bin = (p.out_bam_coord_nbins - 1) as usize;
        } else if bam_output.n_bins > 1 {
            let i_bin_i = servicefuns_l239_binarysearch1a(
                align_g,
                &p.out_bam_sorting_bin_start,
                bam_output.n_bins as i32 - 1,
            );
            if i_bin_i < 0 {
                return Err(format!(
                    "malformed BAM record for coordinate output: genomic coordinate {} is outside sorting bins",
                    align_g
                ));
            }
            i_bin = i_bin_i as usize;
        }
    }
    if i_bin >= bam_output.bin_buffers.len() || i_bin >= bam_output.bin_bytes.len() {
        return Err(format!(
            "malformed BAM coordinate output state: bin {} is outside {} bins",
            i_bin,
            bam_output.bin_buffers.len()
        ));
    }

    let limit = if i_bin > 0 || bam_output.n_bins > 1 {
        bam_output.bin_size
    } else {
        bam_output.bin_size1
    };
    let needed = bam_output.bin_bytes[i_bin]
        .checked_add(bam_size)
        .and_then(|v| v.checked_add(8))
        .ok_or_else(|| {
            "malformed BAM coordinate output state: bin byte count overflow".to_string()
        })?;
    if needed > limit {
        if bam_output.n_bins > 1 || i_bin == (p.out_bam_coord_nbins - 1) as usize {
            bam_output.bin_streams[i_bin].extend_from_slice(&bam_output.bin_buffers[i_bin]);
            bam_output.bin_buffers[i_bin].clear();
            bam_output.bin_bytes[i_bin] = 0;
        } else {
            bamoutput_l118_bamoutput_coordbins(bam_output, p)?;
            bamoutput_l77_bamoutput_coordonealign(bam_output, p, bam_in, bam_size, i_read)?;
            return Ok(());
        }
    }

    bam_output.bin_buffers[i_bin].extend_from_slice(&bam_in[..bam_size as usize]);
    bam_output.bin_bytes[i_bin] = bam_output.bin_bytes[i_bin]
        .checked_add(bam_size)
        .ok_or_else(|| {
            "malformed BAM coordinate output state: bin byte count overflow".to_string()
        })?;
    bam_output.bin_buffers[i_bin].extend_from_slice(&i_read.to_le_bytes());
    bam_output.bin_bytes[i_bin] = bam_output.bin_bytes[i_bin].checked_add(8).ok_or_else(|| {
        "malformed BAM coordinate output state: bin byte count overflow".to_string()
    })?;
    let total_add = bam_size.checked_add(8).ok_or_else(|| {
        "malformed BAM coordinate output state: total byte count overflow".to_string()
    })?;
    bam_output.bin_total_bytes[i_bin] = bam_output.bin_total_bytes[i_bin]
        .checked_add(total_add)
        .ok_or_else(|| {
            "malformed BAM coordinate output state: total byte count overflow".to_string()
        })?;
    bam_output.bin_total_n[i_bin] =
        bam_output.bin_total_n[i_bin]
            .checked_add(1)
            .ok_or_else(|| {
                "malformed BAM coordinate output state: record count overflow".to_string()
            })?;
    Ok(())
}

#[doc = "Original `BAMoutput::coordBins` at STAR/source/BAMoutput.cpp:118. Args: "]
pub fn bamoutput_l118_bamoutput_coordbins(
    bam_output: &mut crate::bam_output::BAMoutput,
    p: &mut crate::parameters_chimeric::Parameters,
) -> Result<(), String> {
    bam_output.n_bins = p.out_bam_coord_nbins;

    if p.out_bam_sorting_bin_start[0] != 0 {
        let mut start_pos = Vec::with_capacity(bam_output.bin_total_n[0] as usize + 1);
        let mut ib = 0usize;
        for _ in 0..bam_output.bin_total_n[0] {
            if ib + 12 > bam_output.bin_buffers[0].len() {
                return Err(
                    "malformed BAM coordinate bin buffer: truncated record header".to_string(),
                );
            }
            let bam_in32_0 =
                u32::from_le_bytes(bam_output.bin_buffers[0][ib..ib + 4].try_into().unwrap());
            let bam_in32_1 = u32::from_le_bytes(
                bam_output.bin_buffers[0][ib + 4..ib + 8]
                    .try_into()
                    .unwrap(),
            );
            let bam_in32_2 = u32::from_le_bytes(
                bam_output.bin_buffers[0][ib + 8..ib + 12]
                    .try_into()
                    .unwrap(),
            );
            start_pos.push(((bam_in32_1 as u64) << 32) | bam_in32_2 as u64);
            let next = ib
                .checked_add(bam_in32_0 as usize)
                .and_then(|v| v.checked_add(4 + 8))
                .ok_or_else(|| {
                    "malformed BAM coordinate bin buffer: record size overflow".to_string()
                })?;
            if next > bam_output.bin_buffers[0].len() {
                return Err(
                    "malformed BAM coordinate bin buffer: truncated record body".to_string()
                );
            }
            ib = next;
        }
        start_pos.sort_unstable();

        p.bam_sorting_log.push_str(&format!(
            "BAM sorting: {} mapped reads\n",
            bam_output.bin_total_n[0]
        ));
        p.bam_sorting_log
            .push_str("BAM sorting bins genomic start loci:\n");

        p.out_bam_sorting_bin_start[0] = 0;
        for ib in 1..(bam_output.n_bins - 1) as usize {
            p.out_bam_sorting_bin_start[ib] = start_pos
                [bam_output.bin_total_n[0] as usize / (bam_output.n_bins as usize - 1) * ib];
            p.bam_sorting_log.push_str(&format!(
                "{}\t{}\t{}\n",
                ib,
                p.out_bam_sorting_bin_start[ib] >> 32,
                (p.out_bam_sorting_bin_start[ib] << 32) >> 32
            ));
        }
    }

    let bin_total_n_old = bam_output.bin_total_n[0];
    let bin_start_old = bam_output.bin_buffers[0].clone();
    bam_output.bin_buffers[0].clear();
    bam_output.bin_bytes[0] = 0;
    bam_output.bin_total_n[0] = 0;
    bam_output.bin_total_bytes[0] = 0;

    let mut ib = 0usize;
    for _ in 0..bin_total_n_old {
        if ib + 4 > bin_start_old.len() {
            return Err(
                "malformed BAM coordinate bin buffer: truncated stored record size".to_string(),
            );
        }
        let bam_in32_0 = u32::from_le_bytes(bin_start_old[ib..ib + 4].try_into().unwrap());
        let ib1 = ib
            .checked_add(bam_in32_0 as usize)
            .and_then(|v| v.checked_add(4))
            .ok_or_else(|| {
                "malformed BAM coordinate bin buffer: stored record size overflow".to_string()
            })?;
        if ib1
            .checked_add(8)
            .is_none_or(|end| end > bin_start_old.len())
        {
            return Err(
                "malformed BAM coordinate bin buffer: truncated stored read order".to_string(),
            );
        }
        let i_read = u64::from_le_bytes(bin_start_old[ib1..ib1 + 8].try_into().unwrap());
        bamoutput_l77_bamoutput_coordonealign(
            bam_output,
            p,
            &bin_start_old[ib..ib1],
            bam_in32_0 as u64 + 4,
            i_read,
        )?;
        ib = ib1 + 8;
    }
    Ok(())
}

#[doc = "Original `BAMoutput::coordFlush` at STAR/source/BAMoutput.cpp:168. Args: "]
pub fn bamoutput_l168_bamoutput_coordflush(
    bam_output: &mut crate::bam_output::BAMoutput,
    p: &mut crate::parameters_chimeric::Parameters,
) -> Result<(), String> {
    if bam_output.n_bins == 1 {
        bamoutput_l118_bamoutput_coordbins(bam_output, p)?;
    }
    for i_bin in 0..bam_output.n_bins as usize {
        bam_output.bin_streams[i_bin].extend_from_slice(&bam_output.bin_buffers[i_bin]);
        bam_output.bin_buffers[i_bin].clear();
        bam_output.bin_bytes[i_bin] = 0;
    }
    Ok(())
}

#[doc = "Original `BAMoutput::coordUnmappedPrepareBySJout` at STAR/source/BAMoutput.cpp:179. Args: "]
pub fn bamoutput_l179_bamoutput_coordunmappedpreparebysjout(
    bam_output: &mut crate::bam_output::BAMoutput,
    p: &crate::parameters_chimeric::Parameters,
) {
    let i_bin = (p.out_bam_coord_nbins - 1) as usize;
    bam_output.bin_streams[i_bin].extend_from_slice(&bam_output.bin_buffers[i_bin]);
    bam_output.bin_buffers[i_bin].clear();
    bam_output.bin_bytes[i_bin] = 0;
    bam_output.bin_stream_by_sjout[i_bin] = true;
}
