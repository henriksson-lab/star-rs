#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `BAMbinSortByCoordinate` at STAR/source/BAMbinSortByCoordinate.cpp:7. Args: iBin: uint32, binN: uint, binS: uint, nThreads: uint, dirBAMsort: string, P: Parameters, genome: Genome, solo: Solo"]
pub fn bambinsortbycoordinate_l7_bambinsortbycoordinate(
    i_bin: u32,
    bin_n: u64,
    bin_s: u64,
    n_threads: u64,
    dir_bamsort: &str,
    p: &crate::parameters_chimeric::Parameters,
    genome: &crate::genome::Genome,
    temp_files: &[Vec<u8>],
) -> Result<(Vec<u8>, Vec<String>), String> {
    if bin_s == 0 {
        return Ok((Vec::new(), Vec::new()));
    }

    let mut bam_in = vec![0u8; bin_s as usize + 1];
    let mut bam_in_bytes = 0usize;
    let mut removed_files = Vec::new();
    for it in 0..n_threads as usize {
        let bam_in_file = format!("{}{}/{}", dir_bamsort, it, i_bin);
        let file_bytes = temp_files.get(it).map(|v| v.as_slice()).unwrap_or(&[]);
        if !file_bytes.is_empty() {
            let end = bam_in_bytes + file_bytes.len();
            bam_in[bam_in_bytes..end].copy_from_slice(file_bytes);
            bam_in_bytes = end;
        }
        removed_files.push(bam_in_file);
    }

    if bam_in_bytes as u64 != bin_s {
        return Err(format!(
            "EXITING because of FATAL ERROR: number of bytes expected from the BAM bin does not agree with the actual size on disk: Expected bin size={} ; size on disk={} ; bin number={}\n",
            bin_s, bam_in_bytes, i_bin
        ));
    }

    let mut start_pos = vec![[0u64; 3]; bin_n as usize];
    let mut ib = 0usize;
    for ia in 0..bin_n as usize {
        if ib.checked_add(12).is_none_or(|end| end > bam_in_bytes) {
            return Err(format!(
                "EXITING because of FATAL ERROR: truncated temporary BAM coordinate bin: {}",
                i_bin
            ));
        }
        let bam_in32_0 = u32::from_le_bytes(bam_in[ib..ib + 4].try_into().unwrap());
        let bam_in32_1 = u32::from_le_bytes(bam_in[ib + 4..ib + 8].try_into().unwrap());
        let bam_in32_2 = u32::from_le_bytes(bam_in[ib + 8..ib + 12].try_into().unwrap());
        start_pos[ia][0] = ((bam_in32_1 as u64) << 32) | bam_in32_2 as u64;
        start_pos[ia][2] = ib as u64;
        ib = ib
            .checked_add(bam_in32_0 as usize)
            .and_then(|v| v.checked_add(4))
            .ok_or_else(|| {
                format!(
                    "EXITING because of FATAL ERROR: malformed temporary BAM coordinate bin: {}",
                    i_bin
                )
            })?;
        if ib.checked_add(8).is_none_or(|end| end > bam_in_bytes) {
            return Err(format!(
                "EXITING because of FATAL ERROR: truncated temporary BAM coordinate bin: {}",
                i_bin
            ));
        }
        start_pos[ia][1] = u64::from_le_bytes(bam_in[ib..ib + 8].try_into().unwrap());
        ib += 8;
    }
    if ib != bam_in_bytes {
        return Err(format!(
            "EXITING because of FATAL ERROR: temporary BAM coordinate bin has trailing bytes: {}",
            i_bin
        ));
    }

    start_pos.sort_unstable();

    let mut bgzf_bin = bamfunctions_l77_outbamwriteheader(
        &p.sam_header_sorted_coord,
        &genome.chr_name_all,
        &genome.chr_length_all,
    );

    for ia in 0..bin_n as usize {
        let offset = start_pos[ia][2] as usize;
        if offset.checked_add(4).is_none_or(|end| end > bam_in_bytes) {
            return Err(format!(
                "EXITING because of FATAL ERROR: truncated temporary BAM coordinate bin: {}",
                i_bin
            ));
        }
        let size0 = u32::from_le_bytes(bam_in[offset..offset + 4].try_into().unwrap()) as usize + 4;
        if offset
            .checked_add(size0)
            .is_none_or(|end| end > bam_in_bytes)
        {
            return Err(format!(
                "EXITING because of FATAL ERROR: truncated temporary BAM coordinate bin: {}",
                i_bin
            ));
        }
        bgzf_bin.extend_from_slice(&bam_in[offset..offset + size0]);
    }

    Ok((bgzf_bin, removed_files))
}
