#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `BAMbinSortUnmapped` at STAR/source/BAMbinSortUnmapped.cpp:5. Args: iBin: uint32, nThreads: uint, dirBAMsort: string, P: Parameters, genome: Genome, solo: Solo"]
pub fn bambinsortunmapped_l5_bambinsortunmapped(
    i_bin: u32,
    n_threads: u64,
    dir_bamsort: &str,
    p: &crate::parameters_chimeric::Parameters,
    genome: &crate::genome::Genome,
    temp_files: &[Vec<u8>],
) -> Result<(Vec<u8>, Vec<String>), String> {
    let mut bam_in_file = Vec::new();
    for it in 0..n_threads {
        bam_in_file.push(format!("{}{}/{}", dir_bamsort, it, i_bin));
        bam_in_file.push(format!("{}{}/{}.BySJout", dir_bamsort, it, i_bin));
    }

    let mut cursors = vec![0usize; bam_in_file.len()];
    let mut current_records: Vec<Option<(Vec<u8>, u32)>> = vec![None; bam_in_file.len()];
    let mut start_pos = std::collections::BTreeMap::<u64, usize>::new();

    for it in 0..bam_in_file.len() {
        let file_bytes = temp_files.get(it).map(|v| v.as_slice()).unwrap_or(&[]);
        if cursors[it] + 4 <= file_bytes.len() {
            let block_len =
                u32::from_le_bytes(file_bytes[cursors[it]..cursors[it] + 4].try_into().unwrap());
            let bam_size = block_len.checked_add(4).ok_or_else(|| {
                format!(
                    "EXITING because of FATAL ERROR: malformed temporary bam file: {}",
                    bam_in_file[it]
                )
            })?;
            let record_end = cursors[it].checked_add(bam_size as usize).ok_or_else(|| {
                format!(
                    "EXITING because of FATAL ERROR: malformed temporary bam file: {}",
                    bam_in_file[it]
                )
            })?;
            let trailer_end = record_end.checked_add(8).ok_or_else(|| {
                format!(
                    "EXITING because of FATAL ERROR: malformed temporary bam file: {}",
                    bam_in_file[it]
                )
            })?;
            if trailer_end > file_bytes.len() {
                return Err(format!(
                    "EXITING because of FATAL ERROR: truncated temporary bam file: {}",
                    bam_in_file[it]
                ));
            }
            let record = file_bytes[cursors[it]..record_end].to_vec();
            let i_read =
                u64::from_le_bytes(file_bytes[record_end..trailer_end].try_into().unwrap());
            cursors[it] = trailer_end;
            current_records[it] = Some((record, bam_size));
            start_pos.insert(i_read >> 32, it);
        } else if cursors[it] != file_bytes.len() {
            return Err(format!(
                "EXITING because of FATAL ERROR: truncated temporary bam file: {}",
                bam_in_file[it]
            ));
        }
    }

    let mut bgzf_bin = bamfunctions_l77_outbamwriteheader(
        &p.sam_header_sorted_coord,
        &genome.chr_name_all,
        &genome.chr_length_all,
    );

    while !start_pos.is_empty() {
        let (start_key, it) = {
            let (key, value) = start_pos.iter().next().unwrap();
            (*key, *value)
        };
        let start_next = start_pos
            .range((
                std::ops::Bound::Excluded(start_key),
                std::ops::Bound::Unbounded,
            ))
            .next()
            .map(|(key, _)| *key)
            .unwrap_or(u64::MAX);

        loop {
            let (record, bam_size) = current_records[it]
                .take()
                .expect("BAMbinSortUnmapped current record missing");
            bgzf_bin.extend_from_slice(&record[..bam_size as usize]);

            let file_bytes = temp_files.get(it).map(|v| v.as_slice()).unwrap_or(&[]);
            if cursors[it] + 4 <= file_bytes.len() {
                let block_len = u32::from_le_bytes(
                    file_bytes[cursors[it]..cursors[it] + 4].try_into().unwrap(),
                );
                let bam_size = block_len.checked_add(4).ok_or_else(|| {
                    format!(
                        "EXITING because of FATAL ERROR: malformed temporary bam file: {}",
                        bam_in_file[it]
                    )
                })?;
                let record_end = cursors[it].checked_add(bam_size as usize).ok_or_else(|| {
                    format!(
                        "EXITING because of FATAL ERROR: malformed temporary bam file: {}",
                        bam_in_file[it]
                    )
                })?;
                let trailer_end = record_end.checked_add(8).ok_or_else(|| {
                    format!(
                        "EXITING because of FATAL ERROR: malformed temporary bam file: {}",
                        bam_in_file[it]
                    )
                })?;
                if trailer_end > file_bytes.len() {
                    return Err(format!(
                        "EXITING because of FATAL ERROR: truncated temporary bam file: {}",
                        bam_in_file[it]
                    ));
                }
                let record = file_bytes[cursors[it]..record_end].to_vec();
                let i_read =
                    u64::from_le_bytes(file_bytes[record_end..trailer_end].try_into().unwrap());
                cursors[it] = trailer_end;
                current_records[it] = Some((record, bam_size));
                let i_read_top = i_read >> 32;
                if i_read_top > start_next {
                    start_pos.insert(i_read_top, it);
                    break;
                }
            } else {
                if cursors[it] != file_bytes.len() {
                    return Err(format!(
                        "EXITING because of FATAL ERROR: truncated temporary bam file: {}",
                        bam_in_file[it]
                    ));
                }
                break;
            }
        }

        start_pos.remove(&start_key);
    }

    Ok((bgzf_bin, bam_in_file))
}
