#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `bamSortByCoordinate` at STAR/source/bamSortByCoordinate.cpp:8. Args: P: Parameters, RAchunk: ReadAlignChunk, genome: Genome, solo: Solo"]
pub fn bamsortbycoordinate_l8_bamsortbycoordinate(
    p: &mut crate::parameters_chimeric::Parameters,
    ra_chunk: &[crate::read_align_chunk::ReadAlignChunk],
    genome: &crate::genome::Genome,
    temp_files_by_bin: &[Vec<Vec<u8>>],
) -> Result<crate::read_align_chunk::BamSortByCoordinateResult, String> {
    let mut result = crate::read_align_chunk::BamSortByCoordinateResult::default();
    if !p.out_bam_coord {
        return Ok(result);
    }

    p.bam_sorting_log.push_str(&format!(
        "{} ..... started sorting BAM\n",
        timefunctions_l4_timemonthdaytime()
    ));
    let n_bins = p.out_bam_coord_nbins as usize;

    let mut max_mem = 0u64;
    for ibin in 0..n_bins - 1 {
        let mut bin_s = 0u64;
        for it in 0..p.run_thread_n as usize {
            let bam = &ra_chunk[it].chunk_out_bam_coord;
            bin_s += bam.bin_total_bytes[ibin] + 24 * bam.bin_total_n[ibin];
        }
        if bin_s > max_mem {
            max_mem = bin_s;
        }
    }

    let mut unmapped_reads_n = 0u64;
    for it in 0..p.run_thread_n as usize {
        unmapped_reads_n += ra_chunk[it].chunk_out_bam_coord.bin_total_n[n_bins - 1];
    }
    result.max_mem = max_mem;
    result.unmapped_reads_n = unmapped_reads_n;
    p.bam_sorting_log
        .push_str(&format!("Max memory needed for sorting = {}\n", max_mem));

    if max_mem > p.limit_bam_sort_ram {
        return Err(format!(
            "EXITING because of fatal ERROR: not enough memory for BAM sorting: \nSOLUTION: re-run STAR with at least --limitBAMsortRAM {}",
            max_mem + 1_000_000_000
        ));
    } else if max_mem == 0 && unmapped_reads_n == 0 {
        p.bam_sorting_log
            .push_str("WARNING: nothing to sort - no output alignments\n");
        result.output_bam = bamfunctions_l77_outbamwriteheader(
            &p.sam_header_sorted_coord,
            &genome.chr_name_all,
            &genome.chr_length_all,
        );
        return Ok(result);
    }

    result.bin_outputs = vec![Vec::new(); n_bins];
    let mut total_mem = 0u64;
    for ibin1 in 0..n_bins {
        let ibin = n_bins - 1 - ibin1;
        let mut bin_n = 0u64;
        let mut bin_s = 0u64;
        for it in 0..p.run_thread_n as usize {
            let bam = &ra_chunk[it].chunk_out_bam_coord;
            bin_n += bam.bin_total_n[ibin];
            bin_s += bam.bin_total_bytes[ibin];
        }

        if bin_s == 0 {
            continue;
        }

        let files = temp_files_by_bin
            .get(ibin)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        if ibin == n_bins - 1 {
            let (bin_out, removed) = bambinsortunmapped_l5_bambinsortunmapped(
                ibin as u32,
                p.run_thread_n as u64,
                &p.out_bam_sort_tmp_dir,
                p,
                genome,
                files,
            )?;
            result.bin_outputs[ibin] = bin_out;
            result.removed_files.extend(removed);
        } else {
            let new_mem = bin_s + bin_n * 24;
            if total_mem + new_mem >= p.limit_bam_sort_ram {
                return Err(format!(
                    "EXITING because of fatal ERROR: not enough memory for BAM sorting: \nSOLUTION: re-run STAR with at least --limitBAMsortRAM {}",
                    total_mem + new_mem + 1_000_000_000
                ));
            }
            total_mem += new_mem;
            let (bin_out, removed) = bambinsortbycoordinate_l7_bambinsortbycoordinate(
                ibin as u32,
                bin_n,
                bin_s,
                p.run_thread_n as u64,
                &p.out_bam_sort_tmp_dir,
                p,
                genome,
                files,
            )?;
            result.bin_outputs[ibin] = bin_out;
            result.removed_files.extend(removed);
            total_mem -= new_mem;
        }
    }

    for ibin in 0..n_bins {
        let bam_bin_name = format!("{}/b{}", p.out_bam_sort_tmp_dir, ibin);
        if !result.bin_outputs[ibin].is_empty() {
            result.bin_names.push(bam_bin_name);
            result
                .output_bam
                .extend_from_slice(&result.bin_outputs[ibin]);
        }
    }

    Ok(result)
}
