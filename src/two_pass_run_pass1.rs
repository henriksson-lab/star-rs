#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `twoPassRunPass1` at STAR/source/twoPassRunPass1.cpp:9. Args: P: Parameters, genomeMain: Genome, transcriptomeMain: Transcriptome, sjdbLoci: SjdbClass"]
pub fn twopassrunpass1_l9_twopassrunpass1(
    p: &mut crate::parameters_chimeric::Parameters,
    genome_main: &mut crate::genome::Genome,
    transcriptome_main: Option<&crate::transcriptome::Transcriptome>,
    sjdb_loci: &mut crate::sjdb_class::SjdbClass,
    pass1_chunks_in: Option<Vec<crate::read_align_chunk::ReadAlignChunk>>,
    existing_read_files: &std::collections::BTreeSet<String>,
) -> Result<crate::parameters_chimeric::TwoPassRunPass1Result, String> {
    let mut result = crate::parameters_chimeric::TwoPassRunPass1Result::default();
    if !p.two_pass_yes {
        return Ok(result);
    }

    let genome_main1 = genome_sjdb_insert_snapshot(genome_main);
    let mut p1 = p.clone();
    if p1.out_sam_type.is_empty() {
        p1.out_sam_type.push("None".to_string());
    } else {
        p1.out_sam_type[0] = "None".to_string();
    }
    p1.out_sam_bool = false;
    p1.out_bam_unsorted = false;
    p1.out_bam_coord = false;
    p1.p_ch.segment_min = 0;
    p1.quant_yes = false;
    p1.quant_tr_sam_yes = false;
    p1.quant_tr_sam_bam_yes = false;
    p1.quant_gene_full_yes = false;
    p1.quant_ge_count_yes = false;
    p1.quant_gene_yes = false;
    p1.out_sam_unmapped_within = false;
    p1.out_filter_by_sjout_stage = 0;
    p1.out_reads_unmapped = "None".to_string();
    p1.out_file_name_prefix = p.two_pass_dir.clone();
    if p.two_pass_pass1reads_n > 0 {
        p1.read_map_number = std::cmp::min(p.two_pass_pass1reads_n, p.read_map_number);
    }
    p1.wasp_output_mode = "None".to_string();
    p1.wasp_yes = false;
    p1.wasp_sam_tag = false;
    p1.p_solo.type_str = "None".to_string();
    p1.p_solo.solo_type = 0;
    p1.p_ge.transform.out_yes = false;
    p1.p_ge.transform.out_quant = false;
    p1.p_ge.transform.out_sam = false;
    p1.p_ge.transform.out_sj = false;
    result.pass1_parameters = p1.clone();

    let mut stats_all = crate::stats::Stats::default();
    stats_l4_stats_resetn(&mut stats_all);
    let raw_time_start = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    stats_all.time_start = raw_time_start;
    stats_all.time_start_map = raw_time_start;
    result.log_progress.push_str(&format!(
        "{}\tStarted 1st pass mapping\n",
        timefunctions_l14_timemonthdaytime(raw_time_start)
    ));
    result.log_stdout.push_str(&format!(
        "{} ..... started 1st pass mapping\n",
        timefunctions_l14_timemonthdaytime(raw_time_start)
    ));

    let conv_yes = genome_main.genome_out.conv_yes;
    genome_main.genome_out.conv_yes = false;

    let mut pass1_chunks = if let Some(chunks) = pass1_chunks_in {
        chunks
    } else {
        let mut chunks = Vec::with_capacity(p1.run_thread_n.max(0) as usize);
        for ii in 0..p1.run_thread_n {
            chunks.push(readalignchunk_l5_readalignchunk_readalignchunk(
                &p1,
                genome_main,
                transcriptome_main,
                ii,
            )?);
        }
        chunks
    };

    result
        .log_main
        .push_str(&mapthreadsspawn_l6_mapthreadsspawn(
            p1.run_thread_n,
            &vec![0; p1.run_thread_n.max(0) as usize],
            &vec![0; p1.run_thread_n.max(0) as usize],
            || Ok(String::new()),
        )?);
    let output_sj = outputsj_l20_outputsj(&pass1_chunks, &mut p1, genome_main)?;
    result.log_main.push_str(&output_sj.log_main);

    genome_main.genome_out.conv_yes = conv_yes;

    let raw_time_finish = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    result.log_progress.push_str(&format!(
        "{}\tFinished 1st pass mapping\n",
        timefunctions_l14_timemonthdaytime(raw_time_finish)
    ));
    result.log_stdout.push_str(&format!(
        "{} ..... finished 1st pass mapping\n",
        timefunctions_l14_timemonthdaytime(raw_time_finish)
    ));
    result.log_final_out = stats_l99_stats_reportfinal(&mut stats_all, raw_time_finish);

    std::fs::create_dir_all(&p.two_pass_dir).map_err(|e| e.to_string())?;
    p.two_pass_pass2 = true;
    p.two_pass_pass1sj_file = format!("{}/SJ.out.tab", p.two_pass_dir);
    std::fs::write(&p.two_pass_pass1sj_file, &output_sj.sj_out_tab).map_err(|e| e.to_string())?;
    result.pass1_sj_file = p.two_pass_pass1sj_file.clone();

    let sjdb_insert =
        sjdbinsertjunctions_l11_sjdbinsertjunctions(p, genome_main, &genome_main1, sjdb_loci)?;
    result.log_main.push_str(&sjdb_insert.log_main);

    result.killed_read_command_pids = parameters_closereadsfiles_l5_parameters_closereadsfiles(p);
    result.reopened_reads =
        parameters_openreadsfiles_l5_parameters_openreadsfiles(p, existing_read_files)?;
    result.log_main.push_str(&result.reopened_reads.log_main);

    result.pass1_chunks = std::mem::take(&mut pass1_chunks);
    result.output_sj = output_sj;
    result.sjdb_insert = Some(sjdb_insert);
    Ok(result)
}
