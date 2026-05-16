use clap::{Arg, Command};
use std::collections::BTreeSet;
use std::io::Write;
use std::path::Path;

use crate::{
    genome_genomeload_l18_genome_genomeload, genome_l15_genome_genome,
    parameters_l310_parameters_inputparameters, star_l58_main,
    transcriptome_l7_transcriptome_transcriptome, transcriptome_l156_transcriptome_quantsoutput,
};
use crate::{Parameters, StarMainResult, Transcriptome};

pub const PARAMETERS_DEFAULT: &str = include_str!("parametersDefault");

fn write_bgzf_block(
    out: &mut Vec<u8>,
    payload: &[u8],
    compression: flate2::Compression,
) -> Result<(), String> {
    let mut encoder = flate2::write::DeflateEncoder::new(Vec::new(), compression);
    encoder.write_all(payload).map_err(|e| e.to_string())?;
    let compressed = encoder.finish().map_err(|e| e.to_string())?;
    let block_size = 18 + compressed.len() + 8;
    if block_size > 65_536 {
        return Err("BGZF block exceeds 64KB".to_string());
    }
    out.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 255]);
    out.extend_from_slice(&6u16.to_le_bytes());
    out.extend_from_slice(b"BC");
    out.extend_from_slice(&2u16.to_le_bytes());
    out.extend_from_slice(&((block_size - 1) as u16).to_le_bytes());
    out.extend_from_slice(&compressed);
    out.extend_from_slice(&crc32fast::hash(payload).to_le_bytes());
    out.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    Ok(())
}

fn bgzf_compress_bam_bytes(bam: &[u8], compression_level: i32) -> Result<Vec<u8>, String> {
    let compression = if (0..=9).contains(&compression_level) {
        flate2::Compression::new(compression_level as u32)
    } else {
        flate2::Compression::default()
    };
    let mut out = Vec::new();
    let mut start = 0usize;
    while start < bam.len() {
        let mut end = (start + 60_000).min(bam.len());
        loop {
            let before = out.len();
            match write_bgzf_block(&mut out, &bam[start..end], compression) {
                Ok(()) => break,
                Err(_) if end > start + 1 => {
                    out.truncate(before);
                    end = start + (end - start) / 2;
                }
                Err(err) => return Err(err),
            }
        }
        start = end;
    }
    write_bgzf_block(&mut out, &[], compression)?;
    Ok(out)
}

pub fn cli_args() -> Vec<String> {
    let args: Vec<String> = std::env::args().collect();
    let _ = Command::new("STAR")
        .disable_help_flag(true)
        .disable_version_flag(true)
        .arg(
            Arg::new("star_args")
                .num_args(0..)
                .allow_hyphen_values(true)
                .trailing_var_arg(true),
        )
        .try_get_matches_from(args.clone());
    args
}

pub fn run_cli(args: &[String]) -> Result<StarMainResult, String> {
    if args.len() <= 1 || (args.len() == 2 && (args[1] == "-h" || args[1] == "--help")) {
        return star_l58_main(
            args,
            Parameters::default(),
            PARAMETERS_DEFAULT.as_bytes(),
            None,
            None,
            None,
            None,
            &BTreeSet::new(),
            &[],
            &[],
        );
    }

    let mut parameters = Parameters::default();
    let parameter_files = parameter_files_from_args(args)?;
    let mut manifest_name = Option::<String>::None;
    for (_, contents) in &parameter_files {
        for line in contents.lines() {
            let mut fields = line.split_whitespace();
            if fields.next() == Some("readFilesManifest") {
                manifest_name = fields.next().map(ToString::to_string);
            }
        }
    }
    let mut ii_manifest = 1;
    while ii_manifest < args.len() {
        let arg = &args[ii_manifest];
        if arg == "--readFilesManifest" {
            ii_manifest += 1;
            if ii_manifest < args.len() && !args[ii_manifest].starts_with("--") {
                manifest_name = Some(args[ii_manifest].clone());
            }
            while ii_manifest < args.len() && !args[ii_manifest].starts_with("--") {
                ii_manifest += 1;
            }
            continue;
        } else if let Some(value) = arg.strip_prefix("--readFilesManifest=") {
            manifest_name = Some(value.to_string());
        }
        ii_manifest += 1;
    }
    let manifest_contents = match manifest_name.as_deref() {
        Some("-") | None => None,
        Some(file_name) => Some(std::fs::read_to_string(file_name).map_err(|_| {
            format!(
                "EXITING because of FATAL INPUT FILE error: could not open readFileManifest file {}",
                file_name
            )
        })?),
    };
    let mut whitelist_names = Vec::<String>::new();
    let mut ii = 1;
    while ii < args.len() {
        let arg = &args[ii];
        if arg == "--soloCBwhitelist" {
            ii += 1;
            while ii < args.len() && !args[ii].starts_with("--") {
                whitelist_names.extend(args[ii].split_whitespace().map(ToString::to_string));
                ii += 1;
            }
            continue;
        } else if let Some(value) = arg.strip_prefix("--soloCBwhitelist=") {
            whitelist_names.extend(value.split_whitespace().map(ToString::to_string));
        }
        ii += 1;
    }
    for (_, contents) in &parameter_files {
        for line in contents.lines() {
            let mut fields = line.split_whitespace();
            if fields.next() == Some("soloCBwhitelist") {
                whitelist_names.extend(fields.map(ToString::to_string));
            }
        }
    }
    whitelist_names.sort();
    whitelist_names.dedup();
    let mut whitelist_contents = Vec::new();
    for file_name in whitelist_names {
        if file_name != "-" && file_name != "None" {
            let contents = std::fs::read_to_string(&file_name).map_err(|_| {
                format!(
                    "EXITING because of FATAL ERROR: could not open CB whitelist file {}",
                    file_name
                )
            })?;
            whitelist_contents.push((file_name, contents));
        }
    }
    let _scan_state = parameters_l310_parameters_inputparameters(
        &mut parameters,
        args,
        PARAMETERS_DEFAULT,
        &parameter_files,
        manifest_contents.as_deref(),
        &whitelist_contents,
    )?;

    let genome = if parameters
        .run_mode_in
        .first()
        .map(|s| s == "alignReads")
        .unwrap_or(true)
    {
        Some(load_genome_from_parameters(&mut parameters)?)
    } else {
        None
    };

    let existing_read_files = existing_read_files_from_args(args, Some(&parameters));
    if parameters.run_mode_in.first().map(String::as_str) == Some("alignReads")
        || parameters
            .run_mode_in
            .first()
            .map(|s| s.is_empty())
            .unwrap_or(true)
    {
        if !parameters.out_file_tmp.is_empty() {
            std::fs::create_dir_all(&parameters.out_file_tmp).map_err(|e| e.to_string())?;
        }
        if parameters.out_bam_coord && !parameters.out_bam_sort_tmp_dir.is_empty() {
            std::fs::create_dir_all(&parameters.out_bam_sort_tmp_dir).map_err(|e| e.to_string())?;
        }
    }

    let align_reads_mode = parameters.run_mode_in.first().map(String::as_str) == Some("alignReads")
        || parameters
            .run_mode_in
            .first()
            .map(|s| s.is_empty())
            .unwrap_or(true);
    let transcriptome = if parameters.quant_yes && align_reads_mode {
        let tr_info_dir = if parameters.p_ge.sjdb_gtf_file == "-" {
            parameters.p_ge.g_dir.clone()
        } else {
            parameters.sjdb_insert_out_dir.clone()
        };
        let gene_info = std::fs::read_to_string(Path::new(&tr_info_dir).join("geneInfo.tab"))
            .map_err(|e| e.to_string())?;
        let transcript_info =
            std::fs::read_to_string(Path::new(&tr_info_dir).join("transcriptInfo.tab")).ok();
        let exon_info = std::fs::read_to_string(Path::new(&tr_info_dir).join("exonInfo.tab")).ok();
        let exon_ge_tr_info =
            std::fs::read_to_string(Path::new(&tr_info_dir).join("exonGeTrInfo.tab")).ok();
        let (mut transcriptome, _, _log) = transcriptome_l7_transcriptome_transcriptome(
            parameters.quant_yes,
            parameters.p_ge.transform.out_quant,
            &parameters.p_ge.sjdb_gtf_file,
            &parameters.p_ge.g_dir,
            &parameters.sjdb_insert_out_dir,
            &parameters.p_ge.g_dir,
            parameters.quant_tr_sam_yes,
            parameters.quant_ge_count_yes || parameters.quant_gene_yes,
            parameters.quant_gene_full_yes
                || parameters.quant_gene_full_exon_over_intron_yes
                || parameters.quant_gene_full_ex50p_as_yes,
            &gene_info,
            transcript_info.as_deref(),
            exon_info.as_deref(),
            exon_ge_tr_info.as_deref(),
        )?;
        if parameters.quant_ge_count_yes {
            crate::transcriptome_l150_transcriptome_quantsallocate(
                &mut transcriptome,
                true,
            );
        }
        Some(transcriptome)
    } else {
        Some(Transcriptome::default())
    };

    let result = star_l58_main(
        args,
        parameters,
        PARAMETERS_DEFAULT.as_bytes(),
        genome,
        transcriptome,
        None,
        None,
        &existing_read_files,
        &[],
        &[],
    )?;

    if result.exit_code == 0
        && !result.parameters.out_file_name_prefix.is_empty()
        && (args.first().map(|arg0| arg0.as_str()) != Some("STAR")
            || result.parameters.out_file_name_prefix != "./")
    {
        let prefix = &result.parameters.out_file_name_prefix;
        if !result.log_main.is_empty() {
            let path = format!("{}Log.out", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            std::fs::write(path, &result.log_main).map_err(|e| e.to_string())?;
        }
        if !result.log_progress.is_empty() {
            let path = format!("{}Log.progress.out", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            std::fs::write(path, &result.log_progress).map_err(|e| e.to_string())?;
        }
        if result.parameters.out_std != "Log" && !result.log_stdout.is_empty() {
            let path = format!("{}Log.std.out", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            std::fs::write(path, &result.log_stdout).map_err(|e| e.to_string())?;
        }
        if !result.log_final_out.is_empty() {
            let path = format!("{}Log.final.out", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            std::fs::write(path, &result.log_final_out).map_err(|e| e.to_string())?;
        }
        if let Some(solo_result) = &result.solo_cell_filtering
            && let Some(cell_filtering) = &solo_result.cell_filtering
            && let Some(output_results) = &cell_filtering.output_results
        {
            for (path, contents) in &output_results.files {
                if let Some(parent) = Path::new(path).parent()
                    && !parent.as_os_str().is_empty()
                {
                    std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
                }
                std::fs::write(path, contents).map_err(|e| e.to_string())?;
            }
        }
        if let Some(solo_result) = &result.solo_process_and_output {
            for dir in &solo_result.created_directories {
                if !dir.is_empty() {
                    std::fs::create_dir_all(dir)
                        .map_err(|e| format!("could not create STARsolo directory {dir}: {e}"))?;
                }
            }
            for (path, contents) in &solo_result.files {
                if let Some(parent) = Path::new(path).parent()
                    && !parent.as_os_str().is_empty()
                {
                    std::fs::create_dir_all(parent).map_err(|e| {
                        format!("could not create parent directory for STARsolo output {path}: {e}")
                    })?;
                }
                std::fs::write(path, contents)
                    .map_err(|e| format!("could not write STARsolo output {path}: {e}"))?;
            }
            #[cfg(unix)]
            {
                for (target, link) in &solo_result.symlinks {
                    if target.is_empty() || link.is_empty() {
                        continue;
                    }
                    if let Some(parent) = Path::new(link).parent()
                        && !parent.as_os_str().is_empty()
                    {
                        std::fs::create_dir_all(parent).map_err(|e| {
                            format!(
                                "could not create parent directory for STARsolo symlink {link}: {e}"
                            )
                        })?;
                    }
                    match std::os::unix::fs::symlink(target, link) {
                        Ok(()) => {}
                        Err(err) if err.kind() == std::io::ErrorKind::AlreadyExists => {}
                        Err(err) => {
                            return Err(format!(
                                "could not create STARsolo symlink {link} -> {target}: {err}"
                            ));
                        }
                    }
                }
            }
        }
        if result.parameters.out_sam_bool
            && result.parameters.out_std != "SAM"
            && !result.parameters.out_sam_contents.is_empty()
        {
            let path = format!("{}Aligned.out.sam", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            std::fs::write(path, &result.parameters.out_sam_contents).map_err(|e| e.to_string())?;
        }
        if result.parameters.out_bam_unsorted && result.parameters.out_std != "BAM_Unsorted" {
            let path = format!("{}Aligned.out.bam", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            let mut bam = result.parameters.out_bam_unsorted_header.clone();
            for chunk in &result.read_chunks {
                if let Some(out_bam) = &chunk.chunk_out_bam_unsorted {
                    bam.extend_from_slice(&out_bam.bgzf_bam);
                }
            }
            std::fs::write(
                path,
                bgzf_compress_bam_bytes(&bam, result.parameters.out_bam_compression)?,
            )
            .map_err(|e| e.to_string())?;
        }
        if result.parameters.out_bam_coord
            && result.parameters.out_std != "BAM_SortedByCoordinate"
            && let Some(bam_sort) = &result.bam_sort
            && !bam_sort.output_bam.is_empty()
        {
            let path = format!("{}Aligned.sortedByCoord.out.bam", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            std::fs::write(
                path,
                bgzf_compress_bam_bytes(
                    &bam_sort.output_bam,
                    result.parameters.out_bam_compression,
                )?,
            )
            .map_err(|e| e.to_string())?;
        }
        if result.parameters.quant_tr_sam_bam_yes && result.parameters.out_std != "BAM_Quant" {
            let path = format!("{}Aligned.toTranscriptome.out.bam", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            let mut bam = result.parameters.out_quant_bam_header.clone();
            for process_chunk in &result.process_chunks {
                for map_chunk in &process_chunk.map_chunks {
                    bam.extend_from_slice(&map_chunk.quant_bam_output);
                }
            }
            std::fs::write(
                path,
                bgzf_compress_bam_bytes(&bam, result.parameters.out_bam_compression)?,
            )
            .map_err(|e| e.to_string())?;
        }
        if result.parameters.p_ch.out_chim_sam_opened {
            let path = format!("{}Chimeric.out.sam", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            std::fs::write(path, &result.parameters.p_ch.out_chim_sam_contents)
                .map_err(|e| e.to_string())?;
        }
        if result.parameters.p_ch.out_chim_junction_opened {
            let path = format!("{}Chimeric.out.junction", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            std::fs::write(path, &result.parameters.p_ch.out_chim_junction_contents)
                .map_err(|e| e.to_string())?;
        }
        if let Some(output_sj) = &result.output_sj {
            let path = format!("{}SJ.out.tab", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            std::fs::write(path, &output_sj.sj_out_tab).map_err(|e| e.to_string())?;
        }
        if result.parameters.out_reads_unmapped == "Fastx" {
            for imate in 0..result.parameters.read_nends as usize {
                let path = format!("{}Unmapped.out.mate{}", prefix, imate + 1);
                if let Some(parent) = Path::new(&path).parent()
                    && !parent.as_os_str().is_empty()
                {
                    std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
                }
                let mut contents = String::new();
                for process_chunk in &result.process_chunks {
                    for map_chunk in &process_chunk.map_chunks {
                        if let Some(unmapped) = map_chunk.unmapped_fastx_outputs.get(imate) {
                            contents.push_str(unmapped);
                        }
                    }
                }
                std::fs::write(path, contents).map_err(|e| e.to_string())?;
            }
        }
        if let Some(signal) = &result.signal {
            for (path, contents) in &signal.files {
                if let Some(parent) = Path::new(path).parent()
                    && !parent.as_os_str().is_empty()
                {
                    std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
                }
                std::fs::write(path, contents).map_err(|e| e.to_string())?;
            }
        }
        if !result.processed_bam_output.is_empty() {
            let path = format!("{}Processed.out.bam", prefix);
            if let Some(parent) = Path::new(&path).parent()
                && !parent.as_os_str().is_empty()
            {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            std::fs::write(
                path,
                bgzf_compress_bam_bytes(
                    &result.processed_bam_output,
                    result.parameters.out_bam_compression,
                )?,
            )
            .map_err(|e| e.to_string())?;
        }
        if result.parameters.quant_ge_count_yes
            && let Some(transcriptome) = &result.transcriptome
        {
            transcriptome_l156_transcriptome_quantsoutput(
                transcriptome,
                &format!("{}ReadsPerGene.out.tab", prefix),
                &result.stats_all,
            )?;
        }
    }

    Ok(result)
}

pub fn print_result(result: &StarMainResult) {
    if !result.usage.is_empty() {
        print!("{}", result.usage);
    }
    if result.parameters.out_std == "SAM" && !result.parameters.out_sam_contents.is_empty() {
        print!("{}", result.parameters.out_sam_contents);
    } else if result.parameters.out_std == "BAM_Unsorted" && result.parameters.out_bam_unsorted {
        let mut bam = result.parameters.out_bam_unsorted_header.clone();
        for chunk in &result.read_chunks {
            if let Some(out_bam) = &chunk.chunk_out_bam_unsorted {
                bam.extend_from_slice(&out_bam.bgzf_bam);
            }
        }
        if let Ok(bam) = bgzf_compress_bam_bytes(&bam, result.parameters.out_bam_compression) {
            let _ = std::io::stdout().write_all(&bam);
        }
    } else if result.parameters.out_std == "BAM_SortedByCoordinate"
        && let Some(bam_sort) = &result.bam_sort
    {
        if let Ok(bam) =
            bgzf_compress_bam_bytes(&bam_sort.output_bam, result.parameters.out_bam_compression)
        {
            let _ = std::io::stdout().write_all(&bam);
        }
    } else if result.parameters.out_std == "BAM_Quant" && result.parameters.quant_tr_sam_bam_yes {
        let mut bam = result.parameters.out_quant_bam_header.clone();
        for process_chunk in &result.process_chunks {
            for map_chunk in &process_chunk.map_chunks {
                bam.extend_from_slice(&map_chunk.quant_bam_output);
            }
        }
        if let Ok(bam) = bgzf_compress_bam_bytes(&bam, result.parameters.out_bam_compression) {
            let _ = std::io::stdout().write_all(&bam);
        }
    } else if !result.log_stdout.is_empty() {
        print!("{}", result.log_stdout);
    }
    if !result.log_main.is_empty() {
        eprint!("{}", result.log_main);
    }
    if !result.log_progress.is_empty() {
        eprint!("{}", result.log_progress);
    }
    if !result.log_final_out.is_empty() {
        eprint!("{}", result.log_final_out);
    }
}

pub fn parameter_files_from_args(args: &[String]) -> Result<Vec<(String, String)>, String> {
    let mut files = Vec::new();
    let mut ii = 1;
    while ii < args.len() {
        let arg = &args[ii];
        if arg == "--parametersFiles" {
            ii += 1;
            while ii < args.len() && !args[ii].starts_with("--") {
                let value = &args[ii];
                if value != "-" && !value.is_empty() {
                    for file_name in value.split_whitespace() {
                        let contents = std::fs::read_to_string(file_name).map_err(|_| {
                            format!(
                                "EXITING: FATAL INPUT ERROR: could not open parameters file {}\n",
                                file_name
                            )
                        })?;
                        files.push((file_name.to_string(), contents));
                    }
                }
                ii += 1;
            }
            continue;
        } else if let Some(value) = arg.strip_prefix("--parametersFiles=")
            && value != "-"
            && !value.is_empty()
        {
            for file_name in value.split_whitespace() {
                let contents = std::fs::read_to_string(file_name).map_err(|_| {
                    format!(
                        "EXITING: FATAL INPUT ERROR: could not open parameters file {}\n",
                        file_name
                    )
                })?;
                files.push((file_name.to_string(), contents));
            }
        }
        ii += 1;
    }
    Ok(files)
}

pub fn existing_read_files_from_args(
    args: &[String],
    parameters: Option<&Parameters>,
) -> BTreeSet<String> {
    let mut files = BTreeSet::new();
    if let Some(parameters) = parameters {
        for names in &parameters.read_files_names {
            for file_name in names {
                if !file_name.is_empty() && file_name != "-" {
                    files.insert(file_name.clone());
                }
            }
        }
        if !files.is_empty() {
            return files;
        }

        let prefix = &parameters.read_files_prefix_final;
        for file_name in &parameters.read_files_in {
            for split_name in file_name.split(',') {
                if !split_name.is_empty() && split_name != "-" {
                    files.insert(format!("{}{}", prefix, split_name));
                }
            }
        }
        if !files.is_empty() {
            return files;
        }
    }

    let mut ii = 1;
    while ii < args.len() {
        let arg = &args[ii];
        if arg == "--readFilesIn" {
            ii += 1;
            while ii < args.len() && !args[ii].starts_with("--") {
                for file_name in args[ii].split_whitespace() {
                    files.insert(file_name.to_string());
                }
                ii += 1;
            }
            continue;
        } else if let Some(value) = arg.strip_prefix("--readFilesIn=") {
            for file_name in value.split_whitespace() {
                files.insert(file_name.to_string());
            }
        }
        ii += 1;
    }
    files
}

pub fn load_genome_from_parameters(
    parameters: &mut Parameters,
) -> Result<crate::Genome, String> {
    let mut genome = genome_l15_genome_genome(parameters.p_ge.clone());
    let genome_dir = Path::new(&parameters.p_ge.g_dir);
    let genome_parameters =
        crate::io_utils::read_to_string_auto_gzip(genome_dir.join("genomeParameters.txt")).ok();
    let chr_name =
        crate::io_utils::read_to_string_auto_gzip(genome_dir.join("chrName.txt")).map_err(
            |_| {
                format!(
                    "EXITING because of FATAL ERROR: could not open genome file {}/chrName.txt\nSOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permsissions\n",
                    parameters.p_ge.g_dir
                )
            },
        )?;
    let chr_length = crate::io_utils::read_to_string_auto_gzip(genome_dir.join("chrLength.txt"))
        .map_err(|_| {
            format!(
                "EXITING because of FATAL error, could not open file {}/chrLength.txt\nSOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n",
                parameters.p_ge.g_dir
            )
        })?;
    let chr_start = crate::io_utils::read_to_string_auto_gzip(genome_dir.join("chrStart.txt"))
        .map_err(|_| {
            format!(
                "EXITING because of FATAL error, could not open file {}/chrStart.txt\nSOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n",
                parameters.p_ge.g_dir
            )
        })?;
    let genome_contents = crate::io_utils::read_bytes_auto_gzip(genome_dir.join("Genome"))
        .map_err(|_| {
            format!(
                "EXITING because of FATAL error, could not open file {}/Genome\nSOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n",
                parameters.p_ge.g_dir
            )
        })?;
    let sa = crate::io_utils::read_bytes_auto_gzip(genome_dir.join("SA")).map_err(|_| {
        format!(
            "EXITING because of FATAL error, could not open file {}/SA\nSOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n",
            parameters.p_ge.g_dir
        )
    })?;
    let sa_index = crate::io_utils::read_bytes_auto_gzip(genome_dir.join("SAindex")).map_err(|_| {
        format!(
            "EXITING because of FATAL error, could not open file {}/SAindex\nSOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n",
            parameters.p_ge.g_dir
        )
    })?;
    let sjdb_info =
        crate::io_utils::read_to_string_auto_gzip(genome_dir.join("sjdbInfo.txt")).ok();
    let raw_time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    let _ = genome_genomeload_l18_genome_genomeload(
        &mut genome,
        parameters,
        genome_parameters.as_deref(),
        &chr_name,
        &chr_length,
        &chr_start,
        genome_contents,
        sa,
        sa_index,
        sjdb_info.as_deref(),
        raw_time,
    )?;
    Ok(genome)
}
