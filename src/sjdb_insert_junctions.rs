#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `sjdbInsertJunctions` at STAR/source/sjdbInsertJunctions.cpp:11. Args: P: Parameters, mapGen: Genome, mapGen1: Genome, sjdbLoci: SjdbClass"]
pub fn sjdbinsertjunctions_l11_sjdbinsertjunctions(
    p: &mut crate::parameters_chimeric::Parameters,
    map_gen: &mut crate::genome::Genome,
    map_gen1: &crate::genome::Genome,
    sjdb_loci: &mut crate::sjdb_class::SjdbClass,
) -> Result<crate::parameters_chimeric::SjdbInsertJunctionsResult, String> {
    use std::io::Write;

    let mut result = crate::parameters_chimeric::SjdbInsertJunctionsResult::default();

    if map_gen.sjdb_n > 0 && sjdb_loci.chr.is_empty() {
        let file_name = format!("{}/sjdbList.out.tab", p.p_ge.g_dir);
        let contents = std::fs::read_to_string(&file_name).map_err(|_| {
            format!(
                "ERROR_OUT: exiting because of *INPUT FILE* error: could not open input file {}\nSolution: check that the file exists and you have read permission for this file\n          SOLUTION: re-generate the genome in pGe.gDir={}\n",
                file_name, p.p_ge.g_dir
            )
        })?;
        sjdbloadfromstream_l2_sjdbloadfromstream(&contents, sjdb_loci);
        sjdb_loci.priority.resize(sjdb_loci.chr.len(), 30);
        let raw_time = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|duration| duration.as_secs() as libc::time_t)
            .unwrap_or(0);
        result.log_main.push_str(&format!(
            "{}   Loaded database junctions from the generated genome {}: {} total junctions\n\n",
            timefunctions_l14_timemonthdaytime(raw_time),
            file_name,
            sjdb_loci.chr.len()
        ));
    }

    if p.two_pass_pass2 {
        let contents = std::fs::read_to_string(&p.two_pass_pass1sj_file).map_err(|_| {
            format!(
                "FATAL INPUT error, could not open input file with junctions from the 1st pass={}\n",
                p.two_pass_pass1sj_file
            )
        })?;
        sjdbloadfromstream_l2_sjdbloadfromstream(&contents, sjdb_loci);
        sjdb_loci.priority.resize(sjdb_loci.chr.len(), 0);
        let raw_time = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|duration| duration.as_secs() as libc::time_t)
            .unwrap_or(0);
        result.log_main.push_str(&format!(
            "{}   Loaded database junctions from the 1st pass file: {}: {} total junctions\n\n",
            timefunctions_l14_timemonthdaytime(raw_time),
            p.two_pass_pass1sj_file,
            sjdb_loci.chr.len()
        ));
    } else if p.run_mode_in.first().map(|s| s.as_str()) != Some("genomeGenerate") {
        result
            .log_main
            .push_str(&sjdbloadfromfiles_l6_sjdbloadfromfiles(p, sjdb_loci)?);

        let gtf_contents = if map_gen.p_ge.sjdb_gtf_file != "-" && map_gen.sjdb_overhang > 0 {
            Some(
                crate::io_utils::read_to_string_auto_gzip(&map_gen.p_ge.sjdb_gtf_file).map_err(
                    |_| {
                        format!(
                            "FATAL error, could not open file pGe.sjdbGTFfile={}\n",
                            map_gen.p_ge.sjdb_gtf_file
                        )
                    },
                )?,
            )
        } else {
            None
        };
        let (mut gtf, gtf_log) =
            gtf_l7_gtf_gtf(map_gen, p, &p.sjdb_insert_out_dir, gtf_contents.as_deref())?;
        result.log_main.push_str(&gtf_log);
        let gtf_out = gtf_transcriptgenesj_l23_gtf_transcriptgenesj(
            &mut gtf,
            map_gen,
            sjdb_loci,
            &p.sjdb_insert_out_dir,
            &mut result.log_main,
        );
        if gtf.gtf_yes {
            std::fs::create_dir_all(&p.sjdb_insert_out_dir).map_err(|e| e.to_string())?;
            for (name, contents) in [
                ("geneInfo.tab", gtf_out.gene_info_tab.as_str()),
                ("transcriptInfo.tab", gtf_out.transcript_info_tab.as_str()),
                ("exonInfo.tab", gtf_out.exon_info_tab.as_str()),
                ("exonGeTrInfo.tab", gtf_out.exon_ge_tr_info_tab.as_str()),
                (
                    "sjdbList.fromGTF.out.tab",
                    gtf_out.sjdb_list_from_gtf_out_tab.as_str(),
                ),
            ] {
                let path = format!("{}/{}", p.sjdb_insert_out_dir, name);
                std::fs::write(&path, contents).map_err(|e| e.to_string())?;
                result.files_written.push(path);
            }
        }
        result.gtf = Some(gtf_out);
    }

    let n_genome_real = map_gen
        .chr_start
        .get(map_gen.n_chr_real as usize)
        .copied()
        .unwrap_or(map_gen.n_genome);
    let prepare =
        sjdbprepare_l5_sjdbprepare(sjdb_loci, p, n_genome_real, &p.sjdb_insert_out_dir, map_gen)?;
    result.log_main.push_str(&prepare.log_main);
    if !p.sjdb_insert_out_dir.is_empty() {
        std::fs::create_dir_all(&p.sjdb_insert_out_dir).map_err(|err| err.to_string())?;
        for (name, contents) in [
            ("sjdbInfo.txt", prepare.sjdb_info_txt.as_str()),
            ("sjdbList.out.tab", prepare.sjdb_list_out_tab.as_str()),
        ] {
            let path = format!("{}/{}", p.sjdb_insert_out_dir, name);
            std::fs::write(&path, contents).map_err(|err| err.to_string())?;
            result.files_written.push(path);
        }
    }
    let raw_time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    result.log_main.push_str(&format!(
        "{}   Finished preparing junctions\n",
        timefunctions_l14_timemonthdaytime(raw_time)
    ));

    if map_gen.sjdb_n > p.limit_sjdb_insert_nsj {
        return Err(format!(
            "Fatal LIMIT error: the number of junctions to be inserted on the fly ={} is larger than the limitSjdbInsertNsj={}\nFatal LIMIT error: the number of junctions to be inserted on the fly ={} is larger than the limitSjdbInsertNsj={}\nSOLUTION: re-run with at least --limitSjdbInsertNsj {}\n",
            map_gen.sjdb_n,
            p.limit_sjdb_insert_nsj,
            map_gen.sjdb_n,
            p.limit_sjdb_insert_nsj,
            map_gen.sjdb_n
        ));
    }

    let mut gsj = prepare.gsj.clone();
    let mut g = map_gen.g.clone();
    let build = sjdbbuildindex_l16_sjdbbuildindex(p, &mut gsj, &mut g, map_gen, map_gen1)?;

    result.log_main.push_str(&build.log_main);
    let raw_time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as libc::time_t)
        .unwrap_or(0);
    result.log_main.push_str(&format!(
        "{} ..... finished inserting junctions into genome\n",
        timefunctions_l14_timemonthdaytime(raw_time)
    ));

    if map_gen.sa_packed.word_length > 0 && !map_gen.sa_packed.char_array.is_empty() {
        let b = map_gen.n_sa as u64 * map_gen.sa_packed.word_length;
        let need = (b / 8) as usize + 8;
        if map_gen.sa_packed.char_array.len() < need {
            map_gen.sa_packed.char_array.resize(need, 0);
        }
        packedarray_l17_packedarray_writepacked(&mut map_gen.sa_packed, map_gen.n_sa as u64, 0);
    }

    if map_gen.p_ge.sjdb_insert_save == "All" {
        std::fs::create_dir_all(&p.sjdb_insert_out_dir).map_err(|e| e.to_string())?;
        if p.p_ge.g_dir != p.sjdb_insert_out_dir {
            for name in [
                "chrName.txt",
                "chrStart.txt",
                "chrNameLength.txt",
                "chrLength.txt",
            ] {
                let src = format!("{}/{}", p.p_ge.g_dir, name);
                let dst = format!("{}/{}", p.sjdb_insert_out_dir, name);
                streamfuns_l144_copyfile(&src, &dst).map_err(|e| e.to_string())?;
                result.files_written.push(dst);
            }
        }

        map_gen.p_ge.g_file_sizes.clear();
        map_gen.p_ge.g_file_sizes.push(map_gen.n_genome);
        map_gen
            .p_ge
            .g_file_sizes
            .push(map_gen.sa_packed.length_byte);

        let genome_parameters = format!("{}/genomeParameters.txt", p.sjdb_insert_out_dir);
        genomeparameterswrite_l4_genomeparameterswrite(
            &genome_parameters,
            p,
            "ERROR_OUT",
            map_gen,
        )?;
        result.files_written.push(genome_parameters);

        let genome_file = format!("{}/Genome", p.sjdb_insert_out_dir);
        let mut genome_out = streamfuns_l91_ofstropen(&genome_file, "ERROR_OUT")?;
        streamfuns_l51_fstreamwritebig(&mut genome_out, &map_gen.g, map_gen.n_genome)
            .map_err(|e| e.to_string())?;
        result.files_written.push(genome_file);

        let sa_file = format!("{}/SA", p.sjdb_insert_out_dir);
        let mut sa_out = streamfuns_l91_ofstropen(&sa_file, "ERROR_OUT")?;
        streamfuns_l51_fstreamwritebig(
            &mut sa_out,
            &map_gen.sa_packed.char_array,
            map_gen.sa_packed.length_byte,
        )
        .map_err(|e| e.to_string())?;
        result.files_written.push(sa_file);

        let sa_index_file = format!("{}/SAindex", p.sjdb_insert_out_dir);
        let mut sa_index_out = streamfuns_l91_ofstropen(&sa_index_file, "ERROR_OUT")?;
        sa_index_out
            .write_all(&(map_gen.p_ge.g_saindex_nbases as u64).to_ne_bytes())
            .map_err(|e| e.to_string())?;
        for value in map_gen
            .genome_sa_index_start
            .iter()
            .take(map_gen.p_ge.g_saindex_nbases as usize + 1)
        {
            sa_index_out
                .write_all(&(*value as u64).to_ne_bytes())
                .map_err(|e| e.to_string())?;
        }
        streamfuns_l51_fstreamwritebig(
            &mut sa_index_out,
            &map_gen.sai_packed.char_array,
            map_gen.sai_packed.length_byte,
        )
        .map_err(|e| e.to_string())?;
        result.files_written.push(sa_index_file);
    }

    p.win_bin_n = (map_gen.n_genome / (1_u64 << p.win_bin_nbits)) as u32 + 1;

    result.sjdb_prepare = prepare;
    result.sjdb_build_index = build;
    Ok(result)
}

pub fn genome_sjdb_insert_snapshot(
    genome: &crate::genome::Genome,
) -> crate::genome::Genome {
    crate::genome::Genome {
        sjdb_n: genome.sjdb_n,
        sjdb_start: genome.sjdb_start.clone(),
        sjdb_end: genome.sjdb_end.clone(),
        n_sa: genome.n_sa,
        n_genome: genome.n_genome,
        ..Default::default()
    }
}
