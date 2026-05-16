#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Parameters::readFilesInit` at STAR/source/Parameters_readFilesInit.cpp:8. Args: "]
pub fn parameters_readfilesinit_l8_parameters_readfilesinit(
    p: &mut crate::parameters_chimeric::Parameters,
    manifest_contents: Option<&str>,
) -> Result<String, String> {
    let mut log_main = String::new();

    match p.read_files_type.first().map(|s| s.as_str()) {
        Some("Fastx") => {
            p.read_files_type_n = 1;
        }
        Some("SAM") => {
            p.read_files_type_n = 10;
            p.read_files_sam_attr_keep_all = false;
            p.read_files_sam_attr_keep_none = false;
            if p.read_files_sam_attr_keep_in.first().map(|s| s.as_str()) == Some("All") {
                p.read_files_sam_attr_keep_all = true;
            } else if p.read_files_sam_attr_keep_in.first().map(|s| s.as_str()) == Some("None") {
                p.read_files_sam_attr_keep_none = true;
            } else {
                p.read_files_sam_attr_keep.clear();
                for tag in &p.read_files_sam_attr_keep_in {
                    if tag.len() != 2 {
                        return Err("EXITING because of FATAL PARAMETER ERROR: each SAM tags in --readFilesSAMtagsKeep should contain two letters\n                                  SOLUTION: specify only two-letter tags in --readFilesSAMtagsKeep.".to_string());
                    }
                    let bytes = tag.as_bytes();
                    p.read_files_sam_attr_keep
                        .insert(u16::from_ne_bytes([bytes[0], bytes[1]]));
                }
            }
        }
        Some(other) => {
            return Err(format!(
                "EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --readFilesType: {}\nSOLUTION: specify one of the allowed values: Fastx or SAM\n",
                other
            ));
        }
        None => {
            return Err("EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --readFilesType: \nSOLUTION: specify one of the allowed values: Fastx or SAM\n".to_string());
        }
    }

    p.read_files_prefix_final = if p.read_files_prefix == "-" {
        String::new()
    } else {
        p.read_files_prefix.clone()
    };

    if p.read_files_manifest.first().map(|s| s.as_str()) == Some("-") {
        p.read_files_names.resize(p.read_files_in.len(), Vec::new());

        for imate in 0..p.read_files_names.len() {
            p.read_files_names[imate] = p.read_files_in[imate]
                .split(',')
                .map(|s| s.to_string())
                .collect();
            if p.read_files_names[imate]
                .last()
                .map(|s| s.is_empty())
                .unwrap_or(false)
            {
                p.read_files_names[imate].pop();
            }

            if imate > 0 && p.read_files_names[imate].len() != p.read_files_names[imate - 1].len() {
                return Err(format!(
                    "EXITING: because of fatal INPUT ERROR: number of input files for mate{}={} is not equal to that for mate{}={}\nMake sure that the number of files in --readFilesIn is the same for both mates\n",
                    imate + 1,
                    p.read_files_names[imate].len(),
                    imate - 1,
                    p.read_files_names[imate - 1].len()
                ));
            }

            for fn1 in &mut p.read_files_names[imate] {
                *fn1 = format!("{}{}", p.read_files_prefix_final, fn1);
            }
        }

        p.read_files_n = p.read_files_names.first().map(|v| v.len()).unwrap_or(0) as u32;

        if p.out_sam_attr_rgline.first().map(|s| s.as_str()) != Some("-") {
            p.out_sam_attr_rgline_split.clear();
            p.out_sam_attr_rg.clear();
            let mut ii = 0usize;
            while ii < p.out_sam_attr_rgline.len() {
                if ii == 0 || p.out_sam_attr_rgline[ii] == "," {
                    if ii > 0 {
                        ii += 1;
                    }
                    let Some(first_field) = p.out_sam_attr_rgline.get(ii) else {
                        break;
                    };
                    if !first_field.starts_with("ID:") {
                        return Err(format!(
                            "EXITING because of FATAL INPUT ERROR: the first word of a line from --outSAMattrRGline={} does not start with ID:xxx read group identifier\nSOLUTION: re-run STAR with all lines in --outSAMattrRGline starting with ID:xxx\n",
                            first_field
                        ));
                    }
                    p.out_sam_attr_rgline_split.push(first_field.clone());
                    p.out_sam_attr_rg.push(first_field[3..].to_string());
                } else if let Some(last) = p.out_sam_attr_rgline_split.last_mut() {
                    last.push('\t');
                    last.push_str(&p.out_sam_attr_rgline[ii]);
                }
                ii += 1;
            }
        }

        if p.out_sam_attr_rg.len() > 1 && p.out_sam_attr_rg.len() != p.read_files_n as usize {
            return Err(format!(
                "EXITING: because of fatal INPUT ERROR: number of input read files: {} does not agree with number of read group RG entries: {}\nMake sure that the number of RG lines in --outSAMattrRGline is equal to either 1, or the number of input read files in --readFilesIn\n",
                p.read_files_n,
                p.out_sam_attr_rg.len()
            ));
        } else if p.out_sam_attr_rg.len() == 1 {
            for _ in 1..p.read_files_n {
                let rg = p.out_sam_attr_rg[0].clone();
                p.out_sam_attr_rg.push(rg);
            }
        }
    } else {
        let manifest_name = p.read_files_manifest.first().cloned().unwrap_or_default();
        let manifest = manifest_contents.ok_or_else(|| {
            format!(
                "EXITING because of FATAL INPUT FILE error: could not open readFileManifest file {}",
                manifest_name
            )
        })?;
        log_main.push_str(&format!(
            "Reading input file names and read groups from readFileManifest {}\n",
            manifest_name
        ));

        p.read_files_names = vec![Vec::new(), Vec::new()];
        p.out_sam_attr_rgline_split.clear();
        p.out_sam_attr_rg.clear();
        for line in manifest.lines() {
            if line.find(|c| c != ' ' && c != '\t').is_none() {
                continue;
            }
            let mut fields = line.split('\t');
            let f0 = fields.next().unwrap_or("");
            let f1 = fields.next().ok_or_else(|| {
                format!(
                    "EXITING because of FATAL INPUT FILE error: readFileManifest file {} has to contain at least 3 tab separated columns\nSOLUTION: fix the formatting of the readFileManifest file: Read1 <tab> Read2 <tab> ReadGroup. For single-end reads, use - in the 2nd column.\n",
                    manifest_name
                )
            })?;
            let rg_rest = fields.collect::<Vec<_>>().join("\t");
            if rg_rest.is_empty() {
                return Err(format!(
                    "EXITING because of FATAL INPUT FILE error: readFileManifest file {} has to contain at least 3 tab separated columns\nSOLUTION: fix the formatting of the readFileManifest file: Read1 <tab> Read2 <tab> ReadGroup. For single-end reads, use - in the 2nd column.\n",
                    manifest_name
                ));
            }

            p.read_files_names[0].push(format!("{}{}", p.read_files_prefix_final, f0));
            p.read_files_names[1].push(format!("{}{}", p.read_files_prefix_final, f1));
            log_main.push_str(&format!(
                "{}\t{}\t",
                p.read_files_names[0].last().unwrap(),
                p.read_files_names[1].last().unwrap()
            ));

            let mut rg_line = rg_rest;
            if !rg_line.starts_with("ID:") {
                rg_line.insert_str(0, "ID:");
            }
            let tab_pos = rg_line.find('\t').unwrap_or(rg_line.len());
            p.out_sam_attr_rg.push(rg_line[3..tab_pos].to_string());
            log_main.push_str(&rg_line);
            log_main.push('\n');
            p.out_sam_attr_rgline_split.push(rg_line);
        }

        p.read_nends = if p
            .read_files_names
            .get(1)
            .and_then(|v| v.first())
            .and_then(|s| s.as_bytes().last().copied())
            == Some(b'-')
        {
            1
        } else {
            2
        };
        p.read_files_names.truncate(p.read_nends as usize);
        p.read_files_n = p.read_files_names.first().map(|v| v.len()).unwrap_or(0) as u32;
    }

    log_main.push_str(&format!(
        "Number of fastq files for each mate = {}\n",
        p.read_files_n
    ));

    p.read_files_command_string.clear();
    if p.read_files_command.first().map(|s| s.as_str()) == Some("-") {
        if p.read_files_n > 1 {
            p.read_files_command_string = "cat   ".to_string();
        }
    } else {
        for cmd in &p.read_files_command {
            p.read_files_command_string.push_str(cmd);
            p.read_files_command_string.push_str("   ");
        }
    }

    if p.read_files_type_n == 1 {
        p.read_nends = p.read_files_names.len() as u32;
    } else if p.read_files_type_n == 10 {
        if p.read_files_type.len() == 2 && p.read_files_type[1] == "SE" {
            p.read_nends = 1;
        } else if p.read_files_type.len() == 2 && p.read_files_type[1] == "PE" {
            p.read_nends = 2;
        } else {
            return Err("EXITING because of FATAL INPUT ERROR: --readFilesType SAM requires specifying SE or PE reads\nSOLUTION: specify --readFilesType SAM SE for single-end reads or --readFilesType SAM PE for paired-end reads\n".to_string());
        }
    }

    p.read_nmates = p.read_nends;
    Ok(log_main)
}
