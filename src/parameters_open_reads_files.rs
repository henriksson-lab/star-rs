#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Parameters::openReadsFiles` at STAR/source/Parameters_openReadsFiles.cpp:5. Args: "]
pub fn parameters_openreadsfiles_l5_parameters_openreadsfiles(
    p: &mut crate::parameters_chimeric::Parameters,
    existing_files: &std::collections::BTreeSet<String>,
) -> Result<crate::parameters_chimeric::OpenReadsFilesResult, String> {
    let mut result = crate::parameters_chimeric::OpenReadsFilesResult::default();

    if p.read_files_command_string.is_empty() {
        p.read_files_command_pid.resize(p.read_files_in.len(), 0);
        p.read_in_open.resize(p.read_files_in.len(), false);
        for ii in 0..p.read_files_in.len() {
            p.read_files_command_pid[ii] = 0;
            if p.read_in_open[ii] {
                p.read_in_open[ii] = false;
            }
            let rf_name = format!("{}{}", p.read_files_prefix_final, p.read_files_in[ii]);
            if !existing_files.contains(&rf_name) {
                return Err(format!(
                    "EXITING because of fatal input ERROR: could not open readFilesIn={}\n",
                    rf_name
                ));
            }
            p.read_in_open[ii] = true;
            result.opened_inputs.push(rf_name);
        }
    } else {
        p.read_files_in_tmp.clear();
        p.read_files_command_pid.resize(p.read_files_names.len(), 0);
        p.read_in_open.resize(p.read_files_names.len(), false);

        for imate in 0..p.read_files_names.len() {
            let fifo = format!("{}tmp.fifo.read{}", p.out_file_tmp, imate + 1);
            p.read_files_in_tmp.push(fifo.clone());
            result.read_files_in_tmp.push(fifo.clone());
            result
                .log_main
                .push_str(&format!("\n   Input read files for mate {} :\n", imate + 1));

            let command_file_name = format!("{}/readsCommand_read{}", p.out_file_tmp, imate + 1);
            result
                .reads_command_file_names
                .push(command_file_name.clone());
            let mut command_file = String::new();
            if p.sys_shell != "-" && !p.sys_shell.is_empty() {
                command_file.push_str(&format!("#!{}\n", p.sys_shell));
            }
            command_file.push_str(&format!("exec > \"{}\"\n", fifo));

            for ifile in 0..p.read_files_n as usize {
                let Some(file_name) = p.read_files_names.get(imate).and_then(|v| v.get(ifile))
                else {
                    return Err(format!(
                        "EXITING: because of fatal INPUT file error: could not open read file: \nSOLUTION: check that this file exists and has read permision.\n"
                    ));
                };

                if existing_files.contains(file_name) {
                    result
                        .log_main
                        .push_str(&format!("-rw-r--r-- 1 user group 0 {}\n", file_name));
                } else {
                    result
                        .log_main
                        .push_str(&format!(" Could not ls {}\n", file_name));
                    return Err(format!(
                        "EXITING: because of fatal INPUT file error: could not open read file: {}\nSOLUTION: check that this file exists and has read permision.\n",
                        file_name
                    ));
                }

                command_file.push_str(&format!("echo FILE {}\n", ifile));
                command_file.push_str(&format!(
                    "{}   \"{}\"\n",
                    p.read_files_command_string, file_name
                ));
            }

            result
                .log_main
                .push_str(&format!("\n   readsCommandsFile:\n{}\n", command_file));
            result.reads_command_file_contents.push(command_file);
            p.read_files_command_pid[imate] = (imate + 1) as i32;
            result.command_pids.push(p.read_files_command_pid[imate]);
            p.read_in_open[imate] = true;
            result.opened_inputs.push(fifo);
        }
    }

    p.read_files_index = 0;
    if p.read_files_type_n == 10 {
        let names = p.read_files_names.first().cloned().unwrap_or_default();
        let commands = parameters_readsamheader_l6_parameters_readsamheader(
            p,
            &p.read_files_command_string.clone(),
            &names,
            &[],
        )?;
        result.sam_header_commands = commands;
    }

    Ok(result)
}
