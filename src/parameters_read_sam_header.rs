#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Parameters::readSAMheader` at STAR/source/Parameters_readSAMheader.cpp:6. Args: readFilesCommandString: string, readFilesNames: vector<string>"]
pub fn parameters_readsamheader_l6_parameters_readsamheader(
    p: &mut crate::parameters_chimeric::Parameters,
    read_files_command_string: &str,
    read_files_names: &[String],
    command_outputs: &[String],
) -> Result<Vec<String>, String> {
    if read_files_command_string.is_empty() {
        while p.read_in_0.as_bytes().get(p.read_in_0_pos).copied() == Some(b'@') {
            let line_start = p.read_in_0_pos;
            let line_end = p.read_in_0[line_start..]
                .find('\n')
                .map(|offset| line_start + offset)
                .unwrap_or(p.read_in_0.len());
            let str1 = &p.read_in_0[line_start..line_end];
            p.read_in_0_pos = if line_end < p.read_in_0.len() {
                line_end + 1
            } else {
                line_end
            };
            if str1.get(1..3) != Some("HD") && str1.get(1..3) != Some("SQ") {
                p.sam_header_extra.push_str(str1);
                p.sam_header_extra.push('\n');
            }
        }
        return Ok(Vec::new());
    }

    let tmp_fifo = format!("{}tmp.fifo.header", p.out_file_tmp);
    let mut commands = Vec::new();
    for ii in 0..read_files_names.len() {
        let com1 = format!(
            "{}   {} > {}&",
            read_files_command_string, read_files_names[ii], tmp_fifo
        );
        commands.push(com1);

        let tmp_fifo_in = command_outputs.get(ii).map(|s| s.as_str()).unwrap_or("");
        let mut pos = 0usize;
        while tmp_fifo_in.as_bytes().get(pos).copied() == Some(b'@') {
            let line_start = pos;
            let line_end = tmp_fifo_in[line_start..]
                .find('\n')
                .map(|offset| line_start + offset)
                .unwrap_or(tmp_fifo_in.len());
            let str1 = &tmp_fifo_in[line_start..line_end];
            pos = if line_end < tmp_fifo_in.len() {
                line_end + 1
            } else {
                line_end
            };
            if str1.get(1..3) != Some("HD") && str1.get(1..3) != Some("SQ") && !p.two_pass_pass2 {
                p.sam_header_extra.push_str(str1);
                p.sam_header_extra.push('\n');
            }
        }
    }
    Ok(commands)
}
