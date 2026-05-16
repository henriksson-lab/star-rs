#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `sjdbLoadFromFiles` at STAR/source/sjdbLoadFromFiles.cpp:6. Args: P: Parameters, sjdbLoci: SjdbClass"]
pub fn sjdbloadfromfiles_l6_sjdbloadfromfiles(
    p: &crate::parameters_chimeric::Parameters,
    sjdb_loci: &mut crate::sjdb_class::SjdbClass,
) -> Result<String, String> {
    let mut log_main = String::new();

    if p.p_ge.sjdb_file_chr_start_end.first().map(|v| v.as_str()) != Some("-") {
        for file_name in &p.p_ge.sjdb_file_chr_start_end {
            let contents = std::fs::read_to_string(file_name).map_err(|_| {
                format!(
                    "FATAL INPUT error, could not open input file pGe.sjdbFileChrStartEnd={}\n",
                    file_name
                )
            })?;

            sjdbloadfromstream_l2_sjdbloadfromstream(&contents, sjdb_loci);
            sjdb_loci.priority.resize(sjdb_loci.chr.len(), 10);

            let raw_time = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|duration| duration.as_secs() as libc::time_t)
                .unwrap_or(0);
            log_main.push_str(&format!(
                "{}   Loaded database junctions from the pGe.sjdbFileChrStartEnd file(s), total number of junctions:{}\n\n",
                timefunctions_l14_timemonthdaytime(raw_time),
                sjdb_loci.chr.len()
            ));
        }
    }

    Ok(log_main)
}
