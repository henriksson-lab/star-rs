#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ParametersClip::initialize` at STAR/source/ParametersClip_initialize.cpp:6. Args: pPin: Parameters"]
pub fn parametersclip_initialize_l6_parametersclip_initialize(
    clip: &mut crate::parameters_clip::ParametersClip,
    p: &crate::parameters_chimeric::Parameters,
) -> Result<(), String> {
    clip.read_nmates = p.read_nmates;
    clip.read_nends = p.read_nends;

    if clip.adapter_type[0] != "Hamming" && clip.adapter_type[0] != "CellRanger4" {
        return Err(format!(
            "EXITING because of fatal PARAMETER error: --clipAdapterType = {} is not a valid option\nSOLUTION: use valid --clipAdapterType options: Hamming OR CellRanger4\n",
            clip.adapter_type[0]
        ));
    }

    if clip.adapter_type[0] == "CellRanger4" {
        if clip.in_[1].ad_seq.len() > 1 || clip.in_[1].ad_seq[0] != "-" {
            return Err(
                "EXITING because of fatal PARAMETER error: --clipAdapterType CellRanger4 uses fixed sequences for 3' polyA adapters.\nSOLUTION: Do not use --clip5pAdapter* options.\n"
                    .to_string(),
            );
        }
        clip.in_[1].ad_seq[0] = "A".to_string();

        if clip.in_[0].ad_seq.len() > 1 && clip.in_[0].ad_seq[1] != "-" {
            return Err(
                "EXITING because of fatal PARAMETER error: when using --clipAdapterType CellRanger4, only 5' adapter for the 1st mate can be specified.\nSOLUTION: Use only one sequence in --clip5pAdapterSeq (or '-' instead of 2nd sequence).\n"
                    .to_string(),
            );
        }

        if clip.in_[0].ad_seq[0] == "-" {
            clip.in_[0].ad_seq[0] = "AAGCAGTGGTATCAACGCAGAGTACATGGG".to_string();
        }
    } else {
        for im in 0..clip.in_[0].ad_seq.len() {
            if clip.in_[0].ad_seq[im] != "-" {
                return Err(
                    "EXITING because of fatal PARAMETER error: --clip5pAdapterSeq is not supported yet, except for --clipAdapterType CellRanger4.\nSOLUTION: Do not use --clip5pAdapter* options without --clipAdapterType CellRanger4.\n"
                        .to_string(),
                );
            }
        }
    }

    for ip in 0..2 {
        if clip.in_[ip].ad_seq[0] == "-" {
            clip.in_[ip]
                .ad_seq
                .resize(p.read_nmates as usize, "-".to_string());
            clip.in_[ip].ad_mmp.resize(p.read_nmates as usize, 0.0);
        }

        if clip.in_[ip].n[0] == 0 {
            clip.in_[ip].n.resize(p.read_nmates as usize, 0);
        }

        if clip.in_[ip].n_after_ad[0] == 0 {
            clip.in_[ip].n_after_ad.resize(p.read_nmates as usize, 0);
        }
    }

    let p53 = ["5", "3"];
    for ip in 0..2 {
        if clip.in_[ip].ad_seq.len() != p.read_nmates as usize {
            return Err(format!(
                "EXITING because of fatal PARAMETER error: --clip{}pAdapterSeq has to contain {} values to match the number of mates.\nSOLUTION: specify {}values in --clip{}pAdapterSeq , for no clipping use -",
                p53[ip], p.read_nmates, p.read_nmates, p53[ip]
            ));
        }

        if clip.in_[ip].ad_mmp.len() != p.read_nmates as usize {
            return Err(format!(
                "EXITING because of fatal PARAMETER error: --clip{}pAdapterMMp has to contain {} values to match the number of mates.\nSOLUTION: specify {}values in --clip{}pAdapterMMp",
                p53[ip], p.read_nmates, p.read_nmates, p53[ip]
            ));
        }

        if clip.in_[ip].n_after_ad.len() != p.read_nmates as usize {
            return Err(format!(
                "EXITING because of fatal PARAMETER error: --clip{}pAfterAdapterNbases has to contain {} values to match the number of mates.\nSOLUTION: specify {}values in --clip{}pAfterAdapterNbases , for no clipping use 0",
                p53[ip], p.read_nmates, p.read_nmates, p53[ip]
            ));
        }

        if clip.in_[ip].n.len() != p.read_nmates as usize {
            return Err(format!(
                "EXITING because of fatal PARAMETER error: --clip{}pNbases has to contain {} values to match the number of mates.\nSOLUTION: specify {}values in --clip{}pNbases , for no clipping use 0",
                p53[ip], p.read_nmates, p.read_nmates, p53[ip]
            ));
        }
    }
    Ok(())
}

#[doc = "Original `ParametersClip::initializeClipMates` at STAR/source/ParametersClip_initialize.cpp:84. Args: clipMates: vector<vector<ClipMate>>"]
pub fn parametersclip_initialize_l84_parametersclip_initializeclipmates(
    clip: &crate::parameters_clip::ParametersClip,
    clip_mates: &mut Vec<Vec<crate::clip_mate::ClipMate>>,
) {
    clip_mates.resize(clip.read_nends as usize, Vec::new());

    for im in 0..clip_mates.len() {
        clip_mates[im].resize(2, crate::clip_mate::ClipMate::default());

        for ip in 0..2 {
            clip_mates[im][ip].type_ = ip as i32;

            if clip.adapter_type[0] == "CellRanger4" {
                clip_mates[im][ip].type_ += 10;
            }

            if im < clip.read_nmates as usize {
                clipmate_initialize_l5_clipmate_initialize(
                    &mut clip_mates[im][ip],
                    clip.in_[ip].n[im],
                    &clip.in_[ip].ad_seq[im],
                    clip.in_[ip].n_after_ad[im],
                    clip.in_[ip].ad_mmp[im],
                );
            } else {
                clipmate_initialize_l5_clipmate_initialize(&mut clip_mates[im][ip], 0, "-", 0, 0.0);
            }
        }
    }
}
