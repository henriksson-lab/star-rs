#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloReadBarcode::matchCBtoWL` at STAR/source/SoloReadBarcode_getCBandUMI.cpp:9. Args: cbSeq1: string, cbQual1: string, cbWL: vector<uint64>, cbMatch1: int32, cbMatchInd1: vector<uint64>, cbMatchString1: string"]
pub fn soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(
    p_solo: &crate::parameters_solo::ParametersSolo,
    cb_seq1: &str,
    cb_qual1: &str,
    cb_wl: &[u64],
) -> (i32, Vec<u64>, String) {
    let mut cb_match1 = -1;
    let mut cb_match_string1 = String::new();
    let mut cb_match_ind1 = Vec::new();
    let mut cb_b1 = 0u64;
    let pos_n = sequencefuns_l249_convertnuclstrtoint64(cb_seq1, &mut cb_b1);

    if !p_solo.cb_wl_yes {
        if pos_n != -1 {
            cb_match1 = -2;
        } else {
            cb_match_ind1.push(cb_b1);
            cb_match_string1 = cb_b1.to_string();
            cb_match1 = 0;
        }
        return (cb_match1, cb_match_ind1, cb_match_string1);
    }

    if pos_n == -2 {
        cb_match1 = -2;
        return (cb_match1, cb_match_ind1, cb_match_string1);
    } else if pos_n == -1 && !cb_wl.is_empty() {
        let cb_i = servicefuns_l294_binarysearchexact(cb_b1, cb_wl, cb_wl.len() as u64);
        if cb_i >= 0 {
            cb_match_ind1.push(cb_i as u64);
            cb_match_string1 = cb_match_ind1[0].to_string();
            cb_match1 = 0;
            return (cb_match1, cb_match_ind1, cb_match_string1);
        }
    }

    if !p_solo.cb_match_wl.mm1 {
        return (cb_match1, cb_match_ind1, cb_match_string1);
    }

    cb_match1 = 0;
    if pos_n >= 0 {
        let pos_n_usize = pos_n as usize;
        let pos_n_shift = 2 * (cb_seq1.len() - 1 - pos_n_usize);
        let mut matched = false;
        for jj in 0..4u64 {
            let cb_b11 = cb_b1 ^ (jj << pos_n_shift);
            let cb_i1 = if cb_wl.is_empty() {
                -1
            } else {
                servicefuns_l294_binarysearchexact(cb_b11, cb_wl, cb_wl.len() as u64)
            };
            if cb_i1 >= 0 {
                if !p_solo.cb_match_wl.mm1_multi_nbase && matched {
                    cb_match_ind1.clear();
                    cb_match1 = -3;
                    break;
                }
                matched = true;
                cb_match_ind1.push(cb_i1 as u64);
                cb_match1 += 1;
                cb_match_string1.push(' ');
                cb_match_string1.push_str(&cb_i1.to_string());
                cb_match_string1.push(' ');
                cb_match_string1.push(cb_qual1.as_bytes()[pos_n_usize] as char);
            }
        }
    } else {
        for ii in 0..cb_seq1.len() {
            for jj in 1..4u64 {
                let cb_i1 = if cb_wl.is_empty() {
                    -1
                } else {
                    servicefuns_l294_binarysearchexact(
                        cb_b1 ^ (jj << (ii * 2)),
                        cb_wl,
                        cb_wl.len() as u64,
                    )
                };
                if cb_i1 >= 0 {
                    cb_match_ind1.push(cb_i1 as u64);
                    cb_match1 += 1;
                    cb_match_string1.push(' ');
                    cb_match_string1.push_str(&cb_i1.to_string());
                    cb_match_string1.push(' ');
                    cb_match_string1.push(cb_qual1.as_bytes()[cb_seq1.len() - 1 - ii] as char);
                }
            }
        }
    }

    if cb_match1 == 0 {
        cb_match1 = -1;
    } else if cb_match1 == 1 {
        cb_match_string1 = cb_match_ind1[0].to_string();
    } else if !p_solo.cb_match_wl.mm1_multi {
        cb_match1 = -3;
        cb_match_ind1.clear();
        cb_match_string1.clear();
    }

    (cb_match1, cb_match_ind1, cb_match_string1)
}

#[doc = "Original `SoloReadBarcode::addStats` at STAR/source/SoloReadBarcode_getCBandUMI.cpp:93. Args: cbMatch1: int32"]
pub fn soloreadbarcode_getcbandumi_l93_soloreadbarcode_addstats(
    solo_read_barcode: &mut crate::solo_read_barcode::SoloReadBarcode,
    cb_match1: i32,
) {
    if !solo_read_barcode.cb_wl_yes {
        return;
    }

    match cb_match1 {
        0 => {
            solo_read_barcode.cb_read_count_exact[solo_read_barcode.cb_match_ind[0] as usize] += 1;
            solo_read_barcode.stats.v[9] += 1;
        }
        1 => solo_read_barcode.stats.v[10] += 1,
        -1 => solo_read_barcode.stats.v[6] += 1,
        -2 => solo_read_barcode.stats.v[3] += 1,
        -3 => solo_read_barcode.stats.v[8] += 1,
        -11 => solo_read_barcode.stats.v[2] += 1,
        -12 => solo_read_barcode.stats.v[7] += 1,
        -23 => solo_read_barcode.stats.v[4] += 1,
        -24 => solo_read_barcode.stats.v[5] += 1,
        _ => solo_read_barcode.stats.v[11] += 1,
    }
}

#[doc = "Original `SoloReadBarcode::convertCheckUMI` at STAR/source/SoloReadBarcode_getCBandUMI.cpp:133. Args: "]
pub fn soloreadbarcode_getcbandumi_l133_soloreadbarcode_convertcheckumi(
    solo_read_barcode: &mut crate::solo_read_barcode::SoloReadBarcode,
) -> bool {
    if sequencefuns_l249_convertnuclstrtoint64(
        &solo_read_barcode.umi_seq,
        &mut solo_read_barcode.umi_b,
    ) != -1
    {
        solo_read_barcode.umi_check = -23;
        return false;
    }
    if solo_read_barcode.umi_b == solo_read_barcode.homo_polymer[0]
        || solo_read_barcode.umi_b == solo_read_barcode.homo_polymer[1]
        || solo_read_barcode.umi_b == solo_read_barcode.homo_polymer[2]
        || solo_read_barcode.umi_b == solo_read_barcode.homo_polymer[3]
    {
        solo_read_barcode.umi_check = -24;
        return false;
    }
    true
}

#[doc = "Original `SoloReadBarcode::getCBandUMI` at STAR/source/SoloReadBarcode_getCBandUMI.cpp:147. Args: readSeq: char, readQual: char, readLen: uint64, readNameExtraIn: string, readFilesIndex: uint32, readName: char"]
pub fn soloreadbarcode_getcbandumi_l147_soloreadbarcode_getcbandumi(
    solo_read_barcode: &mut crate::solo_read_barcode::SoloReadBarcode,
    p: &mut crate::parameters_chimeric::Parameters,
    read_seq: &[String],
    read_qual: &[String],
    read_len: &[u64],
    read_name_extra_in: &str,
    read_files_index: u32,
    read_name: &str,
) -> Result<(), String> {
    if p.p_solo.solo_type == SOLO_TYPE_NONE {
        return Ok(());
    }

    if p.p_solo.solo_type == SOLO_TYPE_SMART_SEQ {
        solo_read_barcode.cb_seq.clear();
        solo_read_barcode.cb_qual.clear();
        solo_read_barcode.cb_seq_corrected.clear();
        solo_read_barcode.cb_match = 0;
        solo_read_barcode.cb_match_ind = vec![read_files_index as u64];
        solo_read_barcode.cb_match_string = solo_read_barcode.cb_match_ind[0].to_string();
        soloreadbarcode_getcbandumi_l93_soloreadbarcode_addstats(
            solo_read_barcode,
            solo_read_barcode.cb_match,
        );
        return Ok(());
    }

    solo_read_barcode.cb_match = -1;
    solo_read_barcode.cb_match_string.clear();
    solo_read_barcode.cb_match_ind.clear();

    if p.read_files_type_n != 10 {
        let barcode_read = p.p_solo.barcode_read as usize;
        let len = read_len[barcode_read] as usize;
        solo_read_barcode.b_seq = read_seq[barcode_read][..len].to_string();
        solo_read_barcode.b_qual = read_qual[barcode_read][..len].to_string();
    } else {
        let pos1 = read_name_extra_in.find(|c: char| c != ' ' && c != '\t');
        let mut read_name_extra_t = read_name_extra_in.to_string();
        if pos1 == Some(0) {
            read_name_extra_t.insert(0, '\t');
        } else if let Some(pos1) = pos1 {
            read_name_extra_t.replace_range(pos1..pos1 + 1, "\t");
        }

        solo_read_barcode.b_seq.clear();
        solo_read_barcode.b_strings.clear();
        for tag in &p.p_solo.sam_attr_barcode_seq {
            let Some(pos_tag) = read_name_extra_t.find(tag) else {
                return Err(format!(
                    "EXITING because of FATAL ERROR in input read file: could not find barcode sequence SAM attribute {} in read {}\nwith SAM attributes: {}\nSOLUTION: make sure that all reads in the input SAM/BAM have all attributes from --soloInputSAMattrBarcodeSeq\n",
                    tag, read_name, read_name_extra_t
                ));
            };
            let pos_start = pos_tag + 6;
            let pos_end = read_name_extra_t[pos_start..]
                .find('\t')
                .map(|offset| pos_start + offset)
                .unwrap_or(read_name_extra_t.len());
            solo_read_barcode
                .b_strings
                .push(read_name_extra_t[pos_start..pos_end].to_string());
            solo_read_barcode
                .b_seq
                .push_str(solo_read_barcode.b_strings.last().unwrap());
        }

        solo_read_barcode.b_qual.clear();
        if p.p_solo.sam_attr_barcode_qual.is_empty() {
            solo_read_barcode
                .b_qual
                .extend(std::iter::repeat('H').take(solo_read_barcode.b_seq.len()));
        } else {
            for tag in &p.p_solo.sam_attr_barcode_qual {
                let Some(pos_tag) = read_name_extra_t.find(tag) else {
                    return Err(format!(
                        "EXITING because of FATAL ERROR in input read file: could not find barcode qualities SAM attribute {} in read {}\nwith SAM attributes: {}\nSOLUTION: make sure that all reads in the input SAM/BAM have all attributes from --soloInputSAMattrBarcodeQual\n",
                        tag, read_name, read_name_extra_t
                    ));
                };
                let pos_start = pos_tag + 6;
                let pos_end = read_name_extra_t[pos_start..]
                    .find('\t')
                    .map(|offset| pos_start + offset)
                    .unwrap_or(read_name_extra_t.len());
                solo_read_barcode
                    .b_qual
                    .push_str(&read_name_extra_t[pos_start..pos_end]);
            }
        }

        if solo_read_barcode.b_qual.len() != solo_read_barcode.b_seq.len() {
            return Err(format!(
                "EXITING because of FATAL ERROR in input read file: the total length of barcode qualities is {} not equal to the sequence length {}\nRead ID={} ;  Qualities={} ;  Sequence={} ;  Read SAM attributes: {}\nSOLUTION: make sure correct attributes are listed in --soloInputSAMattrBarcodeQual\n",
                solo_read_barcode.b_qual.len(),
                solo_read_barcode.b_seq.len(),
                read_name,
                solo_read_barcode.b_qual,
                solo_read_barcode.b_seq,
                read_name_extra_t
            ));
        }
    }

    if solo_read_barcode.b_seq.len() != p.p_solo.b_l as usize {
        if p.p_solo.b_l > 0 {
            return Err(format!(
                "EXITING because of FATAL ERROR in input read file: the total length of barcode sequence is {} not equal to expected {}\nRead ID={} ;  Sequence={}\nSOLUTION: check the formatting of input read files.\nIf UMI+CB length is not equal to the barcode read length, specify barcode read length with --soloBarcodeReadLength\nTo avoid checking of barcode read length, specify --soloBarcodeReadLength 0",
                solo_read_barcode.b_seq.len(),
                p.p_solo.b_l,
                read_name,
                solo_read_barcode.b_seq
            ));
        } else if solo_read_barcode.b_seq.len() < p.p_solo.cbumi_l as usize {
            let missing = p.p_solo.cbumi_l as usize - solo_read_barcode.b_seq.len();
            solo_read_barcode
                .b_seq
                .extend(std::iter::repeat('N').take(missing));
            solo_read_barcode
                .b_qual
                .extend(std::iter::repeat('H').take(missing));
        }
    }

    if p.p_solo.solo_type != SOLO_TYPE_CB_UMI_SIMPLE {
        let hist_len = if p.p_solo.b_l > 0 {
            p.p_solo.b_l as usize
        } else {
            solo_read_barcode.b_qual.len()
        };
        for ix in 0..hist_len {
            let q = solo_read_barcode.b_qual.as_bytes()[ix] as usize;
            solo_read_barcode.qual_hist[q] += 1;
        }
    }

    if p.p_solo.solo_type == SOLO_TYPE_CB_UMI_SIMPLE {
        if p.p_solo.cb_type_type == 1 {
            let cb_s = p.p_solo.cb_s as usize - 1;
            let cb_e = cb_s + p.p_solo.cb_l as usize;
            let umi_s = p.p_solo.umi_s as usize - 1;
            let umi_e = umi_s + p.p_solo.umi_l as usize;
            solo_read_barcode.cb_seq = solo_read_barcode.b_seq[cb_s..cb_e].to_string();
            solo_read_barcode.umi_seq = solo_read_barcode.b_seq[umi_s..umi_e].to_string();
            solo_read_barcode.cb_qual = solo_read_barcode.b_qual[cb_s..cb_e].to_string();
            solo_read_barcode.umi_qual = solo_read_barcode.b_qual[umi_s..umi_e].to_string();

            for q in solo_read_barcode.cb_qual.as_bytes() {
                solo_read_barcode.qual_hist[*q as usize] += 1;
            }
            for q in solo_read_barcode.umi_qual.as_bytes() {
                solo_read_barcode.qual_hist[*q as usize] += 1;
            }

            let (cb_match, cb_match_ind, cb_match_string) =
                soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(
                    &p.p_solo,
                    &solo_read_barcode.cb_seq,
                    &solo_read_barcode.cb_qual,
                    &p.p_solo.cb_wl,
                );
            solo_read_barcode.cb_match = cb_match;
            solo_read_barcode.cb_match_ind = cb_match_ind;
            solo_read_barcode.cb_match_string = cb_match_string;
        } else if p.p_solo.cb_type_type == 2 {
            solo_read_barcode.cb_seq = solo_read_barcode.b_strings[0].clone();
            solo_read_barcode.umi_seq = solo_read_barcode.b_strings[1].clone();
            let cb1 = if let Some(cb1) = p.p_solo.cb_type_str_map.get(&solo_read_barcode.cb_seq) {
                *cb1
            } else {
                let cb1 = p.p_solo.cb_type_str_map.len() as u32;
                p.p_solo
                    .cb_type_str_map
                    .insert(solo_read_barcode.cb_seq.clone(), cb1);
                cb1
            };
            solo_read_barcode.cb_match_ind.push(cb1 as u64);
            solo_read_barcode.cb_match_string = cb1.to_string();
            solo_read_barcode.cb_match = 0;
        }

        if !soloreadbarcode_getcbandumi_l133_soloreadbarcode_convertcheckumi(solo_read_barcode) {
            solo_read_barcode.cb_match = solo_read_barcode.umi_check;
            solo_read_barcode.cb_match_string.clear();
            solo_read_barcode.cb_match_ind.clear();
            soloreadbarcode_getcbandumi_l93_soloreadbarcode_addstats(
                solo_read_barcode,
                solo_read_barcode.cb_match,
            );
            return Ok(());
        }
    } else if p.p_solo.solo_type == SOLO_TYPE_CB_SAM_TAG_OUT {
        let cb_s = p.p_solo.cb_s as usize - 1;
        let cb_e = cb_s + p.p_solo.cb_l as usize;
        let umi_s = p.p_solo.umi_s as usize - 1;
        let umi_e = umi_s + p.p_solo.umi_l as usize;
        solo_read_barcode.cb_seq = solo_read_barcode.b_seq[cb_s..cb_e].to_string();
        solo_read_barcode.umi_seq = solo_read_barcode.b_seq[umi_s..umi_e].to_string();
        solo_read_barcode.cb_qual = solo_read_barcode.b_qual[cb_s..cb_e].to_string();
        solo_read_barcode.umi_qual = solo_read_barcode.b_qual[umi_s..umi_e].to_string();

        let (cb_match, cb_match_ind, cb_match_string) =
            soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(
                &p.p_solo,
                &solo_read_barcode.cb_seq,
                &solo_read_barcode.cb_qual,
                &p.p_solo.cb_wl,
            );
        solo_read_barcode.cb_match = cb_match;
        solo_read_barcode.cb_match_ind = cb_match_ind;
        solo_read_barcode.cb_match_string = cb_match_string;

        if solo_read_barcode.cb_match == 0 || solo_read_barcode.cb_match == 1 {
            if p.p_solo.cb_wl_yes {
                solo_read_barcode.cb_seq_corrected =
                    p.p_solo.cb_wl_str[solo_read_barcode.cb_match_ind[0] as usize].clone();
            } else {
                solo_read_barcode.cb_seq_corrected = solo_read_barcode.cb_seq.clone();
            }
        } else {
            solo_read_barcode.cb_seq_corrected = "-".to_string();
        }
    } else if p.p_solo.solo_type == SOLO_TYPE_CB_UMI_COMPLEX {
        solo_read_barcode.cb_seq.clear();
        solo_read_barcode.cb_qual.clear();
        solo_read_barcode.umi_seq.clear();
        solo_read_barcode.umi_qual.clear();

        let mut adapter_start = 0u32;
        if p.p_solo.adapter_yes
            && sequencefuns_l341_localalignhammingdist(
                &solo_read_barcode.b_seq,
                &p.p_solo.adapter_seq,
                &mut adapter_start,
            ) > p.p_solo.adapter_mismatches_nmax
        {
            solo_read_barcode.stats.v[0] += 1;
            solo_read_barcode.cb_match = -21;
            return Ok(());
        }

        let Some((umi_seq, umi_qual)) = solobarcode_extractbarcode_l4_solobarcode_extractbarcode(
            &p.p_solo.umi_v,
            &solo_read_barcode.b_seq,
            &solo_read_barcode.b_qual,
            adapter_start,
        ) else {
            solo_read_barcode.stats.v[1] += 1;
            solo_read_barcode.cb_match = -22;
            return Ok(());
        };
        solo_read_barcode.umi_seq = umi_seq;
        solo_read_barcode.umi_qual = umi_qual;

        if p.p_solo.umi_l == 0 {
            p.p_solo.umi_l = solo_read_barcode.umi_seq.len() as u32;
        }

        let mut cb_match_good = true;
        if !soloreadbarcode_getcbandumi_l133_soloreadbarcode_convertcheckumi(solo_read_barcode) {
            cb_match_good = false;
            solo_read_barcode.cb_match = solo_read_barcode.umi_check;
        }

        solo_read_barcode.cb_match_ind = vec![0];
        for cb in &p.p_solo.cb_v {
            let (cb_seq1, cb_qual1) = if let Some((cb_seq1, cb_qual1)) =
                solobarcode_extractbarcode_l4_solobarcode_extractbarcode(
                    cb,
                    &solo_read_barcode.b_seq,
                    &solo_read_barcode.b_qual,
                    adapter_start,
                ) {
                if cb_seq1.len() < cb.min_len as usize
                    || cb_seq1.len() >= cb.wl.len()
                    || cb.wl[cb_seq1.len()].is_empty()
                {
                    if cb_match_good {
                        solo_read_barcode.cb_match = -11;
                        cb_match_good = false;
                    }
                }
                (cb_seq1, cb_qual1)
            } else {
                if cb_match_good {
                    solo_read_barcode.cb_match = -11;
                    cb_match_good = false;
                }
                (String::new(), String::new())
            };
            solo_read_barcode.cb_seq.push_str(&cb_seq1);
            solo_read_barcode.cb_seq.push('_');
            solo_read_barcode.cb_qual.push_str(&cb_qual1);
            solo_read_barcode.cb_qual.push('_');

            if !cb_match_good {
                continue;
            }

            let cb_len1 = cb_seq1.len();
            if p.p_solo.cb_match_wl.edit_dist_2 {
                solo_read_barcode.cb_match = 0;
                let mut cb_b1 = 0u64;
                let pos_n = sequencefuns_l249_convertnuclstrtoint64(&cb_seq1, &mut cb_b1);
                if pos_n != -1 {
                    solo_read_barcode.cb_match = -2;
                    cb_match_good = false;
                } else {
                    let cb_i = servicefuns_l294_binarysearchexact(
                        cb_b1,
                        &cb.wl[cb_len1],
                        cb.wl[cb_len1].len() as u64,
                    );
                    if cb_i >= 0 {
                        solo_read_barcode.cb_match_ind[0] +=
                            cb.wl_factor * (cb_i as u64 + cb.wl_add[cb_len1] as u64);
                    } else {
                        let cb_i = servicefuns_l294_binarysearchexact(
                            cb_b1,
                            &cb.wl_ed[cb_len1],
                            cb.wl_ed[cb_len1].len() as u64,
                        );
                        if cb_i >= 0 {
                            solo_read_barcode.cb_match = 1;
                            let cb_i = cb.wl_ed_ind[cb_len1][cb_i as usize] as u64;
                            solo_read_barcode.cb_match_ind[0] +=
                                cb.wl_factor * (cb_i + cb.wl_add[cb_len1] as u64);
                        } else {
                            solo_read_barcode.cb_match = -1;
                            cb_match_good = false;
                        }
                    }
                }
            } else {
                let (cb_match1, cb_match_ind1, cb_match_string) =
                    soloreadbarcode_getcbandumi_l9_soloreadbarcode_matchcbtowl(
                        &p.p_solo,
                        &cb_seq1,
                        &cb_qual1,
                        &cb.wl[cb_len1],
                    );
                solo_read_barcode.cb_match_string = cb_match_string;
                if cb_match1 < 0 {
                    cb_match_good = false;
                    solo_read_barcode.cb_match = cb_match1;
                } else if cb_match1 > 0 && solo_read_barcode.cb_match > 0 {
                    cb_match_good = false;
                    solo_read_barcode.cb_match = -12;
                } else {
                    solo_read_barcode.cb_match_ind[0] +=
                        cb.wl_factor * (cb_match_ind1[0] + cb.wl_add[cb_len1] as u64);
                    solo_read_barcode.cb_match = solo_read_barcode.cb_match.max(cb_match1);
                }
            }
        }
        solo_read_barcode.cb_seq.pop();
        solo_read_barcode.cb_qual.pop();

        if cb_match_good {
            solo_read_barcode.cb_match_string = solo_read_barcode.cb_match_ind[0].to_string();
        }
    }

    soloreadbarcode_getcbandumi_l93_soloreadbarcode_addstats(
        solo_read_barcode,
        solo_read_barcode.cb_match,
    );
    Ok(())
}
