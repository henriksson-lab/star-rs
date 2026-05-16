#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `signalFromBAM` at STAR/source/signalFromBAM.cpp:5. Args: bamFileName: string, sigFileName: string, P: Parameters"]
pub fn signalfrombam_l5_signalfrombam(
    sig_file_name: &str,
    p: &crate::parameters_chimeric::Parameters,
    target_names: &[String],
    target_lens: &[u32],
    records: &[crate::parameters_chimeric::SignalFromBamRecord],
) -> Result<crate::parameters_chimeric::SignalFromBamResult, String> {
    let mut n_mult = 0.0_f64;
    let mut n_uniq = 0.0_f64;
    let prefix = if p.out_wig_references_prefix.is_empty() {
        "-"
    } else {
        p.out_wig_references_prefix.as_str()
    };

    if p.out_wig_flags.norm == 1 {
        for rec in records {
            if rec.tid < 0 {
                continue;
            }
            let tid = rec.tid as usize;
            if prefix != "-" && !target_names[tid].starts_with(prefix) {
                continue;
            }
            if let Some(nh) = rec.nh {
                if nh == 1 {
                    n_uniq += 1.0;
                } else if nh > 1 {
                    n_mult += 1.0 / nh as f64;
                }
            }
        }
    }

    let sig_n = if p.out_wig_flags.strand { 4 } else { 2 };
    let mut sig_out_file_names = Vec::with_capacity(sig_n as usize);
    sig_out_file_names.push(format!("{sig_file_name}.Unique.str1.out"));
    sig_out_file_names.push(format!("{sig_file_name}.UniqueMultiple.str1.out"));
    if p.out_wig_flags.strand {
        sig_out_file_names.push(format!("{sig_file_name}.Unique.str2.out"));
        sig_out_file_names.push(format!("{sig_file_name}.UniqueMultiple.str2.out"));
    }
    for name in &mut sig_out_file_names {
        name.push_str(if p.out_wig_flags.format == 0 {
            ".bg"
        } else {
            ".wig"
        });
    }

    let mut norm_factor = vec![1.0_f64; sig_n as usize];
    if p.out_wig_flags.norm == 1 {
        norm_factor[0] = 1.0e6 / n_uniq;
        norm_factor[1] = 1.0e6 / (n_uniq + n_mult);
    }
    if p.out_wig_flags.strand {
        norm_factor[2] = norm_factor[0];
        norm_factor[3] = norm_factor[1];
    }

    let mut files = std::collections::BTreeMap::<String, String>::new();
    for name in &sig_out_file_names {
        files.insert(name.clone(), String::new());
    }

    let mut i_chr = -999_i32;
    let mut sig_all = Vec::<f64>::new();
    let mut chr_len = 0_u32;

    for rec_opt in records.iter().map(Some).chain(std::iter::once(None)) {
        let bam_bytes_negative = rec_opt.is_none();
        let rec_tid = rec_opt.map(|rec| rec.tid).unwrap_or(i_chr);
        if rec_tid != i_chr || bam_bytes_negative {
            if i_chr != -999 {
                let chr_name = &target_names[i_chr as usize];
                for is in 0..sig_n as usize {
                    let out = files.get_mut(&sig_out_file_names[is]).unwrap();
                    if p.out_wig_flags.format == 1 {
                        out.push_str(&format!("variableStep chrom={chr_name}\n"));
                    }
                    let mut prev_sig = 0.0_f64;
                    for ig in 0..chr_len as usize {
                        let new_sig = sig_all[sig_n as usize * ig + is];
                        if p.out_wig_flags.format == 0 {
                            if new_sig != prev_sig {
                                if prev_sig != 0.0 {
                                    if p.out_wig_flags.norm == 1 {
                                        out.push_str(&format!(
                                            "{}\t{:.5}\n",
                                            ig,
                                            prev_sig * norm_factor[is]
                                        ));
                                    } else {
                                        out.push_str(&format!(
                                            "{}\t{}\n",
                                            ig,
                                            prev_sig * norm_factor[is]
                                        ));
                                    }
                                }
                                if new_sig != 0.0 {
                                    out.push_str(&format!("{chr_name}\t{ig}\t"));
                                }
                                prev_sig = new_sig;
                            }
                        } else if p.out_wig_flags.format == 1 && new_sig != 0.0 {
                            if p.out_wig_flags.norm == 1 {
                                out.push_str(&format!(
                                    "{}\t{:.5}\n",
                                    ig + 1,
                                    new_sig * norm_factor[is]
                                ));
                            } else {
                                out.push_str(&format!(
                                    "{}\t{}\n",
                                    ig + 1,
                                    new_sig * norm_factor[is]
                                ));
                            }
                        }
                    }
                }
            }
            if bam_bytes_negative {
                break;
            }

            let rec = rec_opt.unwrap();
            i_chr = rec.tid;
            if i_chr == -1 || (prefix != "-" && !target_names[rec.tid as usize].starts_with(prefix))
            {
                i_chr = -999;
                continue;
            }
            chr_len = target_lens[rec.tid as usize] + 1;
            sig_all = vec![0.0; sig_n as usize * chr_len as usize];
        }

        let rec = rec_opt.unwrap();
        if (rec.flag & 0x400) > 0 {
            continue;
        }
        let a_nh = rec.nh.unwrap_or(1);
        if a_nh == 0 {
            continue;
        }
        let mut a_g = rec.pos;
        let mut i_strand = 0_u32;
        if p.out_wig_flags.strand {
            i_strand = u32::from(((rec.flag & 0x10) > 0) == ((rec.flag & 0x80) == 0));
        }

        if p.out_wig_flags.type_ == 1 {
            if (rec.flag & 0x80) > 0 {
                continue;
            }
            if i_strand == 0 {
                if a_nh == 1 {
                    sig_all[(a_g * sig_n as u32 + 2 * i_strand) as usize] += 1.0;
                }
                sig_all[(a_g * sig_n as u32 + 1 + 2 * i_strand) as usize] += 1.0 / a_nh as f64;
                continue;
            }
        }

        for &cig in &rec.cigar {
            let cig_op = cig & 0xf;
            let cig_l = cig >> 4;
            match cig_op {
                BAM_CIGAR_D | BAM_CIGAR_N => {
                    a_g += cig_l;
                }
                BAM_CIGAR_M => {
                    if p.out_wig_flags.type_ == 0
                        || (p.out_wig_flags.type_ == 2 && (rec.flag & 0x80) > 0)
                    {
                        for _ in 0..cig_l {
                            if a_g >= chr_len {
                                return Err(
                                    "BUG: alignment extends past chromosome in signalFromBAM.cpp\n"
                                        .to_string(),
                                );
                            }
                            if a_nh == 1 {
                                sig_all[(a_g * sig_n as u32 + 2 * i_strand) as usize] += 1.0;
                            }
                            sig_all[(a_g * sig_n as u32 + 1 + 2 * i_strand) as usize] +=
                                1.0 / a_nh as f64;
                            a_g += 1;
                        }
                    } else {
                        a_g += cig_l;
                    }
                }
                _ => {}
            }
        }
        if p.out_wig_flags.type_ == 1 {
            a_g = a_g.wrapping_sub(1);
            if a_nh == 1 {
                sig_all[(a_g * sig_n as u32 + 2 * i_strand) as usize] += 1.0;
            }
            sig_all[(a_g * sig_n as u32 + 1 + 2 * i_strand) as usize] += 1.0 / a_nh as f64;
        }
    }

    Ok(crate::parameters_chimeric::SignalFromBamResult {
        files,
        n_unique: n_uniq,
        n_multiple: n_mult,
    })
}
