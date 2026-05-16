#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SNP` at STAR/source/Variation.h:15."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SNP {
    pub n: u32,
    pub loci: Vec<u32>,
    pub loci_v: Vec<u32>,
    pub nt: Vec<[u8; 3]>,
}

#[doc = "Original class `Variation` at STAR/source/Variation.h:30."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Variation {
    pub yes: bool,
    pub snp: SNP,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct VariantInfo {
    pub pos: u64,
    pub len: i32,
    pub seq: [String; 2],
}

#[doc = "Original `Variation::Variation` at STAR/source/Variation.cpp:8. Args: Pin: Parameters, chrStartIn: vector <uint>, chrNameIndexIn: map <string,uint>, yesVar: bool"]
pub fn variation_l8_variation_variation(yes_var: bool) -> crate::variation::Variation {
    crate::variation::Variation {
        yes: yes_var,
        ..Default::default()
    }
}

#[doc = "Original `scanVCF` at STAR/source/Variation.cpp:23. Args: vcf: ifstream, P: Parameters, snp: SNP, chrStart: vector <uint>, chrNameIndex: map <string,uint>"]
pub fn variation_l23_scanvcf(
    vcf: &str,
    snp: &mut crate::variation::SNP,
    chr_start: &[u32],
    chr_name_index: &std::collections::BTreeMap<String, u32>,
    hetero_only: bool,
) -> Result<u32, String> {
    snp.n = 0;
    snp.loci_v.clear();
    snp.nt.clear();
    let mut n_homoz = 0u32;

    for (line_no0, line) in vcf.lines().enumerate() {
        let line_no = line_no0 + 1;
        let mut fields = line.split_whitespace();
        let Some(chr) = fields.next() else {
            continue;
        };
        if chr.starts_with('#') {
            continue;
        }

        let pos = fields
            .next()
            .ok_or_else(|| format!("Malformed VCF line {line_no}: missing POS field"))?
            .parse::<u32>()
            .map_err(|_| format!("Malformed VCF line {line_no}: invalid POS field"))?;
        let _id = fields
            .next()
            .ok_or_else(|| format!("Malformed VCF line {line_no}: missing ID field"))?;
        let ref_ = fields
            .next()
            .ok_or_else(|| format!("Malformed VCF line {line_no}: missing REF field"))?;
        let alt = fields
            .next()
            .ok_or_else(|| format!("Malformed VCF line {line_no}: missing ALT field"))?;
        let _qual = fields
            .next()
            .ok_or_else(|| format!("Malformed VCF line {line_no}: missing QUAL field"))?;
        let _filter = fields
            .next()
            .ok_or_else(|| format!("Malformed VCF line {line_no}: missing FILTER field"))?;
        let _info = fields
            .next()
            .ok_or_else(|| format!("Malformed VCF line {line_no}: missing INFO field"))?;
        let _format = fields
            .next()
            .ok_or_else(|| format!("Malformed VCF line {line_no}: missing FORMAT field"))?;
        let sample = fields.next().unwrap_or("");

        let mut alt_v: Vec<&str> = alt.split(',').collect();
        let max_alt_len = alt_v.iter().map(|value| value.len()).max().unwrap_or(0);
        if ref_.len() == 1 && max_alt_len == 1 {
            alt_v.insert(0, ref_);

            if let Some(&chr_index) = chr_name_index.get(chr) {
                let Some(chr_start_value) = chr_start.get(chr_index as usize) else {
                    return Err(format!(
                        "Malformed VCF line {line_no}: chromosome index is outside chrStart"
                    ));
                };
                if sample.len() < 3 {
                    continue;
                } else if sample.len() > 3 && sample.as_bytes()[3] != b':' {
                    continue;
                } else if sample.as_bytes()[0] == b'0' && sample.as_bytes()[2] == b'0' {
                    continue;
                }
                if sample.as_bytes()[1] != b'/' && sample.as_bytes()[1] != b'|' {
                    continue;
                }

                let Some(allele0) = char::from(sample.as_bytes()[0])
                    .to_digit(10)
                    .map(|value| value as usize)
                else {
                    continue;
                };
                let Some(allele1) = char::from(sample.as_bytes()[2])
                    .to_digit(10)
                    .map(|value| value as usize)
                else {
                    continue;
                };
                if allele0 >= alt_v.len() || allele1 >= alt_v.len() {
                    return Err(format!(
                        "Malformed VCF line {line_no}: genotype allele index is outside REF/ALT alleles"
                    ));
                }
                if alt_v[allele0].as_bytes()[0] == ref_.as_bytes()[0]
                    && alt_v[allele1].as_bytes()[0] == ref_.as_bytes()[0]
                {
                    continue;
                } else if hetero_only && sample.as_bytes()[0] == sample.as_bytes()[2] {
                    n_homoz += 1;
                    continue;
                }

                let nt1 = [
                    sequencefuns_l195_convertnt01234(ref_.as_bytes()[0]),
                    sequencefuns_l195_convertnt01234(alt_v[allele0].as_bytes()[0]),
                    sequencefuns_l195_convertnt01234(alt_v[allele1].as_bytes()[0]),
                ];
                if nt1[0] < 4 && nt1[1] < 4 && nt1[2] < 4 {
                    let Some(pos0) = pos.checked_sub(1) else {
                        return Err(format!(
                            "Malformed VCF line {line_no}: POS field must be 1-based"
                        ));
                    };
                    snp.loci_v.push(pos0 + chr_start_value);
                    snp.nt.push(nt1);
                    snp.n += 1;
                }
            }
        }
    }

    Ok(n_homoz)
}

#[doc = "Original `Variation::loadVCF` at STAR/source/Variation.cpp:81. Args: fileIn: string"]
pub fn variation_l81_variation_loadvcf(
    variation: &mut crate::variation::Variation,
    vcf: &str,
    chr_start: &[u32],
    chr_name_index: &std::collections::BTreeMap<String, u32>,
    hetero_only: bool,
) -> Result<u32, String> {
    let n_homoz = variation_l23_scanvcf(
        vcf,
        &mut variation.snp,
        chr_start,
        chr_name_index,
        hetero_only,
    )?;

    variation.snp.loci = variation.snp.loci_v.clone();
    variation.snp.loci_v.clear();

    if variation.snp.n == 0 {
        return Err(
            "EXITING because of FATAL INPUT FILE ERROR: could not find any SNPs in VCF file"
                .to_string(),
        );
    }

    let nt1 = variation.snp.nt.clone();
    let mut s1: Vec<[u32; 2]> = (0..variation.snp.n as usize)
        .map(|ii| [variation.snp.loci[ii], ii as u32])
        .collect();
    s1.sort_by(|a, b| {
        let cmp = servicefuns_l53_funcompareuint1(&a[0], &b[0]);
        if cmp < 0 {
            std::cmp::Ordering::Less
        } else if cmp > 0 {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    });

    for (ii, value) in s1.iter().enumerate() {
        variation.snp.loci[ii] = value[0];
        variation.snp.nt[ii] = nt1[value[1] as usize];
    }

    Ok(n_homoz)
}

#[doc = "Original `SNP::snpOnBlocks` at STAR/source/Variation.cpp:124. Args: blockStart: uint, blockL: uint, blockShift: int, snpV: vector<vector<array<int,2>>>"]
pub fn variation_l124_snp_snponblocks(
    snp: &crate::variation::SNP,
    block_start: u32,
    block_l: u32,
    block_shift: i32,
    snp_v: &mut Vec<Vec<[i32; 2]>>,
) {
    let mut isnp = servicefuns_l266_binarysearch1b(block_start, &snp.loci, snp.n as i32);
    while (isnp as u32) < snp.n && snp.loci[isnp as usize] < block_start + block_l {
        for ii in 0..2 {
            if snp.nt[isnp as usize][ii + 1] != snp.nt[isnp as usize][0] {
                let snp1 = [
                    (snp.loci[isnp as usize] - block_start) as i32 + block_shift,
                    snp.nt[isnp as usize][ii + 1] as i32,
                ];
                snp_v[ii].push(snp1);
            }
        }
        isnp += 1;
    }
}

#[doc = "Original `Variation::sjdbSnp` at STAR/source/Variation.cpp:139. Args: sjStart: uint, sjEnd: uint, sjdbOverhang1: uint"]
pub fn variation_l139_variation_sjdbsnp(
    variation: &crate::variation::Variation,
    sj_start: u32,
    sj_end: u32,
    sjdb_overhang1: u32,
) -> Vec<Vec<[i32; 2]>> {
    let mut snp_v = vec![Vec::new(), Vec::new()];

    if !variation.yes {
        return vec![Vec::new()];
    }

    variation_l124_snp_snponblocks(
        &variation.snp,
        sj_start - sjdb_overhang1,
        sjdb_overhang1,
        0,
        &mut snp_v,
    );
    variation_l124_snp_snponblocks(
        &variation.snp,
        sj_end + 1,
        sjdb_overhang1,
        sjdb_overhang1 as i32,
        &mut snp_v,
    );

    if (snp_v[0].is_empty() && snp_v[1].is_empty()) || snp_v[0] == snp_v[1] {
        snp_v.pop();
    }

    snp_v
}
