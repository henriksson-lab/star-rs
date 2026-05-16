#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `Parameters::samAttributes` at STAR/source/Parameters_samAttributes.cpp:4. Args: "]
pub fn parameters_samattributes_l4_parameters_samattributes(
    parameters: &mut crate::parameters_chimeric::Parameters,
) -> Result<String, String> {
    let mut log_main = String::new();

    parameters.out_sam_attr_present = Default::default();
    parameters.out_sam_attr_nm_present = false;
    parameters.out_sam_attr_present_quant = parameters.out_sam_attr_present.clone();
    parameters.out_sam_attr_present_quant.nh = true;
    parameters.out_sam_attr_present_quant.hi = true;
    parameters.out_sam_attr_order.clear();
    parameters.out_sam_attr_order_quant.clear();
    parameters.out_sam_attr_order_quant.push(ATTR_NH);
    parameters.out_sam_attr_order_quant.push(ATTR_HI);

    let v_attr1 = if parameters.out_sam_attributes.first().map(|v| v.as_str()) == Some("None") {
        Vec::new()
    } else if parameters.out_sam_attributes.first().map(|v| v.as_str()) == Some("All") {
        ["NH", "HI", "AS", "nM", "NM", "MD", "jM", "jI", "MC", "ch"]
            .into_iter()
            .map(String::from)
            .collect()
    } else if parameters.out_sam_attributes.first().map(|v| v.as_str()) == Some("Standard") {
        ["NH", "HI", "AS", "nM"]
            .into_iter()
            .map(String::from)
            .collect()
    } else {
        parameters.out_sam_attributes.clone()
    };

    for attr in &v_attr1 {
        match attr.as_str() {
            "NH" => {
                parameters.out_sam_attr_order.push(ATTR_NH);
                parameters.out_sam_attr_present.nh = true;
            }
            "HI" => {
                parameters.out_sam_attr_order.push(ATTR_HI);
                parameters.out_sam_attr_present.hi = true;
            }
            "AS" => {
                parameters.out_sam_attr_order.push(ATTR_AS);
                parameters.out_sam_attr_present.as_ = true;
            }
            "NM" => {
                parameters.out_sam_attr_order.push(ATTR_NM);
                parameters.out_sam_attr_present.nm = true;
                parameters.out_sam_attr_nm_present = true;
            }
            "MD" => {
                parameters.out_sam_attr_order.push(ATTR_MD);
                parameters.out_sam_attr_present.md = true;
            }
            "nM" => {
                parameters.out_sam_attr_order.push(ATTR_NM_LOWER);
                parameters.out_sam_attr_present.n_m = true;
            }
            "jM" => {
                parameters.out_sam_attr_order.push(ATTR_JM);
                parameters.out_sam_attr_present.j_m = true;
            }
            "jI" => {
                parameters.out_sam_attr_order.push(ATTR_JI);
                parameters.out_sam_attr_present.j_i = true;
            }
            "vA" => {
                parameters.out_sam_attr_order.push(ATTR_VA);
                parameters.out_sam_attr_present.v_a = true;
            }
            "vG" => {
                parameters.out_sam_attr_order.push(ATTR_VG);
                parameters.out_sam_attr_present.v_g = true;
            }
            "vW" => {
                parameters.out_sam_attr_order.push(ATTR_VW);
                parameters.out_sam_attr_present.v_w = true;
            }
            "ha" => {
                parameters.out_sam_attr_order.push(ATTR_HA);
                parameters.out_sam_attr_present.ha = true;
            }
            "RG" => {
                parameters.out_sam_attr_order.push(ATTR_RG);
                parameters.out_sam_attr_order_quant.push(ATTR_RG);
                parameters.out_sam_attr_present.rg = true;
            }
            "rB" => {
                parameters.out_sam_attr_order.push(ATTR_RB);
                parameters.out_sam_attr_order_quant.push(ATTR_RB);
                parameters.out_sam_attr_present.r_b = true;
            }
            "ch" => {
                parameters.out_sam_attr_order.push(ATTR_CH);
                parameters.out_sam_attr_order_quant.push(ATTR_CH);
                parameters.out_sam_attr_present.ch = true;
            }
            "MC" => {
                parameters.out_sam_attr_order.push(ATTR_MC);
                parameters.out_sam_attr_order_quant.push(ATTR_MC);
                parameters.out_sam_attr_present.mc = true;
            }
            "CR" => {
                parameters.out_sam_attr_order.push(ATTR_CR);
                parameters.out_sam_attr_order_quant.push(ATTR_CR);
                parameters.out_sam_attr_present.cr = true;
            }
            "CY" => {
                parameters.out_sam_attr_order.push(ATTR_CY);
                parameters.out_sam_attr_order_quant.push(ATTR_CY);
                parameters.out_sam_attr_present.cy = true;
            }
            "UR" => {
                parameters.out_sam_attr_order.push(ATTR_UR);
                parameters.out_sam_attr_order_quant.push(ATTR_UR);
                parameters.out_sam_attr_present.ur = true;
            }
            "UY" => {
                parameters.out_sam_attr_order.push(ATTR_UY);
                parameters.out_sam_attr_order_quant.push(ATTR_UY);
                parameters.out_sam_attr_present.uy = true;
            }
            "CB" => {
                parameters.out_sam_attr_order.push(ATTR_CB);
                parameters.out_sam_attr_order_quant.push(ATTR_CB);
                parameters.out_sam_attr_present.cb = true;
            }
            "UB" => {
                parameters.out_sam_attr_order.push(ATTR_UB);
                parameters.out_sam_attr_order_quant.push(ATTR_UB);
                parameters.out_sam_attr_present.ub = true;
            }
            "GX" => {
                parameters.out_sam_attr_order.push(ATTR_GX);
                parameters.out_sam_attr_order_quant.push(ATTR_GX);
                parameters.out_sam_attr_present.gx = true;
            }
            "GN" => {
                parameters.out_sam_attr_order.push(ATTR_GN);
                parameters.out_sam_attr_order_quant.push(ATTR_GN);
                parameters.out_sam_attr_present.gn = true;
            }
            "gx" => {
                parameters.out_sam_attr_order.push(ATTR_GX_LOWER);
                parameters.out_sam_attr_order_quant.push(ATTR_GX_LOWER);
                parameters.out_sam_attr_present.gx_lower = true;
            }
            "gn" => {
                parameters.out_sam_attr_order.push(ATTR_GN_LOWER);
                parameters.out_sam_attr_order_quant.push(ATTR_GN_LOWER);
                parameters.out_sam_attr_present.gn_lower = true;
            }
            "sM" => {
                parameters.out_sam_attr_order.push(ATTR_SM);
                parameters.out_sam_attr_order_quant.push(ATTR_SM);
                parameters.out_sam_attr_present.s_m = true;
            }
            "sS" => {
                parameters.out_sam_attr_order.push(ATTR_SS);
                parameters.out_sam_attr_order_quant.push(ATTR_SS);
                parameters.out_sam_attr_present.s_s = true;
            }
            "sF" => {
                parameters.out_sam_attr_order.push(ATTR_SF);
                parameters.out_sam_attr_order_quant.push(ATTR_SF);
                parameters.out_sam_attr_present.s_f = true;
            }
            "sQ" => {
                parameters.out_sam_attr_order.push(ATTR_SQ);
                parameters.out_sam_attr_order_quant.push(ATTR_SQ);
                parameters.out_sam_attr_present.s_q = true;
            }
            "cN" => {
                parameters.out_sam_attr_order.push(ATTR_CN);
                parameters.out_sam_attr_order_quant.push(ATTR_CN);
                parameters.out_sam_attr_present.c_n = true;
            }
            "XS" => {
                parameters.out_sam_attr_order.push(ATTR_XS);
                parameters.out_sam_attr_present.xs = true;
                if parameters.out_sam_strand_field_type != 1 {
                    log_main.push_str("WARNING --outSAMattributes contains XS, therefore STAR will use --outSAMstrandField intronMotif\n");
                    parameters.out_sam_strand_field_type = 1;
                }
            }
            _ => {
                return Err(format!(
                    "EXITING because of FATAL INPUT ERROR: unknown/unimplemented SAM atrribute (tag): {}\nSOLUTION: re-run STAR with --outSAMattributes that contains only implemented attributes\n",
                    attr
                ));
            }
        }
    }

    if !parameters.var_yes
        && (parameters.out_sam_attr_present.v_a || parameters.out_sam_attr_present.v_g)
    {
        return Err("EXITING because of fatal PARAMETER error: --outSAMattributes contains vA and/or vG tag(s), but --varVCFfile is not set\nSOLUTION: re-run STAR with a --varVCFfile option, or without vA/vG tags in --outSAMattributes\n".to_string());
    }
    if !parameters.wasp_yes && parameters.out_sam_attr_present.v_w {
        return Err("EXITING because of fatal PARAMETER error: --outSAMattributes contains vW tag, but --waspOutputMode is not set\nSOLUTION: re-run STAR with a --waspOutputMode option, or without vW tags in --outSAMattributes\n".to_string());
    }

    if parameters
        .out_sam_attr_rgline
        .first()
        .map(|v| v.as_str())
        .is_some_and(|v| v != "-")
        && !parameters.out_sam_attr_present.rg
    {
        parameters.out_sam_attr_order.push(ATTR_RG);
        parameters.out_sam_attr_order_quant.push(ATTR_RG);
        parameters.out_sam_attr_present.rg = true;
        log_main.push_str("WARNING --outSAMattrRG defines a read group, therefore STAR will output RG attribute\n");
    } else if parameters.out_sam_attr_rg.is_empty() && parameters.out_sam_attr_present.rg {
        return Err("EXITING because of FATAL INPUT ERROR: --outSAMattributes contains RG tag, but --outSAMattrRGline is not set\nSOLUTION: re-run STAR with a valid read group parameter --outSAMattrRGline\n".to_string());
    }

    if parameters.out_sam_strand_field_type == 1 && !parameters.out_sam_attr_present.xs {
        parameters.out_sam_attr_order.push(ATTR_XS);
        log_main.push_str(
            "WARNING --outSAMstrandField=intronMotif, therefore STAR will output XS attribute\n",
        );
    }

    if parameters.wasp_yes && !parameters.out_sam_attr_present.v_w {
        parameters.out_sam_attr_order.push(ATTR_VW);
        parameters.out_sam_attr_order_quant.push(ATTR_VW);
        parameters.out_sam_attr_present.v_w = true;
        log_main
            .push_str("WARNING --waspOutputMode is set, therefore STAR will output vW attribute\n");
    }

    for (yes, tag) in [
        (parameters.out_sam_attr_present.ch, "ch"),
        (parameters.out_sam_attr_present.cr, "CR"),
        (parameters.out_sam_attr_present.cy, "CY"),
        (parameters.out_sam_attr_present.ur, "UR"),
        (parameters.out_sam_attr_present.uy, "UY"),
        (parameters.out_sam_attr_present.cb, "CB"),
        (parameters.out_sam_attr_present.ub, "UB"),
        (parameters.out_sam_attr_present.s_m, "sM"),
        (parameters.out_sam_attr_present.s_s, "sS"),
        (parameters.out_sam_attr_present.s_s, "sF"),
        (parameters.out_sam_attr_present.s_q, "sQ"),
        (parameters.out_sam_attr_present.r_b, "rB"),
        (parameters.out_sam_attr_present.v_g, "vG"),
        (parameters.out_sam_attr_present.v_a, "vA"),
        (parameters.out_sam_attr_present.v_w, "vW"),
        (parameters.out_sam_attr_present.gx, "GX"),
        (parameters.out_sam_attr_present.gn, "GN"),
    ] {
        if let Some(warning) =
            parameters_samattributes_l251_parameters_samattrrequiresbam(parameters, yes, tag)?
        {
            log_main.push_str(&warning);
        }
    }

    if parameters.out_sam_attr_present.gx || parameters.out_sam_attr_present.gn {
        parameters.quant_gene_yes = true;
        parameters.quant_yes = true;
        parameters.p_solo.sam_attr_feature = SOLO_FEATURE_GENE;
    }

    Ok(log_main)
}

#[doc = "Original `Parameters:: samAttrRequiresBAM` at STAR/source/Parameters_samAttributes.cpp:251. Args: attrYes: bool, attrTag: string"]
pub fn parameters_samattributes_l251_parameters_samattrrequiresbam(
    parameters: &crate::parameters_chimeric::Parameters,
    attr_yes: bool,
    attr_tag: &str,
) -> Result<Option<String>, String> {
    if !attr_yes {
        return Ok(None);
    }
    if !parameters.out_bam_unsorted && !parameters.out_bam_coord {
        return Err(format!(
            "EXITING because of fatal PARAMETER error: --outSAMattributes contains {} tag, which requires BAM output.\nSOLUTION: re-run STAR with --outSAMtype BAM Unsorted (and/or) SortedByCoordinate option, or without {} tag in --outSAMattributes\n",
            attr_tag, attr_tag
        ));
    }

    if parameters.out_sam_bool {
        Ok(Some(format!(
            "WARNING: --outSAMattributes contains {} tag. It will be output into BAM file(s), but not SAM file.\n",
            attr_tag
        )))
    } else {
        Ok(None)
    }
}
