#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `ParametersChimeric::initialize` at STAR/source/ParametersChimeric_initialize.cpp:6. Args: pPin: Parameters"]
pub fn parameterschimeric_initialize_l6_parameterschimeric_initialize(
    pc: &mut crate::parameters_chimeric::ParametersChimeric,
    p: &mut crate::parameters_chimeric::Parameters,
    sam_header: &str,
) -> Result<(), String> {
    pc.out_bam = false;
    pc.out_junctions = false;
    pc.out_sam_old = false;
    pc.out_bam_hard_clip = true;

    if pc.segment_min == 0 {
        return Ok(());
    }

    for type1 in pc.out_type.iter() {
        match type1.as_str() {
            "WithinBAM" => pc.out_bam = true,
            "SeparateSAMold" => pc.out_sam_old = true,
            "Junctions" => pc.out_junctions = true,
            "HardClip" => pc.out_bam_hard_clip = true,
            "SoftClip" => pc.out_bam_hard_clip = false,
            _ => {
                return Err(format!(
                    "EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --chimOutType: {}\nSOLUTION: re-run STAR with --chimOutType Junctions , SeparateSAMold  , WithinBAM , HardClip \n",
                    type1
                ));
            }
        }
    }

    if pc.out_sam_old {
        pc.out_chim_sam_opened = true;
        pc.out_chim_sam_contents = sam_header.to_string();
    }

    if pc.out_junctions {
        pc.out_chim_junction_opened = true;
        if pc.multimap_nmax > 0 {
            pc.out_chim_junction_contents =
                "chr_donorA\tbrkpt_donorA\tstrand_donorA\tchr_acceptorB\tbrkpt_acceptorB\tstrand_acceptorB\tjunction_type\trepeat_left_lenA\trepeat_right_lenB\tread_name\tstart_alnA\tcigar_alnA\tstart_alnB\tcigar_alnB\tnum_chim_aln\tmax_poss_aln_score\tnon_chim_aln_score\tthis_chim_aln_score\tbestall_chim_aln_score\tPEmerged_bool\treadgrp\n"
                    .to_string();
        }
    }

    if pc.out_bam && !p.out_bam_unsorted && !p.out_bam_coord {
        return Err("EXITING because of fatal PARAMETERS error: --chimOutType WithinBAM requires BAM output\nSOLUTION: re-run with --outSAMtype BAM Unsorted/SortedByCoordinate\n".to_string());
    }

    if pc.multimap_nmax > 0 && pc.out_sam_old {
        return Err("EXITING because of fatal PARAMETERS error: --chimMultimapNmax > 0 (new chimeric detection) presently only works with --chimOutType Junctions/WithinBAM\nSOLUTION: re-run with --chimOutType Junctions/WithinBAM\n".to_string());
    }

    if p.pe_overlap_nbases_min > 0 && pc.multimap_nmax == 0 && (pc.out_junctions || pc.out_sam_old)
    {
        return Err("EXITING because of fatal PARAMETERS error: --chimMultimapNmax 0 (default old chimeric detection) and --peOverlapNbasesMin > 0 (merging ovelrapping mates) presently only works with --chimOutType WithinBAM\nSOLUTION: re-run with --chimOutType WithinBAM\n".to_string());
    }

    if pc.out_bam && !p.out_sam_attr_nm_present {
        p.out_sam_attr_order.push(ATTR_NM);
        pc.log_main
            .push_str("WARNING --chimOutType=WithinBAM, therefore STAR will output NM attribute\n");
    }

    pc.filter_genomic_n = false;
    for filter in pc.filter_string_in.iter() {
        match filter.as_str() {
            "banGenomicN" => pc.filter_genomic_n = true,
            "None" => {}
            _ => {
                return Err(format!(
                    "EXITING because of fatal PARAMETERS error: unrecognized value of --chimFilter={}\nSOLUTION: use allowed values: banGenomicN || None",
                    filter
                ));
            }
        }
    }
    Ok(())
}
