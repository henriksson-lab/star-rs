#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloRead::record` at STAR/source/SoloRead_record.cpp:3. Args: nTr: uint64, alignOut: Transcript, iRead: uint64, readAnnot: ReadAnnotations"]
pub fn soloread_record_l3_soloread_record(
    solo_read: &mut crate::solo_read::SoloRead,
    p: &crate::parameters_chimeric::Parameters,
    n_tr: u64,
    align_out: &[crate::transcript::Transcript],
    i_read: u64,
    read_annot: &crate::read_annotations::ReadAnnotations,
) {
    if p.p_solo.solo_type == 0 {
        return;
    }
    if p.p_solo.solo_type == 3 {
        return;
    }

    if p.p_solo.read_stats_yes.iter().any(|yes| *yes) {
        soloread_l18_soloread_readflagreset(solo_read);
    }

    if let Some(mut read_bar) = solo_read.read_bar.take() {
        for ii in 0..p.p_solo.n_features as usize {
            if let Some(read_feat) = solo_read.read_feat.get_mut(ii) {
                soloreadfeature_record_l20_soloreadfeature_record(
                    read_feat,
                    p,
                    &mut read_bar,
                    n_tr as u32,
                    align_out,
                    i_read,
                    read_annot,
                );
            }
        }
        solo_read.read_bar = Some(read_bar);
    }
}
