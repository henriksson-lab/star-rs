#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SoloRead` at STAR/source/SoloRead.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloRead {
    pub i_chunk: i32,
    pub read_bar: Option<SoloReadBarcode>,
    pub read_feat: Vec<SoloReadFeature>,
}

#[doc = "Original `SoloRead::SoloRead` at STAR/source/SoloRead.cpp:3. Args: Pin: Parameters, iChunkIn: int32"]
pub fn soloread_l3_soloread_soloread(
    p: &crate::parameters_chimeric::Parameters,
    i_chunk_in: i32,
) -> crate::solo_read::SoloRead {
    let p_solo = &p.p_solo;
    let read_bar = Some(soloreadbarcode_l4_soloreadbarcode_soloreadbarcode(
        p_solo.solo_type,
        p_solo.cb_wl_yes,
        p_solo.cb_wl_size,
        p_solo.umi_l,
    ));

    let mut solo_read = crate::solo_read::SoloRead {
        i_chunk: i_chunk_in,
        read_bar,
        read_feat: Vec::new(),
    };

    if p_solo.solo_type == 0 {
        return solo_read;
    }
    if p_solo.solo_type == 3 {
        return solo_read;
    }

    let n_features = if p_solo.n_features == 0 {
        p_solo.features.len() as u32
    } else {
        p_solo.n_features
    };
    for ii in 0..n_features {
        let feature_type = p_solo.features.get(ii as usize).copied().unwrap_or(ii) as i32;
        solo_read
            .read_feat
            .push(soloreadfeature_l5_soloreadfeature_soloreadfeature(
                feature_type,
                p,
                i_chunk_in,
            ));
    }
    solo_read
}

#[doc = "Original `SoloRead::readFlagReset` at STAR/source/SoloRead.cpp:18. Args: "]
pub fn soloread_l18_soloread_readflagreset(solo_read: &mut crate::solo_read::SoloRead) {
    for ii in 0..solo_read.read_feat.len() {
        solo_read.read_feat[ii].read_flag.flag = 0;
    }
}
