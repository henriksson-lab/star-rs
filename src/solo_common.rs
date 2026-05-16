#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SoloReadFlagClass` at STAR/source/SoloCommon.h:26."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloReadFlagClass {
    pub flag: u32,
    pub flag_counts:
        std::collections::BTreeMap<u64, [u64; crate::include_define::SOLO_READ_FLAG_N_BITS]>,
    pub flag_counts_no_cb: [u64; crate::include_define::SOLO_READ_FLAG_N_BITS],
}
