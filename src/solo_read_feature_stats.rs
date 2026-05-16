#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SoloReadFeatureStats` at STAR/source/SoloReadFeatureStats.h:5."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloReadFeatureStats {
    pub names: Vec<String>,
    pub v: Vec<u64>,
}
