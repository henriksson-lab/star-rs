#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `ParameterInfoBase` at STAR/source/ParameterInfo.h:4."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParameterInfoBase {}

#[doc = "Original class `ParameterInfoScalar` at STAR/source/ParameterInfo.h:58."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParameterInfoScalar {}

#[doc = "Original class `ParameterInfoVector` at STAR/source/ParameterInfo.h:83."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParameterInfoVector {}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParameterScanEntry {
    pub name_string: String,
    pub input_level_allowed: i32,
    pub input_level: i32,
    pub value_line: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParametersScanState {
    pub par_array: Vec<ParameterScanEntry>,
    pub parameter_input_name: Vec<String>,
    pub log_main: String,
}
