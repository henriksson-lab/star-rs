#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SjdbClass` at STAR/source/SjdbClass.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SjdbClass {
    pub chr: Vec<String>,
    pub start: Vec<u64>,
    pub end: Vec<u64>,
    pub str_: Vec<char>,
    pub priority: Vec<u8>,
    pub gene: Vec<std::collections::BTreeSet<u64>>,
}
