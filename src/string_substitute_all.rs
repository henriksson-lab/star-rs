#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `stringSubstituteAll` at STAR/source/stringSubstituteAll.cpp:3. Args: str: std::string, from: std::string, to: std::string"]
pub fn stringsubstituteall_l3_stringsubstituteall(str_: &mut String, from: &str, to: &str) {
    if from.is_empty() {
        return;
    }
    let mut start_pos = 0;
    while let Some(found) = str_[start_pos..].find(from) {
        start_pos += found;
        str_.replace_range(start_pos..start_pos + from.len(), to);
        start_pos += to.len();
    }
}
