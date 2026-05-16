#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `sjdbLoadFromStream` at STAR/source/sjdbLoadFromStream.cpp:2. Args: sjdbStreamIn: ifstream, sjdbLoci: SjdbClass"]
pub fn sjdbloadfromstream_l2_sjdbloadfromstream(
    sjdb_stream_in: &str,
    sjdb_loci: &mut crate::sjdb_class::SjdbClass,
) {
    for one_line in sjdb_stream_in.lines() {
        let mut fields = one_line.split_whitespace();
        let chr1 = fields.next().unwrap_or("");
        let u1 = fields
            .next()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or(0);
        let u2 = fields
            .next()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or(0);
        let mut str1 = fields.next().and_then(|v| v.chars().next()).unwrap_or('\0');

        if !chr1.is_empty() {
            sjdb_loci.chr.push(chr1.to_string());
            sjdb_loci.start.push(u1);
            sjdb_loci.end.push(u2);
            str1 = match str1 {
                '1' | '+' => '+',
                '2' | '-' => '-',
                _ => '.',
            };
            sjdb_loci.str_.push(str1);
        }
    }
}
