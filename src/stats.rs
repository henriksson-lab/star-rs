#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `Stats` at STAR/source/Stats.h:9."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Stats {
    pub read_n: u32,
    pub read_bases: u32,
    pub mapped_reads_u: u32,
    pub mapped_reads_m: u32,
    pub mapped_bases: u32,
    pub mapped_mismatches_n: u32,
    pub mapped_ins_n: u32,
    pub mapped_del_n: u32,
    pub mapped_ins_l: u32,
    pub mapped_del_l: u32,
    pub mapped_portion: f64,
    pub splices_n: [u32; crate::include_define::SJ_MOTIF_SIZE],
    pub splices_nsjdb: u32,
    pub unmapped_other: u32,
    pub unmapped_short: u32,
    pub unmapped_mismatch: u32,
    pub unmapped_multi: u32,
    pub unmapped_all: u32,
    pub chimeric_all: u32,
    pub time_start: libc::time_t,
    pub time_start_map: libc::time_t,
    pub time_finish_map: libc::time_t,
    pub time_last_report: libc::time_t,
    pub time_finish: libc::time_t,
}

#[doc = "Original `Stats::resetN` at STAR/source/Stats.cpp:4. Args: "]
pub fn stats_l4_stats_resetn(stats: &mut crate::stats::Stats) {
    stats.read_n = 0;
    stats.read_bases = 0;
    stats.mapped_mismatches_n = 0;
    stats.mapped_ins_n = 0;
    stats.mapped_del_n = 0;
    stats.mapped_ins_l = 0;
    stats.mapped_del_l = 0;
    stats.mapped_bases = 0;
    stats.mapped_portion = 0.0;
    stats.mapped_reads_u = 0;
    stats.mapped_reads_m = 0;
    stats.unmapped_other = 0;
    stats.unmapped_short = 0;
    stats.unmapped_mismatch = 0;
    stats.unmapped_multi = 0;
    stats.unmapped_all = 0;
    stats.chimeric_all = 0;
    stats.splices_nsjdb = 0;
    for ii in 0..SJ_MOTIF_SIZE {
        stats.splices_n[ii] = 0;
    }
}

#[doc = "Original `Stats::Stats` at STAR/source/Stats.cpp:16. Args: "]
pub fn stats_l16_stats_stats() -> crate::stats::Stats {
    let mut stats = crate::stats::Stats::default();
    stats_l4_stats_resetn(&mut stats);
    stats.time_last_report = 0;
    stats
}

#[doc = "Original `Stats::addStats` at STAR/source/Stats.cpp:21. Args: S: Stats"]
pub fn stats_l21_stats_addstats(
    stats: &mut crate::stats::Stats,
    s: &crate::stats::Stats,
) {
    stats.read_n += s.read_n;
    stats.read_bases += s.read_bases;
    stats.mapped_mismatches_n += s.mapped_mismatches_n;
    stats.mapped_ins_n += s.mapped_ins_n;
    stats.mapped_del_n += s.mapped_del_n;
    stats.mapped_ins_l += s.mapped_ins_l;
    stats.mapped_del_l += s.mapped_del_l;
    stats.mapped_bases += s.mapped_bases;
    stats.mapped_portion += s.mapped_portion;
    stats.mapped_reads_u += s.mapped_reads_u;
    stats.mapped_reads_m += s.mapped_reads_m;
    stats.unmapped_other += s.unmapped_other;
    stats.unmapped_short += s.unmapped_short;
    stats.unmapped_mismatch += s.unmapped_mismatch;
    stats.unmapped_multi += s.unmapped_multi;
    stats.unmapped_all += s.unmapped_all;
    stats.chimeric_all += s.chimeric_all;
    stats.splices_nsjdb += s.splices_nsjdb;
    for ii in 0..SJ_MOTIF_SIZE {
        stats.splices_n[ii] += s.splices_n[ii];
    }
}

#[doc = "Original `Stats::transcriptStats` at STAR/source/Stats.cpp:35. Args: T: Transcript, Lread: uint"]
pub fn stats_l35_stats_transcriptstats(
    stats: &mut crate::stats::Stats,
    t: &crate::transcript::Transcript,
    l_read: u64,
) {
    stats.mapped_mismatches_n += t.n_mm as u32;
    stats.mapped_ins_n += t.n_ins as u32;
    stats.mapped_del_n += t.n_del as u32;
    stats.mapped_ins_l += t.l_ins as u32;
    stats.mapped_del_l += t.l_del as u32;

    if t.n_exons == 0 {
        return;
    }

    let mut mapped_l = 0_u64;
    for ii in 0..t.n_exons as usize {
        mapped_l += t.exons[ii][EX_L];
    }
    for ii in 0..t.n_exons as usize - 1 {
        if t.canon_sj[ii] >= 0 {
            stats.splices_n[t.canon_sj[ii] as usize] += 1;
        }
        if t.sj_annot[ii] == 1 {
            stats.splices_nsjdb += 1;
        }
    }

    stats.mapped_bases += mapped_l as u32;
    stats.mapped_portion += mapped_l as f64 / l_read as f64;
}

#[doc = "Original `Stats::progressReportHeader` at STAR/source/Stats.cpp:62. Args: progressStream: ofstream"]
pub fn stats_l62_stats_progressreportheader() -> String {
    use std::fmt::Write;

    let mut out = String::new();
    write!(
        out,
        "{:>15}{:>9}{:>12}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}\n",
        "Time",
        "Speed",
        "Read",
        "Read",
        "Mapped",
        "Mapped",
        "Mapped",
        "Mapped",
        "Unmapped",
        "Unmapped",
        "Unmapped",
        "Unmapped"
    )
    .unwrap();
    write!(
        out,
        "{:>15}{:>9}{:>12}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}\n",
        " ",
        "M/hr",
        "number",
        "length",
        "unique",
        "length",
        "MMrate",
        "multi",
        "multi+",
        "MM",
        "short",
        "other"
    )
    .unwrap();
    out
}

#[doc = "Original `Stats::progressReport` at STAR/source/Stats.cpp:73. Args: progressStream: ofstream"]
pub fn stats_l73_stats_progressreport(
    stats: &mut crate::stats::Stats,
    time_current: libc::time_t,
) -> Option<String> {
    if (time_current - stats.time_last_report) as f64 >= 60.0 && stats.read_n > 0 {
        let read_n = stats.read_n as f64;
        let elapsed = (time_current - stats.time_start_map) as f64;
        let out = format!(
            "{:>15}{:>9.1}{:>12}{:>9}{:>8.1}%{:>9.1}{:>8.1}%{:>8.1}%{:>8.1}%{:>8.1}%{:>8.1}%{:>8.1}%\n",
            timefunctions_l14_timemonthdaytime(time_current),
            stats.read_n as f64 / 1e6 / elapsed * 3600.0,
            stats.read_n,
            if stats.read_n > 0 {
                stats.read_bases / stats.read_n
            } else {
                0
            },
            if stats.read_n > 0 {
                stats.mapped_reads_u as f64 / read_n * 100.0
            } else {
                0.0
            },
            if stats.read_n > 0 {
                stats.mapped_bases as f64 / stats.mapped_reads_u as f64
            } else {
                0.0
            },
            if stats.read_n > 0 {
                stats.mapped_mismatches_n as f64 / stats.mapped_bases as f64 * 100.0
            } else {
                0.0
            },
            if stats.read_n > 0 {
                stats.mapped_reads_m as f64 / read_n * 100.0
            } else {
                0.0
            },
            if stats.read_n > 0 {
                stats.unmapped_multi as f64 / read_n * 100.0
            } else {
                0.0
            },
            if stats.read_n > 0 {
                stats.unmapped_mismatch as f64 / read_n * 100.0
            } else {
                0.0
            },
            if stats.read_n > 0 {
                stats.unmapped_short as f64 / read_n * 100.0
            } else {
                0.0
            },
            if stats.read_n > 0 {
                stats.unmapped_other as f64 / read_n * 100.0
            } else {
                0.0
            },
        );
        stats.time_last_report = time_current;
        Some(out)
    } else {
        None
    }
}

#[doc = "Original `Stats::reportFinal` at STAR/source/Stats.cpp:99. Args: streamOut: ofstream"]
pub fn stats_l99_stats_reportfinal(
    stats: &mut crate::stats::Stats,
    time_finish: libc::time_t,
) -> String {
    use std::fmt::Write;

    stats.time_finish = time_finish;
    let w1 = 50;
    let read_n = stats.read_n as f64;
    let mut out = String::new();
    writeln!(
        out,
        "{:>w1$}{}",
        "Started job on |\t",
        timefunctions_l14_timemonthdaytime(stats.time_start),
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Started mapping on |\t",
        timefunctions_l14_timemonthdaytime(stats.time_start_map),
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Finished on |\t",
        timefunctions_l14_timemonthdaytime(stats.time_finish),
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}",
        "Mapping speed, Million of reads per hour |\t",
        stats.read_n as f64 / 1e6 / (stats.time_finish - stats.time_start_map) as f64 * 3600.0,
        w1 = w1
    )
    .unwrap();
    out.push('\n');
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of input reads |\t",
        stats.read_n,
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Average input read length |\t",
        if stats.read_n > 0 {
            stats.read_bases / stats.read_n
        } else {
            0
        },
        w1 = w1
    )
    .unwrap();
    writeln!(out, "{:>w1$}", "UNIQUE READS:", w1 = w1).unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Uniquely mapped reads number |\t",
        stats.mapped_reads_u,
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}%",
        "Uniquely mapped reads % |\t",
        if stats.read_n > 0 {
            stats.mapped_reads_u as f64 / read_n * 100.0
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}",
        "Average mapped length |\t",
        if stats.mapped_reads_u > 0 {
            stats.mapped_bases as f64 / stats.mapped_reads_u as f64
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();

    writeln!(
        out,
        "{:>w1$}{}",
        "Number of splices: Total |\t",
        stats.splices_n.iter().sum::<u32>(),
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of splices: Annotated (sjdb) |\t",
        stats.splices_nsjdb,
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of splices: GT/AG |\t",
        stats.splices_n[1] + stats.splices_n[2],
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of splices: GC/AG |\t",
        stats.splices_n[3] + stats.splices_n[4],
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of splices: AT/AC |\t",
        stats.splices_n[5] + stats.splices_n[6],
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of splices: Non-canonical |\t",
        stats.splices_n[0],
        w1 = w1
    )
    .unwrap();

    writeln!(
        out,
        "{:>w1$}{:.2}%",
        "Mismatch rate per base, % |\t",
        stats.mapped_mismatches_n as f64 / stats.mapped_bases as f64 * 100.0,
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}%",
        "Deletion rate per base |\t",
        if stats.mapped_bases > 0 {
            stats.mapped_del_l as f64 / stats.mapped_bases as f64 * 100.0
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}",
        "Deletion average length |\t",
        if stats.mapped_del_n > 0 {
            stats.mapped_del_l as f64 / stats.mapped_del_n as f64
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}%",
        "Insertion rate per base |\t",
        if stats.mapped_bases > 0 {
            stats.mapped_ins_l as f64 / stats.mapped_bases as f64 * 100.0
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}",
        "Insertion average length |\t",
        if stats.mapped_ins_n > 0 {
            stats.mapped_ins_l as f64 / stats.mapped_ins_n as f64
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();
    writeln!(out, "{:>w1$}", "MULTI-MAPPING READS:", w1 = w1).unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of reads mapped to multiple loci |\t",
        stats.mapped_reads_m,
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}%",
        "% of reads mapped to multiple loci |\t",
        if stats.read_n > 0 {
            stats.mapped_reads_m as f64 / read_n * 100.0
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of reads mapped to too many loci |\t",
        stats.unmapped_multi,
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}%",
        "% of reads mapped to too many loci |\t",
        if stats.read_n > 0 {
            stats.unmapped_multi as f64 / read_n * 100.0
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();
    writeln!(out, "{:>w1$}", "UNMAPPED READS:", w1 = w1).unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of reads unmapped: too many mismatches |\t",
        stats.unmapped_mismatch,
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}%",
        "% of reads unmapped: too many mismatches |\t",
        if stats.read_n > 0 {
            stats.unmapped_mismatch as f64 / read_n * 100.0
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of reads unmapped: too short |\t",
        stats.unmapped_short,
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}%",
        "% of reads unmapped: too short |\t",
        if stats.read_n > 0 {
            stats.unmapped_short as f64 / read_n * 100.0
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of reads unmapped: other |\t",
        stats.unmapped_other,
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}%",
        "% of reads unmapped: other |\t",
        if stats.read_n > 0 {
            stats.unmapped_other as f64 / read_n * 100.0
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();
    writeln!(out, "{:>w1$}", "CHIMERIC READS:", w1 = w1).unwrap();
    writeln!(
        out,
        "{:>w1$}{}",
        "Number of chimeric reads |\t",
        stats.chimeric_all,
        w1 = w1
    )
    .unwrap();
    writeln!(
        out,
        "{:>w1$}{:.2}%",
        "% of chimeric reads |\t",
        if stats.read_n > 0 {
            stats.chimeric_all as f64 / read_n * 100.0
        } else {
            0.0
        },
        w1 = w1
    )
    .unwrap();

    out
}

#[doc = "Original `Stats::writeLines` at STAR/source/Stats.cpp:147. Args: streamOut: ofstream, outType: vector<int>, commStr: string, outStr: string"]
pub fn stats_l147_stats_writelines(
    stats: &crate::stats::Stats,
    out_type: &[i32],
    comm_str: &str,
    out_str: &str,
) -> String {
    use std::fmt::Write;

    let mut out = String::new();
    for tt in out_type {
        if *tt == 1 {
            if !out_str.is_empty() {
                writeln!(out, "{comm_str} {out_str}").unwrap();
            }
            writeln!(
                out,
                "{comm_str} Nreads {}\tNreadsUnique {}\tNreadsMulti {}",
                stats.read_n, stats.mapped_reads_u, stats.mapped_reads_m
            )
            .unwrap();
        }
    }
    out
}
