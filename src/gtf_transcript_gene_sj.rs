#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `GTF::transcriptGeneSJ` at STAR/source/GTF_transcriptGeneSJ.cpp:23. Args: dirOut: string"]
pub fn gtf_transcriptgenesj_l23_gtf_transcriptgenesj(
    gtf: &mut crate::gtf::GTF,
    genome: &crate::genome::Genome,
    sjdb_loci: &mut crate::sjdb_class::SjdbClass,
    _dir_out: &str,
    log_main: &mut String,
) -> crate::gtf::GtfTranscriptGeneSjOutput {
    let mut out = crate::gtf::GtfTranscriptGeneSjOutput::default();
    if !gtf.gtf_yes {
        return out;
    }

    gtf.exon_n = gtf.exon_loci.len() as u64;
    gtf.exon_loci.sort_by_key(|ex| (ex[GTF_EX_T], ex[GTF_EX_S]));

    let exon_n = gtf.exon_n as usize;
    let mut exge_loci = Vec::<[u64; 5]>::with_capacity(exon_n);
    for exon in &gtf.exon_loci {
        exge_loci.push([
            exon[GTF_EX_S],
            exon[GTF_EX_E],
            gtf.transcript_strand[exon[GTF_EX_T] as usize] as u64,
            exon[GTF_EX_G],
            exon[GTF_EX_T],
        ]);
    }
    exge_loci.sort();

    out.exon_ge_tr_info_tab
        .push_str(&format!("{}\n", gtf.exon_n));
    for exge in &exge_loci {
        out.exon_ge_tr_info_tab.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            exge[0], exge[1], exge[2], exge[3], exge[4]
        ));
    }

    out.gene_info_tab
        .push_str(&format!("{}\n", gtf.gene_id.len()));
    for ig in 0..gtf.gene_id.len() {
        out.gene_info_tab.push_str(&format!(
            "{}\t{}\t{}\n",
            gtf.gene_id[ig], gtf.gene_attr[ig][0], gtf.gene_attr[ig][1]
        ));
    }

    let mut extr_loci = Vec::<[u64; 6]>::with_capacity(exon_n);
    if exon_n > 0 {
        extr_loci.resize(exon_n, [0; 6]);
        let mut trex1 = 0usize;
        for iex in 0..=exon_n {
            if iex == exon_n || gtf.exon_loci[iex][GTF_EX_T] != gtf.exon_loci[trex1][GTF_EX_T] {
                for iex1 in trex1..iex {
                    extr_loci[iex1][1] = gtf.exon_loci[iex - 1][GTF_EX_E];
                }
                if iex == exon_n {
                    break;
                }
                trex1 = iex;
            }
            extr_loci[iex][0] = gtf.exon_loci[trex1][GTF_EX_S];
            extr_loci[iex][2] = gtf.exon_loci[iex][GTF_EX_T];
            extr_loci[iex][3] = gtf.exon_loci[iex][GTF_EX_S];
            extr_loci[iex][4] = gtf.exon_loci[iex][GTF_EX_E];
            extr_loci[iex][5] = gtf.exon_loci[iex][GTF_EX_G];
        }
    }
    extr_loci.sort_by(|a, b| a[..5].cmp(&b[..5]));

    out.transcript_info_tab
        .push_str(&format!("{}\n", gtf.transcript_id.len()));
    out.exon_info_tab.push_str(&format!("{}\n", gtf.exon_n));
    if exon_n > 0 {
        let mut trid = extr_loci[0][2];
        let mut trex = 0usize;
        let mut trstart = extr_loci[0][0];
        let mut trend = extr_loci[0][1];
        let mut exlen = 0u64;
        for iex in 0..=exon_n {
            if iex == exon_n || extr_loci[iex][2] != trid {
                out.transcript_info_tab.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    gtf.transcript_id[trid as usize],
                    extr_loci[iex - 1][0],
                    extr_loci[iex - 1][1],
                    trend,
                    gtf.transcript_strand[trid as usize],
                    iex - trex,
                    trex,
                    extr_loci[iex - 1][5]
                ));
                if iex == exon_n {
                    break;
                }
                trid = extr_loci[iex][2];
                trstart = extr_loci[iex][0];
                trex = iex;
                trend = std::cmp::max(trend, extr_loci[iex - 1][1]);
                exlen = 0;
            }
            out.exon_info_tab.push_str(&format!(
                "{}\t{}\t{}\n",
                extr_loci[iex][3] - trstart,
                extr_loci[iex][4] - trstart,
                exlen
            ));
            exlen += extr_loci[iex][4] - extr_loci[iex][3] + 1;
        }
    }

    let mut sj_loci = Vec::<[u64; 4]>::new();
    if exon_n > 0 {
        let mut tr_id_n = gtf.exon_loci[0][GTF_EX_T];
        for iex in 1..exon_n {
            if tr_id_n == gtf.exon_loci[iex][GTF_EX_T] {
                if gtf.exon_loci[iex][GTF_EX_S] <= gtf.exon_loci[iex - 1][GTF_EX_E] + 1 {
                } else if gtf.exon_loci[iex][GTF_EX_S] <= gtf.exon_loci[iex - 1][GTF_EX_E] {
                    let chr1 = genome.chr_bin
                        [(gtf.exon_loci[iex][GTF_EX_S] >> genome.p_ge.g_chr_bin_nbits) as usize]
                        as usize;
                    log_main.push_str(&format!(
                        "WARNING: while processing pGe.sjdbGTFfile={}: overlapping exons:\n",
                        genome.p_ge.sjdb_gtf_file
                    ));
                    log_main.push_str(&format!(
                        "{}\t{}\t{}\n",
                        genome.chr_name[chr1],
                        gtf.exon_loci[iex - 1][GTF_EX_S] + 1 - genome.chr_start[chr1],
                        gtf.exon_loci[iex - 1][GTF_EX_E] + 1 - genome.chr_start[chr1]
                    ));
                    log_main.push_str(&format!(
                        "{}\t{}\t{}\n",
                        genome.chr_name[chr1],
                        gtf.exon_loci[iex][GTF_EX_S] + 1 - genome.chr_start[chr1],
                        gtf.exon_loci[iex][GTF_EX_E] + 1 - genome.chr_start[chr1]
                    ));
                } else {
                    sj_loci.push([
                        gtf.exon_loci[iex - 1][GTF_EX_E] + 1,
                        gtf.exon_loci[iex][GTF_EX_S] - 1,
                        gtf.transcript_strand[tr_id_n as usize] as u64,
                        gtf.exon_loci[iex][GTF_EX_G] + 1,
                    ]);
                }
            } else {
                tr_id_n = gtf.exon_loci[iex][GTF_EX_T];
            }
        }
    }
    sj_loci.sort_by_key(|sj| (sj[0], sj[1]));

    let strand_char = ['.', '+', '-'];
    let sjdb_n1 = sjdb_loci.chr.len();
    sjdb_loci
        .gene
        .resize(sjdb_n1, std::collections::BTreeSet::new());
    for ii in 0..sj_loci.len() {
        if ii == 0
            || sj_loci[ii][0] != sj_loci[ii - 1][0]
            || sj_loci[ii][1] != sj_loci[ii - 1][1]
            || sj_loci[ii][2] != sj_loci[ii - 1][2]
        {
            let chr1 =
                genome.chr_bin[(sj_loci[ii][0] >> genome.p_ge.g_chr_bin_nbits) as usize] as usize;
            sjdb_loci.chr.push(genome.chr_name[chr1].clone());
            sjdb_loci
                .start
                .push(sj_loci[ii][0] + 1 - genome.chr_start[chr1]);
            sjdb_loci
                .end
                .push(sj_loci[ii][1] + 1 - genome.chr_start[chr1]);
            sjdb_loci.str_.push(strand_char[sj_loci[ii][2] as usize]);
            sjdb_loci
                .gene
                .push(std::iter::once(sj_loci[ii][3]).collect());
        } else if let Some(gene) = sjdb_loci.gene.last_mut() {
            gene.insert(sj_loci[ii][3]);
        }
    }

    for ii in sjdb_n1..sjdb_loci.chr.len() {
        out.sjdb_list_from_gtf_out_tab.push_str(&format!(
            "{}\t{}\t{}\t{}",
            sjdb_loci.chr[ii], sjdb_loci.start[ii], sjdb_loci.end[ii], sjdb_loci.str_[ii]
        ));
        let mut genes = sjdb_loci.gene[ii].iter();
        if let Some(first) = genes.next() {
            out.sjdb_list_from_gtf_out_tab
                .push_str(&format!("\t{}", first));
            for gene in genes {
                out.sjdb_list_from_gtf_out_tab
                    .push_str(&format!(",{}", gene));
            }
        }
        out.sjdb_list_from_gtf_out_tab.push('\n');
    }

    sjdb_loci.priority.resize(sjdb_loci.chr.len(), 20);
    let added = sjdb_loci.chr.len() - sjdb_n1;
    log_main.push_str(&format!(
        "Processing pGe.sjdbGTFfile={}, found:\n\t\t{} transcripts\n\t\t{} exons (non-collapsed)\n\t\t{} collapsed junctions\nTotal junctions: {}\n",
        genome.p_ge.sjdb_gtf_file,
        gtf.transcript_id.len(),
        gtf.exon_n,
        added,
        sjdb_loci.chr.len()
    ));
    log_main.push_str(" ..... finished GTF processing\n\n");

    out.n_junctions_added = added as u64;
    out
}
