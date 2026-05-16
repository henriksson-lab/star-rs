#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `GTF::superTranscript` at STAR/source/GTF_superTranscript.cpp:9. Args: "]
pub fn gtf_supertranscript_l9_gtf_supertranscript(
    gtf: &mut crate::gtf::GTF,
    genome: &mut crate::genome::Genome,
    p: &mut crate::parameters_chimeric::Parameters,
) -> crate::gtf::GtfSuperTranscriptOutput {
    let mut out = crate::gtf::GtfSuperTranscriptOutput::default();
    if p.p_ge.g_type_string != "Transcriptome" && p.p_ge.g_type_string != "SuperTranscriptome" {
        return out;
    }

    for ii in 0..genome.n_chr_real as usize {
        out.full_genome_chr_name_txt
            .push_str(&format!("{}\n", genome.chr_name[ii]));
        out.full_genome_chr_start_txt
            .push_str(&format!("{}\n", genome.chr_start[ii]));
        out.full_genome_chr_length_txt
            .push_str(&format!("{}\n", genome.chr_length[ii]));
        out.full_genome_chr_name_length_txt.push_str(&format!(
            "{}\t{}\n",
            genome.chr_name[ii], genome.chr_length[ii]
        ));
    }
    if genome.chr_start.len() > genome.n_chr_real as usize {
        out.full_genome_chr_start_txt.push_str(&format!(
            "{}\n",
            genome.chr_start[genome.n_chr_real as usize]
        ));
    }
    out.full_genome_sequence = genome.g[..genome.n_genome as usize].to_vec();

    let n_minus_strand_offset = genome.n_genome;
    for exon in gtf.exon_loci.iter_mut() {
        let trans_id = exon[GTF_EX_T] as usize;
        if gtf.transcript_strand[trans_id] == 2 {
            let temp = exon[GTF_EX_S];
            exon[GTF_EX_S] = 2 * n_minus_strand_offset - 1 - exon[GTF_EX_E];
            exon[GTF_EX_E] = 2 * n_minus_strand_offset - 1 - temp;
        }
    }

    let mut merged_intervals = Vec::<[u64; 2]>::new();
    if gtf.exon_loci.is_empty() {
        return out;
    }

    gtf.exon_loci.sort_by_key(|exon| exon[GTF_EX_S]);
    let mut gap_value = gtf.exon_loci[0][GTF_EX_S];
    let mut curr = [gtf.exon_loci[0][GTF_EX_S], gtf.exon_loci[0][GTF_EX_E]];
    gtf.exon_loci[0][GTF_EX_S] = 0;
    gtf.exon_loci[0][GTF_EX_E] -= gap_value;

    for ii in 1..gtf.exon_loci.len() {
        if gtf.exon_loci[ii][GTF_EX_S] <= curr[1] + 1 {
            curr[1] = std::cmp::max(curr[1], gtf.exon_loci[ii][GTF_EX_E]);
        } else {
            gap_value += gtf.exon_loci[ii][GTF_EX_S] - curr[1] - 1;
            merged_intervals.push(curr);
            curr = [gtf.exon_loci[ii][GTF_EX_S], gtf.exon_loci[ii][GTF_EX_E]];
        }
        gtf.exon_loci[ii][GTF_EX_S] -= gap_value;
        gtf.exon_loci[ii][GTF_EX_E] -= gap_value;
    }
    merged_intervals.push(curr);

    gtf.super_trome.seq_concat.clear();
    for interval in &merged_intervals {
        for pos in interval[0]..=interval[1] {
            gtf.super_trome.seq_concat.push(genome.g[pos as usize]);
        }
    }
    out.log_main.push_str(&format!(
        "SuperTranscriptome (condensed) genome length = {}\n",
        gtf.super_trome.seq_concat.len()
    ));

    gtf.transcript_start_end
        .resize(gtf.transcript_id.len(), [u64::MAX, 0]);
    for exon in &gtf.exon_loci {
        let trans_id = exon[GTF_EX_T] as usize;
        gtf.transcript_start_end[trans_id][0] =
            std::cmp::min(gtf.transcript_start_end[trans_id][0], exon[GTF_EX_S]);
        gtf.transcript_start_end[trans_id][1] =
            std::cmp::max(gtf.transcript_start_end[trans_id][1], exon[GTF_EX_E]);
    }

    let mut super_tr_start_end = Vec::<[u64; 2]>::new();
    let mut merged_intervals_super_tr_index = vec![0u64; merged_intervals.len()];
    let mut mi_i = 0usize;
    let mut mi_len = 0u64;
    let mut transcript_start_end_sorted = gtf.transcript_start_end.clone();
    transcript_start_end_sorted.sort_by_key(|tr| tr[0]);
    let mut curr_tr = transcript_start_end_sorted[0];
    for tr in &transcript_start_end_sorted {
        while mi_i < merged_intervals.len() && mi_len < curr_tr[1] {
            mi_len += merged_intervals[mi_i][1] - merged_intervals[mi_i][0] + 1;
            merged_intervals_super_tr_index[mi_i] = super_tr_start_end.len() as u64;
            mi_i += 1;
        }
        curr_tr[1] = std::cmp::max(curr_tr[1], mi_len.saturating_sub(1));

        if tr[0] <= curr_tr[1] {
            curr_tr[1] = std::cmp::max(curr_tr[1], tr[1]);
        } else {
            super_tr_start_end.push(curr_tr);
            curr_tr = *tr;
        }
    }
    super_tr_start_end.push(curr_tr);

    gtf.super_trome.seq.clear();
    let mut max_st_len = 0u64;
    for tr in &super_tr_start_end {
        gtf.super_trome
            .seq
            .push(gtf.super_trome.seq_concat[tr[0] as usize..=tr[1] as usize].to_vec());
        max_st_len = std::cmp::max(max_st_len, tr[1] - tr[0]);
    }
    out.log_main.push_str(&format!(
        "Number of superTranscripts = {};   max length = {}\n",
        super_tr_start_end.len(),
        max_st_len
    ));

    gtf.super_trome.tr_index.resize(gtf.transcript_id.len(), 0);
    let mut ist = 0usize;
    for exon in &gtf.exon_loci {
        if exon[GTF_EX_S] > super_tr_start_end[ist][1] {
            ist += 1;
        }
        gtf.super_trome.tr_index[exon[GTF_EX_T] as usize] = ist as u64;
    }

    gtf.super_trome
        .tr_start_end
        .resize(gtf.transcript_start_end.len(), [0, 0]);
    for ii in 0..gtf.transcript_start_end.len() {
        let st_i = gtf.super_trome.tr_index[ii] as usize;
        gtf.super_trome.tr_start_end[ii][0] =
            gtf.transcript_start_end[ii][0] - super_tr_start_end[st_i][0];
        gtf.super_trome.tr_start_end[ii][1] =
            gtf.transcript_start_end[ii][1] - super_tr_start_end[st_i][0];
    }

    gtf.exon_loci
        .sort_by_key(|exon| (exon[GTF_EX_T], exon[GTF_EX_S]));
    gtf.transcript_seq.clear();
    gtf.transcript_seq
        .resize(gtf.transcript_id.len(), Vec::new());
    for exon in &gtf.exon_loci {
        let trans_id = exon[GTF_EX_T] as usize;
        gtf.transcript_seq[trans_id].extend_from_slice(
            &gtf.super_trome.seq_concat[exon[GTF_EX_S] as usize..=exon[GTF_EX_E] as usize],
        );
    }

    let num_to_char = ['A', 'C', 'G', 'T', 'N'];
    for ii in 0..gtf.transcript_seq.len() {
        out.transcript_sequences_fasta
            .push_str(&format!(">{}\n", gtf.transcript_id[ii]));
        for &base in &gtf.transcript_seq[ii] {
            out.transcript_sequences_fasta
                .push(num_to_char[base as usize]);
        }
        out.transcript_sequences_fasta.push('\n');
    }
    for ii in 0..gtf.super_trome.seq.len() {
        out.super_transcript_sequences_fasta
            .push_str(&format!(">st{}\n", ii));
        for &base in &gtf.super_trome.seq[ii] {
            out.super_transcript_sequences_fasta
                .push(num_to_char[base as usize]);
        }
        out.super_transcript_sequences_fasta.push('\n');
    }

    gtf.super_trome.sj.clear();
    for ii in 1..gtf.exon_loci.len() {
        if gtf.exon_loci[ii][GTF_EX_T] == gtf.exon_loci[ii - 1][GTF_EX_T]
            && gtf.exon_loci[ii][GTF_EX_S] > gtf.exon_loci[ii - 1][GTF_EX_E] + 1
        {
            let sti = gtf.super_trome.tr_index[gtf.exon_loci[ii][GTF_EX_T] as usize];
            let sts = super_tr_start_end[sti as usize][0];
            gtf.super_trome.sj.push(crate::super_transcriptome::sjInfo {
                start: (gtf.exon_loci[ii - 1][GTF_EX_E] - sts) as u32,
                end: (gtf.exon_loci[ii][GTF_EX_S] - sts) as u32,
                tr: gtf.exon_loci[ii][GTF_EX_T] as u32,
                super_: sti as u32,
            });
        }
    }
    let (sj_tsv, sj_log) =
        supertranscriptome_l4_supertranscriptome_sjcollapse(&mut gtf.super_trome);
    out.super_transcript_sj_tsv = sj_tsv;
    out.log_main.push_str(&sj_log);

    if p.p_ge.g_type_string == "Transcriptome" {
        gtf_supertranscript_l232_genome_concatenatechromosomes(
            genome,
            &gtf.transcript_seq,
            &gtf.transcript_id,
            genome.genome_chr_bin_nbases as u64,
        );
        gtf.gtf_yes = false;
        p.sjdb_insert_yes = false;
    } else if p.p_ge.g_type_string == "SuperTranscriptome" {
        let super_transcript_id: Vec<String> = (0..super_tr_start_end.len())
            .map(|ii| format!("st{}", ii))
            .collect();
        gtf_supertranscript_l232_genome_concatenatechromosomes(
            genome,
            &gtf.super_trome.seq,
            &super_transcript_id,
            genome.genome_chr_bin_nbases as u64,
        );

        gtf.transcript_strand.resize(super_tr_start_end.len(), 1);
        gtf.exon_loci.sort_by_key(|exon| exon[GTF_EX_S]);
        let mut ist = 0usize;
        for exon in gtf.exon_loci.iter_mut() {
            if exon[GTF_EX_S] > super_tr_start_end[ist][1] {
                ist += 1;
            }
            exon[GTF_EX_S] += genome.chr_start[ist] - super_tr_start_end[ist][0];
            exon[GTF_EX_E] += genome.chr_start[ist] - super_tr_start_end[ist][0];
        }
    }

    out.conversion_to_full_genome_tsv.push_str(&format!(
        "{}\t{}\n",
        merged_intervals.len(),
        n_minus_strand_offset
    ));
    let mut cond_gstart = 0u64;
    for ib in 0..merged_intervals.len() {
        let block = merged_intervals[ib];
        let len1 = block[1] - block[0] + 1;
        let i_sutr = merged_intervals_super_tr_index[ib] as usize;
        out.conversion_to_full_genome_tsv.push_str(&format!(
            "{}\t{}\t{}\n",
            genome.chr_start[i_sutr] + cond_gstart - super_tr_start_end[i_sutr][0],
            len1,
            block[0]
        ));
        cond_gstart += len1;
    }

    out
}

#[doc = "Original `Genome::concatenateChromosomes` at STAR/source/GTF_superTranscript.cpp:232. Args: vecSeq: vector<vector<uint8>>, vecName: vector<string>, padBin: uint64"]
pub fn gtf_supertranscript_l232_genome_concatenatechromosomes(
    genome: &mut crate::genome::Genome,
    vec_seq: &[Vec<u8>],
    vec_name: &[String],
    pad_bin: u64,
) {
    genome.n_chr_real = vec_seq.len() as u32;
    genome.chr_length.resize(genome.n_chr_real as usize, 0);
    genome.chr_start.resize(genome.n_chr_real as usize + 1, 0);
    genome.chr_name = vec_name.to_vec();
    genome.chr_name_index.clear();

    for ii in 0..genome.n_chr_real as usize {
        genome.chr_length[ii] = vec_seq[ii].len() as u64;
        genome.chr_start[ii + 1] =
            genome.chr_start[ii] + ((genome.chr_length[ii] + 1) / pad_bin + 1) * pad_bin;
        genome
            .chr_name_index
            .insert(genome.chr_name[ii].clone(), ii as u64);
    }

    genome.n_genome = *genome.chr_start.last().unwrap_or(&0);
    genome.g = vec![GENOME_SPACING_CHAR; (2 * genome.n_genome) as usize];
    for ii in 0..genome.n_chr_real as usize {
        let start = genome.chr_start[ii] as usize;
        genome.g[start..start + vec_seq[ii].len()].copy_from_slice(&vec_seq[ii]);
    }

    for ii in 0..genome.n_genome as usize {
        genome.g[2 * genome.n_genome as usize - 1 - ii] = if genome.g[ii] < 4 {
            3 - genome.g[ii]
        } else {
            genome.g[ii]
        };
    }
}
