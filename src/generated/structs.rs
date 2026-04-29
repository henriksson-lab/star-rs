#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

#[doc = "Original class `BAMoutput` at STAR/source/BAMoutput.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BAMoutput {
    pub n_bins: u32,
    pub bin_total_n: Vec<u64>,
    pub bin_total_bytes: Vec<u64>,
    pub bam_array_size: u64,
    pub bam_array: Vec<u8>,
    pub bin_size: u64,
    pub bin_size1: u64,
    pub bin_bytes: Vec<u64>,
    pub bin_bytes1: u64,
    pub bin_buffers: Vec<Vec<u8>>,
    pub bin_streams: Vec<Vec<u8>>,
    pub bin_stream_by_sjout: Vec<bool>,
    pub bgzf_bam: Vec<u8>,
    pub bam_dir: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BamCore {
    pub tid: i32,
    pub pos: i32,
    pub bin: u32,
    pub qual: u32,
    pub l_qname: u32,
    pub flag: u32,
    pub n_cigar: u32,
    pub l_qseq: i32,
    pub mtid: i32,
    pub mpos: i32,
    pub isize: i32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Bam1 {
    pub core: BamCore,
    pub l_data: i32,
    pub m_data: i32,
    pub data_offset: usize,
}

#[doc = "Original class `Cell` at STAR/source/opal/opal.cpp:1206."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Cell {
    pub h: i64,
    pub e: i64,
    pub f: i64,
}

#[doc = "Original struct `CellEH` at STAR/source/opal/opal.cpp:156."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CellEH {
    pub h: (),
    pub e: (),
}

#[doc = "Original class `Chain` at STAR/source/Chain.h:16."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Chain {
    pub chain_file_name: String,
    pub chr_chains: std::collections::BTreeMap<String, OneChain>,
}

#[doc = "Original class `ChimericAlign` at STAR/source/ChimericAlign.h:14."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericAlign {
    pub seg1: ChimericSegment,
    pub seg2: ChimericSegment,
    pub chim_score: i32,
    pub chim_j1: u32,
    pub chim_j2: u32,
    pub chim_repeat1: u32,
    pub chim_repeat2: u32,
    pub chim_motif: i32,
    pub chim_str: i32,
    pub stitching_done: bool,
    pub al1: Transcript,
    pub al2: Transcript,
    pub ex1: u32,
    pub ex2: u32,
    pub junction_overhang_min: u32,
}

#[doc = "Original class `ChimericDetection` at STAR/source/ChimericDetection.h:12."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericDetection {
    pub p: Parameters,
    pub ra: Option<ReadAlign>,
    pub tr_all: Vec<Vec<Transcript>>,
    pub n_w: u32,
    pub n_win_tr: Vec<u32>,
    pub read1: [Vec<u8>; 2],
    pub out_gen: Genome,
    pub chim_aligns: Vec<ChimericAlign>,
    pub ostream_chim_junction_attached: bool,
}

#[doc = "Original class `ChimericSegment` at STAR/source/ChimericSegment.h:9."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericSegment {
    pub align: Transcript,
    pub ro_s: u32,
    pub ro_e: u32,
    pub str_: u32,
    pub segment_min: u32,
    pub segment_read_gap_max: u32,
}

#[doc = "Original class `ChimericTranscript` at STAR/source/ChimericTranscript.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericTranscript {}

#[doc = "Original class `ClipCR4` at STAR/source/ClipCR4.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ClipCR4 {
    pub db_n: u32,
    pub score_matrix: Vec<i32>,
    pub read_len: u32,
    pub alphabet_length: i32,
    pub gap_open: i32,
    pub gap_ext: i32,
    pub db_seq_arr: Vec<u8>,
    pub db_seqs_len: Vec<i32>,
    pub store_clip: Vec<u32>,
    pub opal_res: Vec<OpalSearchResult>,
}

#[doc = "Original class `ClipMate` at STAR/source/ClipMate.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ClipMate {
    pub type_: i32,
    pub n: u32,
    pub n_after_ad: u32,
    pub ad_seq: String,
    pub ad_seq_num: Vec<u8>,
    pub ad_mmp: f64,
    pub clipped_info: u8,
    pub clipped_ad_n: u32,
    pub clipped_ad_mm: u32,
    pub clipped_n: u32,
    pub cr4: Option<ClipCR4>,
}

#[doc = "Original struct `Data` at STAR/source/SimpleGoodTuring/sgt.h:73."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Data {}

#[doc = "Original class `GTF` at STAR/source/GTF.h:10."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GTF {
    pub gtf_yes: bool,
    pub exon_n: u64,
    pub exon_loci: Vec<[u64; crate::generated::functions::GTF_EX_L]>,
    pub transcript_strand: Vec<u32>,
    pub transcript_id: Vec<String>,
    pub gene_id: Vec<String>,
    pub gene_attr: Vec<[String; 2]>,
    pub transcript_seq: Vec<Vec<u8>>,
    pub transcript_start_end: Vec<[u64; 2]>,
    pub super_trome: SuperTranscriptome,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GtfTranscriptGeneSjOutput {
    pub n_junctions_added: u64,
    pub exon_ge_tr_info_tab: String,
    pub gene_info_tab: String,
    pub transcript_info_tab: String,
    pub exon_info_tab: String,
    pub sjdb_list_from_gtf_out_tab: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GtfSuperTranscriptOutput {
    pub transcript_sequences_fasta: String,
    pub super_transcript_sequences_fasta: String,
    pub super_transcript_sj_tsv: String,
    pub conversion_to_full_genome_tsv: String,
    pub full_genome_chr_name_txt: String,
    pub full_genome_chr_start_txt: String,
    pub full_genome_chr_length_txt: String,
    pub full_genome_chr_name_length_txt: String,
    pub full_genome_sequence: Vec<u8>,
    pub log_main: String,
}

#[doc = "Original struct `GeneInfo1` at STAR/source/Transcriptome_alignExonOverlap.cpp:25."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GeneInfo1 {
    pub g: u32,
    pub ia: u32,
    pub ot: [bool; 6],
}

#[doc = "Original struct `GeneStrOverlapAlign` at STAR/source/Transcriptome_alignExonOverlap.cpp:12."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GeneStrOverlapAlign {
    pub g: u32,
    pub ov: i64,
    pub ia: u32,
    pub exl: u32,
    pub str: bool,
    pub sjc: bool,
}

#[doc = "Original class `Genome` at STAR/source/Genome.h:13."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genome {
    pub g: Vec<u8>,
    pub sa: Vec<u32>,
    pub sa_packed: PackedArray,
    pub sa_insert: PackedArray,
    pub sa_pass2: Vec<u32>,
    pub sai: Vec<u32>,
    pub sai_packed: PackedArray,
    pub genome_sa_index_start: Vec<u32>,
    pub n_genome: u64,
    pub n_sa: u32,
    pub n_sa_byte: u32,
    pub n_chr_real: u32,
    pub genome_chr_bin_nbases: u32,
    pub chr_bin_n: u32,
    pub gstrand_bit: u32,
    pub gstrand_mask: u32,
    pub sai_mark_absent_mask_c: u32,
    pub sai_mark_nmask: u32,
    pub sai_mark_nmask_c: u32,
    pub sjdb_overhang: u32,
    pub sjdb_length: u32,
    pub sj_gstart: u64,
    pub sj_dstart: Vec<u64>,
    pub sj_astart: Vec<u64>,
    pub sjdb_start: Vec<u32>,
    pub sjdb_end: Vec<u32>,
    pub sjdb_motif: Vec<u8>,
    pub sjdb_shift_left: Vec<u8>,
    pub sjdb_shift_right: Vec<u8>,
    pub sjdb_strand: Vec<u8>,
    pub sjdb_n: u32,
    pub p_ge: ParametersGenome,
    pub chr_start: Vec<u64>,
    pub chr_length: Vec<u64>,
    pub chr_bin: Vec<u32>,
    pub chr_name: Vec<String>,
    pub chr_name_all: Vec<String>,
    pub chr_name_index: std::collections::BTreeMap<String, u64>,
    pub chr_length_all: Vec<u32>,
    pub genome_insert_l: u64,
    pub n_sa_insert: u64,
    pub var: Variation,
    pub genome_out: GenomeOut,
    pub align_sjdb_overhang_min: u32,
    pub align_intron_min: u32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GenomeOut {
    pub conv_yes: bool,
    pub gaps_are_junctions: bool,
    pub conv_file: String,
    pub conv_blocks: Vec<[u64; 3]>,
    pub n_minus_strand_offset: u64,
}

#[doc = "Original class `InOutStreams` at STAR/source/InOutStreams.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InOutStreams {
    pub log_stdout_attached: bool,
    pub out_sam_attached: bool,
    pub log_stdout_flushed: bool,
    pub out_sam_flushed: bool,
    pub log_stdout_file_flushed: bool,
    pub out_sam_file_flushed: bool,
    pub out_chim_sam_flushed: bool,
    pub out_chim_junction_flushed: bool,
    pub log_progress_flushed: bool,
    pub log_main_flushed: bool,
    pub log_final_flushed: bool,
    pub out_local_chains_flushed: bool,
    pub out_sam_file_closed: bool,
    pub out_chim_sam_closed: bool,
    pub out_chim_junction_closed: bool,
    pub log_progress_closed: bool,
    pub log_final_closed: bool,
    pub out_local_chains_closed: bool,
    pub out_unmapped_reads_open: [bool; 2],
    pub out_unmapped_reads_flushed: [bool; 2],
    pub out_unmapped_reads_closed: [bool; 2],
}

#[doc = "Original class `Junction` at STAR/source/OutSJ.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Junction {
    pub gen_out: Genome,
    pub record: Option<JunctionRecord>,
}

#[doc = "Original class `MultiMappers` at STAR/source/ParametersSolo.h:45."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MultiMappers {
    pub types_in: Vec<String>,
    pub types: Vec<i32>,
    pub type_main: i32,
    pub yes_multi: bool,
    pub yes_n: u32,
    pub yes_b: [bool; 5],
    pub count_ind_i: [u32; 5],
    pub count_ind_main: u32,
}

#[doc = "Original class `OneChain` at STAR/source/Chain.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OneChain {
    pub b_n: u32,
    pub chr1: String,
    pub chr2: String,
    pub b_start1: Vec<u32>,
    pub b_start2: Vec<u32>,
    pub b_len: Vec<u32>,
}

#[doc = "Original struct `OpalSearchResult` at STAR/source/opal/opal.h:47."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OpalSearchResult {
    pub score_set: i32,
    pub score: i32,
    pub end_location_target: i32,
    pub end_location_query: i32,
    pub start_location_target: i32,
    pub start_location_query: i32,
    pub alignment: Option<Vec<u8>>,
    pub alignment_length: i32,
}

#[doc = "Original class `OutSJ` at STAR/source/OutSJ.h:36."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OutSJ {
    pub n: u64,
    pub n_store: u64,
    pub junctions: Vec<JunctionRecord>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct JunctionRecord {
    pub start: u32,
    pub gap: u32,
    pub strand: i8,
    pub motif: i32,
    pub annot: u8,
    pub count_unique: u32,
    pub count_multiple: u32,
    pub overhang_left: u16,
    pub overhang_right: u16,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ExitWithErrorResult {
    pub stream_out1: String,
    pub stream_out2: String,
    pub error_int: i32,
    pub thread_mutex_locked: bool,
    pub in_out_deleted: bool,
}

#[doc = "Original class `PackedArray` at STAR/source/PackedArray.h:6."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PackedArray {
    pub bit_rec_mask: u64,
    pub word_comp_length: u64,
    pub array_allocated: bool,
    pub word_length: u64,
    pub length: u64,
    pub length_byte: u64,
    pub char_array: Vec<u8>,
}

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

#[doc = "Original class `Parameters` at STAR/source/ParametersChimeric.h:6."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Parameters {
    pub command_line_full: String,
    pub command_line: String,
    pub version_genome: String,
    pub out_file_tmp: String,
    pub run_thread_n: i32,
    pub run_rng_seed: i32,
    pub run_dir_perm: u32,
    pub out_file_name_prefix: String,
    pub out_tmp_keep: String,
    pub run_restart_type: i32,
    pub read_nmates: u32,
    pub read_nends: u32,
    pub read_quality_score_base: u8,
    pub out_std: String,
    pub p_ge: ParametersGenome,
    pub p_ch: ParametersChimeric,
    pub input_bam_file: String,
    pub bam_remove_duplicates_type: String,
    pub bam_remove_duplicates_yes: bool,
    pub bam_remove_duplicates_mark_multi: bool,
    pub bam_remove_duplicates_mate2bases_n: u32,
    pub out_bam_unsorted: bool,
    pub out_bam_coord: bool,
    pub out_bam_compression: i32,
    pub out_sam_bool: bool,
    pub out_wig_references_prefix: String,
    pub out_wig_flags: OutWigFlags,
    pub out_sam_type: Vec<String>,
    pub out_sam_header_hd: Vec<String>,
    pub out_sam_header_pg: Vec<String>,
    pub out_sam_header_comment_file: String,
    pub out_sam_attr_rgline_split: Vec<String>,
    pub quant_tr_sam_yes: bool,
    pub quant_tr_sam_bam_yes: bool,
    pub sam_header_hd: String,
    pub sam_header: String,
    pub sam_header_sorted_coord: String,
    pub out_sam_contents: String,
    pub out_bam_unsorted_header: Vec<u8>,
    pub out_quant_bam_header: Vec<u8>,
    pub p_solo: ParametersSolo,
    pub sj_all: [Vec<u64>; 2],
    pub genome_dir: String,
    pub out_sam_attr_nm_present: bool,
    pub out_sam_attr_order: Vec<i32>,
    pub out_sam_attr_order_quant: Vec<i32>,
    pub out_sam_attributes: Vec<String>,
    pub out_sam_attr_present: SamAttrPresent,
    pub out_sam_attr_present_quant: SamAttrPresent,
    pub out_sam_attr_rgline: Vec<String>,
    pub out_sam_attr_rg: Vec<String>,
    pub out_sam_strand_field_type: i32,
    pub out_sam_mapq_unique: i32,
    pub out_sam_tlen: i32,
    pub out_sam_flag_or: u16,
    pub out_sam_flag_and: u16,
    pub out_sam_attr_ih_start: u32,
    pub out_sam_mode: String,
    pub out_sam_order: String,
    pub out_sam_read_id: String,
    pub out_sam_read_id_number: bool,
    pub out_sam_unmapped_keep_pairs: bool,
    pub out_sam_unmapped_within: bool,
    pub out_sam_filter_yes: bool,
    pub out_sam_filter_keep_only_added_references: bool,
    pub out_sam_filter_keep_all_added_references: bool,
    pub out_sam_mult_nmax: u64,
    pub out_filter_multimap_score_range: i32,
    pub out_filter_multimap_nmax: u64,
    pub out_filter_mismatch_nmax: u32,
    pub out_filter_mismatch_nover_read_lmax: f64,
    pub out_filter_score_min: i32,
    pub out_filter_score_min_over_lread: f64,
    pub out_filter_match_nmin: u32,
    pub out_filter_match_nmin_over_lread: f64,
    pub seed_search_lmax: u32,
    pub seed_search_start_lmax: u32,
    pub seed_search_start_lmax_over_lread: f64,
    pub seed_split_min: u32,
    pub seed_map_min: u32,
    pub out_multimapper_order_random: bool,
    pub out_sam_mult_nmax_is_limited: bool,
    pub out_sam_primary_flag: String,
    pub align_ends_protrude_concordant_pair: bool,
    pub align_ends_type: AlignEndsType,
    pub align_ends_protrude: AlignEndsProtrude,
    pub align_insertion_flush: AlignInsertionFlush,
    pub align_soft_clip_at_reference_ends_yes: bool,
    pub align_sj_overhang_min: u32,
    pub align_sjdb_overhang_min: u32,
    pub align_spliced_mate_map_lmin: u32,
    pub align_spliced_mate_map_lmin_over_lmate: f64,
    pub align_intron_max: u32,
    pub align_intron_min: u32,
    pub align_mates_gap_max: u32,
    pub align_sj_stitch_mismatch_nmax: Vec<i32>,
    pub out_filter_intron_strands: String,
    pub out_filter_intron_motifs: String,
    pub score_genomic_length_log2scale: f64,
    pub score_gap: i32,
    pub score_gap_noncan: i32,
    pub score_gap_gcag: i32,
    pub score_gap_atac: i32,
    pub score_del_base: i32,
    pub score_del_open: i32,
    pub score_ins_base: i32,
    pub score_ins_open: i32,
    pub score_stitch_sj_shift: i32,
    pub read_files_type: Vec<String>,
    pub read_files_type_n: i32,
    pub read_files_prefix: String,
    pub read_files_prefix_final: String,
    pub read_files_names: Vec<Vec<String>>,
    pub read_files_n: u32,
    pub read_files_manifest: Vec<String>,
    pub read_files_command: Vec<String>,
    pub read_files_command_string: String,
    pub read_files_sam_attr_keep_in: Vec<String>,
    pub read_files_sam_attr_keep_all: bool,
    pub read_files_sam_attr_keep_none: bool,
    pub read_files_sam_attr_keep: std::collections::BTreeSet<u16>,
    pub read_files_in_tmp: Vec<String>,
    pub read_files_index: u32,
    pub read_mates_lengths_in: String,
    pub sys_shell: String,
    pub read_name_separator_char: Vec<char>,
    pub out_qs_conversion_add: i32,
    pub var_yes: bool,
    pub wasp_yes: bool,
    pub quant_ge_count_yes: bool,
    pub quant_gene_yes: bool,
    pub quant_gene_full_yes: bool,
    pub quant_gene_full_exon_over_intron_yes: bool,
    pub quant_gene_full_ex50p_as_yes: bool,
    pub quant_yes: bool,
    pub quant_tr_sam_indel: bool,
    pub quant_tr_sam_single_end: bool,
    pub quant_tr_sam_soft_clip: bool,
    pub pe_overlap_nbases_min: u32,
    pub pe_overlap_mmp: f64,
    pub out_filter_mismatch_nover_lmax: f64,
    pub seed_multimap_nmax: u32,
    pub seed_per_read_nmax: u32,
    pub win_bin_chr_nbits: u32,
    pub win_bin_nbits: u32,
    pub win_anchor_dist_nbins: u32,
    pub win_anchor_multimap_nmax: u32,
    pub win_flank_nbins: u32,
    pub win_bin_n: u32,
    pub align_windows_per_read_nmax: u32,
    pub genome_insert_chr_ind_first: u32,
    pub out_bam_sorting_thread_n: u32,
    pub out_bam_sorting_thread_nactual: u32,
    pub out_bam_sorting_bins_n: u32,
    pub out_bam_coord_nbins: u32,
    pub chunk_out_bam_size_bytes: u64,
    pub limit_out_sam_one_read_bytes: u64,
    pub limit_out_sj_one_read: u64,
    pub limit_sjdb_insert_nsj: u32,
    pub out_bam_sorting_bin_start: Vec<u64>,
    pub bam_sorting_log: String,
    pub limit_bam_sort_ram: u64,
    pub limit_genome_generate_ram: u64,
    pub limit_nreads_soft: u64,
    pub out_bam_sort_tmp_dir: String,
    pub out_bam_file_coord_name: String,
    pub chunk_in_size_bytes_array: u64,
    pub chunk_in_size_bytes: u64,
    pub i_read_all: u64,
    pub read_map_number: u64,
    pub align_transcripts_per_read_nmax: u32,
    pub align_transcripts_per_window_nmax: u32,
    pub max_nsplit: u32,
    pub seed_per_window_nmax: u32,
    pub seed_none_loci_per_window: u32,
    pub win_read_coverage_relative_min: f64,
    pub win_read_coverage_bases_min: u32,
    pub p_clip: ParametersClip,
    pub sam_header_extra: String,
    pub read_in_0: String,
    pub read_in_0_pos: usize,
    pub read_files_in: Vec<String>,
    pub read_in_open: Vec<bool>,
    pub read_files_command_pid: Vec<i32>,
    pub two_pass_yes: bool,
    pub two_pass_pass2: bool,
    pub two_pass_pass1sj_file: String,
    pub two_pass_dir: String,
    pub two_pass_pass1reads_n: u64,
    pub run_mode_in: Vec<String>,
    pub genome_num_to_nt: Vec<u8>,
    pub limit_out_sj_collapsed: u64,
    pub out_filter_by_sjout_stage: i32,
    pub out_filter_type: String,
    pub out_sj: bool,
    pub out_sj_filter_reads: String,
    pub out_reads_unmapped: String,
    pub out_sjfilter_count_unique_min: Vec<i32>,
    pub out_sjfilter_count_total_min: Vec<i32>,
    pub out_sjfilter_overhang_min: Vec<i32>,
    pub out_sjfilter_intron_max_vs_read_n: Vec<u32>,
    pub out_sjfilter_dist_to_other_sj_min: Vec<i32>,
    pub sj_novel_n: u32,
    pub sj_novel_start: Vec<u32>,
    pub sj_novel_end: Vec<u32>,
    pub sjdb_insert_yes: bool,
    pub sjdb_insert_out_dir: String,
    pub wasp_output_mode: String,
    pub wasp_sam_tag: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct OutputSjResult {
    pub sj_out_tab: String,
    pub sj_start_gap_tsv: String,
    pub log_main: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AlignEndsType {
    pub in_: String,
    pub ext: [[bool; 2]; 2],
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AlignEndsProtrude {
    pub in_: Vec<String>,
    pub n_bases_max: i32,
    pub concordant_pair: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AlignInsertionFlush {
    pub in_: String,
    pub flush_right: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct OpenReadsFilesResult {
    pub opened_inputs: Vec<String>,
    pub read_files_in_tmp: Vec<String>,
    pub reads_command_file_names: Vec<String>,
    pub reads_command_file_contents: Vec<String>,
    pub command_pids: Vec<i32>,
    pub log_main: String,
    pub sam_header_commands: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SjdbPrepareResult {
    pub gsj: Vec<u8>,
    pub sjdb_info_txt: String,
    pub sjdb_list_out_tab: String,
    pub log_main: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SjdbBuildIndexResult {
    pub log_main: String,
    pub n_gsj: u32,
    pub sj_new: u64,
    pub n_ind: u32,
    pub ind_array: Vec<[u64; 2]>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SjdbInsertJunctionsResult {
    pub log_main: String,
    pub sjdb_prepare: SjdbPrepareResult,
    pub sjdb_build_index: SjdbBuildIndexResult,
    pub gtf: Option<GtfTranscriptGeneSjOutput>,
    pub files_written: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GenomeTransformGenomeResult {
    pub log_main: String,
    pub transform_blocks_tsv: String,
    pub debug: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GenomeGenerateResult {
    pub log_main: String,
    pub log_stdout: String,
    pub files_written: Vec<String>,
    pub gtf: Option<GtfTranscriptGeneSjOutput>,
    pub super_transcriptome: Option<GtfSuperTranscriptOutput>,
    pub sjdb_insert: Option<SjdbInsertJunctionsResult>,
    pub n_genome_true: u64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GenomeLoadResult {
    pub log_stdout: String,
    pub log_main: String,
    pub load_and_exit: bool,
    pub removed_shared_memory: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TwoPassRunPass1Result {
    pub log_progress: String,
    pub log_stdout: String,
    pub log_main: String,
    pub pass1_parameters: Parameters,
    pub pass1_chunks: Vec<ReadAlignChunk>,
    pub output_sj: OutputSjResult,
    pub sjdb_insert: Option<SjdbInsertJunctionsResult>,
    pub log_final_out: String,
    pub killed_read_command_pids: Vec<i32>,
    pub reopened_reads: OpenReadsFilesResult,
    pub pass1_sj_file: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StarMainResult {
    pub exit_code: i32,
    pub usage: String,
    pub log_stdout: String,
    pub log_main: String,
    pub log_progress: String,
    pub log_final_out: String,
    pub parameters: Parameters,
    pub genome: Option<Genome>,
    pub genome_generate: Vec<GenomeGenerateResult>,
    pub sjdb_insert: Option<SjdbInsertJunctionsResult>,
    pub two_pass: Option<TwoPassRunPass1Result>,
    pub output_sj: Option<OutputSjResult>,
    pub bam_sort: Option<BamSortByCoordinateResult>,
    pub signal: Option<SignalFromBamResult>,
    pub processed_bam_output: Vec<u8>,
    pub transcriptome: Option<Transcriptome>,
    pub stats_all: Stats,
    pub read_chunks: Vec<ReadAlignChunk>,
    pub process_chunks: Vec<ReadAlignChunkProcessChunksResult>,
    pub killed_read_command_pids: Vec<i32>,
    pub removed_tmp: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct OutWigFlags {
    pub yes: bool,
    pub strand: bool,
    pub type_: i32,
    pub format: i32,
    pub norm: i32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SignalFromBamRecord {
    pub tid: i32,
    pub pos: u32,
    pub flag: u16,
    pub cigar: Vec<u32>,
    pub nh: Option<u32>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SignalFromBamResult {
    pub files: std::collections::BTreeMap<String, String>,
    pub n_unique: f64,
    pub n_multiple: f64,
}

#[doc = "Original class `ParametersChimeric` at STAR/source/ParametersChimeric.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParametersChimeric {
    pub segment_min: u32,
    pub junction_overhang_min: u32,
    pub segment_read_gap_max: u32,
    pub score_min: i32,
    pub score_drop_max: i32,
    pub score_separation: i32,
    pub score_junction_non_gtag: i32,
    pub main_segment_mult_nmax: u32,
    pub multimap_score_range: u32,
    pub multimap_nmax: u32,
    pub nonchim_score_drop_min: u32,
    pub out_type: Vec<String>,
    pub out_bam: bool,
    pub out_junctions: bool,
    pub out_sam_old: bool,
    pub out_bam_hard_clip: bool,
    pub out_junction_format: Vec<i32>,
    pub filter_string_in: Vec<String>,
    pub filter_genomic_n: bool,
    pub out_chim_sam_opened: bool,
    pub out_chim_sam_contents: String,
    pub out_chim_junction_opened: bool,
    pub out_chim_junction_contents: String,
    pub log_main: String,
}

#[doc = "Original class `ParametersClip` at STAR/source/ParametersClip.h:21."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParametersClip {
    pub adapter_type: Vec<String>,
    pub in_: [ReadClipInput; 2],
    pub read_nmates: u32,
    pub read_nends: u32,
}

#[doc = "Original class `ParametersGenome` at STAR/source/ParametersGenome.h:9."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParametersGenome {
    pub g_dir: String,
    pub transform: ParametersGenomeTransform,
    pub g_type_string: String,
    pub g_fasta_files: Vec<String>,
    pub g_chain_files: Vec<String>,
    pub g_load: String,
    pub g_chr_bin_nbits: u32,
    pub g_saindex_nbases: u32,
    pub g_sasparse_d: u32,
    pub g_suffix_length_max: u32,
    pub sjdb_overhang: u32,
    pub sjdb_file_chr_start_end: Vec<String>,
    pub sjdb_gtf_file: String,
    pub sjdb_gtf_chr_prefix: String,
    pub sjdb_gtf_feature_exon: String,
    pub sjdb_gtf_tag_exon_parent_transcript: String,
    pub sjdb_gtf_tag_exon_parent_gene: String,
    pub sjdb_gtf_tag_exon_parent_gene_name: Vec<String>,
    pub sjdb_gtf_tag_exon_parent_gene_type: Vec<String>,
    pub sjdb_insert_save: String,
    pub g_file_sizes: Vec<u64>,
    pub sjdb_score: i32,
    pub chr_set_mito_strings: Vec<String>,
    pub chr_set_mito: std::collections::BTreeSet<u64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParametersGenomeTransform {
    pub type_string: String,
    pub type_: i32,
    pub vcf_file: String,
    pub output: Vec<String>,
    pub out_yes: bool,
    pub out_sam: bool,
    pub out_sj: bool,
    pub out_quant: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SamAttrPresent {
    pub nh: bool,
    pub hi: bool,
    pub as_: bool,
    pub nm: bool,
    pub md: bool,
    pub n_m: bool,
    pub j_m: bool,
    pub j_i: bool,
    pub rg: bool,
    pub mc: bool,
    pub xs: bool,
    pub ch: bool,
    pub v_a: bool,
    pub v_g: bool,
    pub v_w: bool,
    pub r_b: bool,
    pub ha: bool,
    pub cr: bool,
    pub cy: bool,
    pub ur: bool,
    pub uy: bool,
    pub cb: bool,
    pub ub: bool,
    pub gx: bool,
    pub gn: bool,
    pub gx_lower: bool,
    pub gn_lower: bool,
    pub s_m: bool,
    pub s_s: bool,
    pub s_q: bool,
    pub s_f: bool,
    pub c_n: bool,
}

#[doc = "Original class `ParametersSolo` at STAR/source/ParametersSolo.h:14."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParametersSolo {
    pub yes: bool,
    pub type_str: String,
    pub solo_type: i32,
    pub barcode_read: u32,
    pub barcode_read_separate: bool,
    pub barcode_start: u32,
    pub barcode_end: u32,
    pub sam_attr_yes: bool,
    pub read_info_yes: Vec<bool>,
    pub read_index_yes: Vec<bool>,
    pub read_stats_yes: Vec<bool>,
    pub read_stats_type: String,
    pub cb_wl_yes: bool,
    pub cb_wl_size: u32,
    pub cb_wl: Vec<u64>,
    pub cb_wl_str: Vec<String>,
    pub solo_cb_whitelist: Vec<String>,
    pub cb_position_str: Vec<String>,
    pub umi_position_str: String,
    pub cb_v: Vec<SoloBarcode>,
    pub umi_v: SoloBarcode,
    pub adapter_yes: bool,
    pub adapter_seq: String,
    pub adapter_mismatches_nmax: u32,
    pub cb_l: u32,
    pub cb_s: u32,
    pub umi_s: u32,
    pub b_l: u32,
    pub cbumi_l: u32,
    pub cb_type_type: i32,
    pub cb_type_str_map: std::collections::BTreeMap<String, u32>,
    pub sam_attr_barcode_seq: Vec<String>,
    pub sam_attr_barcode_qual: Vec<String>,
    pub features: Vec<u32>,
    pub feature_yes: Vec<bool>,
    pub feature_ind: Vec<i32>,
    pub feature_first: i32,
    pub n_features: u32,
    pub out_file_names: Vec<String>,
    pub out_format_features_gene_field3: String,
    pub umi_l: u32,
    pub umi_mask_low: u32,
    pub umi_dedup: UMIdedup,
    pub umi_filtering: SoloUmiFiltering,
    pub multi_map: MultiMappers,
    pub cb_match_wl: CBMatchWL,
    pub sam_attr_feature: i32,
    pub qs_base: i8,
    pub qs_max: i8,
    pub cb_min_p: f64,
    pub redistr_reads_nfiles: u32,
    pub cluster_cb_file: String,
    pub cell_filter: SoloCellFilter,
    pub strand: i32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloCellFilter {
    pub type_: Vec<String>,
    pub top_cells: u32,
    pub knee: SoloCellFilterKnee,
    pub ed_cr: SoloCellFilterEmptyDropsCr,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloCellFilterKnee {
    pub n_expected_cells: f64,
    pub max_percentile: f64,
    pub max_min_ratio: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloCellFilterEmptyDropsCr {
    pub ind_min: u32,
    pub ind_max: u32,
    pub umi_min: u32,
    pub umi_min_frac_median: f64,
    pub cand_max_n: u32,
    pub fdr: f64,
    pub sim_n: u32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CBMatchWL {
    pub type_: String,
    pub mm1: bool,
    pub mm1_multi: bool,
    pub one_exact: bool,
    pub mm1_multi_pc: bool,
    pub mm1_multi_nbase: bool,
    pub edit_dist_2: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloUmiFiltering {
    pub type_: Vec<String>,
    pub yes: bool,
    pub multi_gene_umi: bool,
    pub multi_gene_umi_all: bool,
    pub multi_gene_umi_cr: bool,
}

#[doc = "Original class `Quantifications` at STAR/source/Quantifications.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Quantifications {
    pub gene_counts: QuantificationGeneCounts,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct QuantificationGeneCounts {
    pub n_ge: u32,
    pub n_type: i32,
    pub c_multi: u64,
    pub c_ambig: Vec<u64>,
    pub c_none: Vec<u64>,
    pub g_count: Vec<Vec<u64>>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct QuantTranscriptomeBamRequest {
    pub transcript: Transcript,
    pub n_align_t: u32,
    pub i_align_t: u32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct QuantTranscriptomeResult {
    pub n_align_t: u32,
    pub align_t: Vec<Transcript>,
    pub bam_requests: Vec<QuantTranscriptomeBamRequest>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericBamOutputRequest {
    pub al1: Transcript,
    pub al2: Transcript,
    pub i_tr: u32,
    pub chim_n: u32,
    pub is_best_chim_align: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericDetectionOldOutputResult {
    pub chim_n: u32,
    pub chim_sam: String,
    pub chim_junction: String,
    pub bam_requests: Vec<ChimericBamOutputRequest>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignGenomeTransformResult {
    pub al_best: Transcript,
    pub al_mult: Vec<Transcript>,
    pub al_n: u32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AlignBamRequest {
    pub transcript: Transcript,
    pub n_tr_out: u32,
    pub i_tr_out: u32,
    pub tr_chr_start: u64,
    pub mate_chr: u32,
    pub mate_start: u32,
    pub mate_strand: i8,
    pub align_type: i32,
    pub mate_map: Option<[bool; 2]>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AlignBamResult {
    pub n_lines: u32,
    pub records: [Vec<u8>; 2],
    pub record_sizes: [u32; 2],
    pub sam_flags: [u16; 2],
    pub mapq: [i32; 2],
    pub n_cigar: [u32; 2],
    pub signal_records: [Option<SignalFromBamRecord>; 2],
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericBamOutputResult {
    pub chim_represent: i32,
    pub chim_type: i32,
    pub representative_request_index: i32,
    pub supplementary_request_index: i32,
    pub bam_requests: Vec<AlignBamRequest>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericDetectionRequest {
    pub detector: String,
    pub n_w: u32,
    pub read_length: Vec<u32>,
    pub max_non_chim_align_score: i32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericDetectionResult {
    pub request: Option<ChimericDetectionRequest>,
    pub old_output_requested: bool,
    pub old_output: Option<ChimericDetectionOldOutputResult>,
    pub mult_output: Option<ChimericDetectionMultResult>,
    pub mult_map_select_requested: bool,
    pub mapped_filter_requested: bool,
    pub pe_tr_chim: Vec<Transcript>,
    pub chim_record: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChimericDetectionMultResult {
    pub chim_record: bool,
    pub chim_n: u32,
    pub chim_score_best: i32,
    pub min_score_to_consider: i32,
    pub chim_junction: String,
    pub bam_outputs: Vec<ChimericBamOutputResult>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeOverlapMergeMapResult {
    pub map_one_read_requested: bool,
    pub chimeric_detection: Option<ChimericDetectionResult>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WriteSamResult {
    pub out_bam_bytes: u64,
    pub unmap_type: i32,
    pub mate_mapped: [bool; 2],
    pub sam: String,
    pub bam_requests: Vec<AlignBamRequest>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct OutputAlignmentsResult {
    pub out_bam_bytes: u64,
    pub unmap_type: i32,
    pub out_filter_by_sjout_pass: bool,
    pub mate_mapped: [bool; 2],
    pub sam: String,
    pub write_sam: WriteSamResult,
    pub quant_transcriptome: Option<QuantTranscriptomeResult>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignOneReadResult {
    pub status: i32,
    pub map_one_read_requested: bool,
    pub splice_graph_log: String,
    pub pe_overlap: PeOverlapMergeMapResult,
    pub chimeric_detection: Option<ChimericDetectionResult>,
    pub wasp: WaspMapResult,
    pub output_alignments: Option<OutputAlignmentsResult>,
    pub tr_mult: Vec<Transcript>,
    pub aligns_gen_out: ReadAlignGenomeTransformResult,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WaspMapRequest {
    pub alleles: Vec<u8>,
    pub read1: [Vec<u8>; 3],
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WaspMapOutcome {
    pub unmap_type: i32,
    pub n_tr: u64,
    pub align: Transcript,
    pub transformed: Option<ReadAlignGenomeTransformResult>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WaspMapResult {
    pub wasp_type: i32,
    pub requests: Vec<WaspMapRequest>,
}

#[doc = "Original class `ReadAlign` at STAR/source/ReadAlign.h:22."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlign {
    pub chim_record: bool,
    pub wasp_type: i32,
    pub unmap_type: i32,
    pub i_read: u64,
    pub stats_ra: Stats,
    pub tr_best: Transcript,
    pub l_read: u32,
    pub out_filter_mismatch_nmax_total: u32,
    pub read_nmates: u32,
    pub map_marker: u32,
    pub n_a: u32,
    pub n_p: u32,
    pub n_split: u32,
    pub n_w: u32,
    pub n_wall: u32,
    pub n_tr: u64,
    pub revert_strand: bool,
    pub tr_mult: Vec<Transcript>,
    pub n_win_tr: Vec<u64>,
    pub tr_all: Vec<Vec<Transcript>>,
    pub n_um: [u32; 2],
    pub stored_lmin: u32,
    pub uniq_lmax: u32,
    pub uniq_lmax_ind: u32,
    pub mult_lmax: u32,
    pub mult_lmax_n: u32,
    pub mult_nmin_l: u32,
    pub mult_nmin: u32,
    pub mult_nmax: u32,
    pub mult_nmax_l: u32,
    pub chim_n: u32,
    pub chim_str: u32,
    pub chim_repeat0: u32,
    pub chim_repeat1: u32,
    pub chim_j0: u32,
    pub chim_j1: u32,
    pub chim_motif: i32,
    pub tr_chim: Vec<Transcript>,
    pub max_score_mate: Vec<i32>,
    pub read_length: Vec<u32>,
    pub read_length_original: Vec<u32>,
    pub read_length_pair_original: u32,
    pub read_name: String,
    pub i_read_all: u64,
    pub read_filter: i32,
    pub read_files_index: u32,
    pub read_file_type: i32,
    pub out_bam_bytes: u64,
    pub rng_mult_order_seed: u64,
    pub rng_uniform_real_0_to_1: [f64; 2],
    pub align_tr_all: Vec<Transcript>,
    pub split_r: [Vec<u32>; 3],
    pub tr_array: Vec<Transcript>,
    pub tr_array_pointer_n: usize,
    pub tr_init: Box<Transcript>,
    pub aligns_gen_out_al_mult: Vec<Transcript>,
    pub read0: Vec<Vec<u8>>,
    pub qual0: Vec<Vec<u8>>,
    pub read_name_mates: Vec<Vec<u8>>,
    pub read_name_extra: Vec<String>,
    pub out_bam_one_align_nbytes: Vec<u32>,
    pub out_bam_one_align: Vec<Vec<u8>>,
    pub chunk_out_chim_junction_opened: bool,
    pub solo_read: SoloRead,
    pub spl_graph_present: bool,
    pub qual_hist: Vec<Vec<u64>>,
    pub pe_ov: ReadAlignPeOverlap,
    pub read1: [Vec<u8>; 3],
    pub clip_mates: Vec<Vec<ClipMate>>,
    pub mates_cigar: Vec<String>,
    pub win_bin: [Vec<u32>; 2],
    pub wc: Vec<[u32; crate::generated::functions::WC_SIZE]>,
    pub wa: Vec<Vec<[u32; crate::generated::functions::WA_SIZE]>>,
    pub n_wa: Vec<u32>,
    pub n_wap: Vec<u32>,
    pub wal_rec: Vec<u32>,
    pub w_last_anchor: Vec<u32>,
    pub wa_incl: Vec<bool>,
    pub score_seed_best: Vec<i32>,
    pub score_seed_best_mm: Vec<u32>,
    pub score_seed_best_ind: Vec<u32>,
    pub seed_chain: Vec<u32>,
    pub pc: Vec<[u32; crate::generated::functions::PC_SIZE]>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignPeOverlap {
    pub yes: bool,
    pub n_ov: u32,
    pub ov_s: u32,
    pub mate_start: [u32; 2],
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StoredAlign {
    pub i_dir: u32,
    pub shift: u32,
    pub n_rep: u32,
    pub l: u32,
    pub ind_start_end: [u32; 2],
    pub i_frag: u32,
}

#[doc = "Original class `ReadAlignChunk` at STAR/source/ReadAlignChunk.h:12."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignChunk {
    pub i_thread: i32,
    pub i_chunk_in: u32,
    pub no_reads_left: bool,
    pub chunk_tr: Option<Transcriptome>,
    pub ra: ReadAlign,
    pub chunk_in: Vec<Vec<u8>>,
    pub chunk_in_size_bytes_total: Vec<u64>,
    pub read_in_stream_n: usize,
    pub chunk_out_bam: Vec<u8>,
    pub chunk_out_bam_total: u64,
    pub chunk_out_bam_file_name: String,
    pub chunk_out_bam_unsorted: Option<BAMoutput>,
    pub chunk_out_bam_coord: BAMoutput,
    pub chunk_out_bam_quant: Option<BAMoutput>,
    pub chunk_out_sj: OutSJ,
    pub chunk_out_sj1: OutSJ,
    pub chunk_out_chim_sam_path: Option<String>,
    pub chunk_out_chim_junction_path: Option<String>,
    pub chunk_out_unmapped_reads_paths: Vec<String>,
    pub chunk_out_filter_by_sjout_files: Vec<String>,
    pub wasp_ra_present: bool,
    pub pe_merge_ra_present: bool,
    pub log_main: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignChunkMapChunkResult {
    pub direct_sam_output: Vec<u8>,
    pub paired_keep_input_order_tmp: Vec<u8>,
    pub chimeric_sam_output: String,
    pub chimeric_junction_output: String,
    pub unmapped_fastx_outputs: Vec<String>,
    pub signal_records: Vec<SignalFromBamRecord>,
    pub quant_bam_output: Vec<u8>,
    pub paired_keep_input_order_tmp_name: Option<String>,
    pub paired_keep_input_order_final_name: Option<String>,
    pub progress_report: Option<String>,
    pub log_main: String,
    pub reads_processed: u64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAlignChunkProcessChunksResult {
    pub map_chunks: Vec<ReadAlignChunkMapChunkResult>,
    pub chunk_inputs: Vec<Vec<Vec<u8>>>,
    pub log_main: String,
    pub chunks_read: u32,
    pub paired_keep_input_order_cat_after_chunks: Vec<u32>,
    pub flushed_bam_unsorted: bool,
    pub flushed_bam_coord: bool,
    pub flushed_bam_quant: bool,
    pub chim_sam_cat_path: Option<String>,
    pub chim_junction_cat_path: Option<String>,
    pub unmapped_fastx_cat_paths: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BamSortByCoordinateResult {
    pub output_bam: Vec<u8>,
    pub bin_outputs: Vec<Vec<u8>>,
    pub bin_names: Vec<String>,
    pub removed_files: Vec<String>,
    pub max_mem: u64,
    pub unmapped_reads_n: u64,
}

#[doc = "Original class `ReadAnnotFeature` at STAR/source/ReadAnnotations.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAnnotFeature {
    pub f_set: std::collections::BTreeSet<u32>,
    pub f_align: Vec<std::collections::BTreeSet<u32>>,
    pub ov_type: u32,
}

#[doc = "Original class `ReadAnnotations` at STAR/source/ReadAnnotations.h:20."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadAnnotations {
    pub annot_features: Vec<ReadAnnotFeature>,
    pub gene_exon_overlap: Vec<i32>,
    pub transcript_concordant: Vec<[u32; 2]>,
    pub gene_velocyto_simple: [u32; 2],
    pub tr_velocyto_type: Vec<TrTypeStruct>,
}

#[doc = "Original class `ReadClipInput` at STAR/source/ParametersClip.h:12."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadClipInput {
    pub n: Vec<u32>,
    pub n_after_ad: Vec<u32>,
    pub ad_seq: Vec<String>,
    pub ad_mmp: Vec<f64>,
}

#[doc = "Original class `ReadSoloFeatures` at STAR/source/SoloReadFeature_record.cpp:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadSoloFeatures {
    pub gene: u32,
    pub gene_mult: Vec<u32>,
    pub sj: Vec<[u64; 2]>,
    pub sj_annot: bool,
    pub ind_annot_tr: usize,
    pub align_out: Vec<Transcript>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TrTypeStruct {
    pub tr: u32,
    pub type_: u8,
}

#[doc = "Original class `SNP` at STAR/source/Variation.h:15."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SNP {
    pub n: u32,
    pub loci: Vec<u32>,
    pub loci_v: Vec<u32>,
    pub nt: Vec<[u8; 3]>,
}

#[doc = "Original class `SharedMemory` at STAR/source/SharedMemory.h:85."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SharedMemory {
    pub key: i32,
    pub counter_key: i32,
    pub unload_last: bool,
    pub shm_id: i32,
    pub shared_counter_id: i32,
    pub counter_mem_attached: bool,
    pub mapped: bool,
    pub length: usize,
    pub is_allocator: bool,
    pub needs_allocation: bool,
    pub exception: SharedMemoryException,
    pub shared_objects_use_count_value: i32,
    pub clean_count: u32,
}

#[doc = "Original class `SharedMemoryException` at STAR/source/SharedMemory.h:33."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SharedMemoryException {
    pub error_code: i32,
    pub error_detail: i32,
}

#[doc = "Original class `Simd` at STAR/source/opal/opal.cpp:551."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Simd {}

#[doc = "Original struct `Simd<char>` at STAR/source/opal/opal.cpp:554."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Simd_char_ {
    pub numseqs: i64,
    pub satarthm: bool,
}

#[doc = "Original struct `Simd<int>` at STAR/source/opal/opal.cpp:580."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Simd_int_ {
    pub numseqs: i64,
    pub satarthm: bool,
}

#[doc = "Original struct `Simd<short>` at STAR/source/opal/opal.cpp:567."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Simd_short_ {
    pub numseqs: i64,
    pub satarthm: bool,
}

#[doc = "Original struct `SimdSW` at STAR/source/opal/opal.cpp:95."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SimdSW {}

#[doc = "Original struct `SimdSW<char>` at STAR/source/opal/opal.cpp:98."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SimdSW_char_ {
    pub numseqs: i64,
    pub satarthm: bool,
    pub negrange: bool,
}

#[doc = "Original struct `SimdSW<int>` at STAR/source/opal/opal.cpp:126."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SimdSW_int_ {
    pub numseqs: i64,
    pub satarthm: bool,
    pub negrange: bool,
}

#[doc = "Original struct `SimdSW<short>` at STAR/source/opal/opal.cpp:112."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SimdSW_short_ {
    pub numseqs: i64,
    pub satarthm: bool,
    pub negrange: bool,
}

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

#[doc = "Original class `Solo` at STAR/source/Solo.h:11."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Solo {
    pub p_solo: ParametersSolo,
    pub solo_feat: Vec<SoloFeature>,
    pub read_bar_sum: Option<SoloReadBarcode>,
}

#[doc = "Original class `SoloBarcode` at STAR/source/SoloBarcode.h:9."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloBarcode {
    pub anchor_type: [i32; 2],
    pub anchor_dist: [i32; 2],
    pub adapter_length: i32,
    pub wl: Vec<Vec<u64>>,
    pub wl_ed: Vec<Vec<u64>>,
    pub wl_ed_ind: Vec<Vec<u32>>,
    pub wl_factor: u64,
    pub wl_add: Vec<u32>,
    pub min_len: u32,
    pub total_size: u32,
    pub i_cb: u32,
    pub i_len: u32,
}

#[doc = "Original class `SoloFeature` at STAR/source/SoloFeature.h:18."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloFeature {
    pub feature_type: i32,
    pub features_number: i64,
    pub read_bar_sum: Option<SoloReadBarcode>,
    pub read_feat_sum: Option<SoloReadFeature>,
    pub read_feat_all: Vec<SoloReadFeature>,
    pub read_feat_all_len: usize,
    pub n_reads_mapped: u64,
    pub n_reads_input: u64,
    pub n_cb: u32,
    pub rgu_stride: u32,
    pub r_cbp: Vec<Vec<u32>>,
    pub cb_feature_umi_map: Vec<()>,
    pub count_cell_gene_umi: Vec<u32>,
    pub count_cell_gene_umi_index: Vec<u32>,
    pub count_mat_mult_i: Vec<u32>,
    pub count_mat_mult_m: Vec<f64>,
    pub count_mat_mult_s: u32,
    pub count_mat_stride: u32,
    pub output_prefix: String,
    pub output_prefix_filtered: String,
    pub copied_features_tsv: String,
    pub ind_cb_wl: Vec<u32>,
    pub ind_cb: Vec<u32>,
    pub n_gene_per_cb: Vec<u32>,
    pub n_gene_per_cb_multi: Vec<u32>,
    pub n_read_per_cb: Vec<u32>,
    pub n_read_per_cb_total: Vec<u32>,
    pub n_read_per_cb_unique: Vec<u32>,
    pub n_umi_per_cb: Vec<u32>,
    pub n_umi_per_cb_multi: Vec<u32>,
    pub n_umi_per_cb_sorted: Vec<u32>,
    pub read_info: Vec<SoloFeatureReadInfo>,
    pub read_flag_counts: SoloReadFlagClass,
    pub sj_all: [Vec<u64>; 2],
    pub filtered_cells: SoloFilteredCells,
    pub redistr_files_cb_first: Vec<u32>,
    pub redistr_files_cb_index: Vec<u32>,
    pub redistr_files_nreads: Vec<u64>,
    pub redistr_files_streams: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloFeatureReadInfo {
    pub cb: i64,
    pub umi: u64,
}

#[doc = "Original class `SoloFilteredCells` at STAR/source/SoloFilteredCells.h:4."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloFilteredCells {
    pub filt_vec_bool: Vec<bool>,
    pub n_cells: u32,
    pub n_read_in_cells: u64,
    pub median_read_per_cell_unique: u64,
    pub mean_read_per_cell_unique: u64,
    pub n_umi_in_cells: u64,
    pub median_umi_per_cell: u64,
    pub mean_umi_per_cell: u64,
    pub n_gene_in_cells: u64,
    pub median_gene_per_cell: u64,
    pub mean_gene_per_cell: u64,
    pub n_gene_detected: u64,
    pub n_cells_simple: u64,
    pub n_read_in_cells_unique: u64,
    pub n_read_per_cell_unique: Vec<u32>,
    pub n_gene_per_cell: Vec<u32>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloOutputResults {
    pub created_directory: String,
    pub files: std::collections::BTreeMap<String, String>,
    pub symlinks: Vec<(String, String)>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloCellFilteringResult {
    pub log_main: String,
    pub empty_drops_requested: bool,
    pub empty_drops: Option<SoloEmptyDropsCrResult>,
    pub output_results: Option<SoloOutputResults>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloEmptyDropsCrResult {
    pub log_main: String,
    pub early_return: Option<String>,
    pub feat_det_n: u32,
    pub min_umi: u32,
    pub candidate_first: u32,
    pub candidate_last: u32,
    pub extra_cells: u32,
    pub p_values: Vec<(u32, f64, f64)>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloFeatureProcessRecordsResult {
    pub created_directories: Vec<String>,
    pub files: std::collections::BTreeMap<String, String>,
    pub symlinks: Vec<(String, String)>,
    pub log_main: String,
    pub stats_output: Option<SoloStatsOutput>,
    pub cell_filtering: Option<SoloCellFilteringResult>,
    pub returned_after_quant_transcript: bool,
    pub count_cb_gene_umi_called: bool,
    pub quant_transcript_called: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloProcessAndOutputResult {
    pub files: std::collections::BTreeMap<String, String>,
    pub symlinks: Vec<(String, String)>,
    pub created_directories: Vec<String>,
    pub log_stdout: String,
    pub log_main: String,
    pub feature_results: Vec<SoloFeatureProcessRecordsResult>,
    pub returned_after_barcode_output: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloFeatureQuantTranscriptResult {
    pub files: std::collections::BTreeMap<String, String>,
    pub log_main: String,
    pub stdout: String,
    pub cluster_expression: std::collections::BTreeMap<u32, Vec<f64>>,
    pub returned_no_cluster_file: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloConstructorCellFilteringResult {
    pub solo: Solo,
    pub cell_filtering: Option<SoloCellFilteringResult>,
    pub log_stdout: String,
    pub log_main: String,
    pub exited: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloStatsOutput {
    pub files: std::collections::BTreeMap<String, String>,
}

#[doc = "Original class `SoloRead` at STAR/source/SoloRead.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloRead {
    pub i_chunk: i32,
    pub read_bar: Option<SoloReadBarcode>,
    pub read_feat: Vec<SoloReadFeature>,
}

#[doc = "Original class `SoloReadBarcode` at STAR/source/SoloReadBarcode.h:8."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloReadBarcode {
    pub solo_type: i32,
    pub cb_wl_yes: bool,
    pub cb_wl_size: u32,
    pub umi_l: u32,
    pub homo_polymer: [u64; 4],
    pub umi_seq: String,
    pub cb_seq: String,
    pub cb_qual: String,
    pub umi_qual: String,
    pub b_seq: String,
    pub b_qual: String,
    pub b_strings: Vec<String>,
    pub cb_seq_corrected: String,
    pub umi_b: u64,
    pub umi_check: i32,
    pub cb_match: i32,
    pub cb_match_string: String,
    pub cb_match_ind: Vec<u64>,
    pub cb_read_count_exact: Vec<u32>,
    pub qual_hist: Vec<u64>,
    pub stats: SoloReadBarcodeStats,
}

#[doc = "Original class `SoloReadBarcodeStats` at STAR/source/SoloReadBarcodeStats.h:5."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloReadBarcodeStats {
    pub names: Vec<String>,
    pub v: Vec<u64>,
}

#[doc = "Original class `SoloReadFeature` at STAR/source/SoloReadFeature.h:15."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloReadFeature {
    pub feature_type: i32,
    pub read_info_yes: bool,
    pub read_index_yes: bool,
    pub stream_reads_path: Option<String>,
    pub stream_reads: String,
    pub cb_wl_yes: bool,
    pub cb_wl_size: u32,
    pub cb_read_count: Vec<u32>,
    pub cb_read_count_map: std::collections::BTreeMap<u64, u32>,
    pub transcript_dist_count: Vec<u32>,
    pub read_flag: SoloReadFlagClass,
    pub stats: SoloReadFeatureStats,
}

#[doc = "Original class `SoloReadFeatureStats` at STAR/source/SoloReadFeatureStats.h:5."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloReadFeatureStats {
    pub names: Vec<String>,
    pub v: Vec<u64>,
}

#[doc = "Original class `SoloReadFlagClass` at STAR/source/SoloCommon.h:26."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SoloReadFlagClass {
    pub flag: u32,
    pub flag_counts:
        std::collections::BTreeMap<u64, [u64; crate::generated::functions::SOLO_READ_FLAG_N_BITS]>,
    pub flag_counts_no_cb: [u64; crate::generated::functions::SOLO_READ_FLAG_N_BITS],
}

#[doc = "Original class `SpliceGraph` at STAR/source/SpliceGraph.h:12."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpliceGraph {
    pub super_trome: SuperTranscriptome,
    pub ra: Option<ReadAlign>,
    pub super_tr_seed_count: Vec<u16>,
    pub scoring_matrix: Vec<Vec<i32>>,
    pub score_two_columns: [Vec<i32>; 2],
    pub direction_matrix: Vec<u8>,
    pub sj_dindex: Vec<u32>,
    pub gap_penalty: i8,
    pub match_score: i8,
    pub mismatch_penalty: i8,
    pub align_info: SpliceGraphAlignInfo,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpliceGraphAlignInfo {
    pub n_map: u32,
    pub n_mm: u32,
    pub n_i: u32,
    pub n_d: u32,
    pub n_sj: u32,
    pub a_start: [u32; 2],
    pub a_end: [u32; 2],
}

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
    pub splices_n: [u32; crate::generated::functions::SJ_MOTIF_SIZE],
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

#[doc = "Original class `SuperTranscript` at STAR/source/SuperTranscriptome.h:14."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SuperTranscript {
    pub seq: Vec<u8>,
    pub length: u32,
    pub sj_c: Vec<[u32; 3]>,
    pub sj_donor: Vec<u32>,
}

#[doc = "Original class `SuperTranscriptome` at STAR/source/SuperTranscriptome.h:23."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SuperTranscriptome {
    pub sj_nmax: u32,
    pub sj_donor_nmax: u32,
    pub n: u32,
    pub seq_concat: Vec<u8>,
    pub seq: Vec<Vec<u8>>,
    pub tr_index: Vec<u64>,
    pub tr_start_end: Vec<[u64; 2]>,
    pub sj: Vec<sjInfo>,
    pub super_trs: Vec<SuperTranscript>,
}

#[doc = "Original class `ThreadControl` at STAR/source/ThreadControl.h:9."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ThreadControl {
    pub chunk_in_n: u32,
    pub chunk_out_n: u32,
}

#[doc = "Original class `Transcript` at STAR/source/Transcript.h:10."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Transcript {
    pub exons: Vec<[u32; crate::generated::functions::EX_SIZE]>,
    pub cigar: Vec<[u32; 2]>,
    pub canon_sj: Vec<i32>,
    pub sj_annot: Vec<u8>,
    pub shift_sj: Vec<[u32; 2]>,
    pub sj_str: Vec<u8>,
    pub sj_yes: bool,
    pub intron_motifs: [u32; 3],
    pub sj_motif_strand: u32,
    pub n_exons: u32,
    pub l_read: u32,
    pub read_length: Vec<u32>,
    pub read_nmates: u32,
    pub read_length_original: Vec<u32>,
    pub read_length_pair_original: u32,
    pub c_start: u32,
    pub r_start: u32,
    pub ro_start: u32,
    pub ro_str: u32,
    pub r_length: u32,
    pub mapped_length: u32,
    pub g_start: u32,
    pub g_length: u32,
    pub chr: u32,
    pub str_: u32,
    pub i_frag: i32,
    pub primary_flag: bool,
    pub max_score: i32,
    pub n_match: u32,
    pub n_mm: u32,
    pub n_gap: u32,
    pub l_gap: u32,
    pub l_del: u32,
    pub l_ins: u32,
    pub n_del: u32,
    pub n_ins: u32,
    pub n_unique: u32,
    pub n_anchor: u32,
    pub extend_l: u32,
    pub haplo_type: u32,
    pub var_ind: Vec<i32>,
    pub var_gen_coord: Vec<u32>,
    pub var_read_coord: Vec<u32>,
    pub var_allele: Vec<u8>,
    pub i_read: u64,
    pub read_name: String,
}

#[doc = "Original class `Transcriptome` at STAR/source/Transcriptome.h:13."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Transcriptome {
    pub n_tr: u32,
    pub n_ge: u32,
    pub ge_id: Vec<String>,
    pub ge_name: Vec<String>,
    pub ge_biotype: Vec<String>,
    pub tr_id: Vec<String>,
    pub tr_s: Vec<u32>,
    pub tr_e: Vec<u32>,
    pub tr_e_max: Vec<u32>,
    pub tr_ex_n: Vec<u16>,
    pub tr_ex_i: Vec<u32>,
    pub tr_str: Vec<u8>,
    pub tr_gene: Vec<u32>,
    pub tr_len: Vec<u32>,
    pub ex_se: Vec<u32>,
    pub ex_len_cum: Vec<u32>,
    pub ex_g: TranscriptomeExG,
    pub gene_full: TranscriptomeGeneFull,
    pub quants: Quantifications,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TranscriptomeExG {
    pub n_ex: u64,
    pub s: Vec<u64>,
    pub e: Vec<u64>,
    pub e_max: Vec<u64>,
    pub str_: Vec<u8>,
    pub g: Vec<u32>,
    pub t: Vec<u32>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TranscriptomeGeneFull {
    pub s: Vec<u64>,
    pub e: Vec<u64>,
    pub e_max: Vec<u64>,
    pub str_: Vec<u8>,
    pub g: Vec<u32>,
}

#[doc = "Original class `UMIdedup` at STAR/source/ParametersSolo.h:16."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct UMIdedup {
    pub types_in: Vec<String>,
    pub types: Vec<i32>,
    pub type_main: i32,
    pub yes_n: u32,
    pub yes_b: [bool; 6],
    pub count_ind_i: [u32; 6],
    pub count_ind_main: u32,
}

#[doc = "Original class `Variation` at STAR/source/Variation.h:30."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Variation {
    pub yes: bool,
    pub snp: SNP,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct VariantInfo {
    pub pos: u64,
    pub len: i32,
    pub seq: [String; 2],
}

#[doc = "Original struct `simde_mm_loadu_si128_s` at STAR/source/opal/simde_avx2.h:18343."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct simde_mm_loadu_si128_s {}

#[doc = "Original struct `sjInfo` at STAR/source/SuperTranscriptome.h:7."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct sjInfo {
    pub start: u32,
    pub end: u32,
    pub tr: u32,
    pub super_: u32,
}

#[doc = "Original struct `type_cbMMind` at STAR/source/SoloBarcode.cpp:49."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct type_cbMMind {
    pub cb: i64,
    pub ind: u32,
    pub mm: u32,
}
