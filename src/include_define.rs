#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

pub const OPAL_MODE_NW: i32 = 0;

pub const OPAL_MODE_HW: i32 = 1;

pub const OPAL_MODE_OV: i32 = 2;

pub const OPAL_MODE_SW: i32 = 3;

pub const OPAL_ERR_OVERFLOW: i32 = 1;

pub const OPAL_ERR_NO_SIMD_SUPPORT: i32 = 2;

pub const OPAL_ERR_INVALID_MODE: i32 = 3;

pub const OPAL_OVERFLOW_SIMPLE: i32 = 0;

pub const OPAL_OVERFLOW_BUCKETS: i32 = 1;

pub const OPAL_SEARCH_SCORE: i32 = 0;

pub const OPAL_SEARCH_SCORE_END: i32 = 1;

pub const OPAL_SEARCH_ALIGNMENT: i32 = 2;

pub const OPAL_ALIGN_MATCH: u8 = 0;

pub const OPAL_ALIGN_DEL: u8 = 1;

pub const OPAL_ALIGN_INS: u8 = 2;

pub const OPAL_ALIGN_MISMATCH: u8 = 3;

pub const MARK_FRAG_SPACER_BASE: u8 = 11;

pub const MARKER_READ_TOO_SHORT: u64 = 999_910;

pub const MARKER_NO_GOOD_WINDOW: u64 = 999_902;

pub const MARKER_NO_GOOD_PIECES: u64 = 999_909;

pub const MARKER_ALL_PIECES_EXCEED_SEED_MULTIMAP_NMAX: u64 = 999_908;

pub const EX_R: usize = 0;

pub const EX_G: usize = 1;

pub const EX_L: usize = 2;

pub const EX_IFRAG: usize = 3;

pub const EX_SJA: usize = 4;

pub const EX_SIZE: usize = 5;

pub const GENOME_EX_T: usize = 0;

pub const GENOME_EX_S: usize = 1;

pub const GENOME_EX_E: usize = 2;

pub const GENOME_EX_G: usize = 3;

pub const GENOME_EX_L: usize = 4;

pub const GTF_EX_T: usize = 0;

pub const GTF_EX_S: usize = 1;

pub const GTF_EX_E: usize = 2;

pub const GTF_EX_G: usize = 3;

pub const GTF_EX_L: usize = 4;

/// Per-read max exon count. C++ STAR uses 20 for short reads, 1000 for
/// `COMPILE_FOR_LONG_READS`. We match the short-read default so the
/// fixed-size arrays inside `Transcript` stay compact.
pub const MAX_N_EXONS: usize = 20;
pub const MAX_N_MATES: usize = 3;

pub const MAX_SJ_REPEAT_SEARCH: u32 = 255;

pub const WA_LENGTH: usize = 0;

pub const WA_R_START: usize = 1;

pub const WA_G_START: usize = 2;

pub const WA_N_REP: usize = 3;

pub const WA_ANCHOR: usize = 4;

pub const WA_I_FRAG: usize = 5;

pub const WA_SJ_A: usize = 6;

pub const WA_SIZE: usize = 7;

pub const WC_STR: usize = 0;

pub const WC_CHR: usize = 1;

pub const WC_G_START: usize = 2;

pub const WC_G_END: usize = 3;

pub const WC_SIZE: usize = 4;

pub const EXIT_CREATE_EXTEND_WINDOWS_WITH_ALIGN_TOO_MANY_WINDOWS: i32 = 101;

pub const PC_R_START: usize = 0;

pub const PC_LENGTH: usize = 1;

pub const PC_STR: usize = 2;

pub const PC_DIR: usize = 3;

pub const PC_NREP: usize = 4;

pub const PC_SASTART: usize = 5;

pub const PC_SAEND: usize = 6;

pub const PC_IFRAG: usize = 7;

pub const PC_SIZE: usize = 8;

pub const BAM_CIGAR_M: u32 = 0;

pub const BAM_CIGAR_I: u32 = 1;

pub const BAM_CIGAR_D: u32 = 2;

pub const BAM_CIGAR_N: u32 = 3;

pub const BAM_CIGAR_S: u32 = 4;

pub const BAM_CIGAR_H: u32 = 5;

pub const UINT_WIN_BIN_MAX: u32 = u16::MAX as u32;

pub const MARKER_TOO_MANY_ANCHORS_PER_WINDOW: u64 = 999_905;

pub const SJ_MOTIF_SIZE: usize = 7;

pub const SOLO_READ_BARCODE_N_STATS: usize = 12;

pub const SOLO_READ_FEATURE_N_STATS: usize = 11;

pub const SOLO_READ_FLAG_N_BITS: usize = 15;

pub const SOLO_READ_FLAG_CB_MATCH: u32 = 0;

pub const SOLO_READ_FLAG_CB_PERFECT: u32 = 1;

pub const SOLO_READ_FLAG_CB_MM_UNIQUE: u32 = 2;

pub const SOLO_READ_FLAG_CB_MM_MULTIPLE: u32 = 3;

pub const SOLO_READ_FLAG_GENOME_U: u32 = 4;

pub const SOLO_READ_FLAG_GENOME_M: u32 = 5;

pub const SOLO_READ_FLAG_FEATURE_U: u32 = 6;

pub const SOLO_READ_FLAG_FEATURE_M: u32 = 7;

pub const SOLO_READ_FLAG_EXONIC: u32 = 8;

pub const SOLO_READ_FLAG_INTRONIC: u32 = 9;

pub const SOLO_READ_FLAG_EXONIC_AS: u32 = 10;

pub const SOLO_READ_FLAG_INTRONIC_AS: u32 = 11;

pub const SOLO_READ_FLAG_MITO: u32 = 12;

pub const SOLO_READ_FLAG_COUNTED_U: u32 = 13;

pub const SOLO_READ_FLAG_COUNTED_M: u32 = 14;

pub const SOLO_READ_FEATURE_STAT_NO_UNMAPPED: usize = 0;

pub const SOLO_READ_FEATURE_STAT_NO_NO_FEATURE: usize = 1;

pub const SOLO_READ_FEATURE_STAT_MULTI_FEATURE: usize = 2;

pub const SOLO_READ_FEATURE_STAT_SUB_MULTI_FEATURE_MULTI_GENOMIC: usize = 3;

pub const SOLO_READ_FEATURE_STAT_NO_TOO_MANY_WL_MATCHES: usize = 4;

pub const SOLO_READ_FEATURE_STAT_NO_MM_TO_WL_WITHOUT_EXACT: usize = 5;

pub const SOLO_READ_FEATURE_STAT_YES_WL_MATCH: usize = 6;

pub const SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_EXACT: usize = 7;

pub const SOLO_READ_FEATURE_STAT_YES_SUB_WL_MATCH_UNIQUE_FEATURE: usize = 8;

pub const SOLO_READ_FEATURE_STAT_YES_CELL_BARCODES: usize = 9;

pub const SOLO_READ_FEATURE_STAT_YES_UMIS: usize = 10;

pub const SOLO_FEATURE_SJ: i32 = 0;

pub const SOLO_FEATURE_TRANSCRIPT3P: i32 = 1;

pub const SOLO_FEATURE_GENE_FULL: i32 = 2;

pub const SOLO_FEATURE_GENE_FULL_EXON_OVER_INTRON: i32 = 3;

pub const SOLO_FEATURE_GENE_FULL_EX50P_AS: i32 = 4;

pub const SOLO_TYPE_NONE: i32 = 0;

pub const SOLO_TYPE_CB_UMI_SIMPLE: i32 = 1;

pub const SOLO_TYPE_CB_UMI_COMPLEX: i32 = 2;

pub const SOLO_TYPE_CB_SAM_TAG_OUT: i32 = 3;

pub const SOLO_TYPE_SMART_SEQ: i32 = 4;

pub const SOLO_FEATURE_GENE: i32 = 5;

pub const SOLO_FEATURE_VELOCYTO_SIMPLE: i32 = 6;

pub const SOLO_FEATURE_VELOCYTO: i32 = 7;

pub const GENOME_SPACING_CHAR: u8 = 5;

pub const DEF_READ_NAME_LENGTH_MAX: usize = 50_000;

pub const DEF_READ_SEQ_LENGTH_MAX: usize = 650;

pub const ATTR_NH: i32 = 1;

pub const ATTR_HI: i32 = 2;

pub const ATTR_AS: i32 = 3;

pub const ATTR_NM: i32 = 4;

pub const ATTR_MD: i32 = 5;

pub const ATTR_NM_LOWER: i32 = 6;

pub const ATTR_JM: i32 = 7;

pub const ATTR_JI: i32 = 8;

pub const ATTR_XS: i32 = 9;

pub const ATTR_RG: i32 = 10;

pub const ATTR_VG: i32 = 11;

pub const ATTR_VA: i32 = 12;

pub const ATTR_VW: i32 = 13;

pub const ATTR_CH: i32 = 14;

pub const ATTR_MC: i32 = 15;

pub const ATTR_RB: i32 = 16;

pub const ATTR_CR: i32 = 17;

pub const ATTR_CY: i32 = 18;

pub const ATTR_UR: i32 = 19;

pub const ATTR_UY: i32 = 20;

pub const ATTR_CB: i32 = 21;

pub const ATTR_UB: i32 = 22;

pub const ATTR_GX: i32 = 23;

pub const ATTR_GN: i32 = 24;

pub const ATTR_SM: i32 = 25;

pub const ATTR_SS: i32 = 26;

pub const ATTR_SQ: i32 = 27;

pub const ATTR_HA: i32 = 28;

pub const ATTR_CN: i32 = 29;

pub const ATTR_GX_LOWER: i32 = 30;

pub const ATTR_GN_LOWER: i32 = 31;

pub const ATTR_SF: i32 = 32;

pub const SJ_SAM_ANNOTATED_MOTIF_SHIFT: i32 = 20;

pub const SPLICEGRAPH_MAX_SEQ_LENGTH: usize = 100_000;

pub const ALIGN_VS_TRANSCRIPT_INTRON: i32 = 0;

pub const ALIGN_VS_TRANSCRIPT_EXON_INTRON: i32 = 1;

pub const ALIGN_VS_TRANSCRIPT_EXON_INTRON_SPAN: i32 = 2;

pub const ALIGN_VS_TRANSCRIPT_CONCORDANT: i32 = 3;

pub const EXIT_CODE_INPUT_FILES: i32 = 104;

pub const EXIT_CODE_GENOME_FILES: i32 = 105;

pub const EXIT_CODE_SHM: i32 = 106;

pub const EXIT_CODE_MEMORY_ALLOCATION: i32 = 108;

pub const SHM_EOPENFAILED: i32 = 4;

pub const SHM_EEXISTS: i32 = 5;

pub const SHM_EFTRUNCATE: i32 = 6;

pub const SHM_EMAPFAILED: i32 = 7;

pub const SHM_ECLOSE: i32 = 8;

pub const SHM_EUNLINK: i32 = 9;

pub const SHM_ECOUNTERCREATE: i32 = 10;

pub const SHM_ECOUNTERREMOVE: i32 = 11;

pub const SHM_ECOUNTERUSE: i32 = 12;
