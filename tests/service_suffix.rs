use star_rs::generated::functions::*;
use star_rs::generated::structs::{
    Genome, PackedArray, Parameters, ParametersGenome, ReadAlign, StoredAlign, Transcript,
};

#[test]
fn service_comparators_match_lexicographic_rules() {
    assert_eq!(servicefuns_l7_sum1d(&[1u32, 2, 3, 100], 3), 6);
    assert_eq!(servicefuns_l13_funcomparenumbers(&7, &3), 1);
    assert_eq!(servicefuns_l13_funcomparenumbers(&3, &7), -1);
    assert_eq!(servicefuns_l13_funcomparenumbers(&7, &7), 0);
    assert_eq!(servicefuns_l26_funcomparenumbersreverse(&7, &3), -1);
    assert_eq!(servicefuns_l26_funcomparenumbersreverse(&3, &7), 1);
    assert_eq!(
        servicefuns_l39_funcomparenumbersreverseshift::<u32, 2>(&[0, 1, 9], &[5, 4, 7]),
        -1
    );
    assert_eq!(servicefuns_l53_funcompareuint1(&7, &3), 1);
    assert_eq!(servicefuns_l53_funcompareuint1(&7, &7), 0);
    assert_eq!(servicefuns_l66_funcompareuint2(&[1, 9], &[1, 8]), 1);
    assert_eq!(servicefuns_l66_funcompareuint2(&[1, 8], &[1, 8]), 0);
    assert_eq!(servicefuns_l66_funcompareuint2(&[0, 99], &[1, 0]), -1);

    assert_eq!(
        servicefuns_l84_funcomparearrays::<u32, 3>(&[1, 2, 3], &[1, 2, 4]),
        -1
    );
    assert_eq!(
        servicefuns_l101_funcomparearraysreverse::<u32, 3>(&[1, 2, 3], &[1, 2, 4]),
        1
    );
    assert_eq!(
        servicefuns_l118_funcomparearraysshift::<u32, 2, 1>(&[9, 1, 5], &[8, 1, 4]),
        1
    );
    assert_eq!(
        servicefuns_l135_funcomparetypesecondfirst(&[9, 1], &[8, 1]),
        1
    );
    assert_eq!(
        servicefuns_l135_funcomparetypesecondfirst(&[9, 1], &[9, 1]),
        0
    );
    assert_eq!(
        servicefuns_l135_funcomparetypesecondfirst(&[9, 0], &[0, 1]),
        -1
    );
    assert_eq!(
        servicefuns_l153_funcomparetypeshift::<u32, 1>(&[9, 2], &[8, 3]),
        -1
    );
}

#[test]
fn service_split_string_preserves_empty_fields_like_getline() {
    let mut elems = vec!["old".to_string()];
    let max_l = servicefuns_l167_splitstring("aa,,bbbb,", ',', &mut elems);
    assert_eq!(max_l, 4);
    assert_eq!(elems, vec!["aa", "", "bbbb"]);
}

#[test]
fn service_binary_search_variants_match_boundary_contracts() {
    let arr = [1, 2, 2, 2, 5, 8];
    let mut i1 = 99;
    assert!(servicefuns_l212_binarysearch_leleft(
        2,
        &arr,
        arr.len() as u32,
        &mut i1
    ));
    assert_eq!(i1, 1);
    assert!(!servicefuns_l212_binarysearch_leleft(
        0,
        &arr,
        arr.len() as u32,
        &mut i1
    ));

    assert_eq!(servicefuns_l192_binarysearch1(2, &arr, arr.len() as u32), 3);
    assert_eq!(servicefuns_l192_binarysearch1(4, &arr, arr.len() as u32), 3);
    assert_eq!(
        servicefuns_l192_binarysearch1(0, &arr, arr.len() as u32),
        u32::MAX
    );

    assert_eq!(
        servicefuns_l239_binarysearch1a(2, &arr, arr.len() as i32),
        3
    );
    assert_eq!(
        servicefuns_l239_binarysearch1a(9, &arr, arr.len() as i32),
        5
    );
    assert_eq!(
        servicefuns_l239_binarysearch1a(0, &arr, arr.len() as i32),
        -1
    );

    assert_eq!(
        servicefuns_l266_binarysearch1b(2, &arr, arr.len() as i32),
        1
    );
    assert_eq!(
        servicefuns_l266_binarysearch1b(0, &arr, arr.len() as i32),
        0
    );
    assert_eq!(
        servicefuns_l266_binarysearch1b(9, &arr, arr.len() as i32),
        -1
    );

    assert_eq!(
        servicefuns_l294_binarysearchexact(5, &arr, arr.len() as u64),
        4
    );
    assert_eq!(
        servicefuns_l294_binarysearchexact(4, &arr, arr.len() as u64),
        -1
    );
}

#[test]
fn suffix_leaf_math_matches_original_overflow_safe_logic() {
    assert_eq!(
        suffixarrayfuns_l4_medianuint2(u32::MAX - 1, u32::MAX),
        u32::MAX - 1
    );
    assert_eq!(suffixarrayfuns_l397_funcalcsai(&[0, 1, 2, 3], 3), 27);
    assert_eq!(suffixarrayfuns_l397_funcalcsai(&[0, 1, 4, 3], 3), -1);

    assert_eq!(
        suffixarrayfuns_l210_comparerefends(100, 10, 20, true, true),
        1
    );
    assert_eq!(
        suffixarrayfuns_l210_comparerefends(100, 30, 20, true, true),
        -1
    );
    assert_eq!(
        suffixarrayfuns_l210_comparerefends(100, 10, 20, false, false),
        1
    );
    assert_eq!(
        suffixarrayfuns_l210_comparerefends(100, 90, 20, false, false),
        -1
    );
    assert_eq!(
        suffixarrayfuns_l210_comparerefends(100, 10, u64::MAX, false, false),
        -1
    );
}

#[test]
fn sj_align_split_detects_only_alignments_crossing_sj_spacer() {
    let genome = Genome {
        sj_gstart: 1_000,
        sjdb_length: 21,
        sjdb_overhang: 10,
        sj_dstart: vec![10_000, 20_000],
        sj_astart: vec![11_000, 21_000],
        ..Default::default()
    };
    let mut a1_d = 0;
    let mut a_length_d = 0;
    let mut a1_a = 0;
    let mut a_length_a = 0;
    let mut isj = 0;

    assert!(sjalignsplit_l3_sjalignsplit(
        1_000 + 21 + 8,
        5,
        &genome,
        &mut a1_d,
        &mut a_length_d,
        &mut a1_a,
        &mut a_length_a,
        &mut isj,
    ));
    assert_eq!(isj, 1);
    assert_eq!(a_length_d, 2);
    assert_eq!(a_length_a, 3);
    assert_eq!(a1_d, 20_008);
    assert_eq!(a1_a, 21_000);

    assert!(!sjalignsplit_l3_sjalignsplit(
        1_000 + 3,
        6,
        &genome,
        &mut a1_d,
        &mut a_length_d,
        &mut a1_a,
        &mut a_length_a,
        &mut isj,
    ));
    assert!(!sjalignsplit_l3_sjalignsplit(
        1_000 + 12,
        6,
        &genome,
        &mut a1_d,
        &mut a_length_d,
        &mut a1_a,
        &mut a_length_a,
        &mut isj,
    ));
}

#[test]
fn suffix_compare_seq_to_genome_matches_forward_reverse_and_spacer_rules() {
    let genome = Genome {
        g: vec![0, 1, 2, 3, 4, 5],
        sa: vec![1, 2, 0b1_0000 | 1],
        n_genome: 6,
        gstrand_bit: 4,
        gstrand_mask: 0b1111,
        ..Default::default()
    };
    let fwd = [9u8, 1, 2, 3, 9];
    let rev = [9u8, 4, 3, 2, 9];

    let mut comp = false;
    assert_eq!(
        suffixarrayfuns_l10_compareseqtogenome(&genome, [&fwd, &rev], 1, 3, 0, 0, true, &mut comp),
        3
    );

    let greater = [9u8, 1, 4, 3, 9];
    assert_eq!(
        suffixarrayfuns_l10_compareseqtogenome(
            &genome,
            [&greater, &rev],
            1,
            3,
            0,
            0,
            true,
            &mut comp,
        ),
        1
    );
    assert!(comp);

    assert_eq!(
        suffixarrayfuns_l10_compareseqtogenome(&genome, [&fwd, &rev], 1, 3, 0, 2, true, &mut comp),
        3
    );

    let spacer_genome = Genome {
        g: vec![0, 5, 2, 9],
        sa: vec![0, 3],
        n_genome: 4,
        gstrand_bit: 4,
        gstrand_mask: 0b1111,
        ..Default::default()
    };
    let spaced = [0u8, 5, 3];
    let mut comp1 = 0;
    assert_eq!(
        suffixarrayfuns_l221_compareseqtogenome1(
            &spacer_genome,
            [&spaced, &spaced],
            0,
            3,
            0,
            0,
            true,
            2,
            &mut comp1,
        ),
        1
    );
    assert_eq!(comp1, 1);
}

#[test]
fn suffix_array_search1_returns_exact_or_first_greater_suffix_index() {
    let genome = Genome {
        g: vec![0, 1, 2, 3],
        sa: vec![0, 1, 3],
        n_genome: 4,
        gstrand_bit: 4,
        gstrand_mask: 0b1111,
        ..Default::default()
    };
    let query_exact = [1u8, 2];
    let mut l = 0;
    assert_eq!(
        suffixarrayfuns_l297_suffixarraysearch1(
            &genome,
            [&query_exact, &query_exact],
            0,
            2,
            u64::MAX,
            true,
            0,
            2,
            &mut l,
        ),
        1
    );
    assert_eq!(l, 2);

    let query_between = [2u8, 0];
    l = 0;
    assert_eq!(
        suffixarrayfuns_l297_suffixarraysearch1(
            &genome,
            [&query_between, &query_between],
            0,
            2,
            u64::MAX,
            true,
            0,
            2,
            &mut l,
        ),
        2
    );
    assert_eq!(l, 0);
}

#[test]
fn suffix_find_mult_range_and_max_mappable_length_expand_exact_match_range() {
    let genome = Genome {
        g: vec![1, 2, 3, 1, 2, 3, 1, 2, 4],
        sa: vec![0, 3, 6],
        n_genome: 9,
        gstrand_bit: 4,
        gstrand_mask: 0b1111,
        ..Default::default()
    };
    let query = [1u8, 2, 3];

    assert_eq!(
        suffixarrayfuns_l106_findmultrange(
            &genome,
            1,
            3,
            2,
            2,
            1,
            3,
            2,
            2,
            [&query, &query],
            true,
            0,
        ),
        1
    );

    let mut l = 0;
    let mut range = [u32::MAX; 2];
    assert_eq!(
        suffixarrayfuns_l133_maxmappablelength(
            &genome,
            [&query, &query],
            0,
            3,
            0,
            2,
            true,
            &mut l,
            &mut range,
        ),
        2
    );
    assert_eq!(l, 3);
    assert_eq!(range, [0, 1]);
}

#[test]
fn genome_sa_index_find_next_index_matches_step_and_binary_search() {
    let genome = Genome {
        n_sa: 6,
        n_genome: 7,
        gstrand_bit: 31,
        gstrand_mask: u32::MAX,
        p_ge: ParametersGenome {
            g_saindex_nbases: 2,
            ..Default::default()
        },
        ..Default::default()
    };
    let g = [0u8, 0, 0, 1, 1, 2, 3];
    let sa = [0u64, 1, 2, 3, 4, 5];
    let mut i_l4 = -1;
    let mut ind_full =
        suffixarrayfuns_l353_funcalcsaifromsa(&g, &sa, 7, 31, u32::MAX as u64, 0, 2, &mut i_l4);
    let mut isa = 0;

    genomesaindex_l178_funsaifindnextindex(&g, &sa, 2, &mut isa, &mut ind_full, &mut i_l4, &genome);

    assert_eq!(isa, 2);
    assert_eq!(ind_full, 1);
    assert_eq!(i_l4, -1);
}

#[test]
fn genome_sa_index_find_next_index_preserves_cpp_end_jump_behavior() {
    let genome = Genome {
        n_sa: 6,
        n_genome: 7,
        gstrand_bit: 31,
        gstrand_mask: u32::MAX,
        p_ge: ParametersGenome {
            g_saindex_nbases: 2,
            ..Default::default()
        },
        ..Default::default()
    };
    let g = [0u8, 0, 0, 1, 1, 2, 3];
    let sa = [0u64, 1, 2, 3, 4, 5];
    let mut i_l4 = -1;
    let mut ind_full =
        suffixarrayfuns_l353_funcalcsaifromsa(&g, &sa, 7, 31, u32::MAX as u64, 4, 2, &mut i_l4);
    let mut isa = 4;

    genomesaindex_l178_funsaifindnextindex(&g, &sa, 2, &mut isa, &mut ind_full, &mut i_l4, &genome);

    assert_eq!(isa, 6);
    assert_eq!(ind_full, 11);
    assert_eq!(i_l4, -1);
}

#[test]
fn genome_sa_index_chunk_writes_present_absent_and_tail_entries() {
    let absent = 1 << 6;
    let genome = Genome {
        n_sa: 6,
        n_genome: 7,
        gstrand_bit: 31,
        gstrand_mask: u32::MAX,
        genome_sa_index_start: vec![0, 4, 20],
        sai_mark_absent_mask_c: absent,
        sai_mark_nmask_c: 1 << 5,
        p_ge: ParametersGenome {
            g_saindex_nbases: 2,
            ..Default::default()
        },
        ..Default::default()
    };
    let g = [0u8, 0, 0, 1, 1, 2, 3];
    let sa = [0u64, 1, 2, 3, 4, 5];
    let mut sai = vec![0u32; 20];

    genomesaindex_l117_genomesaindexchunk(&g, &sa, &mut sai, 0, 5, &genome).unwrap();

    assert_eq!(
        sai,
        vec![
            0,
            3,
            5,
            6 | absent,
            0,
            2,
            3 | absent,
            3 | absent,
            3 | absent,
            3,
            4,
            5 | absent,
            5 | absent,
            5 | absent,
            5 | absent,
            5,
            6 | absent,
            6 | absent,
            6 | absent,
            6 | absent,
        ]
    );
}

#[test]
fn genome_sa_index_chunk_marks_previous_entry_when_suffix_contains_n() {
    let absent = 1 << 6;
    let n_mask = 1 << 5;
    let genome = Genome {
        n_sa: 4,
        n_genome: 6,
        gstrand_bit: 31,
        gstrand_mask: u32::MAX,
        genome_sa_index_start: vec![0, 4, 20],
        sai_mark_absent_mask_c: absent,
        sai_mark_nmask_c: n_mask,
        p_ge: ParametersGenome {
            g_saindex_nbases: 2,
            ..Default::default()
        },
        ..Default::default()
    };
    let g = [0u8, 0, 0, 4, 0, 1];
    let sa = [0u64, 1, 2, 4];
    let mut sai = vec![0u32; 20];

    genomesaindex_l117_genomesaindexchunk(&g, &sa, &mut sai, 0, 3, &genome).unwrap();

    assert_eq!(sai[4], n_mask);
    assert_eq!(sai[5], 3);
    assert_eq!(sai[1], 4 | absent);
}

#[test]
fn genome_sa_index_initializes_offsets_masks_and_table() {
    let mut genome = Genome {
        n_genome: 7,
        n_sa: 6,
        gstrand_bit: 4,
        gstrand_mask: 0b1111,
        p_ge: ParametersGenome {
            g_saindex_nbases: 2,
            ..Default::default()
        },
        ..Default::default()
    };
    let g = [0u8, 0, 0, 1, 1, 2, 3];
    let sa = [0u64, 1, 2, 3, 4, 5];

    let sai = genomesaindex_l6_genomesaindex(&g, &sa, &mut genome).unwrap();

    assert_eq!(genome.genome_sa_index_start, vec![0, 4, 20]);
    assert_eq!(genome.sai_mark_nmask_c, 1 << 5);
    assert_eq!(genome.sai_mark_nmask, !(1 << 5));
    assert_eq!(genome.sai_mark_absent_mask_c, 1 << 6);
    assert_eq!(sai, genome.sai);
    assert_eq!(sai[0], 0);
    assert_eq!(sai[1], 3);
    assert_eq!(sai[3], 6 | (1 << 6));
    assert_eq!(sai[15], 5);
}

#[test]
fn suffix_max_mappable_length_reports_longest_partial_range() {
    let genome = Genome {
        g: vec![1, 2, 5, 1, 2, 6, 1, 3, 0],
        sa: vec![0, 3, 6],
        n_genome: 9,
        gstrand_bit: 4,
        gstrand_mask: 0b1111,
        ..Default::default()
    };
    let query = [1u8, 2, 7];
    let mut l = 0;
    let mut range = [u32::MAX; 2];

    assert_eq!(
        suffixarrayfuns_l133_maxmappablelength(
            &genome,
            [&query, &query],
            0,
            3,
            0,
            2,
            true,
            &mut l,
            &mut range,
        ),
        2
    );
    assert_eq!(l, 2);
    assert_eq!(range, [0, 1]);
}

#[test]
fn readalign_max_mappable_length_two_strands_uses_sa_index_and_records_store_aligns() {
    let mut sai = vec![0u32; 20];
    sai[4 + 6] = 0;
    sai[4 + 7] = 2;
    let genome = Genome {
        g: vec![1, 2, 3, 1, 2, 3, 1, 2, 4],
        sa: vec![0, 3, 6],
        sai,
        genome_sa_index_start: vec![0, 4, 20],
        n_genome: 9,
        n_sa: 3,
        sj_gstart: 100,
        gstrand_bit: 4,
        gstrand_mask: 0b1111,
        sai_mark_absent_mask_c: 1 << 6,
        sai_mark_nmask: !(1 << 5),
        sai_mark_nmask_c: 1 << 5,
        ..Default::default()
    };
    let p_ge = ParametersGenome {
        g_saindex_nbases: 2,
        g_sasparse_d: 1,
        ..Default::default()
    };
    let read = [1u8, 2, 3];
    let mut max_lbest = 0;
    let mut stored = Vec::new();

    assert_eq!(
        readalign_maxmappablelength2strands_l5_readalign_maxmappablelength2strands(
            &genome,
            &p_ge,
            [&read, &read],
            0,
            3,
            0,
            0,
            0,
            &mut max_lbest,
            0,
            &mut stored,
        )
        .unwrap(),
        2
    );
    assert_eq!(max_lbest, 3);
    assert_eq!(
        stored,
        vec![StoredAlign {
            i_dir: 0,
            shift: 0,
            n_rep: 2,
            l: 3,
            ind_start_end: [0, 1],
            i_frag: 0,
        }]
    );
}

#[test]
fn readalign_max_mappable_length_two_strands_accepts_unique_index_hit_without_sa_search() {
    let mut sai = vec![0u32; 20];
    sai[4 + 6] = 2;
    sai[4 + 7] = 3;
    let genome = Genome {
        g: vec![1, 2, 3, 1, 2, 3, 1, 2, 4],
        sa: vec![0, 3, 6],
        sai,
        genome_sa_index_start: vec![0, 4, 20],
        n_genome: 9,
        n_sa: 3,
        sj_gstart: 100,
        gstrand_bit: 4,
        gstrand_mask: 0b1111,
        sai_mark_absent_mask_c: 1 << 6,
        sai_mark_nmask: !(1 << 5),
        sai_mark_nmask_c: 1 << 5,
        ..Default::default()
    };
    let p_ge = ParametersGenome {
        g_saindex_nbases: 2,
        g_sasparse_d: 1,
        ..Default::default()
    };
    let read = [1u8, 2, 4];
    let mut max_lbest = 0;
    let mut stored = Vec::new();

    assert_eq!(
        readalign_maxmappablelength2strands_l5_readalign_maxmappablelength2strands(
            &genome,
            &p_ge,
            [&read, &read],
            0,
            3,
            0,
            0,
            0,
            &mut max_lbest,
            1,
            &mut stored,
        )
        .unwrap(),
        1
    );
    assert_eq!(max_lbest, 3);
    assert_eq!(
        stored,
        vec![StoredAlign {
            i_dir: 0,
            shift: 0,
            n_rep: 1,
            l: 3,
            ind_start_end: [2, 2],
            i_frag: 1,
        }]
    );
}

#[test]
fn readalign_map_one_read_records_unique_seed_and_stitches_window() {
    let mut sai = vec![0u32; 20];
    sai[4 + 6] = 0;
    sai[4 + 7] = 2;
    let genome = Genome {
        g: vec![1, 2, 3, 1, 2, 3, 1, 2, 4],
        sa: vec![0, 3, 6],
        sai,
        genome_sa_index_start: vec![0, 4, 20],
        n_genome: 9,
        n_sa: 3,
        gstrand_bit: 4,
        gstrand_mask: 0b1111,
        sai_mark_absent_mask_c: 1 << 6,
        sai_mark_nmask: !(1 << 5),
        sai_mark_nmask_c: 1 << 5,
        sj_gstart: 100,
        chr_bin: vec![0; 16],
        chr_start: vec![0],
        chr_length: vec![9],
        p_ge: ParametersGenome {
            g_saindex_nbases: 2,
            g_sasparse_d: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    let p = Parameters {
        max_nsplit: 4,
        seed_split_min: 1,
        seed_map_min: 0,
        seed_search_start_lmax: 50,
        seed_search_start_lmax_over_lread: 10.0,
        seed_multimap_nmax: 10,
        seed_per_read_nmax: 8,
        seed_per_window_nmax: 8,
        win_anchor_multimap_nmax: 10,
        win_bin_nbits: 0,
        win_bin_n: 16,
        win_bin_chr_nbits: 4,
        win_anchor_dist_nbins: 1,
        align_windows_per_read_nmax: 4,
        align_transcripts_per_window_nmax: 4,
        align_transcripts_per_read_nmax: 16,
        out_filter_match_nmin: 0,
        out_filter_intron_motifs: "None".to_string(),
        out_filter_multimap_score_range: 0,
        align_soft_clip_at_reference_ends_yes: true,
        ..Default::default()
    };
    let mut ra = ReadAlign {
        l_read: 3,
        read1: [vec![1, 2, 3], vec![2, 1, 0], vec![0, 1, 2]],
        read_length: vec![3, 0],
        read_length_original: vec![3, 0],
        read_nmates: 1,
        out_filter_mismatch_nmax_total: 10,
        max_score_mate: vec![i32::MIN],
        split_r: [vec![0; 4], vec![0; 4], vec![0; 4]],
        tr_init: Box::new(Transcript::default()),
        ..Default::default()
    };

    assert_eq!(
        readalign_maponeread_l6_readalign_maponeread(&mut ra, &genome, &p).unwrap(),
        0
    );

    assert_eq!(ra.map_marker, 0);
    assert_eq!(ra.n_split, 1);
    assert_eq!(ra.n_a, 2);
    assert_eq!(ra.n_p, 1);
    assert_eq!(ra.pc[0][PC_R_START], 0);
    assert_eq!(ra.pc[0][PC_LENGTH], 3);
    assert_eq!(ra.pc[0][PC_NREP], 2);
    assert_eq!(ra.n_w, 2);
    assert_eq!(ra.n_tr, 2);
    assert_eq!(ra.tr_best.max_score, 3);
    assert_eq!(ra.tr_all[0][0].exons[0][EX_G], 0);
    assert_eq!(ra.tr_all[0][0].exons[0][EX_L], 3);
}

#[test]
fn readalign_map_one_read_marks_no_good_pieces_after_quality_split() {
    let p = Parameters {
        max_nsplit: 4,
        seed_split_min: 3,
        out_filter_match_nmin: 0,
        ..Default::default()
    };
    let genome = Genome::default();
    let mut ra = ReadAlign {
        l_read: 2,
        read1: [
            vec![MARK_FRAG_SPACER_BASE, MARK_FRAG_SPACER_BASE],
            vec![0; 2],
            vec![0; 2],
        ],
        split_r: [vec![0; 4], vec![0; 4], vec![0; 4]],
        tr_init: Box::new(Transcript::default()),
        ..Default::default()
    };

    readalign_maponeread_l6_readalign_maponeread(&mut ra, &genome, &p).unwrap();

    assert_eq!(ra.n_split, 0);
    assert_eq!(ra.map_marker, MARKER_NO_GOOD_PIECES);
    assert_eq!(ra.tr_best.r_length, 0);
    assert_eq!(ra.n_w, 0);
}

#[test]
fn suffix_uint_comparators_use_suffix_bytes_after_primary_key() {
    let genome = [1u8, 2, 5, 1, 3, 5, 1, 2, 5];

    assert_eq!(
        funcompareuintandsuffixes_l6_funcompareuintandsuffixes(&[2, 0], &[3, 3], &genome),
        -1
    );
    assert_eq!(
        funcompareuintandsuffixes_l6_funcompareuintandsuffixes(&[4, 0], &[3, 3], &genome),
        1
    );
    assert_eq!(
        funcompareuintandsuffixes_l6_funcompareuintandsuffixes(&[2, 0], &[2, 3], &genome),
        -1
    );
    assert_eq!(
        funcompareuintandsuffixes_l6_funcompareuintandsuffixes(&[2, 6], &[2, 0], &genome),
        1
    );

    assert_eq!(
        funcompareuintandsuffixesmemcmp_l7_funcompareuintandsuffixesmemcmp(
            &[2, 0],
            &[2, 3],
            &genome,
            3,
        ),
        -1
    );
    assert_eq!(
        funcompareuintandsuffixesmemcmp_l7_funcompareuintandsuffixesmemcmp(
            &[2, 6],
            &[2, 0],
            &genome,
            3,
        ),
        1
    );
    assert_eq!(
        funcompareuintandsuffixesmemcmp_l7_funcompareuintandsuffixesmemcmp(
            &[2, 0],
            &[2, 6],
            &genome,
            3,
        ),
        -1
    );
}

#[test]
fn insert_seq_sa_shifts_existing_indices_with_64_bit_packed_array() {
    let mut sa = packedarray_l3_packedarray_packedarray();
    packedarray_l8_packedarray_definebits(&mut sa, 33, 4);
    packedarray_l31_packedarray_allocatearray(&mut sa);
    let n2bit = 1u64 << 32;
    packedarray_l17_packedarray_writepacked(&mut sa, 0, 0);
    packedarray_l17_packedarray_writepacked(&mut sa, 1, 2);
    packedarray_l17_packedarray_writepacked(&mut sa, 2, n2bit | 1);
    packedarray_l17_packedarray_writepacked(&mut sa, 3, n2bit | 2);

    let mut sa1 = PackedArray::default();
    let mut sai = PackedArray::default();
    let n_ind = insertseqsa_l18_insertseqsa(
        &mut sa,
        &mut sa1,
        &mut sai,
        &[0, 1],
        &[GENOME_SPACING_CHAR],
        2,
        1,
        2,
        &Parameters::default(),
        &Genome::default(),
    )
    .unwrap();

    assert_eq!(n_ind, 0);
    assert_eq!(packedarray_h18_packedarray_index(&sa, 0), 0);
    assert_eq!(packedarray_h18_packedarray_index(&sa, 1), 3);
    assert_eq!(packedarray_h18_packedarray_index(&sa, 2), n2bit | 1);
    assert_eq!(packedarray_h18_packedarray_index(&sa, 3), n2bit | 3);
    assert_eq!(sa1.word_length, 33);
    assert_eq!(sa1.length, 4);
    assert_eq!(packedarray_h18_packedarray_index(&sa1, 0), 0);
    assert_eq!(packedarray_h18_packedarray_index(&sa1, 1), 3);
    assert_eq!(packedarray_h18_packedarray_index(&sa1, 2), n2bit | 1);
    assert_eq!(packedarray_h18_packedarray_index(&sa1, 3), n2bit | 3);
}

#[test]
fn sjdb_build_index_inserts_new_junction_suffix_and_copies_sequence() {
    let mut sa = packedarray_l3_packedarray_packedarray();
    packedarray_l8_packedarray_definebits(&mut sa, 33, 1);
    packedarray_l31_packedarray_allocatearray(&mut sa);
    packedarray_l17_packedarray_writepacked(&mut sa, 0, 0);
    let mut sa2 = PackedArray::default();
    let mut sai = PackedArray::default();

    let mut gsj = vec![0, 1, 0, 0, 0];
    let mut genome_seq = vec![1, GENOME_SPACING_CHAR, GENOME_SPACING_CHAR];
    genome_seq.resize(12, GENOME_SPACING_CHAR);
    let mut map_gen = Genome {
        g: genome_seq.clone(),
        sa: vec![0],
        n_genome: 3,
        n_sa: 1,
        n_chr_real: 0,
        chr_start: vec![10],
        gstrand_bit: 32,
        gstrand_mask: u32::MAX,
        sjdb_n: 1,
        sjdb_length: 2,
        sjdb_start: vec![100],
        sjdb_end: vec![200],
        ..Default::default()
    };
    let map_gen1 = Genome {
        n_genome: 3,
        n_sa: 1,
        sjdb_n: 0,
        ..Default::default()
    };

    let result = sjdbbuildindex_l16_sjdbbuildindex(
        &Parameters::default(),
        &mut gsj,
        &mut genome_seq,
        &mut sa,
        &mut sa2,
        &mut sai,
        &mut map_gen,
        &map_gen1,
    )
    .unwrap();

    assert_eq!(result.sj_new, 1);
    assert_eq!(result.n_ind, 2);
    assert_eq!(map_gen.n_genome, 12);
    assert_eq!(map_gen.n_sa, 3);
    assert_eq!(&genome_seq[10..12], &[0, GENOME_SPACING_CHAR]);
    let packed_sa: Vec<u64> = (0..map_gen.n_sa as u64)
        .map(|ii| packedarray_h18_packedarray_index(&sa, ii))
        .collect();
    assert!(packed_sa.contains(&10));
    assert!(packed_sa.contains(&0));
    assert_eq!(map_gen.sj_gstart, 10);
    assert!(
        result
            .log_main
            .contains("Finished inserting junction indices")
    );
}

#[test]
fn sjdb_build_index_returns_without_junctions() {
    let mut gsj = Vec::new();
    let mut g = vec![0, 1, 2];
    let mut sa = PackedArray::default();
    let mut sa2 = PackedArray::default();
    let mut sai = PackedArray::default();
    let mut map_gen = Genome::default();
    let map_gen1 = Genome::default();

    let result = sjdbbuildindex_l16_sjdbbuildindex(
        &Parameters::default(),
        &mut gsj,
        &mut g,
        &mut sa,
        &mut sa2,
        &mut sai,
        &mut map_gen,
        &map_gen1,
    )
    .unwrap();

    assert_eq!(result.n_ind, 0);
    assert!(result.log_main.is_empty());
    assert_eq!(g, vec![0, 1, 2]);
}

#[test]
fn genome_suffix_comparator_matches_reverse_word_and_sentinel_rules() {
    let genome = [
        9u8, 9, 9, 9, 9, 9, 9, 9, //
        1, 2, 3, 4, 5, 7, 8, 9, //
        1, 2, 3, 4, 5, 7, 8, 10, //
        1, 2, 3, 4, 6, 7, 8, 9, //
    ];

    assert_eq!(
        genome_genomegenerate_l29_funcomparesuffixes(15, 23, &genome, 1),
        -1
    );
    assert_eq!(
        genome_genomegenerate_l29_funcomparesuffixes(31, 15, &genome, 1),
        1
    );
    assert_eq!(
        genome_genomegenerate_l29_funcomparesuffixes(23, 15, &genome, 1),
        1
    );

    let equal_prefix = [0u8, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1];
    assert_eq!(
        genome_genomegenerate_l29_funcomparesuffixes(15, 7, &equal_prefix, 1),
        -1
    );
}

#[test]
fn fun_calc_sai_from_sa_matches_forward_reverse_and_non_base_stop() {
    let genome = [
        0u8, 1, 2, 3, 0, 1, 2, 3, 5, 1, 2, 3, 0, 1, 2, 3, //
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 5, 3, 0, 1, 2, 3,
    ];
    let gstrand_bit = 5;
    let gstrand_mask = (1u64 << gstrand_bit) - 1;
    let sa = [2u64, (1u64 << gstrand_bit) | 12, 8];

    let mut il4 = 99;
    assert_eq!(
        suffixarrayfuns_l353_funcalcsaifromsa(
            &genome,
            &sa,
            genome.len() as u64,
            gstrand_bit,
            gstrand_mask,
            0,
            4,
            &mut il4,
        ),
        0b10_11_00_01
    );
    assert_eq!(il4, -1);

    assert_eq!(
        suffixarrayfuns_l353_funcalcsaifromsa(
            &genome,
            &sa,
            genome.len() as u64,
            gstrand_bit,
            gstrand_mask,
            1,
            4,
            &mut il4,
        ),
        0b00_01_10_11
    );
    assert_eq!(il4, -1);

    assert_eq!(
        suffixarrayfuns_l353_funcalcsaifromsa(
            &genome,
            &sa,
            genome.len() as u64,
            gstrand_bit,
            gstrand_mask,
            2,
            4,
            &mut il4,
        ),
        0
    );
    assert_eq!(il4, 0);
}

#[test]
fn linux_proc_memory_has_star_suffix_format() {
    let value = systemfunctions_l6_linuxprocmemory();
    assert!(value.ends_with('\n'));
    assert!(value.contains("VmRSS") || value == "\n");
}
