use star_rs::generated::functions::*;
use star_rs::generated::structs::OpalSearchResult;

unsafe extern "C" {
    fn tzset();
}

#[test]
fn time_month_day_time_matches_star_format_for_epoch() {
    // STAR formats "%b %d %H:%M:%SS" and then erases the last byte.
    unsafe {
        std::env::set_var("TZ", "UTC");
        tzset();
    }
    assert_eq!(timefunctions_l14_timemonthdaytime(0), "Jan 01 00:00:00");
}

#[test]
fn opal_array_and_gap_primitives_match_cpp_logic() {
    assert_eq!(opal_l1031_arraymax(&[3, 1, 7, 7, 2], 5), 7);
    assert_eq!(opal_l1031_arraymax(&[3, 1, 7], 99), 7);
    assert_eq!(opal_l1031_arraymax::<i32>(&[], 1), 0);
    assert_eq!(opal_l1031_arraymax(&[3, 1, 7], -1), 0);
    assert_eq!(opal_l1048_gappenalty(0, 5, 2), 0);
    assert_eq!(opal_l1048_gappenalty(1, 5, 2), 5);
    assert_eq!(opal_l1048_gappenalty(4, 5, 2), 11);
}

#[test]
fn opal_simd_sw_lane_primitives_match_intrinsic_semantics() {
    assert_eq!(opal_l87_simdisallzeroes(&[0, 0, 0]), 1);
    assert_eq!(opal_l87_simdisallzeroes(&[0, 1, 0]), 0);
    assert_eq!(opal_l148_print_mmxxxi(&[-3i8, 0, 127]), "-3 0 127 ");

    let mut a8 = [0i8; 32];
    let mut b8 = [0i8; 32];
    a8[0] = 120;
    b8[0] = 20;
    a8[1] = -120;
    b8[1] = 20;
    a8[2] = -1;
    b8[2] = 1;
    assert_eq!(opal_l103_add(a8, b8)[0], 127);
    assert_eq!(opal_l104_sub(a8, b8)[1], -128);
    assert_eq!(opal_l105_min(a8, b8)[2], 1);
    assert_eq!(opal_l106_max(a8, b8)[2], -1);
    assert_eq!(opal_l107_cmpgt(a8, b8)[0], -1);
    assert_eq!(opal_l107_cmpgt(a8, b8)[1], 0);
    assert_eq!(opal_l108_set1(260)[0], 4);

    let mut a16 = [0i16; 16];
    let mut b16 = [0i16; 16];
    a16[0] = 32_000;
    b16[0] = 1_000;
    a16[1] = -32_000;
    b16[1] = 1_000;
    assert_eq!(opal_l117_add(a16, b16)[0], i16::MAX);
    assert_eq!(opal_l118_sub(a16, b16)[1], i16::MIN);
    assert_eq!(opal_l119_min(a16, b16)[1], -32_000);
    assert_eq!(opal_l120_max(a16, b16)[0], 32_000);
    assert_eq!(opal_l121_cmpgt(a16, b16)[0], -1);
    assert_eq!(opal_l122_set1(70_000)[0], 4_464);

    let a32 = [i32::MAX, i32::MIN, 5, 0, 0, 0, 0, 0];
    let b32 = [1, 1, 7, 0, 0, 0, 0, 0];
    assert_eq!(opal_l131_add(a32, b32)[0], i32::MIN);
    assert_eq!(opal_l132_sub(a32, b32)[1], i32::MAX);
    assert_eq!(opal_l133_min(a32, b32)[2], 5);
    assert_eq!(opal_l134_max(a32, b32)[2], 7);
    assert_eq!(opal_l135_cmpgt(a32, b32)[2], 0);
    assert_eq!(opal_l136_set1(-3), [-3; 8]);
}

#[test]
fn opal_simd_lane_primitives_use_signed_byte_minmax() {
    let mut a8 = [0i8; 32];
    let mut b8 = [0i8; 32];
    a8[0] = -1;
    b8[0] = 1;
    a8[1] = 120;
    b8[1] = 20;
    assert_eq!(opal_l558_add(a8, b8)[1], 127);
    assert_eq!(opal_l559_sub(a8, b8)[0], -2);
    assert_eq!(opal_l560_min(a8, b8)[0], -1);
    assert_eq!(opal_l561_max(a8, b8)[0], 1);
    assert_eq!(opal_l562_cmpgt(a8, b8)[1], -1);
    assert_eq!(opal_l563_set1(255)[0], -1);

    let a16 = [10i16, -10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let b16 = [5i16, -20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    assert_eq!(opal_l571_add(a16, b16)[0], 15);
    assert_eq!(opal_l572_sub(a16, b16)[1], 10);
    assert_eq!(opal_l573_min(a16, b16)[1], -20);
    assert_eq!(opal_l574_max(a16, b16)[0], 10);
    assert_eq!(opal_l575_cmpgt(a16, b16)[0], -1);
    assert_eq!(opal_l576_set1(-7)[0], -7);

    let a32 = [10i32, -10, i32::MAX, 0, 0, 0, 0, 0];
    let b32 = [5i32, -20, 1, 0, 0, 0, 0, 0];
    assert_eq!(opal_l584_add(a32, b32)[2], i32::MIN);
    assert_eq!(opal_l585_sub(a32, b32)[1], 10);
    assert_eq!(opal_l586_min(a32, b32)[1], -20);
    assert_eq!(opal_l587_max(a32, b32)[0], 10);
    assert_eq!(opal_l588_cmpgt(a32, b32)[0], -1);
    assert_eq!(opal_l589_set1(11), [11; 8]);
}

#[test]
fn opal_band_borders_cover_modes() {
    assert_eq!(opal_l1059_calculatebottombandborderov(8, 10, 7, 5, 2, 3), 3);
    assert_eq!(opal_l1074_calculatetopbandborderhw(8, 10, 7, 5, 2, 3), 0);
    assert_eq!(opal_l1089_calculatebottombandborderhw(8, 10, 7, 5, 2, 3), 3);
    assert_eq!(opal_l1106_calculatebottombandbordernw(8, 10, 7, 5, 2, 3), 3);
    assert_eq!(
        opal_l1153_calculatebandborders(8, OPAL_MODE_OV, 10, 7, 5, 2, 3),
        (3, 2)
    );
    assert_eq!(
        opal_l1153_calculatebandborders(1000, OPAL_MODE_SW, 10, 7, 5, 2, 3),
        (-1, -1)
    );
    assert_eq!(
        opal_l1153_calculatebandborders(8, OPAL_MODE_HW, 10, 7, 5, 2, 3),
        (3, 0)
    );
    assert_eq!(
        opal_l1153_calculatebandborders(8, OPAL_MODE_NW, 10, 7, 5, 2, 3),
        (3, 0)
    );
}

#[test]
fn opal_reverse_and_revert_match_cpp_logic() {
    assert_eq!(opal_l1188_createreversecopy(b"ACGT", 4), b"TGCA");
    assert_eq!(opal_l1188_createreversecopy(b"ACGT", 6), b"TGCA");
    assert_eq!(opal_l1188_createreversecopy(b"ACGT", -1), b"");

    let mut values = [1, 2, 3, 4, 5];
    opal_l1197_revertarray(&mut values, 5);
    assert_eq!(values, [5, 4, 3, 2, 1]);

    let mut short_values = [1, 2, 3];
    opal_l1197_revertarray(&mut short_values, 5);
    assert_eq!(short_values, [3, 2, 1]);
    opal_l1197_revertarray(&mut short_values, -1);
    assert_eq!(short_values, [3, 2, 1]);
}

#[test]
fn opal_load_next_sequence_skips_calculated_and_sets_null_at_end() {
    let db: [&[u8]; 3] = [b"AA", b"CCC", b"G"];
    let lengths = [2, 3, 1];
    let calculated = [true, false, true];
    let mut next = 0;
    let mut curr_idx = -9;
    let mut curr_pos = None;
    let mut curr_len = -9;
    let mut ended = 0;

    assert!(opal_l474_loadnextsequence(
        &mut next,
        3,
        &mut curr_idx,
        &mut curr_pos,
        &mut curr_len,
        &db,
        &lengths,
        &calculated,
        &mut ended,
    ));
    assert_eq!(next, 2);
    assert_eq!(curr_idx, 1);
    assert_eq!(curr_pos, Some(&b"CCC"[..]));
    assert_eq!(curr_len, 3);
    assert_eq!(ended, 1);

    assert!(!opal_l474_loadnextsequence(
        &mut next,
        3,
        &mut curr_idx,
        &mut curr_pos,
        &mut curr_len,
        &db,
        &lengths,
        &calculated,
        &mut ended,
    ));
    assert_eq!(next, 3);
    assert_eq!(curr_idx, -1);
    assert_eq!(curr_pos, None);
    assert_eq!(curr_len, -1);
    assert_eq!(ended, 2);
}

#[test]
fn opal_search_database_sw_scores_end_locations_and_skip_state() {
    let score_matrix = [2, -2, -2, -2, 2, -2, -2, -2, 2];
    let db: [&[u8]; 3] = [b"\x00\x01\x02", b"\x02\x02", b"\x00\x02\x01\x02"];
    let lengths = [3, 2, 4];
    let mut results = vec![
        OpalSearchResult::default(),
        OpalSearchResult {
            score: 99,
            score_set: 1,
            end_location_query: 7,
            end_location_target: 8,
            ..Default::default()
        },
        OpalSearchResult::default(),
    ];
    let skip = [false, true, false];

    assert_eq!(
        opal_l498_searchdatabasesw(
            b"\x00\x01\x02",
            3,
            &db,
            3,
            &lengths,
            3,
            1,
            &score_matrix,
            3,
            &mut results,
            OPAL_SEARCH_SCORE_END,
            Some(&skip),
            OPAL_OVERFLOW_SIMPLE,
        ),
        0
    );
    assert_eq!(results[0].score, 6);
    assert_eq!(results[0].end_location_query, 2);
    assert_eq!(results[0].end_location_target, 2);
    assert_eq!(results[1].score, 99);
    assert_eq!(results[1].end_location_query, 7);
    assert_eq!(results[2].score, 4);
    assert_eq!(results[2].end_location_query, 2);
    assert_eq!(results[2].end_location_target, 3);

    assert_eq!(
        opal_l498_searchdatabasesw(
            b"\x00\x01\x02",
            3,
            &db[..1],
            1,
            &lengths[..1],
            3,
            1,
            &score_matrix,
            3,
            &mut results[..1],
            OPAL_SEARCH_SCORE,
            None,
            OPAL_OVERFLOW_SIMPLE,
        ),
        0
    );
    assert_eq!(results[0].score, 6);
    assert_eq!(results[0].end_location_query, -1);
    assert_eq!(results[0].end_location_target, -1);

    let mut direct = vec![OpalSearchResult::default(), OpalSearchResult::default()];
    assert_eq!(
        opal_l167_searchdatabasesw(
            b"\x00\x01\x02",
            3,
            &db[..2],
            2,
            &lengths[..2],
            3,
            1,
            &score_matrix,
            3,
            &mut direct,
            OPAL_SEARCH_SCORE_END,
            &[false, true],
            OPAL_OVERFLOW_SIMPLE,
        ),
        0
    );
    assert_eq!(direct[0].score, 6);
    assert_eq!(direct[0].end_location_query, 2);
    assert_eq!(direct[0].end_location_target, 2);
    assert_eq!(direct[1], OpalSearchResult::default());

    let seqs = vec![vec![0u8]; 1025];
    let db_many: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let lengths_many = vec![1; 1025];
    let mut many_results = vec![OpalSearchResult::default(); 1025];
    many_results[1024].score = 77;
    let mut skip_many = vec![false; 1025];
    skip_many[1024] = true;
    assert_eq!(
        opal_l498_searchdatabasesw(
            b"\x00",
            1,
            &db_many,
            1025,
            &lengths_many,
            3,
            1,
            &score_matrix,
            3,
            &mut many_results,
            OPAL_SEARCH_SCORE,
            Some(&skip_many),
            OPAL_OVERFLOW_BUCKETS,
        ),
        0
    );
    assert_eq!(many_results[0].score, 2);
    assert_eq!(many_results[1024].score, 77);
}

#[test]
fn opal_search_database_wrapper_dispatches_sw_and_alignment_modes() {
    let score_matrix = [2, -2, -2, -2, 2, -2, -2, -2, 2];
    let db: [&[u8]; 2] = [b"\x00\x01\x02", b"\x02\x02"];
    let lengths = [3, 2];
    let mut results = vec![
        OpalSearchResult::default(),
        OpalSearchResult {
            score_set: 1,
            score: 44,
            end_location_query: 1,
            end_location_target: 1,
            alignment: Some(vec![9]),
            alignment_length: 1,
            ..Default::default()
        },
    ];

    assert_eq!(
        opal_l1437_opalsearchdatabase(
            b"\x00\x01\x02",
            3,
            &db,
            2,
            &lengths,
            3,
            1,
            &score_matrix,
            3,
            &mut results,
            OPAL_SEARCH_SCORE_END,
            OPAL_MODE_SW,
            OPAL_OVERFLOW_SIMPLE,
        ),
        0
    );
    assert_eq!(results[0].score, 6);
    assert_eq!(results[0].end_location_query, 2);
    assert_eq!(results[0].end_location_target, 2);
    assert_eq!(results[0].alignment, None);
    assert_eq!(results[0].alignment_length, -1);
    assert_eq!(results[0].start_location_query, -1);
    assert_eq!(results[1].score, 44);
    assert_eq!(results[1].alignment, None);
    assert_eq!(results[1].alignment_length, -1);

    let mut aligned = vec![OpalSearchResult::default()];
    assert_eq!(
        opal_l1437_opalsearchdatabase(
            b"\x00\x01\x02",
            3,
            &db[..1],
            1,
            &lengths[..1],
            3,
            1,
            &score_matrix,
            3,
            &mut aligned,
            OPAL_SEARCH_ALIGNMENT,
            OPAL_MODE_SW,
            OPAL_OVERFLOW_SIMPLE,
        ),
        0
    );
    assert_eq!(aligned[0].score, 6);
    assert_eq!(aligned[0].start_location_query, 0);
    assert_eq!(aligned[0].start_location_target, 0);
    assert_eq!(aligned[0].end_location_query, 2);
    assert_eq!(aligned[0].end_location_target, 2);
    assert_eq!(
        aligned[0].alignment,
        Some(vec![OPAL_ALIGN_MATCH, OPAL_ALIGN_MATCH, OPAL_ALIGN_MATCH])
    );
    assert_eq!(aligned[0].alignment_length, 3);

    assert_eq!(
        opal_l1437_opalsearchdatabase(
            b"\x00",
            1,
            &db[..1],
            1,
            &lengths[..1],
            3,
            1,
            &score_matrix,
            3,
            &mut aligned,
            OPAL_SEARCH_SCORE,
            99,
            OPAL_OVERFLOW_SIMPLE,
        ),
        OPAL_ERR_INVALID_MODE
    );
}

#[test]
fn opal_search_database_char_sw_scores_all_sequences() {
    let score_matrix = [2, -2, -2, -2, 2, -2, -2, -2, 2];
    let db: [&[u8]; 2] = [b"\x00\x01\x02", b"\x02\x02"];
    let lengths = [3, 2];
    let mut results = vec![OpalSearchResult::default(), OpalSearchResult::default()];

    assert_eq!(
        opal_l1526_opalsearchdatabasecharsw(
            b"\x00\x01\x02",
            3,
            &db,
            2,
            &lengths,
            3,
            1,
            &score_matrix,
            3,
            &mut results,
        ),
        0
    );

    assert_eq!(results[0].score_set, 1);
    assert_eq!(results[0].score, 6);
    assert_eq!(results[0].end_location_query, -1);
    assert_eq!(results[0].end_location_target, -1);
    assert_eq!(results[1].score_set, 1);
    assert_eq!(results[1].score, 2);
}

#[test]
fn opal_template_search_database_scores_modes_without_alignment_cleanup() {
    let score_matrix = [2, -2, -2, -2, 2, -2, -2, -2, 2];
    let db: [&[u8]; 2] = [b"\x00\x01\x02", b"\x00\x02"];
    let lengths = [3, 2];
    let mut results = vec![
        OpalSearchResult::default(),
        OpalSearchResult {
            score: 77,
            score_set: 1,
            alignment: Some(vec![9]),
            alignment_length: 1,
            ..Default::default()
        },
    ];
    let skip = [false, true];

    assert_eq!(
        opal_l986_searchdatabase(
            b"\x00\x01\x02",
            3,
            &db,
            2,
            &lengths,
            3,
            1,
            &score_matrix,
            3,
            &mut results,
            OPAL_SEARCH_SCORE_END,
            Some(&skip),
            OPAL_OVERFLOW_SIMPLE,
            OPAL_MODE_NW,
        ),
        0
    );
    assert_eq!(results[0].score, 6);
    assert_eq!(results[0].end_location_query, 2);
    assert_eq!(results[0].end_location_target, 2);
    assert_eq!(results[1].score, 77);
    assert_eq!(results[1].alignment, Some(vec![9]));

    assert_eq!(
        opal_l986_searchdatabase(
            b"\x00",
            1,
            &db[..1],
            1,
            &lengths[..1],
            3,
            1,
            &score_matrix,
            3,
            &mut results[..1],
            OPAL_SEARCH_SCORE,
            None,
            OPAL_OVERFLOW_SIMPLE,
            99,
        ),
        OPAL_ERR_INVALID_MODE
    );

    let mut direct = vec![OpalSearchResult::default(), OpalSearchResult::default()];
    assert_eq!(
        opal_l597_searchdatabase(
            b"\x00\x01\x02",
            3,
            &db,
            2,
            &lengths,
            3,
            1,
            &score_matrix,
            3,
            &mut direct,
            OPAL_SEARCH_SCORE_END,
            &[false, true],
            OPAL_OVERFLOW_SIMPLE,
            OPAL_MODE_NW,
        ),
        0
    );
    assert_eq!(direct[0].score, 6);
    assert_eq!(direct[0].end_location_query, 2);
    assert_eq!(direct[0].end_location_target, 2);
    assert_eq!(direct[1], OpalSearchResult::default());

    let seqs = vec![vec![0u8]; 1025];
    let db_many: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let lengths_many = vec![1; 1025];
    let mut many_results = vec![OpalSearchResult::default(); 1025];
    many_results[1024].score = 77;
    let mut skip_many = vec![false; 1025];
    skip_many[1024] = true;
    assert_eq!(
        opal_l986_searchdatabase(
            b"\x00",
            1,
            &db_many,
            1025,
            &lengths_many,
            3,
            1,
            &score_matrix,
            3,
            &mut many_results,
            OPAL_SEARCH_SCORE,
            Some(&skip_many),
            OPAL_OVERFLOW_BUCKETS,
            OPAL_MODE_NW,
        ),
        0
    );
    assert_eq!(many_results[0].score, 2);
    assert_eq!(many_results[1024].score, 77);
}

#[test]
fn opal_find_alignment_builds_nw_alignment_and_gap_trace() {
    let score_matrix = [2, -2, -2, -2, 2, -2, -2, -2, 2];
    let mut exact = OpalSearchResult::default();
    opal_l1238_findalignment(
        &[0, 1, 2],
        3,
        &[0, 1, 2],
        3,
        3,
        1,
        &score_matrix,
        3,
        -100,
        &mut exact,
        OPAL_MODE_NW,
    );
    assert_eq!(exact.score_set, 1);
    assert_eq!(exact.score, 6);
    assert_eq!(exact.start_location_target, 0);
    assert_eq!(exact.start_location_query, 0);
    assert_eq!(exact.end_location_target, 2);
    assert_eq!(exact.end_location_query, 2);
    assert_eq!(
        exact.alignment,
        Some(vec![OPAL_ALIGN_MATCH, OPAL_ALIGN_MATCH, OPAL_ALIGN_MATCH])
    );
    assert_eq!(exact.alignment_length, 3);

    let mut with_gap = OpalSearchResult::default();
    opal_l1238_findalignment(
        &[0, 1],
        2,
        &[0, 2, 1],
        3,
        3,
        1,
        &score_matrix,
        3,
        -100,
        &mut with_gap,
        OPAL_MODE_NW,
    );
    assert_eq!(with_gap.score, 1);
    assert_eq!(
        with_gap.alignment,
        Some(vec![OPAL_ALIGN_MATCH, OPAL_ALIGN_INS, OPAL_ALIGN_MATCH])
    );
    assert_eq!(with_gap.alignment_length, 3);
}

#[test]
fn opal_search_result_state_matches_cpp_api() {
    let mut result = OpalSearchResult {
        score_set: 1,
        score: 42,
        end_location_target: 8,
        end_location_query: 9,
        start_location_target: 2,
        start_location_query: 3,
        alignment: Some(vec![0, 1, 2]),
        alignment_length: 3,
    };

    opal_l1553_opalinitsearchresult(&mut result);
    assert_eq!(opal_l1561_opalsearchresultisempty(&result), 1);
    assert_eq!(result.start_location_target, -1);
    assert_eq!(result.start_location_query, -1);
    assert_eq!(result.end_location_target, -1);
    assert_eq!(result.end_location_query, -1);
    assert_eq!(result.alignment, None);
    assert_eq!(result.alignment_length, 0);

    opal_l1565_opalsearchresultsetscore(&mut result, 17);
    assert_eq!(opal_l1561_opalsearchresultisempty(&result), 0);
    assert_eq!(result.score, 17);
}
