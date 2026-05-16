use star_rs::*;
use star_rs::{
    ChimericSegment, ClipMate, Genome, GenomeOut, Parameters, ParametersChimeric, ParametersGenome,
    ReadAlign, SharedMemory, Stats, Transcript, Variation,
};
use star_rs::{SoloReadFeature, SoloReadFeatureStats};

// --- helpers for filling fixed-size Transcript array fields from slice literals ---
fn tr_exons_arr(items: &[[u64; 5]]) -> [[u64; 5]; 20] {
    let mut a = [[0u64; 5]; 20];
    for (i, v) in items.iter().enumerate() {
        a[i] = *v;
    }
    a
}
fn tr_canon_sj_arr(items: &[i32]) -> [i32; 20] {
    let mut a = [0i32; 20];
    for (i, v) in items.iter().enumerate() {
        a[i] = *v;
    }
    a
}
fn tr_sj_u8_arr(items: &[u8]) -> [u8; 20] {
    let mut a = [0u8; 20];
    for (i, v) in items.iter().enumerate() {
        a[i] = *v;
    }
    a
}
#[allow(dead_code)]
fn tr_shift_sj_arr(items: &[[u64; 2]]) -> [[u64; 2]; 20] {
    let mut a = [[0u64; 2]; 20];
    for (i, v) in items.iter().enumerate() {
        a[i] = *v;
    }
    a
}
#[allow(dead_code)]
fn tr_read_length_arr(items: &[u64]) -> [u64; 3] {
    let mut a = [0u64; 3];
    for (i, v) in items.iter().enumerate() {
        a[i] = *v;
    }
    a
}

fn bam_record(name: &[u8], pos: u32, flag: u32, cigar: &[u32]) -> Vec<u32> {
    bam_record_with_seq(name, pos, flag, cigar, 0, &[])
}

fn bam_record_with_seq(
    name: &[u8],
    pos: u32,
    flag: u32,
    cigar: &[u32],
    seq_len: u32,
    seq: &[u8],
) -> Vec<u32> {
    let bytes_len = 9 * 4 + name.len() + cigar.len() * 4 + 16;
    let mut bytes = vec![0u8; bytes_len];
    let write_u32 = |bytes: &mut [u8], idx: usize, value: u32| {
        bytes[idx * 4..idx * 4 + 4].copy_from_slice(&value.to_ne_bytes());
    };
    write_u32(&mut bytes, 2, pos);
    write_u32(&mut bytes, 3, name.len() as u32);
    write_u32(&mut bytes, 4, (flag << 16) | cigar.len() as u32);
    write_u32(&mut bytes, 5, seq_len);
    bytes[9 * 4..9 * 4 + name.len()].copy_from_slice(name);
    let cigar_start = 9 * 4 + name.len();
    for (i, value) in cigar.iter().enumerate() {
        bytes[cigar_start + i * 4..cigar_start + i * 4 + 4].copy_from_slice(&value.to_ne_bytes());
    }
    let seq_start = cigar_start + cigar.len() * 4;
    bytes[seq_start..seq_start + seq.len()].copy_from_slice(seq);
    let words = (bytes.len() + 3) / 4;
    let mut out = vec![0u32; words];
    for i in 0..words {
        let mut word = [0u8; 4];
        let start = i * 4;
        let end = (start + 4).min(bytes.len());
        word[..end - start].copy_from_slice(&bytes[start..end]);
        out[i] = u32::from_ne_bytes(word);
    }
    out
}

fn bam_record_with_seq_aux(
    name: &[u8],
    pos: u32,
    flag: u32,
    cigar: &[u32],
    seq_len: u32,
    seq: &[u8],
    nh: i32,
    as_score: Option<i32>,
) -> Vec<u32> {
    let mut aux = Vec::new();
    aux.extend_from_slice(b"NHi");
    aux.extend_from_slice(&nh.to_ne_bytes());
    if let Some(as_score) = as_score {
        aux.extend_from_slice(b"ASi");
        aux.extend_from_slice(&as_score.to_ne_bytes());
    }
    let bytes_len = 9 * 4 + name.len() + cigar.len() * 4 + seq.len() + seq_len as usize + aux.len();
    let mut bytes = vec![0u8; bytes_len];
    let write_u32 = |bytes: &mut [u8], idx: usize, value: u32| {
        bytes[idx * 4..idx * 4 + 4].copy_from_slice(&value.to_ne_bytes());
    };
    write_u32(&mut bytes, 2, pos);
    write_u32(&mut bytes, 3, name.len() as u32);
    write_u32(&mut bytes, 4, (flag << 16) | cigar.len() as u32);
    write_u32(&mut bytes, 5, seq_len);
    bytes[9 * 4..9 * 4 + name.len()].copy_from_slice(name);
    let cigar_start = 9 * 4 + name.len();
    for (i, value) in cigar.iter().enumerate() {
        bytes[cigar_start + i * 4..cigar_start + i * 4 + 4].copy_from_slice(&value.to_ne_bytes());
    }
    let seq_start = cigar_start + cigar.len() * 4;
    bytes[seq_start..seq_start + seq.len()].copy_from_slice(seq);
    let qual_start = seq_start + seq_len.div_ceil(2) as usize;
    let aux_start = qual_start + seq_len as usize;
    bytes[aux_start..aux_start + aux.len()].copy_from_slice(&aux);
    let words = (bytes.len() + 3) / 4;
    let mut out = vec![0u32; words];
    for i in 0..words {
        let mut word = [0u8; 4];
        let start = i * 4;
        let end = (start + 4).min(bytes.len());
        word[..end - start].copy_from_slice(&bytes[start..end]);
        out[i] = u32::from_ne_bytes(word);
    }
    out
}

#[test]
fn string_substitute_all_matches_cpp_find_replace_progression() {
    let mut value = "xx".to_string();
    stringsubstituteall_l3_stringsubstituteall(&mut value, "x", "yx");
    assert_eq!(value, "yxyx");

    stringsubstituteall_l3_stringsubstituteall(&mut value, "", "z");
    assert_eq!(value, "yxyx");
}

#[test]
fn warning_message_matches_star_prefix_and_newline() {
    assert_eq!(
        errorwarning_l25_warningmessage("low disk space"),
        "!!!!! WARNING: low disk space\n"
    );
}

#[test]
fn exit_with_error_records_streams_exit_code_and_cleanup_side_effects() {
    let p = Parameters {
        run_thread_n: 2,
        ..Default::default()
    };
    let observed = errorwarning_l8_exitwitherror("bad input", true, false, 74, &p, 0);
    let expected_time = timefunctions_l14_timemonthdaytime(0);

    assert_eq!(
        observed.stream_out1,
        format!(
            "\nbad input\n{} ...... FATAL ERROR, exiting\n",
            expected_time
        )
    );
    assert_eq!(observed.stream_out2, "");
    assert_eq!(observed.error_int, 74);
    assert!(observed.thread_mutex_locked);
    assert!(observed.in_out_deleted);
}

#[test]
fn shared_memory_names_match_posix_key_formatting() {
    let shared = SharedMemory {
        key: 12345,
        ..Default::default()
    };
    assert_eq!(
        sharedmemory_l87_sharedmemory_getposixobjectkey(&shared),
        "/12345"
    );
    assert_eq!(
        sharedmemory_l94_sharedmemory_countername(&shared),
        "/shared_use_counter12345"
    );
}

#[test]
fn shared_memory_state_machine_matches_allocate_close_clean_paths() {
    let mut shared = sharedmemory_l22_sharedmemory_sharedmemory(100, true);
    assert_eq!(shared.key, 100);
    assert_eq!(shared.counter_key, 101);
    assert_eq!(shared.shm_id, -1);
    assert_eq!(shared.shared_counter_id, 101);
    assert!(shared.counter_mem_attached);
    assert!(shared.needs_allocation);

    sharedmemory_l68_sharedmemory_allocate(&mut shared, 4096).unwrap();
    assert_eq!(shared.shm_id, 100);
    assert_eq!(shared.length, 4096 + std::mem::size_of::<usize>());
    assert!(shared.mapped);
    assert!(shared.is_allocator);
    assert!(!shared.needs_allocation);
    assert_eq!(
        sharedmemory_l157_sharedmemory_getsharedobjectinfo(&shared).unwrap(),
        4096 + std::mem::size_of::<usize>()
    );

    let already = sharedmemory_l68_sharedmemory_allocate(&mut shared, 1).unwrap_err();
    assert_eq!(already.error_code, 3);

    sharedmemory_l191_sharedmemory_close(&mut shared).unwrap();
    assert!(!shared.mapped);
    assert_eq!(shared.shm_id, -1);
    assert!(sharedmemory_l219_sharedmemory_unlink(&mut shared).unwrap());
    assert!(shared.needs_allocation);
    assert!(!sharedmemory_l219_sharedmemory_unlink(&mut shared).unwrap());

    shared.shm_id = 100;
    shared.needs_allocation = false;
    shared.shared_counter_id = 101;
    shared.counter_mem_attached = true;
    sharedmemory_l237_sharedmemory_clean(&mut shared).unwrap();
    assert_eq!(shared.clean_count, 1);
    assert_eq!(shared.shared_counter_id, -1);
    assert!(!shared.counter_mem_attached);
}

#[test]
fn shared_memory_destructor_unloads_only_last_attached_job() {
    let mut shared = sharedmemory_l22_sharedmemory_sharedmemory(200, true);
    sharedmemory_l68_sharedmemory_allocate(&mut shared, 10).unwrap();
    shared.shared_objects_use_count_value = 3;
    let log = sharedmemory_l37_sharedmemory_sharedmemory(&mut shared).unwrap();
    assert_eq!(
        log,
        "2 other job(s) are attached to the shared memory segment, will not remove it.\n"
    );
    assert_eq!(shared.clean_count, 0);

    let mut last = sharedmemory_l22_sharedmemory_sharedmemory(201, true);
    sharedmemory_l68_sharedmemory_allocate(&mut last, 10).unwrap();
    last.shared_objects_use_count_value = 1;
    let log = sharedmemory_l37_sharedmemory_sharedmemory(&mut last).unwrap();
    assert_eq!(
        log,
        "No other jobs are attached to the shared memory segment, removing it.\n"
    );
    assert_eq!(last.clean_count, 1);

    let mut keep = sharedmemory_l22_sharedmemory_sharedmemory(202, false);
    sharedmemory_l68_sharedmemory_allocate(&mut keep, 10).unwrap();
    assert_eq!(
        sharedmemory_l37_sharedmemory_sharedmemory(&mut keep).unwrap(),
        ""
    );
    assert_eq!(keep.clean_count, 0);
}

#[test]
fn binary_search2_matches_duplicate_x_scan_behavior() {
    let x = [1, 2, 2, 2, 5, 8];
    let y = [7, 9, 4, 3, 1, 0];
    assert_eq!(binarysearch2_l3_binarysearch2(2, 4, &x, &y, 6), 2);
    assert_eq!(binarysearch2_l3_binarysearch2(2, 8, &x, &y, 6), -1);
    assert_eq!(binarysearch2_l3_binarysearch2(4, 1, &x, &y, 6), -1);
    assert_eq!(binarysearch2_l3_binarysearch2(9, 1, &x, &y, 6), -1);
    assert_eq!(
        binarysearch2_l3_binarysearch2(2, 4, &x[..3], &y[..3], 99),
        2
    );
    assert_eq!(binarysearch2_l3_binarysearch2(2, 4, &x, &y, -1), -1);
}

#[test]
fn blocks_overlap_matches_read_and_genome_diagonal_rules() {
    let t1 = Transcript {
        n_exons: 2,
        exons: tr_exons_arr(&[[0, 100, 10, 0, 0], [20, 120, 10, 0, 0]]),
        ..Default::default()
    };
    let t2 = Transcript {
        n_exons: 2,
        exons: tr_exons_arr(&[[5, 105, 10, 0, 0], [20, 200, 10, 0, 0]]),
        ..Default::default()
    };
    assert_eq!(blocksoverlap_l3_blocksoverlap(&t1, &t2), 5);
}

#[test]
fn bam_duplicate_name_and_cigar_primitives_match_layout_arithmetic() {
    let a = bam_record(b"readA", 100, 0, &[(2 << 4) | 4, (5 << 4), (3 << 4) | 4]);
    let b = bam_record(b"readB", 110, 0x80, &[(10 << 4)]);
    assert_eq!(bamremoveduplicates_l13_funcomparenames(&a, &b), -1);
    assert_eq!(bamremoveduplicates_l34_funstartextends(&a), 98);
    assert_eq!(bamremoveduplicates_l34_funstartextends(&b), 110);

    let mut cout = [0u32; 100];
    let n = bamremoveduplicates_l43_funcigarextends(&a, &mut cout);
    assert_eq!(n, 1);
    assert_eq!(cout[0], 10 << 4);

    assert_eq!(bamremoveduplicates_l61_funcomparecigarsextends(&a, &b), 0);
}

#[test]
fn bam_duplicate_coord_flag_cigar_seq_comparator_matches_original_order() {
    let r1a = bam_record(b"readA", 100, 0, &[(5 << 4)]);
    let r1b = bam_record(b"readB", 100, 0, &[(5 << 4)]);
    let r2a = bam_record_with_seq(b"mateA", 300, 0, &[(6 << 4)], 6, &[0x12, 0x34, 0x56]);
    let r2b = bam_record_with_seq(b"mateB", 300, 0, &[(6 << 4)], 6, &[0x12, 0x35, 0x56]);

    assert_eq!(
        bamremoveduplicates_l72_funcomparecoordflagcigarseq(&r1a, &r2a, &r1b, &r2b, 6),
        -1
    );

    let shifted = bam_record(b"readC", 101, 0, &[(5 << 4)]);
    assert_eq!(
        bamremoveduplicates_l72_funcomparecoordflagcigarseq(&shifted, &r2a, &r1b, &r2b, 6),
        1
    );

    let r2_rev_a = bam_record_with_seq(
        b"mateC",
        300,
        0x10,
        &[(6 << 4)],
        7,
        &[0x91, 0x23, 0x45, 0x67],
    );
    let r2_rev_b = bam_record_with_seq(
        b"mateD",
        300,
        0x10,
        &[(6 << 4)],
        7,
        &[0x92, 0x23, 0x45, 0x67],
    );
    assert_eq!(
        bamremoveduplicates_l72_funcomparecoordflagcigarseq(&r1a, &r2_rev_a, &r1b, &r2_rev_b, 6),
        -1
    );
}

#[test]
fn bam_remove_duplicates_marks_candidates_and_unmarks_best_scoring_pair() {
    let cigar = [(6 << 4)];
    let seq = [0x12, 0x34, 0x56];
    let mut records = vec![
        bam_record_with_seq_aux(b"readA", 100, 0, &cigar, 6, &seq, 1, Some(20)),
        bam_record_with_seq_aux(b"readA", 300, 0x80, &cigar, 6, &seq, 1, Some(20)),
        bam_record_with_seq_aux(b"readB", 100, 0, &cigar, 6, &seq, 1, Some(25)),
        bam_record_with_seq_aux(b"readB", 300, 0x80, &cigar, 6, &seq, 1, Some(25)),
        bam_record_with_seq_aux(b"readC", 500, 0, &cigar, 6, &seq, 2, Some(5)),
    ];

    bamremoveduplicates_l114_bamremoveduplicates(&mut records, 6, true).unwrap();
    let dup_word = 0x400 << 16;
    assert_ne!(records[0][4] & dup_word, 0);
    assert_ne!(records[1][4] & dup_word, 0);
    assert_eq!(records[2][4] & dup_word, 0);
    assert_eq!(records[3][4] & dup_word, 0);
    assert_ne!(records[4][4] & dup_word, 0);

    let mut no_mark_multi = vec![bam_record_with_seq_aux(
        b"readD",
        500,
        0,
        &cigar,
        6,
        &seq,
        2,
        Some(5),
    )];
    bamremoveduplicates_l114_bamremoveduplicates(&mut no_mark_multi, 6, false).unwrap();
    assert_eq!(no_mark_multi[0][4] & dup_word, 0);
}

#[test]
fn bam_remove_duplicates_reports_missing_required_tags() {
    let cigar = [(6 << 4)];
    let seq = [0x12, 0x34, 0x56];
    let mut missing_as = vec![
        bam_record_with_seq_aux(b"readA", 100, 0, &cigar, 6, &seq, 1, None),
        bam_record_with_seq_aux(b"readA", 300, 0x80, &cigar, 6, &seq, 1, None),
    ];
    assert!(
        bamremoveduplicates_l114_bamremoveduplicates(&mut missing_as, 6, false)
            .unwrap_err()
            .contains("SAM tag AS is missing")
    );
}

#[test]
fn bam_remove_duplicates_reports_malformed_records_without_panic() {
    let mut short = vec![vec![0u32; 4]];
    assert!(
        bamremoveduplicates_l114_bamremoveduplicates(&mut short, 6, false)
            .unwrap_err()
            .contains("malformed BAM record")
    );

    let mut no_cigar = vec![bam_record_with_seq_aux(
        b"readA",
        100,
        0,
        &[],
        6,
        &[0x12, 0x34, 0x56],
        1,
        Some(20),
    )];
    assert!(
        bamremoveduplicates_l114_bamremoveduplicates(&mut no_cigar, 6, false)
            .unwrap_err()
            .contains("malformed BAM CIGAR")
    );

    let cigar = [(3 << 4)];
    let seq = [0x12, 0x34];
    let mut too_short = vec![bam_record_with_seq_aux(
        b"readA",
        100,
        0,
        &cigar,
        3,
        &seq,
        1,
        Some(20),
    )];
    assert!(
        bamremoveduplicates_l114_bamremoveduplicates(&mut too_short, 6, false)
            .unwrap_err()
            .contains("shorter than bamRemoveDuplicatesMate2basesN")
    );
}

#[test]
fn transcript_state_primitives_match_original_methods() {
    let mut tr = transcript_l3_transcript_transcript();
    tr.max_score = 7;
    tr.n_match = 11;
    tr.n_mm = 2;
    tr.n_gap = 3;
    tr.l_gap = 13;
    tr.l_del = 5;
    tr.n_del = 1;
    tr.l_ins = 4;
    tr.n_ins = 2;
    tr.n_unique = 6;
    tr.n_anchor = 8;
    tr.primary_flag = true;
    transcript_l8_transcript_reset(&mut tr);
    assert_eq!(tr.max_score, 0);
    assert_eq!(tr.n_match, 0);
    assert!(!tr.primary_flag);

    let tr_in = Transcript {
        max_score: 7,
        n_match: 11,
        n_mm: 2,
        n_gap: 3,
        l_gap: 13,
        l_del: 5,
        n_del: 1,
        l_ins: 4,
        n_ins: 2,
        n_unique: 6,
        n_anchor: 8,
        ..Default::default()
    };
    transcript_l28_transcript_add(&mut tr, &tr_in);
    assert_eq!(tr.max_score, 7);
    assert_eq!(tr.n_unique, 6);
    assert_eq!(tr.n_anchor, 0);

    let tr_sj = Transcript {
        n_exons: 3,
        exons: tr_exons_arr(&[[0, 100, 10, 0, 0], [20, 140, 15, 0, 0], [40, 200, 10, 0, 0]]),
        canon_sj: tr_canon_sj_arr(&[1, -1]),
        sj_annot: tr_sj_u8_arr(&[0, 1]),
        ..Default::default()
    };
    let mut sj = Vec::new();
    let mut annot_yes = false;
    transcript_l38_transcript_extractsplicejunctions(&tr_sj, &mut sj, &mut annot_yes);
    assert_eq!(sj, vec![[110, 30]]);
    assert!(!annot_yes);

    let tr_span = Transcript {
        n_exons: 2,
        exons: tr_exons_arr(&[[5, 100, 10, 0, 0], [25, 160, 20, 0, 0]]),
        c_start: 1000,
        l_read: 60,
        ..Default::default()
    };
    assert_eq!(
        transcript_l53_transcript_chrstartlengthextended(&tr_span),
        (995u64 << 32) | 100
    );
}

#[test]
fn extend_align_matches_scoring_boundaries_and_end_to_end_modes() {
    let read = [0u8, 1, 2, 3, 0, MARK_FRAG_SPACER_BASE, 2];
    let genome = [0u8, 1, 3, 3, 0, 2, 5];
    let mut tr = Transcript::default();

    assert!(extendalign_l6_extendalign(
        &read, &genome, 0, 0, 1, 1, 5, 0, 0, 3, 0.6, false, &mut tr,
    ));
    assert_eq!(tr.extend_l, 5);
    assert_eq!(tr.max_score, 3);
    assert_eq!(tr.n_match, 4);
    assert_eq!(tr.n_mm, 1);

    let mut e2e = Transcript::default();
    assert!(extendalign_l6_extendalign(
        &read, &genome, 0, 0, 1, 1, 4, 0, 0, 3, 0.0, true, &mut e2e,
    ));
    assert_eq!(e2e.extend_l, 4);
    assert_eq!(e2e.max_score, 2);
    assert_eq!(e2e.n_match, 3);
    assert_eq!(e2e.n_mm, 1);

    let mut boundary = Transcript::default();
    assert!(extendalign_l6_extendalign(
        &read,
        &genome,
        0,
        6,
        1,
        1,
        2,
        0,
        0,
        3,
        0.0,
        true,
        &mut boundary,
    ));
    assert_eq!(boundary.extend_l, 0);
    assert_eq!(boundary.max_score, -999_999_999);
    assert_eq!(boundary.n_mm, 4);
}

#[test]
fn stitch_gap_indel_scores_deletion_position_and_mismatches() {
    let read = [9u8, 0, 1, 2, 9];
    let genome = [9u8, 0, 0, 3, 1, 2, 9, 9, 9];
    let mut i_rbest = 99;
    let mut n_mm = 99;

    assert_eq!(
        stitchgapindel_l4_stitchgapindel(
            0,
            0,
            4,
            7,
            5,
            100,
            100,
            &read,
            &genome,
            -2,
            -2,
            &mut i_rbest,
            &mut n_mm,
        ),
        14
    );
    assert_eq!(i_rbest, 1);
    assert_eq!(n_mm, 0);

    let genome_mm = [9u8, 0, 0, 3, 1, 0, 9, 9, 9];
    assert_eq!(
        stitchgapindel_l4_stitchgapindel(
            0,
            0,
            4,
            7,
            5,
            100,
            100,
            &read,
            &genome_mm,
            -2,
            -2,
            &mut i_rbest,
            &mut n_mm,
        ),
        12
    );
    assert_eq!(i_rbest, 1);
    assert_eq!(n_mm, 1);

    assert_eq!(
        stitchgapindel_l4_stitchgapindel(
            0,
            0,
            4,
            4,
            5,
            100,
            100,
            &read,
            &genome,
            -2,
            -2,
            &mut i_rbest,
            &mut n_mm,
        ),
        -1
    );
}

#[test]
fn stitch_align_to_transcript_uses_annotated_sjdb_fast_path() {
    let mut tr = Transcript {
        n_exons: 1,
        exons: tr_exons_arr(&[[0, 10, 3, 0, 0], [0; EX_SIZE]]),
        canon_sj: tr_canon_sj_arr(&[0; 2]),
        sj_annot: tr_sj_u8_arr(&[0; 2]),
        shift_sj: tr_shift_sj_arr(&[[0; 2]; 2]),
        sj_str: tr_sj_u8_arr(&[0; 2]),
        ..Default::default()
    };
    let p = Parameters {
        p_ge: ParametersGenome {
            sjdb_score: 2,
            ..Default::default()
        },
        ..Default::default()
    };
    let genome = Genome {
        g: vec![0; 64],
        sjdb_motif: vec![1],
        sjdb_shift_left: vec![0],
        sjdb_shift_right: vec![0],
        sjdb_strand: vec![1],
        ..Default::default()
    };

    let score = stitchaligntotranscript_l9_stitchaligntotranscript(
        2, 12, 3, 20, 4, 0, 0, &p, &[0; 16], &genome, &mut tr, 10,
    );

    assert_eq!(score, 6);
    assert_eq!(tr.n_exons, 2);
    assert_eq!(tr.exons[1], [3, 20, 4, 0, 0]);
    assert_eq!(tr.canon_sj[0], 1);
    assert_eq!(tr.sj_annot[0], 1);
    assert_eq!(tr.sj_str[0], 1);
    assert_eq!(tr.n_match, 4);
}

#[test]
fn stitch_window_aligns_records_single_window_alignment() {
    let wa = vec![{
        let mut a = [0u32; WA_SIZE];
        a[WA_LENGTH] = 4;
        a[WA_R_START] = 0;
        a[WA_G_START] = 10;
        a[WA_N_REP] = 1;
        a[WA_ANCHOR] = 1;
        a[WA_I_FRAG] = 0;
        a[WA_SJ_A] = u32::MAX;
        a
    }];
    let mut incl = vec![false; 1];
    let r = vec![0, 1, 2, 3];
    let genome = Genome {
        g: vec![0; 128],
        n_genome: 128,
        chr_start: vec![0],
        chr_length: vec![128],
        ..Default::default()
    };
    let p = Parameters {
        out_filter_intron_motifs: "None".to_string(),
        align_transcripts_per_window_nmax: 4,
        out_filter_multimap_score_range: 0,
        align_soft_clip_at_reference_ends_yes: true,
        ..Default::default()
    };
    let mut ra = ReadAlign {
        out_filter_mismatch_nmax_total: 10,
        read_length: vec![4],
        max_score_mate: vec![i32::MIN],
        ..Default::default()
    };
    let mut w_tr = vec![Transcript::default(); 4];
    let mut n_win_tr = 0;

    stitchwindowaligns_l8_stitchwindowaligns(
        0,
        wa.len() as u32,
        0,
        &mut incl,
        0,
        0,
        Transcript::default(),
        4,
        &wa,
        &r,
        &genome,
        &p,
        &mut w_tr,
        &mut n_win_tr,
        &mut ra,
    )
    .unwrap();

    assert_eq!(n_win_tr, 1);
    assert_eq!(w_tr[0].n_exons, 1);
    assert_eq!(w_tr[0].exons[0][EX_R], 0);
    assert_eq!(w_tr[0].exons[0][EX_G], 10);
    assert_eq!(w_tr[0].exons[0][EX_L], 4);
    assert_eq!(w_tr[0].n_match, 4);
    assert_eq!(w_tr[0].n_unique, 1);
    assert_eq!(w_tr[0].n_anchor, 1);
    assert_eq!(w_tr[0].max_score, 4);
    assert_eq!(w_tr[0].i_frag, 0);
    assert_eq!(ra.max_score_mate[0], 4);
}

#[test]
fn stitch_window_aligns_updates_genomic_end_after_right_extension() {
    let wa = vec![{
        let mut a = [0u32; WA_SIZE];
        a[WA_LENGTH] = 3;
        a[WA_R_START] = 0;
        a[WA_G_START] = 10;
        a[WA_N_REP] = 1;
        a[WA_ANCHOR] = 1;
        a[WA_I_FRAG] = 0;
        a[WA_SJ_A] = u32::MAX;
        a
    }];
    let mut incl = vec![false; 1];
    let r = vec![0, 0, 0, 0, 0];
    let genome = Genome {
        g: vec![0; 128],
        n_genome: 128,
        chr_start: vec![0],
        chr_length: vec![128],
        ..Default::default()
    };
    let p = Parameters {
        out_filter_intron_motifs: "None".to_string(),
        align_transcripts_per_window_nmax: 4,
        align_soft_clip_at_reference_ends_yes: true,
        ..Default::default()
    };
    let mut ra = ReadAlign {
        out_filter_mismatch_nmax_total: 10,
        read_length: vec![5],
        max_score_mate: vec![i32::MIN],
        ..Default::default()
    };
    let mut w_tr = vec![Transcript::default(); 4];
    let mut n_win_tr = 0;

    stitchwindowaligns_l8_stitchwindowaligns(
        0,
        wa.len() as u32,
        0,
        &mut incl,
        0,
        0,
        Transcript::default(),
        5,
        &wa,
        &r,
        &genome,
        &p,
        &mut w_tr,
        &mut n_win_tr,
        &mut ra,
    )
    .unwrap();

    assert_eq!(n_win_tr, 1);
    assert_eq!(w_tr[0].exons[0][EX_L], 5);
    assert_eq!(w_tr[0].g_length, 5);
    assert_eq!(w_tr[0].max_score, 5);
}

#[test]
fn stitch_window_aligns_replaces_lower_scoring_subset() {
    let old = Transcript {
        n_exons: 1,
        exons: tr_exons_arr(&[[0, 10, 3, 0, 0]]),
        mapped_length: 3,
        max_score: 3,
        g_length: 3,
        ..Default::default()
    };
    let wa = vec![{
        let mut a = [0u32; WA_SIZE];
        a[WA_LENGTH] = 5;
        a[WA_R_START] = 0;
        a[WA_G_START] = 10;
        a[WA_N_REP] = 1;
        a[WA_ANCHOR] = 1;
        a[WA_I_FRAG] = 0;
        a[WA_SJ_A] = u32::MAX;
        a
    }];
    let mut incl = vec![false; 1];
    let genome = Genome {
        g: vec![0; 128],
        n_genome: 128,
        chr_start: vec![0],
        chr_length: vec![128],
        ..Default::default()
    };
    let p = Parameters {
        out_filter_intron_motifs: "None".to_string(),
        align_transcripts_per_window_nmax: 4,
        align_soft_clip_at_reference_ends_yes: true,
        ..Default::default()
    };
    let mut ra = ReadAlign {
        out_filter_mismatch_nmax_total: 10,
        read_length: vec![5],
        max_score_mate: vec![i32::MIN],
        ..Default::default()
    };
    let mut w_tr = vec![old, Transcript::default(), Transcript::default()];
    let mut n_win_tr = 1;

    stitchwindowaligns_l8_stitchwindowaligns(
        0,
        1,
        0,
        &mut incl,
        0,
        0,
        Transcript::default(),
        5,
        &wa,
        &[0; 5],
        &genome,
        &p,
        &mut w_tr,
        &mut n_win_tr,
        &mut ra,
    )
    .unwrap();

    assert_eq!(n_win_tr, 1);
    assert_eq!(w_tr[0].exons[0][EX_L], 5);
    assert_eq!(w_tr[0].max_score, 5);
}

#[test]
fn stitch_window_aligns_filters_noncanonical_introns_and_reports_bad_mode() {
    let mut tr = Transcript {
        n_exons: 2,
        exons: tr_exons_arr(&[[0, 10, 5, 0, 0], [10, 25, 5, 0, 0]]),
        canon_sj: tr_canon_sj_arr(&[0]),
        sj_annot: tr_sj_u8_arr(&[0]),
        sj_str: tr_sj_u8_arr(&[0]),
        shift_sj: tr_shift_sj_arr(&[[0; 2]]),
        r_start: 0,
        g_start: 10,
        n_match: 10,
        ..Default::default()
    };
    tr.exons[0][EX_IFRAG] = 0;
    tr.exons[1][EX_IFRAG] = 0;
    let genome = Genome {
        g: vec![0; 128],
        n_genome: 128,
        chr_start: vec![0],
        chr_length: vec![128],
        ..Default::default()
    };
    let mut ra = ReadAlign {
        out_filter_mismatch_nmax_total: 10,
        read_length: vec![15],
        max_score_mate: vec![i32::MIN],
        ..Default::default()
    };
    let mut w_tr = vec![Transcript::default(); 2];
    let mut n_win_tr = 0;
    let mut incl = Vec::new();

    stitchwindowaligns_l8_stitchwindowaligns(
        0,
        0,
        10,
        &mut incl,
        14,
        29,
        tr.clone(),
        15,
        &[],
        &[0; 15],
        &genome,
        &Parameters {
            out_filter_intron_motifs: "RemoveNoncanonical".to_string(),
            align_soft_clip_at_reference_ends_yes: true,
            ..Default::default()
        },
        &mut w_tr,
        &mut n_win_tr,
        &mut ra,
    )
    .unwrap();
    assert_eq!(n_win_tr, 0);

    let err = stitchwindowaligns_l8_stitchwindowaligns(
        0,
        0,
        10,
        &mut incl,
        14,
        29,
        tr,
        15,
        &[],
        &[0; 15],
        &genome,
        &Parameters {
            out_filter_intron_motifs: "BadMode".to_string(),
            align_soft_clip_at_reference_ends_yes: true,
            ..Default::default()
        },
        &mut w_tr,
        &mut n_win_tr,
        &mut ra,
    )
    .unwrap_err();
    assert!(err.contains("unrecognized value of --outFilterIntronMotifs=BadMode"));
}

#[test]
fn stitch_align_to_transcript_fills_same_fragment_gap_without_indel() {
    let mut tr = Transcript {
        n_exons: 1,
        exons: tr_exons_arr(&[[0, 0, 3, 0, u32::MAX], [0; EX_SIZE]]),
        canon_sj: tr_canon_sj_arr(&[0; 2]),
        sj_annot: tr_sj_u8_arr(&[0; 2]),
        shift_sj: tr_shift_sj_arr(&[[0; 2]; 2]),
        sj_str: tr_sj_u8_arr(&[0; 2]),
        ..Default::default()
    };
    let p = Parameters {
        align_sj_stitch_mismatch_nmax: vec![10; 4],
        ..Default::default()
    };
    let genome = Genome {
        g: vec![0, 1, 2, 3, 0, 1, 2, 3],
        ..Default::default()
    };
    let read = vec![0, 1, 2, 3, 0, 1, 2, 3];

    let score = stitchaligntotranscript_l9_stitchaligntotranscript(
        2,
        2,
        5,
        5,
        2,
        0,
        u32::MAX,
        &p,
        &read,
        &genome,
        &mut tr,
        10,
    );

    assert_eq!(score, -1_000_007);
    assert_eq!(tr.n_exons, 1);
    assert_eq!(tr.exons[0], [0, 0, 3, 0, u32::MAX]);
    assert_eq!(tr.n_match, 0);
    assert_eq!(tr.n_mm, 0);
}

#[test]
fn stitch_align_to_transcript_records_same_fragment_insertion() {
    let mut tr = Transcript {
        n_exons: 1,
        exons: tr_exons_arr(&[[0, 0, 3, 0, u32::MAX], [0; EX_SIZE]]),
        canon_sj: tr_canon_sj_arr(&[0; 2]),
        sj_annot: tr_sj_u8_arr(&[0; 2]),
        shift_sj: tr_shift_sj_arr(&[[0; 2]; 2]),
        sj_str: tr_sj_u8_arr(&[0; 2]),
        ..Default::default()
    };
    let p = Parameters {
        align_intron_min: 10,
        align_sj_stitch_mismatch_nmax: vec![10; 4],
        score_ins_base: -2,
        score_ins_open: -2,
        ..Default::default()
    };
    let genome = Genome {
        g: vec![0, 1, 2, 3, 0, 1, 2, 3],
        ..Default::default()
    };
    let read = vec![0, 1, 2, 9, 3, 0, 1, 2];

    let score = stitchaligntotranscript_l9_stitchaligntotranscript(
        2,
        2,
        5,
        4,
        2,
        0,
        u32::MAX,
        &p,
        &read,
        &genome,
        &mut tr,
        10,
    );

    assert_eq!(score, -1);
    assert_eq!(tr.n_exons, 2);
    assert_eq!(tr.exons[0][EX_L], 3);
    assert_eq!(tr.exons[1], [4, 3, 3, 0, u32::MAX]);
    assert_eq!(tr.canon_sj[0], -2);
    assert_eq!(tr.n_ins, 1);
    assert_eq!(tr.l_ins, 1);
    assert_eq!(tr.n_match, 3);
}

#[test]
fn read_align_stitch_window_seeds_records_single_seed_transcript() {
    let mut read_align = ReadAlign {
        l_read: 4,
        out_filter_mismatch_nmax_total: 10,
        tr_init: Box::new(Transcript::default()),
        tr_best: Transcript {
            str_: 0,
            ..Default::default()
        },
        n_wa: vec![1],
        wa: vec![vec![[4, 0, 10, 1, 1, 0, u32::MAX]]],
        tr_all: vec![Vec::new()],
        n_win_tr: vec![0],
        max_score_mate: vec![0],
        ..Default::default()
    };
    let p = Parameters {
        align_sj_overhang_min: 1,
        out_filter_intron_motifs: "None".to_string(),
        ..Default::default()
    };
    let genome = Genome {
        g: vec![0; 1000],
        ..Default::default()
    };

    assert!(
        readalign_stitchwindowseeds_l12_readalign_stitchwindowseeds(
            &mut read_align,
            0,
            0,
            None,
            &[0; 8],
            &genome,
            &p,
        )
        .unwrap()
    );

    assert_eq!(read_align.n_win_tr[0], 1);
    let tr = &read_align.tr_all[0][0];
    assert_eq!(tr.n_exons, 1);
    assert_eq!(tr.exons[0], [0, 10, 4, 0, u32::MAX]);
    assert_eq!(tr.max_score, 4);
    assert_eq!(tr.n_match, 4);
    assert_eq!(tr.r_length, 4);
    assert_eq!(tr.i_frag, 0);
    assert_eq!(read_align.max_score_mate[0], 4);
}

#[test]
fn read_align_stitch_window_seeds_records_second_best_with_exclusion_mask() {
    let mut read_align = ReadAlign {
        l_read: 3,
        out_filter_mismatch_nmax_total: 10,
        tr_init: Box::new(Transcript::default()),
        tr_best: Transcript {
            str_: 0,
            ..Default::default()
        },
        n_wa: vec![2],
        wa: vec![vec![
            [4, 0, 10, 1, 1, 0, u32::MAX],
            [3, 0, 20, 1, 1, 0, u32::MAX],
        ]],
        tr_all: vec![vec![Transcript::default()]],
        n_win_tr: vec![1],
        max_score_mate: vec![0, 0],
        ..Default::default()
    };
    let p = Parameters {
        align_sj_overhang_min: 1,
        out_filter_intron_motifs: "None".to_string(),
        ..Default::default()
    };
    let genome = Genome {
        g: vec![0; 1000],
        ..Default::default()
    };

    assert!(
        readalign_stitchwindowseeds_l12_readalign_stitchwindowseeds(
            &mut read_align,
            0,
            0,
            Some(&[true, false]),
            &[0; 1000],
            &genome,
            &p,
        )
        .unwrap()
    );

    assert_eq!(read_align.n_win_tr[0], 2);
    assert_eq!(read_align.tr_all[0][1].exons[0], [0, 20, 3, 0, u32::MAX]);
    assert!(read_align.wa_incl[1]);
    assert!(!read_align.wa_incl[0]);
    assert_eq!(read_align.max_score_mate[0], 3);
}

#[test]
fn transcript_align_score_recomputes_matches_mismatches_and_junction_scores() {
    let read0 = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
    let read2 = [3u8, 2, 1, 0, 3, 2, 1, 0, 3, 2];
    let genome = [0u8, 1, 1, 3, 0, 1, 2, 3, 2, 1, 0, 3, 2, 1, 0, 3];

    let mut tr = Transcript {
        n_exons: 3,
        exons: tr_exons_arr(&[[0, 0, 4, 0, 0], [4, 6, 3, 0, 0], [7, 10, 2, 0, 0]]),
        sj_annot: tr_sj_u8_arr(&[1, 0]),
        canon_sj: tr_canon_sj_arr(&[0, -1]),
        max_score: 99,
        n_match: 99,
        n_mm: 99,
        ..Default::default()
    };
    let score = transcript_alignscore_l4_transcript_alignscore(
        &mut tr, &read0, &read2, &genome, 10, -2, -5, -3, -7, -11, -13, -17, -19, 0.0,
    );
    assert_eq!(tr.n_match, 4);
    assert_eq!(tr.n_mm, 5);
    assert_eq!(score, -1);
    assert_eq!(tr.max_score, -1);

    let mut reverse = Transcript {
        n_exons: 1,
        exons: tr_exons_arr(&[[1, 12, 3, 0, 0]]),
        ro_str: 1,
        ..Default::default()
    };
    assert_eq!(
        transcript_alignscore_l4_transcript_alignscore(
            &mut reverse,
            &read0,
            &read2,
            &genome,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0.0,
        ),
        3
    );
    assert_eq!(reverse.n_match, 3);
    assert_eq!(reverse.n_mm, 0);
}

#[test]
fn primary_align_mark_matches_best_reordering_and_flag_rules() {
    let mut unique = vec![Transcript {
        max_score: 5,
        ..Default::default()
    }];
    funprimaryalignmark_l3_funprimaryalignmark(&mut unique, 1, 5, false, false, "", &[]);
    assert!(unique[0].primary_flag);

    let mut all_best = vec![
        Transcript {
            max_score: 7,
            c_start: 10,
            ..Default::default()
        },
        Transcript {
            max_score: 3,
            c_start: 20,
            ..Default::default()
        },
        Transcript {
            max_score: 7,
            c_start: 30,
            ..Default::default()
        },
    ];
    funprimaryalignmark_l3_funprimaryalignmark(
        &mut all_best,
        3,
        7,
        false,
        false,
        "AllBestScore",
        &[],
    );
    assert_eq!(
        all_best
            .iter()
            .map(|tr| tr.primary_flag)
            .collect::<Vec<_>>(),
        vec![true, false, true]
    );

    let mut limited = vec![
        Transcript {
            max_score: 3,
            c_start: 10,
            ..Default::default()
        },
        Transcript {
            max_score: 7,
            c_start: 20,
            ..Default::default()
        },
        Transcript {
            max_score: 7,
            c_start: 30,
            ..Default::default()
        },
    ];
    funprimaryalignmark_l3_funprimaryalignmark(&mut limited, 3, 7, false, true, "", &[]);
    assert_eq!(limited[0].c_start, 20);
    assert!(limited[0].primary_flag);

    let mut random = vec![
        Transcript {
            max_score: 9,
            c_start: 1,
            ..Default::default()
        },
        Transcript {
            max_score: 9,
            c_start: 2,
            ..Default::default()
        },
        Transcript {
            max_score: 4,
            c_start: 3,
            ..Default::default()
        },
        Transcript {
            max_score: 4,
            c_start: 4,
            ..Default::default()
        },
    ];
    funprimaryalignmark_l3_funprimaryalignmark(&mut random, 4, 9, true, false, "", &[0.0, 0.0]);
    assert_eq!(
        random.iter().map(|tr| tr.c_start).collect::<Vec<_>>(),
        vec![2, 1, 4, 3]
    );
    assert!(random[0].primary_flag);
}

#[test]
fn read_align_mult_map_select_filters_sets_coordinates_and_primary_flags() {
    let mut read_align = ReadAlign {
        n_w: 2,
        n_wap: vec![2, 2],
        n_win_tr: vec![2, 2],
        l_read: 100,
        tr_best: Transcript {
            max_score: 40,
            ..Default::default()
        },
        ..Default::default()
    };
    let genome = Genome {
        chr_start: vec![1000, 5000],
        ..Default::default()
    };
    let tr_all = vec![
        vec![
            Transcript {
                max_score: 40,
                chr: 1,
                str_: 1,
                ro_str: 1,
                r_start: 10,
                r_length: 20,
                g_start: 5020,
                ..Default::default()
            },
            Transcript {
                max_score: 35,
                r_start: 5,
                r_length: 10,
                g_start: 5015,
                ..Default::default()
            },
        ],
        vec![
            Transcript {
                max_score: 39,
                chr: 0,
                str_: 0,
                ro_str: 0,
                r_start: 7,
                r_length: 12,
                g_start: 1033,
                ..Default::default()
            },
            Transcript {
                max_score: 20,
                ..Default::default()
            },
        ],
    ];

    let selected = readalign_multmapselect_l8_readalign_multmapselect(
        &mut read_align,
        &genome,
        &tr_all,
        1,
        10,
        false,
        false,
        "AllBestScore",
        &[],
    )
    .unwrap();

    assert_eq!(read_align.n_tr, 2);
    assert_eq!(
        selected
            .iter()
            .map(|tr| (
                tr.max_score,
                tr.chr,
                tr.str_,
                tr.ro_str,
                tr.ro_start,
                tr.c_start
            ))
            .collect::<Vec<_>>(),
        vec![(40, 1, 1, 1, 70, 20), (39, 0, 0, 0, 7, 33)]
    );
    assert_eq!(
        selected
            .iter()
            .map(|tr| tr.primary_flag)
            .collect::<Vec<_>>(),
        vec![true, false]
    );
}

#[test]
fn read_align_mult_map_select_preserves_old_primary_and_limit_paths() {
    let genome = Genome {
        chr_start: vec![0],
        ..Default::default()
    };
    let tr_all = vec![vec![
        Transcript {
            max_score: 10,
            chr: 0,
            g_start: 3,
            ..Default::default()
        },
        Transcript {
            max_score: 9,
            chr: 0,
            g_start: 5,
            ..Default::default()
        },
    ]];
    let mut old_way = ReadAlign {
        n_w: 1,
        n_wap: vec![2],
        n_win_tr: vec![2],
        l_read: 20,
        tr_best: Transcript {
            max_score: 10,
            ..Default::default()
        },
        ..Default::default()
    };
    let selected = readalign_multmapselect_l8_readalign_multmapselect(
        &mut old_way,
        &genome,
        &tr_all,
        1,
        10,
        false,
        false,
        "",
        &[],
    )
    .unwrap();
    assert_eq!(selected.iter().filter(|tr| tr.primary_flag).count(), 1);
    assert!(old_way.tr_best.primary_flag);

    let mut too_multi = ReadAlign {
        n_w: 1,
        n_wap: vec![2],
        n_win_tr: vec![2],
        l_read: 20,
        tr_best: Transcript {
            max_score: 10,
            ..Default::default()
        },
        ..Default::default()
    };
    let limited = readalign_multmapselect_l8_readalign_multmapselect(
        &mut too_multi,
        &genome,
        &tr_all,
        1,
        1,
        false,
        false,
        "",
        &[],
    )
    .unwrap();
    assert_eq!(too_multi.n_tr, 2);
    assert_eq!(limited.len(), 2);
    assert!(limited.iter().all(|tr| !tr.primary_flag));
    assert!(!too_multi.tr_best.primary_flag);
}

#[test]
fn read_align_mult_map_select_orders_compact_nonbest_repeat_ties() {
    let genome = Genome {
        chr_start: vec![0],
        ..Default::default()
    };
    let tr_all = vec![vec![
        Transcript {
            max_score: 35,
            chr: 0,
            r_start: 0,
            r_length: 35,
            g_start: 223130,
            n_mm: 0,
            ..Default::default()
        },
        Transcript {
            max_score: 33,
            chr: 0,
            r_start: 0,
            r_length: 34,
            g_start: 223128,
            n_mm: 0,
            ..Default::default()
        },
        Transcript {
            max_score: 33,
            chr: 0,
            r_start: 0,
            r_length: 34,
            g_start: 223129,
            n_mm: 0,
            ..Default::default()
        },
    ]];
    let mut read_align = ReadAlign {
        n_w: 1,
        n_win_tr: vec![3],
        l_read: 51,
        tr_best: tr_all[0][0].clone(),
        ..Default::default()
    };

    let selected = readalign_multmapselect_l8_readalign_multmapselect(
        &mut read_align,
        &genome,
        &tr_all,
        2,
        20,
        false,
        false,
        "",
        &[],
    )
    .unwrap();

    assert_eq!(
        selected
            .iter()
            .map(|tr| (tr.max_score, tr.c_start))
            .collect::<Vec<_>>(),
        vec![(35, 223130), (33, 223128), (33, 223129)]
    );

    let tr_all_r_start_tie = vec![vec![
        Transcript {
            max_score: 38,
            chr: 0,
            r_start: 4,
            r_length: 43,
            g_start: 223122,
            n_mm: 2,
            ..Default::default()
        },
        Transcript {
            max_score: 37,
            chr: 0,
            r_start: 0,
            r_length: 46,
            g_start: 223128,
            n_mm: 4,
            ..Default::default()
        },
        Transcript {
            max_score: 37,
            chr: 0,
            r_start: 0,
            r_length: 46,
            g_start: 223129,
            n_mm: 4,
            ..Default::default()
        },
    ]];
    let mut mismatch_read_align = ReadAlign {
        n_w: 1,
        n_win_tr: vec![3],
        l_read: 51,
        tr_best: tr_all_r_start_tie[0][0].clone(),
        ..Default::default()
    };

    let mismatch_selected = readalign_multmapselect_l8_readalign_multmapselect(
        &mut mismatch_read_align,
        &genome,
        &tr_all_r_start_tie,
        1,
        20,
        false,
        false,
        "",
        &[],
    )
    .unwrap();

    assert_eq!(
        mismatch_selected
            .iter()
            .map(|tr| (tr.max_score, tr.n_mm, tr.c_start))
            .collect::<Vec<_>>(),
        vec![(38, 2, 223122), (37, 4, 223128), (37, 4, 223129)]
    );
}

#[test]
fn read_align_mult_map_select_orders_nonbest_softclip_before_insertion_tie() {
    let genome = Genome {
        chr_start: vec![0],
        ..Default::default()
    };
    let tr_all = vec![vec![
        Transcript {
            max_score: 33,
            chr: 0,
            r_start: 11,
            r_length: 38,
            g_start: 223125,
            n_mm: 2,
            ..Default::default()
        },
        Transcript {
            max_score: 32,
            chr: 0,
            r_start: 14,
            r_length: 37,
            g_start: 223127,
            n_mm: 2,
            ..Default::default()
        },
        Transcript {
            max_score: 32,
            chr: 0,
            r_start: 11,
            r_length: 40,
            g_start: 223125,
            n_mm: 1,
            ..Default::default()
        },
    ]];
    let mut read_align = ReadAlign {
        n_w: 1,
        n_win_tr: vec![3],
        l_read: 51,
        tr_best: tr_all[0][0].clone(),
        ..Default::default()
    };

    let selected = readalign_multmapselect_l8_readalign_multmapselect(
        &mut read_align,
        &genome,
        &tr_all,
        1,
        20,
        false,
        false,
        "",
        &[],
    )
    .unwrap();

    assert_eq!(
        selected
            .iter()
            .map(|tr| (tr.max_score, tr.r_start, tr.c_start, tr.n_mm))
            .collect::<Vec<_>>(),
        vec![
            (33, 11, 223125, 2),
            (32, 14, 223127, 2),
            (32, 11, 223125, 1)
        ]
    );
}

#[test]
fn read_align_mult_map_select_marks_nonfirst_tr_best_primary_after_coordinate_fill() {
    let genome = Genome {
        chr_start: vec![0],
        ..Default::default()
    };
    let best = Transcript {
        max_score: 38,
        chr: 0,
        str_: 0,
        ro_str: 0,
        r_start: 1,
        ro_start: 1,
        r_length: 45,
        mapped_length: 45,
        g_start: 189665,
        g_length: 45,
        n_match: 42,
        n_mm: 3,
        n_unique: 1,
        n_anchor: 1,
        n_exons: 1,
        exons: tr_exons_arr(&[[1, 189665, 45, 0, u32::MAX]]),
        canon_sj: tr_canon_sj_arr(&[-1]),
        shift_sj: tr_shift_sj_arr(&[[0, 0]]),
        sj_str: tr_sj_u8_arr(&[0]),
        sj_annot: tr_sj_u8_arr(&[0]),
        read_name: "@SRR10143877.15951".to_string(),
        ..Default::default()
    };
    let tr_all = vec![
        vec![Transcript {
            max_score: 38,
            chr: 0,
            str_: 1,
            ro_str: 1,
            r_start: 1,
            ro_start: 1,
            r_length: 49,
            mapped_length: 49,
            g_start: 165869,
            g_length: 50,
            n_match: 46,
            n_mm: 3,
            n_exons: 2,
            exons: tr_exons_arr(&[[1, 165869, 29, 0, u32::MAX], [30, 165899, 20, 0, u32::MAX]]),
            canon_sj: tr_canon_sj_arr(&[-1, 0]),
            shift_sj: tr_shift_sj_arr(&[[0, 1], [0, 0]]),
            sj_str: tr_sj_u8_arr(&[0, 0]),
            sj_annot: tr_sj_u8_arr(&[0, 0]),
            ..Default::default()
        }],
        vec![best.clone()],
    ];
    let mut read_align = ReadAlign {
        n_w: 2,
        n_win_tr: vec![1, 1],
        l_read: 51,
        tr_best: Transcript {
            c_start: 0,
            primary_flag: false,
            ..best
        },
        ..Default::default()
    };

    let selected = readalign_multmapselect_l8_readalign_multmapselect(
        &mut read_align,
        &genome,
        &tr_all,
        0,
        20,
        false,
        false,
        "",
        &[],
    )
    .unwrap();

    assert_eq!(read_align.tr_best.c_start, 189665);
    assert_eq!(
        selected
            .iter()
            .map(|tr| tr.primary_flag)
            .collect::<Vec<_>>(),
        vec![false, true]
    );
}

#[test]
fn read_align_mult_map_select_uses_recorded_transcripts_not_anchor_counts() {
    let genome = Genome {
        chr_start: vec![0],
        ..Default::default()
    };
    let tr_all = vec![vec![Transcript {
        max_score: 12,
        chr: 0,
        r_length: 10,
        g_start: 4,
        ..Default::default()
    }]];
    let mut read_align = ReadAlign {
        n_w: 1,
        n_wap: vec![0],
        n_win_tr: vec![1],
        l_read: 20,
        tr_best: Transcript {
            max_score: 12,
            ..Default::default()
        },
        ..Default::default()
    };

    let selected = readalign_multmapselect_l8_readalign_multmapselect(
        &mut read_align,
        &genome,
        &tr_all,
        0,
        10,
        false,
        false,
        "",
        &[],
    )
    .unwrap();

    assert_eq!(read_align.n_tr, 1);
    assert_eq!(selected.len(), 1);
    assert!(selected[0].primary_flag);
}

#[test]
fn read_align_mult_map_select_recovers_zero_anchor_singleton_best_window() {
    let genome = Genome {
        chr_start: vec![0],
        ..Default::default()
    };
    let tr_all = vec![vec![
        Transcript {
            max_score: 31,
            chr: 0,
            r_length: 34,
            g_start: 8,
            ..Default::default()
        },
        Transcript {
            max_score: 20,
            chr: 0,
            r_length: 22,
            g_start: 40,
            ..Default::default()
        },
        Transcript {
            max_score: 12,
            chr: 0,
            r_length: 18,
            g_start: 80,
            ..Default::default()
        },
    ]];
    let mut read_align = ReadAlign {
        n_w: 1,
        n_wap: vec![0],
        n_win_tr: vec![3],
        l_read: 36,
        tr_best: Transcript {
            max_score: 31,
            ..Default::default()
        },
        ..Default::default()
    };

    let selected = readalign_multmapselect_l8_readalign_multmapselect(
        &mut read_align,
        &genome,
        &tr_all,
        0,
        20,
        false,
        false,
        "",
        &[],
    )
    .unwrap();

    assert_eq!(read_align.n_tr, 1);
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].max_score, 31);
    assert!(selected[0].primary_flag);
}

#[test]
fn quantifications_constructor_and_add_quants_match_gene_count_layout() {
    let mut q = quantifications_l3_quantifications_quantifications(4);
    assert_eq!(q.gene_counts.n_ge, 4);
    assert_eq!(q.gene_counts.n_type, 3);
    assert_eq!(q.gene_counts.c_multi, 0);
    assert_eq!(q.gene_counts.c_ambig, vec![0, 0, 0]);
    assert_eq!(q.gene_counts.c_none, vec![0, 0, 0]);
    assert_eq!(q.gene_counts.g_count, vec![vec![0; 4]; 3]);

    q.gene_counts.c_multi = 10;
    q.gene_counts.c_ambig = vec![1, 2, 3];
    q.gene_counts.c_none = vec![4, 5, 6];
    q.gene_counts.g_count = vec![vec![1, 0, 2, 0], vec![0, 3, 0, 4], vec![5, 6, 7, 8]];

    let mut incoming = quantifications_l3_quantifications_quantifications(4);
    incoming.gene_counts.c_multi = 7;
    incoming.gene_counts.c_ambig = vec![10, 20, 30];
    incoming.gene_counts.c_none = vec![40, 50, 60];
    incoming.gene_counts.g_count = vec![
        vec![10, 11, 12, 13],
        vec![20, 21, 22, 23],
        vec![30, 31, 32, 33],
    ];

    quantifications_l25_quantifications_addquants(&mut q, &incoming);

    assert_eq!(q.gene_counts.c_multi, 17);
    assert_eq!(q.gene_counts.c_ambig, vec![11, 22, 33]);
    assert_eq!(q.gene_counts.c_none, vec![44, 55, 66]);
    assert_eq!(
        q.gene_counts.g_count,
        vec![
            vec![11, 11, 14, 13],
            vec![20, 24, 22, 27],
            vec![35, 37, 39, 41],
        ]
    );
}

#[test]
fn log_multinomial_pdf_sparse_matches_star_indexing() {
    let amb_profile_log_p = [0.0f64, -0.2, -0.7, -1.4];
    let count_cell_gene_umi = [
        1u32, 2, 99, //
        2, 3, 99, //
        3, 1, 99, //
    ];
    let log_factorial = [
        0.0,
        0.0,
        (2.0f64).ln(),
        (6.0f64).ln(),
        (24.0f64).ln(),
        (120.0f64).ln(),
        (720.0f64).ln(),
    ];

    let observed = solofeature_emptydrops_cr_l219_logmultinomialpdfsparse(
        &amb_profile_log_p,
        &count_cell_gene_umi,
        3,
        1,
        0,
        3,
        &log_factorial,
    );
    let expected = log_factorial[6] - (log_factorial[2] + log_factorial[3] + log_factorial[1])
        + amb_profile_log_p[1] * 2.0
        + amb_profile_log_p[2] * 3.0
        + amb_profile_log_p[3] * 1.0;
    assert!((observed - expected).abs() < 1e-12);
}

#[test]
fn stats_state_primitives_match_original_methods() {
    let mut stats = stats_l16_stats_stats();
    assert_eq!(stats.read_n, 0);
    assert_eq!(stats.time_last_report, 0);

    stats.read_n = 5;
    stats.read_bases = 500;
    stats.mapped_reads_u = 3;
    stats.mapped_reads_m = 2;
    stats.mapped_bases = 100;
    stats.mapped_mismatches_n = 1;
    stats.mapped_ins_n = 2;
    stats.mapped_del_n = 3;
    stats.mapped_ins_l = 4;
    stats.mapped_del_l = 5;
    stats.mapped_portion = 1.5;
    stats.splices_n = [1, 2, 3, 4, 5, 6, 7];
    stats.splices_nsjdb = 8;
    stats.unmapped_other = 9;
    stats.unmapped_short = 10;
    stats.unmapped_mismatch = 11;
    stats.unmapped_multi = 12;
    stats.unmapped_all = 13;
    stats.chimeric_all = 14;

    let add = Stats {
        read_n: 1,
        read_bases: 2,
        mapped_reads_u: 3,
        mapped_reads_m: 4,
        mapped_bases: 5,
        mapped_mismatches_n: 6,
        mapped_ins_n: 7,
        mapped_del_n: 8,
        mapped_ins_l: 9,
        mapped_del_l: 10,
        mapped_portion: 2.5,
        splices_n: [7, 6, 5, 4, 3, 2, 1],
        splices_nsjdb: 2,
        unmapped_other: 3,
        unmapped_short: 4,
        unmapped_mismatch: 5,
        unmapped_multi: 6,
        unmapped_all: 7,
        chimeric_all: 8,
        ..Default::default()
    };
    stats_l21_stats_addstats(&mut stats, &add);
    assert_eq!(stats.read_n, 6);
    assert_eq!(stats.mapped_del_l, 15);
    assert_eq!(stats.mapped_portion, 4.0);
    assert_eq!(stats.splices_n, [8, 8, 8, 8, 8, 8, 8]);

    stats_l4_stats_resetn(&mut stats);
    assert_eq!(stats.read_n, 0);
    assert_eq!(stats.splices_n, [0; SJ_MOTIF_SIZE]);

    let tr = Transcript {
        n_exons: 3,
        exons: tr_exons_arr(&[[0, 100, 10, 0, 0], [20, 140, 15, 0, 0], [40, 200, 5, 0, 0]]),
        canon_sj: tr_canon_sj_arr(&[1, -1]),
        sj_annot: tr_sj_u8_arr(&[1, 0]),
        n_mm: 2,
        n_ins: 1,
        n_del: 3,
        l_ins: 4,
        l_del: 6,
        ..Default::default()
    };
    stats_l35_stats_transcriptstats(&mut stats, &tr, 60);
    assert_eq!(stats.mapped_mismatches_n, 2);
    assert_eq!(stats.mapped_bases, 30);
    assert_eq!(stats.splices_n[1], 1);
    assert_eq!(stats.splices_nsjdb, 1);
    assert_eq!(stats.mapped_portion, 0.5);
}

#[test]
fn stats_text_helpers_match_original_emission_shape() {
    let header = stats_l62_stats_progressreportheader();
    assert!(header.contains("Time"));
    assert!(header.contains("M/hr"));
    assert!(header.ends_with('\n'));

    let stats = Stats {
        read_n: 10,
        mapped_reads_u: 7,
        mapped_reads_m: 2,
        ..Default::default()
    };
    assert_eq!(
        stats_l147_stats_writelines(&stats, &[0, 1, 2, 1], "#", "sample"),
        "# sample\n# Nreads 10\tNreadsUnique 7\tNreadsMulti 2\n# sample\n# Nreads 10\tNreadsUnique 7\tNreadsMulti 2\n"
    );
    assert_eq!(
        stats_l147_stats_writelines(&stats, &[1], "#", ""),
        "# Nreads 10\tNreadsUnique 7\tNreadsMulti 2\n"
    );

    let mut progress = Stats {
        read_n: 100,
        read_bases: 10_000,
        mapped_reads_u: 70,
        mapped_reads_m: 20,
        mapped_bases: 7_000,
        mapped_mismatches_n: 14,
        unmapped_multi: 3,
        unmapped_mismatch: 2,
        unmapped_short: 4,
        unmapped_other: 1,
        time_start_map: 1_000,
        time_last_report: 1_030,
        ..Default::default()
    };
    assert_eq!(stats_l73_stats_progressreport(&mut progress, 1_080), None);
    let report = stats_l73_stats_progressreport(&mut progress, 1_091).unwrap();
    assert!(report.ends_with('\n'));
    assert!(report.contains("100"));
    assert!(report.contains("70.0%"));
    assert_eq!(progress.time_last_report, 1_091);

    let mut final_stats = Stats {
        read_n: 100,
        read_bases: 10_000,
        mapped_reads_u: 70,
        mapped_reads_m: 20,
        mapped_bases: 7_000,
        mapped_mismatches_n: 14,
        mapped_ins_n: 2,
        mapped_del_n: 3,
        mapped_ins_l: 6,
        mapped_del_l: 9,
        splices_n: [1, 2, 3, 4, 5, 6, 7],
        splices_nsjdb: 8,
        unmapped_multi: 3,
        unmapped_mismatch: 2,
        unmapped_short: 4,
        unmapped_other: 1,
        chimeric_all: 5,
        time_start: 1_000,
        time_start_map: 1_100,
        ..Default::default()
    };
    let final_report = stats_l99_stats_reportfinal(&mut final_stats, 4_700);
    assert!(final_report.contains("Number of input reads |\t100"));
    assert!(final_report.contains("Uniquely mapped reads % |\t70.00%"));
    assert!(final_report.contains("Number of splices: Total |\t28"));
    assert!(final_report.contains("% of chimeric reads |\t5.00%"));
    assert_eq!(final_stats.time_finish, 4_700);
}

#[test]
fn chimeric_segment_primitives_match_original_logic() {
    let p_ch = ParametersChimeric {
        segment_min: 10,
        segment_read_gap_max: 5,
        ..Default::default()
    };
    let align = Transcript {
        exons: tr_exons_arr(&[[0, 100, 20, 0, 0]]),
        intron_motifs: [0, 1, 0],
        n_exons: 1,
        l_read: 100,
        read_length: tr_read_length_arr(&[50]),
        r_length: 20,
        max_score: 20,
        str_: 0,
        ..Default::default()
    };
    let seg = chimericsegment_l3_chimericsegment_chimericsegment(&p_ch, &align);
    assert_eq!(seg.str_, 1);
    assert_eq!(seg.ro_s, 0);
    assert_eq!(seg.ro_e, 19);
    assert!(chimericsegment_l19_chimericsegment_segmentcheck(&seg));

    let reverse_align = Transcript {
        exons: tr_exons_arr(&[[20, 100, 20, 0, 0]]),
        intron_motifs: [0, 0, 1],
        n_exons: 1,
        l_read: 100,
        read_length: tr_read_length_arr(&[50]),
        r_length: 20,
        max_score: 25,
        str_: 1,
        ..Default::default()
    };
    let reverse_seg = chimericsegment_l3_chimericsegment_chimericsegment(&p_ch, &reverse_align);
    assert_eq!(reverse_seg.str_, 1);
    assert_eq!(reverse_seg.ro_s, 59);
    assert_eq!(reverse_seg.ro_e, 78);

    let bad_seg = ChimericSegment {
        align: Transcript {
            r_length: 9,
            intron_motifs: [1, 0, 0],
            ..Default::default()
        },
        segment_min: 10,
        ..Default::default()
    };
    assert!(!chimericsegment_l19_chimericsegment_segmentcheck(&bad_seg));

    let seg1 = ChimericSegment {
        align: Transcript {
            read_length: tr_read_length_arr(&[50]),
            max_score: 40,
            ..Default::default()
        },
        ro_s: 0,
        ro_e: 30,
        segment_min: 10,
        segment_read_gap_max: 5,
        ..Default::default()
    };
    let seg2 = ChimericSegment {
        align: Transcript {
            read_length: tr_read_length_arr(&[50]),
            max_score: 35,
            ..Default::default()
        },
        ro_s: 25,
        ro_e: 60,
        segment_min: 10,
        segment_read_gap_max: 5,
        ..Default::default()
    };
    assert_eq!(
        chimericdetection_chimericdetectionmult_l6_chimericalignscore(&seg1, &seg2),
        69
    );
}

#[test]
fn transcript_generate_cigarp_matches_softclip_indel_junction_and_mate_gap_rules() {
    let tr = Transcript {
        exons: tr_exons_arr(&[[2, 100, 8, 0, 0], [12, 118, 5, 0, 0], [20, 130, 4, 0, 0]]),
        canon_sj: tr_canon_sj_arr(&[-1, 1]),
        sj_annot: tr_sj_u8_arr(&[0, 0]),
        n_exons: 3,
        read_nmates: 1,
        read_length_original: tr_read_length_arr(&[30]),
        read_length_pair_original: 30,
        ..Default::default()
    };
    assert_eq!(
        transcript_generatecigarp_l4_transcript_generatecigarp(&tr),
        "2S8M2I10D5M3I7N4M6S"
    );

    let mate_gap = Transcript {
        exons: tr_exons_arr(&[[0, 100, 20, 0, 0], [35, 140, 10, 1, 0]]),
        canon_sj: tr_canon_sj_arr(&[-3]),
        sj_annot: tr_sj_u8_arr(&[0]),
        n_exons: 2,
        read_nmates: 2,
        str_: 0,
        read_length_original: tr_read_length_arr(&[30, 25]),
        read_length_pair_original: 60,
        ..Default::default()
    };
    assert_eq!(
        transcript_generatecigarp_l4_transcript_generatecigarp(&mate_gap),
        "20M10S20p4S10M15S"
    );

    let overlap = Transcript {
        exons: tr_exons_arr(&[[0, 100, 20, 0, 0], [20, 110, 10, 1, 0]]),
        canon_sj: tr_canon_sj_arr(&[-3]),
        sj_annot: tr_sj_u8_arr(&[0]),
        n_exons: 2,
        read_nmates: 2,
        str_: 0,
        read_length_original: tr_read_length_arr(&[30, 25]),
        read_length_pair_original: 60,
        ..Default::default()
    };
    assert_eq!(
        transcript_generatecigarp_l4_transcript_generatecigarp(&overlap),
        "20M-10p10M"
    );
}

#[test]
fn read_align_cigar_helpers_match_transcript_cigarp_and_mate_specific_rules() {
    let read_align = ReadAlign {
        read_length_original: vec![30, 25],
        read_length_pair_original: 60,
        ..Default::default()
    };
    let tr = Transcript {
        exons: tr_exons_arr(&[[0, 100, 20, 0, 0], [35, 140, 10, 1, 0]]),
        canon_sj: tr_canon_sj_arr(&[-3]),
        sj_annot: tr_sj_u8_arr(&[0]),
        n_exons: 2,
        str_: 0,
        ..Default::default()
    };
    assert_eq!(
        readalign_outputtranscriptcigarp_l4_readalign_outputtranscriptcigarp(&read_align, &tr, 2),
        "20M10S20p4S10M15S"
    );

    let mut read_align = ReadAlign {
        read_length: vec![30, 25],
        read_length_original: vec![30, 25],
        clip_mates: vec![
            vec![
                ClipMate {
                    clipped_n: 2,
                    ..Default::default()
                },
                ClipMate {
                    clipped_n: 9,
                    ..Default::default()
                },
            ],
            vec![
                ClipMate {
                    clipped_n: 4,
                    ..Default::default()
                },
                ClipMate {
                    clipped_n: 3,
                    ..Default::default()
                },
            ],
        ],
        ..Default::default()
    };
    let tr = Transcript {
        exons: tr_exons_arr(&[[2, 100, 8, 0, 0], [12, 118, 5, 0, 0], [31, 200, 10, 1, 0]]),
        canon_sj: tr_canon_sj_arr(&[-1, 1]),
        sj_annot: tr_sj_u8_arr(&[0, 0]),
        n_exons: 3,
        str_: 0,
        ..Default::default()
    };
    readalign_calccigar_l3_readalign_calccigar(&mut read_align, &tr, 2, 1, 0);
    assert_eq!(read_align.mates_cigar, vec!["4S8M2I10D5M11S", "3S10M12S"]);
}

#[test]
fn transcript_convert_genome_cigar_splits_blocks_and_trims_short_terminal_junctions() {
    let gen_out = Genome {
        genome_out: GenomeOut {
            conv_blocks: vec![[0, 10, 0], [10, 10, 15], [20, 10, 30], [u64::MAX, 0, 0]],
            n_minus_strand_offset: 1_000,
            ..Default::default()
        },
        align_sjdb_overhang_min: 6,
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 4,
            ..Default::default()
        },
        chr_bin: vec![7, 8, 9],
        ..Default::default()
    };
    let source = Transcript {
        g_start: 5,
        str_: 0,
        cigar: vec![[0, 20]],
        ..Default::default()
    };
    let mut out = Transcript::default();

    assert!(
        transcript_convertgenomecigar_l2_transcript_convertgenomecigar(&source, &gen_out, &mut out)
    );
    assert_eq!(out.g_start, 15);
    assert_eq!(out.cigar, vec![[4, 5], [0, 10], [4, 5]]);
    assert_eq!(out.g_length, 10);
    assert_eq!(out.str_, 0);
    assert_eq!(out.chr, 7);
}

#[test]
fn transcript_convert_genome_cigar_flips_minus_strand_output_coordinates() {
    let gen_out = Genome {
        genome_out: GenomeOut {
            conv_blocks: vec![[0, 10, 0], [10, 10, 15], [20, 10, 40], [u64::MAX, 0, 0]],
            n_minus_strand_offset: 40,
            ..Default::default()
        },
        align_sjdb_overhang_min: 3,
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 4,
            ..Default::default()
        },
        chr_bin: vec![0, 1, 2],
        ..Default::default()
    };
    let source = Transcript {
        g_start: 22,
        str_: 1,
        cigar: vec![[0, 5]],
        ..Default::default()
    };
    let mut out = Transcript::default();

    assert!(
        transcript_convertgenomecigar_l2_transcript_convertgenomecigar(&source, &gen_out, &mut out)
    );
    assert_eq!(out.g_start, 33);
    assert_eq!(out.g_length, 5);
    assert_eq!(out.cigar, vec![[0, 5]]);
    assert_eq!(out.str_, 0);
    assert_eq!(out.chr, 2);
}

#[test]
fn transcript_transform_genome_splits_blocks_and_recomputes_deletion_junctions() {
    let gen_out = Genome {
        genome_out: GenomeOut {
            conv_blocks: vec![[0, 10, 0], [10, 10, 15], [20, 10, 30], [u64::MAX, 0, 0]],
            ..Default::default()
        },
        align_intron_min: 10,
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 4,
            ..Default::default()
        },
        chr_bin: vec![3, 4, 5],
        ..Default::default()
    };
    let source = Transcript {
        exons: tr_exons_arr(&[[0, 5, 20, 0, 0]]),
        n_exons: 1,
        str_: 1,
        ..Default::default()
    };
    let mut out = Transcript::default();

    assert!(transcript_transformgenome_l4_transcript_transformgenome(
        &source, &gen_out, &mut out
    ));
    assert_eq!(
        &out.exons[..3],
        &[[0, 5, 5, 0, 0], [5, 15, 10, 0, 0], [15, 30, 5, 0, 0]]
    );
    assert_eq!(out.n_exons, 3);
    assert_eq!(&out.canon_sj[..2], &[-1, -1]);
    assert_eq!(&out.sj_annot[..2], &[0u8, 0]);
    assert_eq!(out.str_, 1);
    assert_eq!(out.chr, 3);
}

#[test]
fn transcript_transform_genome_rejects_annotated_zero_motif_shift_without_overhang() {
    let gen_out = Genome {
        genome_out: GenomeOut {
            conv_blocks: vec![[0, 1_000, 0], [u64::MAX, 0, 0]],
            ..Default::default()
        },
        sjdb_start: vec![105],
        sjdb_end: vec![114],
        sjdb_motif: vec![0],
        sjdb_shift_left: vec![5],
        sjdb_n: 1,
        p_ge: ParametersGenome {
            g_chr_bin_nbits: 4,
            ..Default::default()
        },
        chr_bin: vec![0; 128],
        ..Default::default()
    };
    let source = Transcript {
        exons: tr_exons_arr(&[[0, 100, 5, 0, 0], [5, 115, 5, 0, 0]]),
        n_exons: 2,
        ..Default::default()
    };
    let mut out = Transcript::default();

    assert!(!transcript_transformgenome_l4_transcript_transformgenome(
        &source, &gen_out, &mut out
    ));
}

#[test]
fn read_align_output_variation_preserves_original_noop_behavior() {
    let var = Variation {
        yes: true,
        ..Default::default()
    };
    let tr = Transcript {
        var_ind: vec![1, 2],
        var_gen_coord: vec![100, 200],
        var_read_coord: vec![3, 4],
        var_allele: vec![0, 1],
        ..Default::default()
    };
    let before = tr.clone();
    readalign_outputvariation_l3_readalign_outputvariation(&var, &tr, 0, 2);
    assert_eq!(tr, before);

    let var_off = Variation {
        yes: false,
        ..Default::default()
    };
    readalign_outputvariation_l3_readalign_outputvariation(&var_off, &tr, 1, 2);
    assert_eq!(tr, before);
}

#[test]
fn read_align_sam_attr_nm_md_matches_mismatches_n_and_indel_gaps() {
    let p = Parameters {
        genome_num_to_nt: b"ACGTN\0".to_vec(),
        ..Default::default()
    };
    let gen_out = Genome {
        g: vec![0, 1, 2, 3, 4, 0, 1, 2, 2, 0, 1, 2],
        ..Default::default()
    };
    let read_align = ReadAlign {
        read1: [
            vec![0, 1, 3, 3, 4, 0, 2, 2],
            Vec::new(),
            vec![3, 3, 3, 3, 3, 3, 3, 3],
        ],
        ..Default::default()
    };
    let tr = Transcript {
        exons: tr_exons_arr(&[[0, 0, 5, 0, 0], [6, 8, 2, 0, 0]]),
        canon_sj: tr_canon_sj_arr(&[-1]),
        ro_str: 0,
        ..Default::default()
    };

    let (tag_nm, tag_md) =
        readalign_alignbam_l9_readalign_samattrnm_md(&read_align, &p, &gen_out, &tr, 0, 1);

    assert_eq!(tag_nm, 7);
    assert_eq!(tag_md, "2G1N0^ACG1A0");
}

#[test]
fn read_align_sam_attr_nm_md_uses_reverse_read_buffer() {
    let p = Parameters {
        genome_num_to_nt: b"ACGTN\0".to_vec(),
        ..Default::default()
    };
    let gen_out = Genome {
        g: vec![3, 2, 1, 0],
        ..Default::default()
    };
    let read_align = ReadAlign {
        read1: [vec![0, 0, 0, 0], Vec::new(), vec![3, 2, 0, 0]],
        ..Default::default()
    };
    let tr = Transcript {
        exons: tr_exons_arr(&[[0, 0, 4, 0, 0]]),
        ro_str: 1,
        ..Default::default()
    };

    let (tag_nm, tag_md) =
        readalign_alignbam_l9_readalign_samattrnm_md(&read_align, &p, &gen_out, &tr, 0, 0);

    assert_eq!(tag_nm, 1);
    assert_eq!(tag_md, "2C1");
}

#[test]
fn read_align_align_bam_writes_mapped_single_end_record_bytes() {
    let mut p = Parameters {
        read_nmates: 1,
        out_sam_mode: "Full".to_string(),
        out_sam_mapq_unique: 255,
        out_sam_flag_and: u16::MAX,
        out_sam_attr_ih_start: 1,
        out_sam_attr_order: vec![
            ATTR_NH,
            ATTR_HI,
            ATTR_AS,
            ATTR_NM_LOWER,
            ATTR_NM,
            ATTR_MD,
            ATTR_RG,
        ],
        out_sam_attr_rg: vec!["rg1".to_string()],
        genome_num_to_nt: b"ACGTN\0".to_vec(),
        ..Default::default()
    };
    p.out_sam_attr_present.mc = false;
    let gen_out = Genome {
        chr_start: vec![0],
        n_chr_real: 1,
        g: vec![0, 1, 2, 3],
        ..Default::default()
    };
    let mut read_align = ReadAlign {
        read_name: "@r1".to_string(),
        read_files_index: 0,
        read_file_type: 2,
        read0: vec![b"ACGT".to_vec()],
        qual0: vec![b"IIII".to_vec()],
        read1: [vec![0, 1, 2, 3], Vec::new(), vec![3, 2, 1, 0]],
        read_length: vec![4],
        read_length_original: vec![4],
        clip_mates: vec![vec![ClipMate::default(), ClipMate::default()]],
        ..Default::default()
    };
    let tr = Transcript {
        exons: tr_exons_arr(&[[0, 0, 4, 0, 0]]),
        n_exons: 1,
        chr: 0,
        str_: 0,
        primary_flag: true,
        max_score: 4,
        ..Default::default()
    };

    let result = readalign_alignbam_l47_readalign_alignbam(
        &mut read_align,
        &p,
        &gen_out,
        &tr,
        1,
        0,
        0,
        u32::MAX,
        u32::MAX,
        0,
        -1,
        None,
        &p.out_sam_attr_order,
    )
    .unwrap();

    let rec = &result.records[0];
    let u32_at = |i: usize| u32::from_ne_bytes(rec[i * 4..i * 4 + 4].try_into().unwrap());
    assert_eq!(result.n_lines, 1);
    assert_eq!(u32_at(0), result.record_sizes[0] - 4);
    assert_eq!(u32_at(1), 0);
    assert_eq!(u32_at(2), 0);
    assert_eq!(u32_at(3) & 0xff, 3);
    assert_eq!((u32_at(3) >> 8) & 0xff, 255);
    assert_eq!(u32_at(4), 1);
    assert_eq!(u32_at(5), 4);
    assert_eq!(&rec[36..39], b"r1\0");
    assert_eq!(
        u32::from_ne_bytes(rec[39..43].try_into().unwrap()),
        4 << 4 | BAM_CIGAR_M
    );
    assert_eq!(&rec[43..45], &[0x12, 0x48]);
    assert_eq!(&rec[45..49], &[40, 40, 40, 40]);
    assert!(rec.windows(3).any(|w| w == b"NH\x69"));
    assert!(rec.windows(3).any(|w| w == b"MDZ"));
    assert!(rec.windows(3).any(|w| w == b"RGZ"));
}

#[test]
fn read_align_align_bam_writes_unmapped_mate_with_mapped_mate_info() {
    let p = Parameters {
        read_nmates: 2,
        out_sam_mode: "Full".to_string(),
        out_sam_flag_and: u16::MAX,
        out_sam_unmapped_keep_pairs: true,
        ..Default::default()
    };
    let gen_out = Genome {
        chr_start: vec![100],
        n_chr_real: 1,
        ..Default::default()
    };
    let mut read_align = ReadAlign {
        read_name: "@pair".to_string(),
        read_filter: b'Y' as i32,
        read_file_type: 1,
        read0: vec![b"AAAA".to_vec(), b"TTTT".to_vec()],
        qual0: vec![b"IIII".to_vec(), b"JJJJ".to_vec()],
        read_length: vec![4, 4],
        read_length_original: vec![4, 4],
        clip_mates: vec![
            vec![ClipMate::default(), ClipMate::default()],
            vec![ClipMate::default(), ClipMate::default()],
        ],
        ..Default::default()
    };
    let tr = Transcript {
        exons: tr_exons_arr(&[[0, 120, 4, 0, 0]]),
        n_exons: 1,
        chr: 0,
        str_: 0,
        primary_flag: false,
        max_score: 4,
        n_mm: 1,
        ..Default::default()
    };

    let result = readalign_alignbam_l47_readalign_alignbam(
        &mut read_align,
        &p,
        &gen_out,
        &tr,
        0,
        0,
        100,
        u32::MAX,
        u32::MAX,
        0,
        4,
        Some([true, false]),
        &[],
    )
    .unwrap();

    let rec = &result.records[1];
    let u32_at = |i: usize| u32::from_ne_bytes(rec[i * 4..i * 4 + 4].try_into().unwrap());
    assert_eq!(result.n_lines, 2);
    assert!(result.records[0].is_empty());
    assert_eq!(result.sam_flags[1], 0x385);
    assert_eq!(u32_at(1), u32::MAX);
    assert_eq!(u32_at(2), u32::MAX);
    assert_eq!(u32_at(4) >> 16, 0x385);
    assert_eq!(u32_at(6), 0);
    assert_eq!(u32_at(7), 20);
    assert_eq!(u32_at(5), 4);
    assert_eq!(&rec[36..41], b"pair\0");
    assert!(rec.windows(3).any(|w| w == b"uTA"));
}

#[test]
fn solo_read_feature_accumulators_match_original_methods() {
    let mut wl = SoloReadFeature {
        cb_wl_yes: true,
        cb_wl_size: 3,
        cb_read_count: vec![1, 2, 3],
        transcript_dist_count: vec![10, 20],
        ..Default::default()
    };
    let wl_in = SoloReadFeature {
        cb_read_count: vec![4, 5, 6],
        transcript_dist_count: vec![1, 2],
        ..Default::default()
    };
    soloreadfeature_l29_soloreadfeature_addcounts(&mut wl, &wl_in);
    assert_eq!(wl.cb_read_count, vec![5, 7, 9]);
    assert_eq!(wl.transcript_dist_count, vec![11, 22]);

    let mut map_feature = SoloReadFeature::default();
    map_feature.cb_read_count_map.insert(10, 1);
    map_feature.cb_read_count_map.insert(20, 2);
    let mut map_in = SoloReadFeature::default();
    map_in.cb_read_count_map.insert(10, 4);
    map_in.cb_read_count_map.insert(30, 5);
    soloreadfeature_l29_soloreadfeature_addcounts(&mut map_feature, &map_in);
    assert_eq!(map_feature.cb_read_count_map.get(&10), Some(&5));
    assert_eq!(map_feature.cb_read_count_map.get(&20), Some(&2));
    assert_eq!(map_feature.cb_read_count_map.get(&30), Some(&5));

    let mut stats_feature = SoloReadFeature {
        stats: SoloReadFeatureStats {
            names: vec!["a".to_string(), "b".to_string()],
            v: vec![1, 2],
        },
        ..Default::default()
    };
    stats_feature.read_flag.flag_counts_no_cb[0] = 1;
    stats_feature.read_flag.flag_counts_no_cb[3] = 2;
    let mut stats_in = SoloReadFeature {
        stats: SoloReadFeatureStats {
            names: vec!["a".to_string(), "b".to_string()],
            v: vec![10, 20],
        },
        ..Default::default()
    };
    stats_in.read_flag.flag_counts_no_cb[0] = 7;
    stats_in.read_flag.flag_counts_no_cb[3] = 8;
    soloreadfeature_l47_soloreadfeature_addstats(&mut stats_feature, &stats_in);
    assert_eq!(stats_feature.stats.v, vec![11, 22]);
    assert_eq!(stats_feature.read_flag.flag_counts_no_cb[0], 8);
    assert_eq!(stats_feature.read_flag.flag_counts_no_cb[3], 10);

    assert_eq!(
        soloreadfeature_l56_soloreadfeature_statsout(&stats_feature),
        "                                                 a             11\n                                                 b             22\n"
    );
}
