use star_rs::*;

#[test]
fn sequence_reverse_complement_and_numeric_complement_match_star_tables() {
    let mut out = [0; 5];
    sequencefuns_l4_complementseqnumbers(&[0, 1, 2, 3, 9], &mut out, 5);
    assert_eq!(out, [3, 2, 1, 0, 9]);
    sequencefuns_l4_complementseqnumbers(&[0, 1], &mut out, 5);
    assert_eq!(out, [3, 2, 1, 0, 9]);

    let mut rc = [0; 6];
    sequencefuns_l16_revcomplementnucleotides(b"ACGTry", &mut rc, 6);
    assert_eq!(&rc, b"ryACGT");
    sequencefuns_l16_revcomplementnucleotides(b"AC", &mut rc, 6);
    assert_eq!(&rc, b"GTACGT");

    let mut seq = "ACGTry".to_string();
    sequencefuns_l56_revcomplementnucleotides(&mut seq);
    assert_eq!(seq, "ryACGT");
}

#[test]
fn sequence_bam_encoding_matches_star_tables() {
    assert_eq!(sequencefuns_l99_nucltonumbam(b'='), 0);
    assert_eq!(sequencefuns_l99_nucltonumbam(b'A'), 1);
    assert_eq!(sequencefuns_l99_nucltonumbam(b'n'), 15);
    assert_eq!(sequencefuns_l99_nucltonumbam(b'?'), 15);

    let mut packed = [0; 3];
    sequencefuns_l122_nuclpackbam(b"ACGTN", &mut packed, 5);
    assert_eq!(packed, [0x12, 0x48, 0xf0]);
    sequencefuns_l122_nuclpackbam(b"ACGTN", &mut packed[..2], 5);
    assert_eq!(packed, [0x12, 0x48, 0xf0]);
}

#[test]
fn sequence_base_number_conversion_round_trips_with_single_n_position() {
    let mut out = [9; 6];
    sequencefuns_l131_convertnucleotidestonumbers(b"AcgTx?", &mut out, 6);
    assert_eq!(out, [0, 1, 2, 3, 4, 4]);
    sequencefuns_l131_convertnucleotidestonumbers(b"AC", &mut out, 6);
    assert_eq!(out, [0, 1, 2, 3, 4, 4]);

    let mut capital = *b"ACGTX";
    sequencefuns_l148_convertcapitalbasestonum(&mut capital, 5);
    assert_eq!(capital, [0, 1, 2, 3, 4]);
    sequencefuns_l148_convertcapitalbasestonum(&mut capital[..2], 5);
    assert_eq!(capital, [4, 4, 2, 3, 4]);

    assert_eq!(sequencefuns_l195_convertnt01234(b't'), 3);
    assert_eq!(sequencefuns_l195_convertnt01234(b'N'), 4);

    let mut int32 = 0;
    assert_eq!(
        sequencefuns_l219_convertnuclstrtoint32("ACNT", &mut int32),
        2
    );
    assert_eq!(int32, 19);
    assert_eq!(sequencefuns_l237_convertnuclint32tostring(int32, 4), "ACAT");
    assert_eq!(
        sequencefuns_l219_convertnuclstrtoint32("ANNT", &mut int32),
        -2
    );

    let mut int64 = 0;
    assert_eq!(
        sequencefuns_l249_convertnuclstrtoint64("ACGT", &mut int64),
        -1
    );
    assert_eq!(sequencefuns_l267_convertnuclint64tostring(int64, 4), "ACGT");
}

#[test]
fn sequence_remove_controls_preserves_star_indexing_behavior() {
    let mut out = [99; 5];
    let kept = sequencefuns_l170_convertnucleotidestonumbersremovecontrols(
        &[b'A', 1, b'C', b'X', b'T'],
        &mut out,
        5,
    );
    assert_eq!(kept, 4);
    assert_eq!(out, [0, 99, 1, 4, 3]);

    let kept =
        sequencefuns_l170_convertnucleotidestonumbersremovecontrols(b"ACGT", &mut out[..2], 4);
    assert_eq!(kept, 2);
    assert_eq!(out, [0, 1, 1, 4, 3]);
}

#[test]
fn sequence_search_and_split_primitives_match_star_logic() {
    assert_eq!(sequencefuns_l280_chrfind(25, 4, &[0, 10, 20, 30]), 2);
    assert_eq!(
        sequencefuns_l293_localsearch(&[4, 0, 1, 2], 4, &[0, 1], 2, 0.0),
        1
    );
    assert_eq!(
        sequencefuns_l293_localsearch(&[4, 0, 1, 2], 8, &[0], 2, 0.0),
        1
    );
    assert_eq!(
        sequencefuns_l317_localsearchnismm(&[4, 0, 1, 2], 4, &[0, 1], 2, 0.0),
        1
    );
    assert_eq!(
        sequencefuns_l317_localsearchnismm(&[4, 0, 1, 2], 8, &[0], 2, 0.0),
        1
    );

    let mut pos = 99;
    assert_eq!(
        sequencefuns_l341_localalignhammingdist("AACCGT", "CNN", &mut pos),
        0
    );
    assert_eq!(pos, 2);

    let mut split = [vec![0; 4], vec![0; 4], vec![0; 4]];
    let n = sequencefuns_l411_qualitysplit(
        &[4, 0, 1, MARK_FRAG_SPACER_BASE, 2, 3, 4, 0],
        8,
        4,
        2,
        &mut split,
    );
    assert_eq!(n, 2);
    assert_eq!(split[0][..2], [1, 4]);
    assert_eq!(split[1][..2], [2, 2]);
    assert_eq!(split[2][..2], [0, 1]);

    let mut split_short = [vec![0; 1], vec![0; 1], vec![0; 1]];
    let n = sequencefuns_l411_qualitysplit(
        &[4, 0, 1, MARK_FRAG_SPACER_BASE],
        8,
        4,
        1,
        &mut split_short,
    );
    assert_eq!(n, 1);
    assert_eq!(split_short[0][0], 1);
    assert_eq!(split_short[1][0], 2);
    assert_eq!(split_short[2][0], 0);
}
