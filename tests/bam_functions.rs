use star_rs::*;
use star_rs::{
    BAMoutput, Bam1, Genome, Parameters, ReadAlignChunk, SoloFeatureReadInfo,
};

fn bam_record(tid: u32, pos: u32, payload: u8) -> Vec<u8> {
    let mut record = Vec::new();
    record.extend_from_slice(&12u32.to_le_bytes());
    record.extend_from_slice(&tid.to_le_bytes());
    record.extend_from_slice(&pos.to_le_bytes());
    record.extend_from_slice(&payload.to_le_bytes());
    record.extend_from_slice(&0u8.to_le_bytes());
    record.extend_from_slice(&0u8.to_le_bytes());
    record.extend_from_slice(&0u8.to_le_bytes());
    record
}

fn bam_record_with_read_order(tid: u32, pos: u32, payload: u8, read_order: u64) -> Vec<u8> {
    let mut record = bam_record(tid, pos, payload);
    record.extend_from_slice(&read_order.to_le_bytes());
    record
}

fn bam_record_with_unmapped_order(tid: u32, pos: u32, payload: u8, read_order_top: u64) -> Vec<u8> {
    bam_record_with_read_order(tid, pos, payload, read_order_top << 32)
}

#[test]
fn bam_cigar_string_and_reg2bin_match_htslib_macros() {
    let cigar = [(10 << 4), (1 << 4) | 1, (5 << 4) | 2, (3 << 4) | 4];
    assert_eq!(
        bamfunctions_l5_bam_cigarstring(&cigar, cigar.len() as i32),
        "10M1I5D3S"
    );
    assert_eq!(bamfunctions_l5_bam_cigarstring(&[], 0), "");
    assert_eq!(bamfunctions_l5_bam_cigarstring(&cigar[..2], 4), "10M1I");

    assert_eq!(bamfunctions_l95_reg2bin(0, 1), 4681);
    assert_eq!(bamfunctions_l95_reg2bin(20_000, 40_000), 585);
    assert_eq!(bamfunctions_l95_reg2bin(0, 600_000_000), 0);
}

#[test]
fn bam_read1_from_array_matches_core_field_unpacking() {
    let block_len = 44i32;
    let x = [
        2u32,
        101,
        (7 << 16) | (30 << 8) | 5,
        (0x41 << 16) | 3,
        12,
        4,
        202,
        99,
    ];
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&block_len.to_ne_bytes());
    for word in x {
        bytes.extend_from_slice(&word.to_ne_bytes());
    }
    bytes.extend_from_slice(&[0u8; 12]);

    let mut bam = Bam1 {
        m_data: 4,
        ..Default::default()
    };
    assert_eq!(bamfunctions_l30_bam_read1_fromarray(&bytes, &mut bam), 48);
    assert_eq!(bam.core.tid, 2);
    assert_eq!(bam.core.pos, 101);
    assert_eq!(bam.core.bin, 7);
    assert_eq!(bam.core.qual, 30);
    assert_eq!(bam.core.l_qname, 5);
    assert_eq!(bam.core.flag, 0x41);
    assert_eq!(bam.core.n_cigar, 3);
    assert_eq!(bam.core.l_qseq, 12);
    assert_eq!(bam.core.mtid, 4);
    assert_eq!(bam.core.mpos, 202);
    assert_eq!(bam.core.isize, 99);
    assert_eq!(bam.l_data, 12);
    assert_eq!(bam.m_data, 16);
    assert_eq!(bam.data_offset, 36);

    let mut bad = Bam1::default();
    let mut short_block = bytes.clone();
    short_block[0..4].copy_from_slice(&16i32.to_ne_bytes());
    assert_eq!(
        bamfunctions_l30_bam_read1_fromarray(&short_block, &mut bad),
        -4
    );

    let truncated_body = bytes[..36].to_vec();
    assert_eq!(
        bamfunctions_l30_bam_read1_fromarray(&truncated_body, &mut bad),
        -4
    );
}

#[test]
fn bamoutput_unsorted_constructor_buffers_and_flushes_like_original() {
    let p = Parameters {
        chunk_out_bam_size_bytes: 12,
        ..Default::default()
    };
    let mut out = bamoutput_l36_bamoutput_bamoutput(vec![1, 2], &p);
    let a1 = [10u8, 11, 12, 13];
    let a2 = [20u8, 21, 22, 23];

    bamoutput_l52_bamoutput_unsortedonealign(&mut out, &a1, 4, 4).unwrap();
    assert_eq!(out.bin_bytes1, 4);
    assert_eq!(out.bgzf_bam, vec![1, 2]);

    bamoutput_l52_bamoutput_unsortedonealign(&mut out, &a2, 4, 16).unwrap();
    assert_eq!(out.bgzf_bam, vec![1, 2, 10, 11, 12, 13]);
    assert_eq!(out.bin_bytes1, 4);

    bamoutput_l70_bamoutput_unsortedflush(&mut out);
    assert_eq!(out.bgzf_bam, vec![1, 2, 10, 11, 12, 13, 20, 21, 22, 23]);
    assert_eq!(out.bin_bytes1, 0);

    let err = bamoutput_l52_bamoutput_unsortedonealign(&mut out, &[1, 2], 4, 4).unwrap_err();
    assert!(err.contains("exceeds input buffer"));

    let err = bamoutput_l52_bamoutput_unsortedonealign(&mut out, &[0; 16], 16, 16).unwrap_err();
    assert!(err.contains("exceeds output buffer"));
}

#[test]
fn bamoutput_coordinate_constructor_rebins_and_flushes_records() {
    let mut p = Parameters {
        out_bam_coord_nbins: 4,
        chunk_out_bam_size_bytes: 120,
        out_bam_sorting_bin_start: vec![1, 0, 0, 0],
        ..Default::default()
    };
    let mut out = bamoutput_l9_bamoutput_bamoutput(2, "/tmp/bam", &p);
    assert_eq!(out.bam_dir, "/tmp/bam2");
    assert_eq!(out.n_bins, 1);
    assert_eq!(out.bin_size, 30);
    assert_eq!(out.bin_size1, 90);

    let r1 = bam_record(0, 100, 1);
    let r2 = bam_record(0, 300, 2);
    let r3 = bam_record(0, 500, 3);
    bamoutput_l77_bamoutput_coordonealign(&mut out, &mut p, &r1, r1.len() as u64, 10).unwrap();
    bamoutput_l77_bamoutput_coordonealign(&mut out, &mut p, &r2, r2.len() as u64, 20).unwrap();
    bamoutput_l77_bamoutput_coordonealign(&mut out, &mut p, &r3, r3.len() as u64, 30).unwrap();

    bamoutput_l168_bamoutput_coordflush(&mut out, &mut p).unwrap();
    assert_eq!(out.n_bins, 4);
    assert_eq!(p.out_bam_sorting_bin_start, vec![0, 300, 500, 0]);
    assert!(p.bam_sorting_log.contains("BAM sorting: 3 mapped reads"));
    assert_eq!(out.bin_total_n[0], 1);
    assert_eq!(out.bin_total_n[1], 1);
    assert_eq!(out.bin_total_n[2], 1);
    assert_eq!(out.bin_bytes, vec![0, 0, 0, 0]);
    assert_eq!(out.bin_streams[0], {
        let mut v = r1.clone();
        v.extend_from_slice(&10u64.to_le_bytes());
        v
    });
    assert_eq!(out.bin_streams[1], {
        let mut v = r2.clone();
        v.extend_from_slice(&20u64.to_le_bytes());
        v
    });
    assert_eq!(out.bin_streams[2], {
        let mut v = r3.clone();
        v.extend_from_slice(&30u64.to_le_bytes());
        v
    });
}

#[test]
fn bamoutput_coordinate_reports_malformed_records_without_panic() {
    let mut p = Parameters {
        out_bam_coord_nbins: 2,
        chunk_out_bam_size_bytes: 120,
        out_bam_sorting_bin_start: vec![0, 0],
        ..Default::default()
    };
    let mut out = bamoutput_l9_bamoutput_bamoutput(0, "/tmp/bam", &p);

    let short = bamoutput_l77_bamoutput_coordonealign(&mut out, &mut p, &[0, 1, 2], 3, 0);
    assert!(short.is_err());

    let oversized = bamoutput_l77_bamoutput_coordonealign(&mut out, &mut p, &[0, 1, 2, 3], 20, 0);
    assert!(oversized.is_err());

    let mut overflow = bamoutput_l9_bamoutput_bamoutput(0, "/tmp/bam", &p);
    overflow.bin_bytes[0] = u64::MAX - 1;
    let record = bam_record(0, 10, 1);
    let err = bamoutput_l77_bamoutput_coordonealign(
        &mut overflow,
        &mut p,
        &record,
        record.len() as u64,
        0,
    )
    .unwrap_err();
    assert!(err.contains("byte count overflow"));
}

#[test]
fn bamoutput_coordinate_unmapped_prepare_flushes_last_bin_to_sjout() {
    let p = Parameters {
        out_bam_coord_nbins: 3,
        chunk_out_bam_size_bytes: 90,
        out_bam_sorting_bin_start: vec![0, 100, 0],
        ..Default::default()
    };
    let mut out = BAMoutput {
        n_bins: 3,
        bin_buffers: vec![Vec::new(), Vec::new(), vec![7, 8, 9]],
        bin_streams: vec![Vec::new(), Vec::new(), vec![1]],
        bin_stream_by_sjout: vec![false; 3],
        bin_bytes: vec![0, 0, 3],
        ..Default::default()
    };

    bamoutput_l179_bamoutput_coordunmappedpreparebysjout(&mut out, &p);

    assert_eq!(out.bin_streams[2], vec![1, 7, 8, 9]);
    assert!(out.bin_buffers[2].is_empty());
    assert_eq!(out.bin_bytes[2], 0);
    assert!(out.bin_stream_by_sjout[2]);
}

#[test]
fn bambinsortbycoordinate_sorts_loaded_thread_bins_and_strips_read_order_trailer() {
    let p = Parameters {
        sam_header_sorted_coord: "@HD\tSO:coordinate\n".to_string(),
        ..Default::default()
    };
    let genome = Genome {
        chr_name_all: vec!["chr1".to_string()],
        chr_length_all: vec![1000],
        ..Default::default()
    };
    let r_late = bam_record(0, 500, 5);
    let r_tie_second = bam_record(0, 100, 2);
    let r_tie_first = bam_record(0, 100, 1);
    let r_mid = bam_record(0, 300, 3);
    let files = vec![
        {
            let mut v = Vec::new();
            v.extend_from_slice(&bam_record_with_read_order(0, 500, 5, 40));
            v.extend_from_slice(&bam_record_with_read_order(0, 100, 2, 20));
            v
        },
        {
            let mut v = Vec::new();
            v.extend_from_slice(&bam_record_with_read_order(0, 100, 1, 10));
            v.extend_from_slice(&bam_record_with_read_order(0, 300, 3, 30));
            v
        },
    ];
    let bin_s = files.iter().map(|v| v.len() as u64).sum();

    let (sorted, removed) = bambinsortbycoordinate_l7_bambinsortbycoordinate(
        2,
        4,
        bin_s,
        2,
        "/tmp/sort/",
        &p,
        &genome,
        &files,
    )
    .unwrap();

    let header =
        bamfunctions_l77_outbamwriteheader("@HD\tSO:coordinate\n", &["chr1".to_string()], &[1000]);
    let mut expected = header;
    expected.extend_from_slice(&r_tie_first);
    expected.extend_from_slice(&r_tie_second);
    expected.extend_from_slice(&r_mid);
    expected.extend_from_slice(&r_late);
    assert_eq!(sorted, expected);
    assert_eq!(removed, vec!["/tmp/sort/0/2", "/tmp/sort/1/2"]);
}

#[test]
fn bambinsortbycoordinate_returns_empty_for_empty_bins_and_errors_on_size_mismatch() {
    let p = Parameters::default();
    let genome = Genome::default();
    assert_eq!(
        bambinsortbycoordinate_l7_bambinsortbycoordinate(1, 0, 0, 2, "d", &p, &genome, &[])
            .unwrap(),
        (Vec::new(), Vec::new())
    );

    let err = bambinsortbycoordinate_l7_bambinsortbycoordinate(
        1,
        1,
        99,
        1,
        "d",
        &p,
        &genome,
        &[bam_record_with_read_order(0, 1, 0, 1)],
    )
    .unwrap_err();
    assert!(err.contains("Expected bin size=99"));

    let err = bambinsortbycoordinate_l7_bambinsortbycoordinate(
        1,
        1,
        2,
        1,
        "d",
        &p,
        &genome,
        &[vec![1, 2]],
    )
    .unwrap_err();
    assert!(err.contains("truncated temporary BAM coordinate bin"));

    let err = bambinsortbycoordinate_l7_bambinsortbycoordinate(
        1,
        1,
        4,
        1,
        "d",
        &p,
        &genome,
        &[u32::MAX.to_le_bytes().to_vec()],
    )
    .unwrap_err();
    assert!(err.contains("truncated temporary BAM coordinate bin"));
}

#[test]
fn bambinsortunmapped_merges_thread_and_bysjout_files_by_top_read_order() {
    let p = Parameters {
        sam_header_sorted_coord: "@HD\tSO:coordinate\n".to_string(),
        ..Default::default()
    };
    let genome = Genome {
        chr_name_all: vec!["chr1".to_string()],
        chr_length_all: vec![1000],
        ..Default::default()
    };
    let r10 = bam_record(0, 10, 10);
    let r20 = bam_record(0, 20, 20);
    let r30 = bam_record(0, 30, 30);
    let r40 = bam_record(0, 40, 40);
    let files = vec![
        {
            let mut v = Vec::new();
            v.extend_from_slice(&bam_record_with_unmapped_order(0, 30, 30, 30));
            v.extend_from_slice(&bam_record_with_unmapped_order(0, 40, 40, 40));
            v
        },
        {
            let mut v = Vec::new();
            v.extend_from_slice(&bam_record_with_unmapped_order(0, 10, 10, 10));
            v
        },
        Vec::new(),
        {
            let mut v = Vec::new();
            v.extend_from_slice(&bam_record_with_unmapped_order(0, 20, 20, 20));
            v
        },
    ];

    let (sorted, removed) =
        bambinsortunmapped_l5_bambinsortunmapped(7, 2, "/tmp/sort/", &p, &genome, &files).unwrap();

    let header =
        bamfunctions_l77_outbamwriteheader("@HD\tSO:coordinate\n", &["chr1".to_string()], &[1000]);
    let mut expected = header;
    expected.extend_from_slice(&r10);
    expected.extend_from_slice(&r20);
    expected.extend_from_slice(&r30);
    expected.extend_from_slice(&r40);
    assert_eq!(sorted, expected);
    assert_eq!(
        removed,
        vec![
            "/tmp/sort/0/7",
            "/tmp/sort/0/7.BySJout",
            "/tmp/sort/1/7",
            "/tmp/sort/1/7.BySJout"
        ]
    );
}

#[test]
fn bambinsortunmapped_reports_truncated_temp_records() {
    let p = Parameters::default();
    let genome = Genome::default();
    let err =
        bambinsortunmapped_l5_bambinsortunmapped(1, 1, "d", &p, &genome, &[vec![12, 0, 0, 0]])
            .unwrap_err();
    assert!(err.contains("truncated temporary bam file"));

    let err = bambinsortunmapped_l5_bambinsortunmapped(1, 1, "d", &p, &genome, &[vec![1, 2]])
        .unwrap_err();
    assert!(err.contains("truncated temporary bam file"));

    let err = bambinsortunmapped_l5_bambinsortunmapped(
        1,
        1,
        "d",
        &p,
        &genome,
        &[u32::MAX.to_le_bytes().to_vec()],
    )
    .unwrap_err();
    assert!(err.contains("malformed temporary bam file"));
}

#[test]
fn bamsortbycoordinate_orchestrates_bins_and_concatenates_existing_bin_outputs() {
    let mut p = Parameters {
        out_bam_coord: true,
        run_thread_n: 2,
        out_bam_coord_nbins: 3,
        limit_bam_sort_ram: 10_000,
        out_bam_sort_tmp_dir: "/tmp/star-sort".to_string(),
        sam_header_sorted_coord: "@HD\tSO:coordinate\n".to_string(),
        ..Default::default()
    };
    let genome = Genome {
        chr_name_all: vec!["chr1".to_string()],
        chr_length_all: vec![1000],
        ..Default::default()
    };
    let c1 = bam_record_with_read_order(0, 200, 2, 20);
    let c2 = bam_record_with_read_order(0, 100, 1, 10);
    let u1 = bam_record_with_unmapped_order(u32::MAX, 0, 9, 30);
    let u2 = bam_record_with_unmapped_order(u32::MAX, 0, 8, 40);
    let ra_chunks = vec![
        ReadAlignChunk {
            chunk_out_bam_coord: BAMoutput {
                bin_total_n: vec![1, 0, 1],
                bin_total_bytes: vec![c1.len() as u64, 0, u1.len() as u64],
                ..Default::default()
            },
            ..Default::default()
        },
        ReadAlignChunk {
            chunk_out_bam_coord: BAMoutput {
                bin_total_n: vec![1, 0, 1],
                bin_total_bytes: vec![c2.len() as u64, 0, u2.len() as u64],
                ..Default::default()
            },
            ..Default::default()
        },
    ];
    let temp_files_by_bin = vec![
        vec![c1.clone(), c2.clone()],
        Vec::new(),
        vec![u1.clone(), Vec::new(), Vec::new(), u2.clone()],
    ];

    let result =
        bamsortbycoordinate_l8_bamsortbycoordinate(&mut p, &ra_chunks, &genome, &temp_files_by_bin)
            .unwrap();

    let header =
        bamfunctions_l77_outbamwriteheader("@HD\tSO:coordinate\n", &["chr1".to_string()], &[1000]);
    let mut bin0 = header.clone();
    bin0.extend_from_slice(&bam_record(0, 100, 1));
    bin0.extend_from_slice(&bam_record(0, 200, 2));
    let mut bin2 = header.clone();
    bin2.extend_from_slice(&bam_record(u32::MAX, 0, 9));
    bin2.extend_from_slice(&bam_record(u32::MAX, 0, 8));
    let mut expected = bin0.clone();
    expected.extend_from_slice(&bin2);

    assert_eq!(result.max_mem, c1.len() as u64 + c2.len() as u64 + 48);
    assert_eq!(result.unmapped_reads_n, 2);
    assert_eq!(
        result.bin_names,
        vec!["/tmp/star-sort/b0", "/tmp/star-sort/b2"]
    );
    assert_eq!(result.bin_outputs[0], bin0);
    assert_eq!(result.bin_outputs[1], Vec::<u8>::new());
    assert_eq!(result.bin_outputs[2], bin2);
    assert_eq!(result.output_bam, expected);
    assert!(p.bam_sorting_log.contains("Max memory needed for sorting"));
}

#[test]
fn bamsortbycoordinate_handles_noop_empty_and_memory_error_paths() {
    let genome = Genome::default();
    let mut p_off = Parameters {
        out_bam_coord: false,
        ..Default::default()
    };
    assert_eq!(
        bamsortbycoordinate_l8_bamsortbycoordinate(&mut p_off, &[], &genome, &[])
            .unwrap()
            .output_bam,
        Vec::<u8>::new()
    );

    let mut p_empty = Parameters {
        out_bam_coord: true,
        run_thread_n: 1,
        out_bam_coord_nbins: 2,
        limit_bam_sort_ram: 10,
        sam_header_sorted_coord: "@HD\n".to_string(),
        ..Default::default()
    };
    let empty_chunk = vec![ReadAlignChunk {
        chunk_out_bam_coord: BAMoutput {
            bin_total_n: vec![0, 0],
            bin_total_bytes: vec![0, 0],
            ..Default::default()
        },
        ..Default::default()
    }];
    let empty = bamsortbycoordinate_l8_bamsortbycoordinate(
        &mut p_empty,
        &empty_chunk,
        &genome,
        &[Vec::new(), Vec::new()],
    )
    .unwrap();
    assert_eq!(
        empty.output_bam,
        bamfunctions_l77_outbamwriteheader("@HD\n", &[], &[])
    );
    assert!(p_empty.bam_sorting_log.contains("WARNING: nothing to sort"));

    let mut p_small = Parameters {
        out_bam_coord: true,
        run_thread_n: 1,
        out_bam_coord_nbins: 2,
        limit_bam_sort_ram: 1,
        ..Default::default()
    };
    let too_large = vec![ReadAlignChunk {
        chunk_out_bam_coord: BAMoutput {
            bin_total_n: vec![1, 0],
            bin_total_bytes: vec![10, 0],
            ..Default::default()
        },
        ..Default::default()
    }];
    let err = bamsortbycoordinate_l8_bamsortbycoordinate(
        &mut p_small,
        &too_large,
        &genome,
        &[Vec::new(), Vec::new()],
    )
    .unwrap_err();
    assert!(err.contains("not enough memory for BAM sorting"));
}

#[test]
fn out_bam_write_header_matches_bam_binary_header_layout() {
    let header = bamfunctions_l77_outbamwriteheader(
        "@HD\tVN:1.6\n",
        &["chr1".to_string(), "chrM".to_string()],
        &[248, 16],
    );
    let mut expected = Vec::new();
    expected.extend_from_slice(b"BAM\x01");
    expected.extend_from_slice(&(11i32).to_ne_bytes());
    expected.extend_from_slice(b"@HD\tVN:1.6\n");
    expected.extend_from_slice(&(2i32).to_ne_bytes());
    expected.extend_from_slice(&(5i32).to_ne_bytes());
    expected.extend_from_slice(b"chr1\0");
    expected.extend_from_slice(&(248i32).to_ne_bytes());
    expected.extend_from_slice(&(5i32).to_ne_bytes());
    expected.extend_from_slice(b"chrM\0");
    expected.extend_from_slice(&(16i32).to_ne_bytes());
    assert_eq!(header, expected);
}

#[test]
fn bam_attr_scalar_writers_match_bam_aux_layout() {
    let mut buf = [0u8; 32];
    assert_eq!(bamfunctions_l106_bamattrarraywrite(-17, b"NM", &mut buf), 7);
    assert_eq!(&buf[..3], b"NMi");
    assert_eq!(i32::from_ne_bytes(buf[3..7].try_into().unwrap()), -17);

    assert_eq!(bamfunctions_l112_bamattrarraywrite(1.5, b"AS", &mut buf), 7);
    assert_eq!(&buf[..3], b"ASf");
    assert_eq!(f32::from_ne_bytes(buf[3..7].try_into().unwrap()), 1.5);

    assert_eq!(
        bamfunctions_l118_bamattrarraywrite(b'Z', b"ch", &mut buf),
        4
    );
    assert_eq!(&buf[..4], b"chAZ");
}

#[test]
fn bam_attr_string_and_array_writers_match_bam_aux_layout() {
    let mut buf = [0u8; 64];
    assert_eq!(
        bamfunctions_l124_bamattrarraywrite("abc", b"ZZ", &mut buf),
        7
    );
    assert_eq!(&buf[..7], b"ZZZabc\0");

    assert_eq!(
        bamfunctions_l130_bamattrarraywrite(&[1, 2, 255], b"BC", &mut buf),
        11
    );
    assert_eq!(&buf[..4], b"BCBc");
    assert_eq!(i32::from_ne_bytes(buf[4..8].try_into().unwrap()), 3);
    assert_eq!(&buf[8..11], &[1, 2, 255]);

    assert_eq!(
        bamfunctions_l138_bamattrarraywrite(&[-1, 7], b"BI", &mut buf),
        16
    );
    assert_eq!(&buf[..4], b"BIBi");
    assert_eq!(i32::from_ne_bytes(buf[4..8].try_into().unwrap()), 2);
    assert_eq!(i32::from_ne_bytes(buf[8..12].try_into().unwrap()), -1);
    assert_eq!(i32::from_ne_bytes(buf[12..16].try_into().unwrap()), 7);
}

#[test]
fn bam_attr_array_write_sam_tags_filters_and_encodes_supported_tags() {
    let mut buf = [0u8; 128];
    let n = bamfunctions_l147_bamattrarraywritesamtags(
        "\tNM:i:-7\tXX:Z:skip\tch:A:Y\tAS:f:1.25\tZZ:Z:abc\t",
        &mut buf,
        false,
        &[
            u16::from_ne_bytes(*b"NM"),
            u16::from_ne_bytes(*b"ch"),
            u16::from_ne_bytes(*b"ZZ"),
        ],
    )
    .unwrap();
    assert_eq!(n, 18);
    assert_eq!(&buf[..3], b"NMi");
    assert_eq!(i32::from_ne_bytes(buf[3..7].try_into().unwrap()), -7);
    assert_eq!(&buf[7..11], b"chAY");
    assert_eq!(&buf[11..18], b"ZZZabc\0");

    let mut all_buf = [0u8; 64];
    let n_all =
        bamfunctions_l147_bamattrarraywritesamtags("AS:f:1.25\tXY:i:2", &mut all_buf, true, &[])
            .unwrap();
    assert_eq!(n_all, 14);
    assert_eq!(&all_buf[..3], b"ASf");
    assert_eq!(f32::from_ne_bytes(all_buf[3..7].try_into().unwrap()), 1.25);
    assert_eq!(&all_buf[7..10], b"XYi");
    assert_eq!(i32::from_ne_bytes(all_buf[10..14].try_into().unwrap()), 2);
}

#[test]
fn bam_attr_array_write_sam_tags_reports_malformed_tags_without_panic() {
    let mut buf = [0u8; 64];
    assert!(
        bamfunctions_l147_bamattrarraywritesamtags("N", &mut buf, true, &[])
            .unwrap_err()
            .contains("malformed SAM attribute")
    );
    assert!(
        bamfunctions_l147_bamattrarraywritesamtags("NM:i:not_an_int", &mut buf, true, &[])
            .unwrap_err()
            .contains("malformed SAM integer attribute")
    );
    assert!(
        bamfunctions_l147_bamattrarraywritesamtags("AS:f:not_a_float", &mut buf, true, &[])
            .unwrap_err()
            .contains("malformed SAM float attribute")
    );
}

#[test]
fn solo_feature_add_bam_tags_appends_cb_ub_and_updates_record_size() {
    let mut bam0 = Vec::new();
    bam0.extend_from_slice(&8u32.to_ne_bytes());
    bam0.extend_from_slice(b"ABCDEFGH");
    let size0 = bam0.len() as u32;
    bam0.extend_from_slice(&(2u64 << 32).to_ne_bytes());

    let read_info = vec![
        SoloFeatureReadInfo {
            cb: -1,
            umi: u64::MAX,
        },
        SoloFeatureReadInfo { cb: 0, umi: 0 },
        SoloFeatureReadInfo {
            cb: 1,
            umi: 0b00_01_10_11,
        },
    ];
    let cb_wlstr = vec!["AAAA".to_string(), "TGCA".to_string()];

    let out =
        solofeature_addbamtags_l5_solofeature_addbamtags(&bam0, size0, &read_info, &cb_wlstr, 4);
    assert_eq!(
        u32::from_ne_bytes(out[..4].try_into().unwrap()),
        out.len() as u32 - 4
    );
    assert_eq!(&out[4..12], b"ABCDEFGH");
    assert_eq!(&out[12..20], b"CBZTGCA\0");
    assert_eq!(&out[20..28], b"UBZACGT\0");

    let missing_index = 0u64 << 32;
    bam0[size0 as usize..size0 as usize + 8].copy_from_slice(&missing_index.to_ne_bytes());
    let missing =
        solofeature_addbamtags_l5_solofeature_addbamtags(&bam0, size0, &read_info, &cb_wlstr, 4);
    assert_eq!(&missing[12..17], b"CBZ-\0");
    assert_eq!(&missing[17..22], b"UBZ-\0");
}
