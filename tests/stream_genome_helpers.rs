use std::io::{Cursor, Seek, Write};

use star_rs::generated::functions::*;
use star_rs::generated::structs::{
    ClipMate, GTF, Genome, PackedArray, Parameters, ParametersGenome, ParametersGenomeTransform,
    VariantInfo,
};

#[test]
fn genome_transform_vector_helpers_match_cpp_value_semantics() {
    let mut v1 = vec![1, 2];
    genome_transformgenome_l10_appendvector(&mut v1, &[3, 4]);
    assert_eq!(v1, vec![1, 2, 3, 4]);

    let joined = genome_transformgenome_l17_concatenatevectors(&v1, &[5, 6]);
    assert_eq!(joined, vec![1, 2, 3, 4, 5, 6]);
    assert_eq!(v1, vec![1, 2, 3, 4]);

    let suffixed = genome_transformgenome_l26_appendstring(
        vec!["chr1".to_string(), "chr2".to_string()],
        "_hap1",
    );
    assert_eq!(suffixed, vec!["chr1_hap1", "chr2_hap1"]);
}

#[test]
fn genome_transform_chr_len_start_sorts_filters_and_repads_chromosomes() {
    let genome = Genome {
        chr_name: vec!["chr1".to_string(), "chr2".to_string()],
        chr_start: vec![0, 32, 64],
        chr_length: vec![10, 20],
        genome_chr_bin_nbases: 8,
        ..Default::default()
    };
    let mut variants = std::collections::BTreeMap::new();
    variants.insert(
        "chr1".to_string(),
        vec![
            VariantInfo {
                pos: 8,
                len: 1,
                seq: ["A".to_string(), "AG".to_string()],
            },
            VariantInfo {
                pos: 4,
                len: -2,
                seq: ["ACG".to_string(), "A".to_string()],
            },
            VariantInfo {
                pos: 5,
                len: 0,
                seq: ["C".to_string(), "T".to_string()],
            },
        ],
    );
    let mut chr_start1 = Vec::new();
    let mut chr_length1 = Vec::new();

    let log = genome_transformgenome_l171_genome_transformchrlenstart(
        &genome,
        &mut variants,
        &mut chr_start1,
        &mut chr_length1,
    );

    assert_eq!(
        variants["chr1"].iter().map(|v| v.pos).collect::<Vec<_>>(),
        vec![4, 8]
    );
    assert_eq!(chr_length1, vec![9, 20]);
    assert_eq!(chr_start1, vec![0, 16, 40]);
    assert_eq!(
        log,
        "chr1: filtered out overlapping variants = 1; remaining variants = 2\n\
Transformed chr length difference: chr1 -1\n\
Transformed chr start difference: chr1 0\n\
Transformed chr start difference: chr2 -16\n"
    );
}

#[test]
fn genome_transform_blocks_write_matches_reversed_tsv_columns() {
    assert_eq!(
        genome_transformgenome_l271_genome_transformblockswrite(&[[10, 5, 100], [20, 7, 200]]),
        "2\t-1\n100\t5\t10\n200\t7\t20\n"
    );
}

#[test]
fn genome_transform_g_and_blocks_rewrites_sequence_and_conversion_blocks() {
    let genome = Genome {
        g: vec![0, 1, 2, 3, 0, 1, 4, 4, 4, 4, 3, 2, 1, 0],
        chr_name: vec!["chr1".to_string(), "chr2".to_string()],
        chr_start: vec![0, 10, 14],
        chr_length: vec![6, 4],
        ..Default::default()
    };
    let mut variants = std::collections::BTreeMap::new();
    variants.insert(
        "chr1".to_string(),
        vec![
            VariantInfo {
                pos: 2,
                len: 1,
                seq: ["C".to_string(), "CT".to_string()],
            },
            VariantInfo {
                pos: 5,
                len: -1,
                seq: ["AC".to_string(), "G".to_string()],
            },
        ],
    );
    let chr_start1 = vec![0, 10, 14];
    let chr_length1 = vec![6, 4];
    let mut blocks = Vec::new();
    let mut g_new = vec![9u8; 14];

    let debug = genome_transformgenome_l215_genome_transformgandblocks(
        &genome,
        &variants,
        &chr_start1,
        &chr_length1,
        &mut blocks,
        &mut g_new,
    );

    assert_eq!(debug, None);
    assert_eq!(&g_new[..6], &[0, 1, 3, 2, 3, 2]);
    assert_eq!(&g_new[10..14], &[3, 2, 1, 0]);
    assert_eq!(blocks, vec![[0, 2, 0], [2, 3, 3], [6, 0, 6], [10, 4, 10]]);
}

#[test]
fn genome_write_chr_info_emits_original_text_files() {
    let dir = std::env::temp_dir().join(format!("star_rs_chr_info_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let genome = Genome {
        n_chr_real: 2,
        chr_name: vec!["chr1".to_string(), "chrM".to_string()],
        chr_start: vec![0, 112, 144],
        chr_length: vec![100, 17],
        ..Default::default()
    };

    genome_genomegenerate_l417_genome_writechrinfo(&genome, dir.to_str().unwrap()).unwrap();

    assert_eq!(
        std::fs::read_to_string(dir.join("chrName.txt")).unwrap(),
        "chr1\nchrM\n"
    );
    assert_eq!(
        std::fs::read_to_string(dir.join("chrStart.txt")).unwrap(),
        "0\n112\n144\n"
    );
    assert_eq!(
        std::fs::read_to_string(dir.join("chrLength.txt")).unwrap(),
        "100\n17\n"
    );
    assert_eq!(
        std::fs::read_to_string(dir.join("chrNameLength.txt")).unwrap(),
        "chr1\t100\nchrM\t17\n"
    );
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn genome_write_genome_sequence_writes_n_genome_prefix() {
    let dir = std::env::temp_dir().join(format!("star_rs_genome_seq_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let genome = Genome {
        g: b"ACGTNextra".to_vec(),
        n_genome: 5,
        ..Default::default()
    };

    genome_genomegenerate_l433_genome_writegenomesequence(&genome, dir.to_str().unwrap()).unwrap();

    assert_eq!(std::fs::read(dir.join("Genome")).unwrap(), b"ACGTN");
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn genome_generate_writes_tiny_fasta_index_files() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_genome_generate_{}_{}",
        std::process::id(),
        unique
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let fasta = dir.join("tiny.fa");
    std::fs::write(&fasta, ">chr1 description\nACGT\n").unwrap();

    let p_ge = ParametersGenome {
        g_dir: dir.to_str().unwrap().to_string(),
        g_type_string: "Full".to_string(),
        g_fasta_files: vec![fasta.to_str().unwrap().to_string()],
        g_chr_bin_nbits: 2,
        g_saindex_nbases: 0,
        g_sasparse_d: 1,
        g_suffix_length_max: 8,
        sjdb_file_chr_start_end: vec!["-".to_string()],
        sjdb_gtf_file: "-".to_string(),
        transform: ParametersGenomeTransform {
            type_string: "None".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    let mut genome = Genome {
        p_ge,
        ..Default::default()
    };
    let mut p = Parameters {
        command_line_full: "STAR --runMode genomeGenerate".to_string(),
        version_genome: "test-version".to_string(),
        run_thread_n: 1,
        limit_genome_generate_ram: 1_000_000,
        limit_sjdb_insert_nsj: 0,
        p_ge: genome.p_ge.clone(),
        ..Default::default()
    };

    let result = genome_genomegenerate_l98_genome_genomegenerate(&mut genome, &mut p, None)
        .expect("tiny genome generation should succeed");

    assert_eq!(result.n_genome_true, 4);
    assert_eq!(genome.n_chr_real, 1);
    assert_eq!(genome.chr_name, vec!["chr1".to_string()]);
    assert_eq!(genome.chr_length, vec![4]);
    assert_eq!(genome.chr_start, vec![0, 8]);
    assert_eq!(genome.n_genome, 8);
    assert_eq!(&genome.g[..4], &[0, 1, 2, 3]);
    assert!(genome.n_sa > 0);
    assert_eq!(genome.sa_packed.length, genome.n_sa as u64);
    assert!(genome.sai.is_empty());
    assert!(
        result
            .log_main
            .contains("starting to generate Genome files")
    );
    assert!(result.log_main.contains("finished successfully"));

    assert_eq!(
        std::fs::read_to_string(dir.join("chrName.txt")).unwrap(),
        "chr1\n"
    );
    assert_eq!(
        std::fs::read_to_string(dir.join("chrLength.txt")).unwrap(),
        "4\n"
    );
    assert_eq!(
        std::fs::read_to_string(dir.join("chrStart.txt")).unwrap(),
        "0\n8\n"
    );
    assert_eq!(
        std::fs::read_to_string(dir.join("chrNameLength.txt")).unwrap(),
        "chr1\t4\n"
    );
    assert_eq!(
        std::fs::read(dir.join("Genome")).unwrap(),
        vec![
            0,
            1,
            2,
            3,
            GENOME_SPACING_CHAR,
            GENOME_SPACING_CHAR,
            GENOME_SPACING_CHAR,
            GENOME_SPACING_CHAR
        ]
    );
    assert!(std::fs::metadata(dir.join("SA")).unwrap().len() > 0);
    assert!(std::fs::metadata(dir.join("SAindex")).unwrap().len() > 0);
    assert!(
        std::fs::read_to_string(dir.join("genomeParameters.txt"))
            .unwrap()
            .contains("genomeFastaFiles")
    );

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn genome_generate_saindex_one_uses_sorted_suffix_array() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_genome_generate_saindex1_{}_{}",
        std::process::id(),
        unique
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let fasta = dir.join("tiny.fa");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();

    let p_ge = ParametersGenome {
        g_dir: dir.to_str().unwrap().to_string(),
        g_type_string: "Full".to_string(),
        g_fasta_files: vec![fasta.to_str().unwrap().to_string()],
        g_chr_bin_nbits: 18,
        g_saindex_nbases: 1,
        g_sasparse_d: 1,
        g_suffix_length_max: 8,
        sjdb_file_chr_start_end: vec!["-".to_string()],
        sjdb_gtf_file: "-".to_string(),
        transform: ParametersGenomeTransform {
            type_string: "None".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    let mut genome = Genome {
        p_ge,
        ..Default::default()
    };
    let mut p = Parameters {
        command_line_full: "STAR --runMode genomeGenerate".to_string(),
        version_genome: "test-version".to_string(),
        run_thread_n: 1,
        limit_genome_generate_ram: 1_000_000,
        limit_sjdb_insert_nsj: 0,
        p_ge: genome.p_ge.clone(),
        ..Default::default()
    };

    genome_genomegenerate_l98_genome_genomegenerate(&mut genome, &mut p, None)
        .expect("genomeSAindexNbases=1 should not produce a decreasing SA index");

    assert_eq!(genome.genome_sa_index_start, vec![0, 4]);
    assert_eq!(genome.sai.len(), 4);
    assert!(genome.sai.iter().any(|&v| v != 0));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn genome_generate_requires_overhang_for_annotation_inputs() {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "star_rs_genome_generate_overhang_{}_{}",
        std::process::id(),
        unique
    ));
    let mut genome = Genome {
        p_ge: ParametersGenome {
            g_dir: dir.to_str().unwrap().to_string(),
            g_fasta_files: vec!["unused.fa".to_string()],
            sjdb_file_chr_start_end: vec!["junctions.tab".to_string()],
            sjdb_gtf_file: "-".to_string(),
            ..Default::default()
        },
        sjdb_overhang: 0,
        ..Default::default()
    };
    let mut p = Parameters {
        limit_genome_generate_ram: 10_000,
        p_ge: genome.p_ge.clone(),
        ..Default::default()
    };

    let err =
        genome_genomegenerate_l98_genome_genomegenerate(&mut genome, &mut p, None).unwrap_err();
    assert!(err.contains("for generating genome with annotations"));

    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn genome_concatenate_chromosomes_pads_indexes_and_adds_reverse_complement() {
    let mut genome = Genome::default();
    let seq = vec![vec![0, 1, 2], vec![3, 4]];
    let names = vec!["chrA".to_string(), "chrB".to_string()];

    gtf_supertranscript_l232_genome_concatenatechromosomes(&mut genome, &seq, &names, 4);

    assert_eq!(genome.n_chr_real, 2);
    assert_eq!(genome.chr_length, vec![3, 2]);
    assert_eq!(genome.chr_start, vec![0, 8, 12]);
    assert_eq!(genome.chr_name, names);
    assert_eq!(genome.chr_name_index.get("chrA"), Some(&0));
    assert_eq!(genome.chr_name_index.get("chrB"), Some(&1));
    assert_eq!(genome.n_genome, 12);
    assert_eq!(&genome.g[0..3], &[0, 1, 2]);
    assert_eq!(&genome.g[3..8], &[GENOME_SPACING_CHAR; 5]);
    assert_eq!(&genome.g[8..10], &[3, 4]);
    assert_eq!(&genome.g[10..12], &[GENOME_SPACING_CHAR; 2]);
    assert_eq!(&genome.g[12..24], &[5, 5, 4, 0, 5, 5, 5, 5, 5, 1, 2, 3]);
}

#[test]
fn genome_scan_fasta_files_counts_metadata_then_loads_numeric_sequence() {
    let dir = std::env::temp_dir().join(format!("star_rs_scan_fasta_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fasta = dir.join("genome.fa");
    std::fs::write(&fasta, ">chr1 description\nACGT\n>chr2\nNN\n").unwrap();

    let mut genome = Genome {
        genome_chr_bin_nbases: 4,
        p_ge: ParametersGenome {
            g_fasta_files: vec![fasta.to_string_lossy().to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    let p = Parameters::default();
    let mut g = vec![GENOME_SPACING_CHAR; 32];

    let n_scan =
        genomescanfastafiles_l5_genomescanfastafiles(&p, &mut g, false, &mut genome).unwrap();
    assert_eq!(n_scan, 12);
    assert_eq!(genome.n_chr_real, 2);
    assert_eq!(
        genome.chr_name,
        vec!["chr1".to_string(), "chr2".to_string()]
    );
    assert_eq!(genome.chr_start, vec![0, 8, 12]);
    assert_eq!(genome.chr_length, vec![4, 2]);
    assert_eq!(genome.chr_name_index.get("chr1"), Some(&0));
    assert_eq!(genome.chr_name_index.get("chr2"), Some(&1));

    let n_run =
        genomescanfastafiles_l5_genomescanfastafiles(&p, &mut g, true, &mut genome).unwrap();
    assert_eq!(n_run, 12);
    assert_eq!(&g[0..4], &[0, 1, 2, 3]);
    assert_eq!(&g[4..8], &[GENOME_SPACING_CHAR; 4]);
    assert_eq!(&g[8..10], &[4, 4]);

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn genome_insert_sequences_moves_sjdb_loads_fasta_and_updates_packed_sa() {
    let dir = std::env::temp_dir().join(format!("star_rs_insert_sequences_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fasta = dir.join("insert.fa");
    std::fs::write(&fasta, ">extra\nN\n").unwrap();

    let mut sa = packedarray_l3_packedarray_packedarray();
    packedarray_l8_packedarray_definebits(&mut sa, 33, 4);
    packedarray_l31_packedarray_allocatearray(&mut sa);
    let n2bit = 1u64 << 32;
    packedarray_l17_packedarray_writepacked(&mut sa, 0, 0);
    packedarray_l17_packedarray_writepacked(&mut sa, 1, 2);
    packedarray_l17_packedarray_writepacked(&mut sa, 2, n2bit | 1);
    packedarray_l17_packedarray_writepacked(&mut sa, 3, n2bit | 2);

    let mut genome = Genome {
        g: vec![0, 1, 7, 8],
        n_genome: 4,
        genome_chr_bin_nbases: 4,
        chr_start: vec![0, 3],
        genome_insert_l: 1,
        sa_packed: sa,
        sa_insert: PackedArray::default(),
        sai_packed: PackedArray::default(),
        p_ge: ParametersGenome {
            g_fasta_files: vec![fasta.to_string_lossy().to_string()],
            ..Default::default()
        },
        ..Default::default()
    };

    let n_ind =
        genome_insertsequences_l9_genome_insertsequences(&mut genome, &Parameters::default())
            .unwrap();

    assert_eq!(n_ind, 0);
    assert_eq!(genome.n_genome, 5);
    assert_eq!(genome.g, vec![0, 1, 4, 7, 8]);
    assert_eq!(packedarray_h18_packedarray_index(&genome.sa_packed, 1), 3);
    assert_eq!(
        packedarray_h18_packedarray_index(&genome.sa_packed, 3),
        n2bit | 3
    );
    assert_eq!(genome.sa_insert.length, 4);
    assert_eq!(packedarray_h18_packedarray_index(&genome.sa_insert, 1), 3);
}

#[test]
fn genome_scan_fasta_files_rejects_non_fasta_input() {
    let dir = std::env::temp_dir().join(format!("star_rs_scan_bad_fasta_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fasta = dir.join("bad.fa");
    std::fs::write(&fasta, "not fasta\n").unwrap();

    let mut genome = Genome {
        genome_chr_bin_nbases: 4,
        p_ge: ParametersGenome {
            g_fasta_files: vec![fasta.to_string_lossy().to_string()],
            ..Default::default()
        },
        ..Default::default()
    };
    let mut g = vec![0; 16];
    let err = genomescanfastafiles_l5_genomescanfastafiles(
        &Parameters::default(),
        &mut g,
        false,
        &mut genome,
    )
    .unwrap_err();
    assert!(err.contains("the file format of the genomeFastaFile"));

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn genome_transform_exon_loci_shifts_edges_and_drops_inverted_exons() {
    let blocks = vec![[0, 5, 100], [10, 5, 200], [20, 5, 300]];
    let mut exons = vec![
        [1, 2, 3, 10, 20],
        [2, 6, 12, 11, 21],
        [3, 14, 18, 12, 22],
        [4, 16, 18, 13, 23],
    ];

    let log = genome_transformgenome_l282_genome_transformexonloci(&mut exons, &blocks);

    assert_eq!(log, "Transform exons: removed 1\n");
    assert_eq!(
        exons,
        vec![
            [1, 102, 103, 10, 20],
            [2, 200, 202, 11, 21],
            [3, 204, 204, 12, 22],
        ]
    );
}

#[test]
fn genome_transform_genome_applies_haploid_vcf_and_updates_annotations() {
    let mut genome = Genome {
        g: vec![0, 1, 2, 3, 0, 4, 4, 4],
        n_genome: 8,
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_name_index: std::collections::BTreeMap::from([("chr1".to_string(), 0)]),
        chr_start: vec![0, 8],
        chr_length: vec![5],
        genome_chr_bin_nbases: 4,
        p_ge: ParametersGenome {
            transform: ParametersGenomeTransform {
                type_: 1,
                vcf_file: "in.vcf".to_string(),
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    let mut gtf = GTF {
        exon_loci: vec![[0, 1, 4, 0]],
        transcript_id: vec!["tr".to_string()],
        gene_id: vec!["gene".to_string()],
        ..Default::default()
    };

    let result = genome_transformgenome_l33_genome_transformgenome(
        &mut genome,
        &mut gtf,
        "chr1\t3\t.\tG\tGT\t.\t.\t.\tGT\t1\n",
    );

    assert_eq!(&genome.g[..6], &[0, 1, 2, 3, 3, 0]);
    assert_eq!(genome.n_genome, 8);
    assert_eq!(genome.chr_start, vec![0, 8]);
    assert_eq!(genome.chr_length, vec![6]);
    assert_eq!(gtf.exon_loci, vec![[0, 1, 5, 0]]);
    assert_eq!(result.transform_blocks_tsv, "2\t-1\n0\t3\t0\n4\t2\t3\n");
    assert!(result.log_main.contains("Old/new genome sizes: 8 8\n"));
}

#[test]
fn genome_transform_genome_builds_diploid_haplotypes_and_suffixes_ids() {
    let mut genome = Genome {
        g: vec![0, 1, 2, 3, 0, 4, 4, 4],
        n_genome: 8,
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_name_index: std::collections::BTreeMap::from([("chr1".to_string(), 0)]),
        chr_start: vec![0, 8],
        chr_length: vec![5],
        genome_chr_bin_nbases: 4,
        p_ge: ParametersGenome {
            transform: ParametersGenomeTransform {
                type_: 2,
                vcf_file: "in.vcf".to_string(),
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    let mut gtf = GTF {
        exon_loci: vec![[0, 1, 3, 0]],
        transcript_strand: vec![1],
        transcript_id: vec!["tr".to_string()],
        gene_id: vec!["gene".to_string()],
        gene_attr: vec![["name".to_string(), "type".to_string()]],
        ..Default::default()
    };

    let result = genome_transformgenome_l33_genome_transformgenome(
        &mut genome,
        &mut gtf,
        "chr1\t2\t.\tC\tT,G\t.\t.\t.\tGT\t1|2\n",
    );

    assert_eq!(&genome.g[0..5], &[0, 3, 2, 3, 0]);
    assert_eq!(&genome.g[8..13], &[0, 2, 2, 3, 0]);
    assert_eq!(genome.n_genome, 16);
    assert_eq!(genome.chr_name, vec!["chr1_h1", "chr1_h2"]);
    assert_eq!(genome.chr_name_index.get("chr1_h1"), Some(&0));
    assert_eq!(genome.chr_name_index.get("chr1_h2"), Some(&1));
    assert_eq!(genome.chr_start, vec![0, 8, 16]);
    assert_eq!(genome.chr_length, vec![5, 5]);
    assert_eq!(gtf.gene_id, vec!["gene_h1", "gene_h2"]);
    assert_eq!(gtf.transcript_id, vec!["tr_h1", "tr_h2"]);
    assert_eq!(gtf.transcript_strand, vec![1, 1]);
    assert_eq!(
        gtf.gene_attr,
        vec![
            ["name".to_string(), "type".to_string()],
            ["name".to_string(), "type".to_string()]
        ]
    );
    assert_eq!(gtf.exon_loci, vec![[0, 1, 3, 0], [1, 9, 11, 1]]);
    assert_eq!(result.transform_blocks_tsv, "2\t-1\n0\t5\t0\n8\t5\t0\n");
}

#[test]
fn genome_load_sjdb_handles_absent_and_present_junction_database() {
    let mut no_sj = Genome {
        n_genome: 100,
        n_chr_real: 1,
        chr_start: vec![0, 100],
        ..Default::default()
    };

    let log = genome_genomeload_l471_genome_loadsjdb(&mut no_sj, None).unwrap();
    assert_eq!(log, "");
    assert_eq!(no_sj.sjdb_n, 0);
    assert_eq!(no_sj.sj_gstart, 101);

    let mut genome = Genome {
        n_genome: 200,
        n_chr_real: 1,
        chr_start: vec![0, 100],
        p_ge: ParametersGenome {
            g_dir: "GenomeDir".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    let log = genome_genomeload_l471_genome_loadsjdb(
        &mut genome,
        Some("2 10\n50 80 1 2 3 1\n120 160 0 4 5 2\n"),
    )
    .unwrap();

    assert_eq!(
        log,
        "Processing splice junctions database sjdbN=2,   pGe.sjdbOverhang=10 \n"
    );
    assert_eq!(genome.sjdb_n, 2);
    assert_eq!(genome.p_ge.sjdb_overhang, 10);
    assert_eq!(genome.sj_gstart, 100);
    assert_eq!(genome.sjdb_start, vec![50, 120]);
    assert_eq!(genome.sjdb_end, vec![80, 160]);
    assert_eq!(genome.sjdb_motif, vec![1, 0]);
    assert_eq!(genome.sjdb_shift_left, vec![2, 4]);
    assert_eq!(genome.sjdb_shift_right, vec![3, 5]);
    assert_eq!(genome.sjdb_strand, vec![1, 2]);
    assert_eq!(genome.sj_dstart, vec![40, 114]);
    assert_eq!(genome.sj_astart, vec![81, 165]);

    let mut missing = Genome {
        n_genome: 200,
        n_chr_real: 1,
        chr_start: vec![0, 100],
        p_ge: ParametersGenome {
            g_dir: "GenomeDir".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    let err = genome_genomeload_l471_genome_loadsjdb(&mut missing, None).unwrap_err();
    assert!(err.contains("could not open file GenomeDir/sjdbInfo.txt"));
}

#[test]
fn genome_load_populates_core_arrays_and_redefines_windows() {
    let mut sa = PackedArray::default();
    packedarray_l8_packedarray_definebits(&mut sa, 33, 3);
    packedarray_l31_packedarray_allocatearray(&mut sa);
    packedarray_l17_packedarray_writepacked(&mut sa, 0, 4);
    packedarray_l17_packedarray_writepacked(&mut sa, 1, 2);
    packedarray_l17_packedarray_writepacked(&mut sa, 2, (1_u64 << 32) | 1);

    let mut sai = PackedArray::default();
    packedarray_l8_packedarray_definebits(&mut sai, 35, 2);
    packedarray_l31_packedarray_allocatearray(&mut sai);
    packedarray_l17_packedarray_writepacked(&mut sai, 0, 7);
    packedarray_l17_packedarray_writepacked(&mut sai, 1, 9);

    let mut sa_index = Vec::new();
    sa_index.extend_from_slice(&1_u32.to_ne_bytes());
    sa_index.extend_from_slice(&0_u32.to_ne_bytes());
    sa_index.extend_from_slice(&2_u32.to_ne_bytes());
    sa_index.extend_from_slice(&sai.char_array);

    let mut genome = Genome {
        p_ge: ParametersGenome {
            g_dir: "GenomeDir".to_string(),
            g_load: "NoSharedMemory".to_string(),
            sjdb_overhang: 10,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut p = Parameters {
        version_genome: "20201".to_string(),
        align_intron_max: 64,
        align_mates_gap_max: 0,
        win_bin_nbits: 1,
        win_anchor_dist_nbins: 3,
        ..Default::default()
    };

    let result = genome_genomeload_l18_genome_genomeload(
        &mut genome,
        &mut p,
        Some(
            "### STAR --runMode genomeGenerate\n\
             ### GstrandBit 32\n\
             versionGenome\t20201\n\
             genomeType\tFull\n\
             genomeSAindexNbases\t14\n\
             genomeChrBinNbits\t6\n\
             genomeSAsparseD\t2\n\
             sjdbOverhang\t10\n\
             genomeTransformType\tHaploid\n\
             genomeFileSizes\t8 16\n",
        ),
        "chr1\n",
        "4\n",
        "0\n8\n",
        &[0, 1, 2, 3, 4, 4, 4, 4],
        &sa.char_array,
        &sa_index,
        None,
        0,
    )
    .unwrap();

    assert!(result.log_stdout.contains("loading genome"));
    assert!(result.log_main.contains("Genome version is compatible"));
    assert_eq!(genome.n_genome, 8);
    assert_eq!(genome.g, vec![0, 1, 2, 3, 4, 4, 4, 4]);
    assert_eq!(genome.n_chr_real, 1);
    assert_eq!(genome.chr_name, vec!["chr1"]);
    assert_eq!(genome.chr_length, vec![4]);
    assert_eq!(genome.chr_start, vec![0, 8]);
    assert_eq!(genome.p_ge.g_saindex_nbases, 1);
    assert_eq!(genome.p_ge.g_chr_bin_nbits, 6);
    assert_eq!(genome.p_ge.g_sasparse_d, 2);
    assert_eq!(genome.gstrand_bit, 32);
    assert_eq!(genome.n_sa, 3);
    assert_eq!(genome.sa, vec![4, 2, (1_u64 << 32) as u32 | 1]);
    assert_eq!(genome.genome_sa_index_start, vec![0, 2]);
    assert_eq!(genome.sai, vec![7, 9]);
    assert_eq!(genome.sjdb_n, 0);
    assert_eq!(genome.sj_gstart, 9);
    assert_eq!(genome.chr_bin, vec![0]);
    assert_eq!(p.win_bin_nbits, 6);
    assert_eq!(p.win_flank_nbins, 2);
    assert_eq!(p.win_anchor_dist_nbins, 4);
    assert_eq!(p.win_bin_chr_nbits, 0);
    assert_eq!(p.win_bin_n, 1);
    assert_eq!(genome.p_ge.transform.type_, 1);
}

#[test]
fn genome_load_rejects_missing_and_incompatible_parameters() {
    let mut genome = Genome {
        p_ge: ParametersGenome {
            g_dir: "GenomeDir".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    let mut p = Parameters {
        version_genome: "current".to_string(),
        ..Default::default()
    };
    let err = genome_genomeload_l18_genome_genomeload(
        &mut genome,
        &mut p,
        None,
        "",
        "",
        "",
        &[],
        &[],
        &[],
        None,
        0,
    )
    .unwrap_err();
    assert!(err.contains("could not open genome file GenomeDir/genomeParameters.txt"));

    let err = genome_genomeload_l18_genome_genomeload(
        &mut genome,
        &mut p,
        Some("versionGenome\told\n"),
        "",
        "",
        "",
        &[],
        &[],
        &[],
        None,
        0,
    )
    .unwrap_err();
    assert!(err.contains("INCOMPATIBLE"));
}

#[test]
fn genome_out_load_populates_output_genome_and_conversion_blocks() {
    let mut genome = Genome {
        p_ge: ParametersGenome {
            g_dir: "OutGenome".to_string(),
            chr_set_mito_strings: vec!["chrM".to_string()],
            ..Default::default()
        },
        genome_out: star_rs::generated::structs::GenomeOut {
            conv_file: "OutGenome/transformGenomeBlocks.tsv".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };

    let log = genome_genomeoutload_l8_genome_genomeoutload(
        &mut genome,
        Some("genomeSAindexNbases 7\ngenomeChrBinNbits 3\ngenomeSAsparseD 4\n"),
        "chr1\nchrM\n",
        "5\n3\n",
        "0\n8\n16\n",
        &[0, 1, 2, 3, 0, 4, 4, 4, 3, 2, 1, 4, 4, 4, 4, 4, 0, 1, 2, 4],
        Some("1 10\n12 15 0 2 3 1\n"),
        "2 -1\n0 5 100\n10 6 200\n",
    )
    .unwrap();

    assert!(log.starts_with("Reading output genome generation parameters:\n"));
    assert!(log.contains("Number of real (reference) chromosomes= 2\n"));
    assert!(log.contains("Processing splice junctions database sjdbN=1"));
    assert_eq!(genome.p_ge.g_saindex_nbases, 7);
    assert_eq!(genome.p_ge.g_chr_bin_nbits, 3);
    assert_eq!(genome.p_ge.g_sasparse_d, 4);
    assert_eq!(genome.genome_chr_bin_nbases, 8);
    assert_eq!(genome.n_genome, 20);
    assert_eq!(genome.g[0..5], [0, 1, 2, 3, 0]);
    assert_eq!(genome.chr_name, vec!["chr1", "chrM"]);
    assert_eq!(genome.chr_length, vec![5, 3]);
    assert_eq!(genome.chr_start, vec![0, 8, 16]);
    assert_eq!(genome.chr_name_index.get("chrM"), Some(&1));
    assert!(genome.p_ge.chr_set_mito.contains(&1));
    assert_eq!(genome.sjdb_n, 1);
    assert_eq!(genome.sj_gstart, 16);
    assert_eq!(genome.sj_dstart, vec![4]);
    assert_eq!(genome.sj_astart, vec![18]);
    assert_eq!(genome.chr_bin, vec![0, 1, 2]);
    assert_eq!(genome.genome_out.n_minus_strand_offset, u64::MAX);
    assert_eq!(
        genome.genome_out.conv_blocks,
        vec![[0, 5, 100], [10, 7, 200], [u64::MAX, 0, 0]]
    );

    let err =
        genome_genomeoutload_l8_genome_genomeoutload(&mut genome, None, "", "", "", &[], None, "")
            .unwrap_err();
    assert!(err.contains("could not open genome file OutGenome/genomeParameters.txt"));
}

#[test]
fn genome_g2str_locus_matches_strand_bit_logic() {
    assert_eq!(genome_genomegenerate_l91_fung2strlocus(7, 100, 4, 0x0f), 7);
    assert_eq!(
        genome_genomegenerate_l91_fung2strlocus(0x13, 100, 4, 0x0f),
        103
    );
}

#[test]
fn stream_read_big_reads_at_most_requested_bytes() {
    let mut cursor = Cursor::new(b"abcdef".to_vec());
    let mut out = [0u8; 8];
    assert_eq!(streamfuns_l39_fstreamreadbig(&mut cursor, &mut out, 4), 4);
    assert_eq!(&out[..4], b"abcd");
    assert_eq!(
        streamfuns_l39_fstreamreadbig(&mut cursor, &mut out[4..], 4),
        2
    );
    assert_eq!(&out[..6], b"abcdef");
}

#[test]
fn stream_write_big_writes_requested_prefix() {
    let mut out = Vec::new();
    streamfuns_l51_fstreamwritebig(&mut out, b"abcdef", 4).unwrap();
    assert_eq!(out, b"abcd");

    streamfuns_l51_fstreamwritebig(&mut out, b"xyz", 0).unwrap();
    assert_eq!(out, b"abcd");
}

#[test]
fn copy_file_copies_bytes() {
    let dir = std::env::temp_dir().join(format!("star-rs-copy-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let src = dir.join("in.bin");
    let dst = dir.join("out.bin");
    std::fs::write(&src, b"a\0bc").unwrap();

    assert_eq!(
        streamfuns_l144_copyfile(src.to_str().unwrap(), dst.to_str().unwrap()).unwrap(),
        4
    );
    assert_eq!(std::fs::read(&dst).unwrap(), b"a\0bc");

    let _ = std::fs::remove_file(src);
    let _ = std::fs::remove_file(dst);
    let _ = std::fs::remove_dir(dir);
}

#[test]
fn read_align_chunk_line_helpers_trim_one_trailing_control_byte() {
    let mut cursor = Cursor::new(b"ACGT\r\nNNNN \nplain\n".to_vec());
    let mut out = Vec::new();

    assert_eq!(
        readalignchunk_processchunks_l284_fastqreadoneline(&mut cursor, &mut out).unwrap(),
        5
    );
    assert_eq!(out, b"ACGT\n");

    assert_eq!(
        readalignchunk_processchunks_l284_fastqreadoneline(&mut cursor, &mut out).unwrap(),
        5
    );
    assert_eq!(out, b"NNNN\n");

    assert_eq!(
        readalignchunk_processchunks_l284_fastqreadoneline(&mut cursor, &mut out).unwrap(),
        6
    );
    assert_eq!(out, b"plain\n");

    let mut text = "read\t".to_string();
    readalignchunk_processchunks_l298_removestringendcontrol(&mut text);
    assert_eq!(text, "read");

    let mut text_without_control = "read".to_string();
    readalignchunk_processchunks_l298_removestringendcontrol(&mut text_without_control);
    assert_eq!(text_without_control, "read");
}

#[test]
fn read_align_chunk_file_helpers_open_cat_and_remove_thread_files() {
    let dir = std::env::temp_dir().join(format!("star-rs-read-align-chunk-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let prefix = dir.join("chunk.thread");
    let prefix = prefix.to_str().unwrap();

    std::fs::write(format!("{prefix}3"), b"old").unwrap();
    let (mut chunk_file, file_name, log) =
        readalignchunk_l116_readalignchunk_chunkfstreamopen(prefix, 3).unwrap();
    assert_eq!(file_name, format!("{prefix}3"));
    assert_eq!(
        log,
        format!("Opening the file: {} ... ok\n", format!("{prefix}3"))
    );
    assert_eq!(std::fs::read(&file_name).unwrap(), b"");

    chunk_file.write_all(b"abc").unwrap();
    chunk_file.seek(std::io::SeekFrom::Start(0)).unwrap();
    chunk_file.write_all(b"XY").unwrap();
    let mut all = Vec::new();
    readalignchunk_l137_readalignchunk_chunkfstreamcat(&mut chunk_file, &mut all, true).unwrap();
    assert_eq!(all, b"XYc");
    assert_eq!(chunk_file.stream_position().unwrap(), 0);

    std::fs::write(format!("{prefix}4"), b"one").unwrap();
    std::fs::write(format!("{prefix}5"), b"two").unwrap();
    let mut i_c = 4;
    readalignchunk_l151_readalignchunk_chunkfilescat(&mut all, prefix, &mut i_c).unwrap();
    assert_eq!(all, b"XYconetwo");
    assert_eq!(i_c, 6);
    assert!(!std::path::Path::new(&format!("{prefix}4")).exists());
    assert!(!std::path::Path::new(&format!("{prefix}5")).exists());

    drop(chunk_file);
    let _ = std::fs::remove_file(file_name);
    let _ = std::fs::remove_dir(dir);
}

#[test]
fn read_load_parses_fastq_applies_quality_conversion_and_name_trimming() {
    let mut cursor = Cursor::new(b"@read/1 42 Y 3 CB:Z:abc\nACGTN\n+clip\n!!!!!\n".to_vec());
    let p = Parameters {
        read_name_separator_char: vec!['/'],
        out_qs_conversion_add: 100,
        ..Default::default()
    };
    let mut l_read = 0;
    let mut l_read_original = 0;
    let mut read_name = String::new();
    let mut seq = String::new();
    let mut seq_num = Vec::new();
    let mut qual = String::new();
    let mut clip = vec![
        ClipMate {
            type_: -1,
            ..Default::default()
        },
        ClipMate {
            type_: -1,
            ..Default::default()
        },
    ];
    let mut i_read_all = 0;
    let mut read_files_index = 0;
    let mut read_filter = 0;
    let mut read_name_extra = String::new();

    assert_eq!(
        readload_l4_readload(
            &mut cursor,
            &p,
            &mut l_read,
            &mut l_read_original,
            &mut read_name,
            &mut seq,
            &mut seq_num,
            &mut qual,
            &mut clip,
            &mut i_read_all,
            &mut read_files_index,
            &mut read_filter,
            &mut read_name_extra,
        )
        .unwrap(),
        2
    );

    assert_eq!(read_name, "@read");
    assert_eq!(i_read_all, 42);
    assert_eq!(read_filter, b'Y');
    assert_eq!(read_files_index, 3);
    assert_eq!(read_name_extra, "CB:Z:abc");
    assert_eq!(seq, "ACGTN");
    assert_eq!(seq_num, vec![0, 1, 2, 3, 4]);
    assert_eq!(l_read, 5);
    assert_eq!(l_read_original, 5);
    assert_eq!(clip[0].clipped_info, b'+');
    assert_eq!(qual, "~~~~~");

    assert_eq!(
        readload_l4_readload(
            &mut cursor,
            &p,
            &mut l_read,
            &mut l_read_original,
            &mut read_name,
            &mut seq,
            &mut seq_num,
            &mut qual,
            &mut clip,
            &mut i_read_all,
            &mut read_files_index,
            &mut read_filter,
            &mut read_name_extra,
        )
        .unwrap(),
        -1
    );
}

#[test]
fn read_load_parses_fasta_and_reports_quality_length_mismatch() {
    let p = Parameters {
        read_name_separator_char: vec![' '],
        ..Default::default()
    };
    let mut cursor = Cursor::new(b">read description\nTGCA\n".to_vec());
    let mut l_read = 0;
    let mut l_read_original = 0;
    let mut read_name = String::new();
    let mut seq = String::new();
    let mut seq_num = Vec::new();
    let mut qual = String::new();
    let mut clip = vec![
        ClipMate {
            type_: -1,
            ..Default::default()
        },
        ClipMate {
            type_: -1,
            ..Default::default()
        },
    ];
    let mut i_read_all = 9;
    let mut read_files_index = 8;
    let mut read_filter = b'N';
    let mut read_name_extra = String::new();

    assert_eq!(
        readload_l4_readload(
            &mut cursor,
            &p,
            &mut l_read,
            &mut l_read_original,
            &mut read_name,
            &mut seq,
            &mut seq_num,
            &mut qual,
            &mut clip,
            &mut i_read_all,
            &mut read_files_index,
            &mut read_filter,
            &mut read_name_extra,
        )
        .unwrap(),
        1
    );
    assert_eq!(read_name, ">read");
    assert_eq!(seq_num, vec![3, 2, 1, 0]);
    assert_eq!(qual, "AAAA");

    let mut bad = Cursor::new(b"@bad\nACGT\n+\nIII\n".to_vec());
    let err = readload_l4_readload(
        &mut bad,
        &p,
        &mut l_read,
        &mut l_read_original,
        &mut read_name,
        &mut seq,
        &mut seq_num,
        &mut qual,
        &mut clip,
        &mut i_read_all,
        &mut read_files_index,
        &mut read_filter,
        &mut read_name_extra,
    )
    .unwrap_err();
    assert!(err.contains("quality string length is not equal to sequence length"));
}
