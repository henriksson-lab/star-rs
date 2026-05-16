#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `appendVector` at STAR/source/Genome_transformGenome.cpp:10. Args: v1: vector<T>, v2: vector<T>"]
pub fn genome_transformgenome_l10_appendvector<T: Clone>(v1: &mut Vec<T>, v2: &[T]) {
    v1.reserve(v2.len());
    v1.extend_from_slice(v2);
}

#[doc = "Original `concatenateVectors` at STAR/source/Genome_transformGenome.cpp:17. Args: v1: vector<T>, v2: vector<T>"]
pub fn genome_transformgenome_l17_concatenatevectors<T: Clone>(v1: &[T], v2: &[T]) -> Vec<T> {
    let mut v_out = Vec::with_capacity(v1.len() + v2.len());
    v_out.extend_from_slice(v1);
    v_out.extend_from_slice(v2);
    v_out
}

#[doc = "Original `appendString` at STAR/source/Genome_transformGenome.cpp:26. Args: vString: vector<string>, strAdd: string"]
pub fn genome_transformgenome_l26_appendstring(
    mut v_string: Vec<String>,
    str_add: &str,
) -> Vec<String> {
    for s in &mut v_string {
        s.push_str(str_add);
    }
    v_string
}

#[doc = "Original `Genome::transformGenome` at STAR/source/Genome_transformGenome.cpp:33. Args: gtf: GTF"]
pub fn genome_transformgenome_l33_genome_transformgenome(
    genome: &mut crate::genome::Genome,
    gtf: &mut crate::gtf::GTF,
    vcf_contents: &str,
) -> crate::parameters_chimeric::GenomeTransformGenomeResult {
    let mut result = crate::parameters_chimeric::GenomeTransformGenomeResult::default();
    if genome.p_ge.transform.type_ == 0 {
        return result;
    }

    result
        .log_main
        .push_str("transformGenome: processing VCF\n");

    let mut vcf_variants =
        vec![
            std::collections::BTreeMap::<String, Vec<crate::variation::VariantInfo>>::new(
            );
            genome.p_ge.transform.type_ as usize
        ];

    for vcf_line in vcf_contents.lines() {
        let mut fields = vcf_line.split_whitespace();
        let Some(chr) = fields.next() else {
            continue;
        };
        if chr.starts_with('#') {
            continue;
        }
        if !genome.chr_name_index.contains_key(chr) {
            result.log_main.push_str(&format!(
                "WARNING: while processing varVCFfile file={}: chromosome '{}' not found in Genome fasta file\n",
                genome.p_ge.transform.vcf_file, chr
            ));
            continue;
        }

        let pos = fields
            .next()
            .and_then(|pos| pos.parse::<u64>().ok())
            .unwrap_or(0);
        let _id = fields.next();
        let ref_seq = fields.next().unwrap_or("").to_string();
        let alt = fields.next().unwrap_or("");
        let _qual = fields.next();
        let _filter = fields.next();
        let _info = fields.next();
        let _format = fields.next();
        let sample = fields.next().unwrap_or("");
        let alt_v: Vec<&str> = alt.split(',').collect();

        if genome.p_ge.transform.type_ == 1 {
            let alt_seq = alt_v.first().copied().unwrap_or("").to_string();
            vcf_variants[0].entry(chr.to_string()).or_default().push(
                crate::variation::VariantInfo {
                    pos,
                    len: alt_seq.len() as i32 - ref_seq.len() as i32,
                    seq: [ref_seq, alt_seq],
                },
            );
        } else if genome.p_ge.transform.type_ == 2 {
            for ih in 0..2usize {
                let gt = sample
                    .as_bytes()
                    .get(ih * 2)
                    .and_then(|b| char::from(*b).to_digit(10))
                    .unwrap_or(0) as usize;
                if gt == 0 {
                    continue;
                }
                let alt_seq = alt_v.get(gt - 1).copied().unwrap_or("").to_string();
                vcf_variants[ih].entry(chr.to_string()).or_default().push(
                    crate::variation::VariantInfo {
                        pos,
                        len: alt_seq.len() as i32 - ref_seq.len() as i32,
                        seq: [ref_seq.clone(), alt_seq],
                    },
                );
            }
        }
    }

    let mut n_genome1 = 0_u64;
    let mut g_new = Vec::<u8>::new();

    if genome.p_ge.transform.type_ == 1 {
        let mut chr_start1 = Vec::new();
        let mut chr_length1 = Vec::new();
        result
            .log_main
            .push_str(&genome_transformgenome_l171_genome_transformchrlenstart(
                genome,
                &mut vcf_variants[0],
                &mut chr_start1,
                &mut chr_length1,
            ));

        n_genome1 = *chr_start1.last().unwrap_or(&0);
        result.log_main.push_str(&format!(
            "Old/new genome sizes: {} {}\n",
            genome.n_genome, n_genome1
        ));
        g_new = vec![4; n_genome1 as usize];

        let mut transform_blocks = Vec::new();
        if let Some(debug) = genome_transformgenome_l215_genome_transformgandblocks(
            genome,
            &vcf_variants[0],
            &chr_start1,
            &chr_length1,
            &mut transform_blocks,
            &mut g_new,
        ) {
            result.debug.push_str(&debug);
        }
        let mut exon_loci = gtf
            .exon_loci
            .iter()
            .map(|exon| [exon[0], exon[1], exon[2], exon[3], 0])
            .collect::<Vec<_>>();
        result
            .log_main
            .push_str(&genome_transformgenome_l282_genome_transformexonloci(
                &mut exon_loci,
                &transform_blocks,
            ));
        gtf.exon_loci = exon_loci
            .into_iter()
            .map(|exon| [exon[0], exon[1], exon[2], exon[3]])
            .collect();

        genome.chr_start = chr_start1;
        genome.chr_length = chr_length1;
        result.transform_blocks_tsv =
            genome_transformgenome_l271_genome_transformblockswrite(&transform_blocks);
    } else if genome.p_ge.transform.type_ == 2 {
        let mut chr_start1 = [Vec::new(), Vec::new()];
        let mut chr_length1 = [Vec::new(), Vec::new()];
        for ih in 0..2usize {
            result
                .log_main
                .push_str(&genome_transformgenome_l171_genome_transformchrlenstart(
                    genome,
                    &mut vcf_variants[ih],
                    &mut chr_start1[ih],
                    &mut chr_length1[ih],
                ));
        }

        let hap0_end = *chr_start1[0].last().unwrap_or(&0);
        for cs in &mut chr_start1[1] {
            *cs += hap0_end;
        }
        n_genome1 = *chr_start1[1].last().unwrap_or(&0);
        result.log_main.push_str(&format!(
            "Old/new genome sizes: {} {}\n",
            genome.n_genome, n_genome1
        ));
        g_new = vec![4; n_genome1 as usize];

        let mut transform_blocks = [Vec::new(), Vec::new()];
        for ih in 0..2usize {
            if let Some(debug) = genome_transformgenome_l215_genome_transformgandblocks(
                genome,
                &vcf_variants[ih],
                &chr_start1[ih],
                &chr_length1[ih],
                &mut transform_blocks[ih],
                &mut g_new,
            ) {
                result.debug.push_str(&debug);
            }
        }

        let mut exon_loci1 = [
            gtf.exon_loci
                .iter()
                .map(|exon| [exon[0], exon[1], exon[2], exon[3], 0])
                .collect::<Vec<_>>(),
            gtf.exon_loci
                .iter()
                .map(|exon| [exon[0], exon[1], exon[2], exon[3], 0])
                .collect::<Vec<_>>(),
        ];
        for exon in &mut exon_loci1[1] {
            exon[GTF_EX_T] += gtf.transcript_id.len() as u64;
            exon[GTF_EX_G] += gtf.gene_id.len() as u64;
        }

        let gene_attr_copy = gtf.gene_attr.clone();
        genome_transformgenome_l10_appendvector(&mut gtf.gene_attr, &gene_attr_copy);
        let transcript_strand_copy = gtf.transcript_strand.clone();
        genome_transformgenome_l10_appendvector(
            &mut gtf.transcript_strand,
            &transcript_strand_copy,
        );
        gtf.gene_id = genome_transformgenome_l17_concatenatevectors(
            &genome_transformgenome_l26_appendstring(gtf.gene_id.clone(), "_h1"),
            &genome_transformgenome_l26_appendstring(gtf.gene_id.clone(), "_h2"),
        );
        gtf.transcript_id = genome_transformgenome_l17_concatenatevectors(
            &genome_transformgenome_l26_appendstring(gtf.transcript_id.clone(), "_h1"),
            &genome_transformgenome_l26_appendstring(gtf.transcript_id.clone(), "_h2"),
        );

        result
            .log_main
            .push_str(&genome_transformgenome_l282_genome_transformexonloci(
                &mut exon_loci1[0],
                &transform_blocks[0],
            ));
        result
            .log_main
            .push_str(&genome_transformgenome_l282_genome_transformexonloci(
                &mut exon_loci1[1],
                &transform_blocks[1],
            ));

        genome.chr_name = genome_transformgenome_l17_concatenatevectors(
            &genome_transformgenome_l26_appendstring(genome.chr_name.clone(), "_h1"),
            &genome_transformgenome_l26_appendstring(genome.chr_name.clone(), "_h2"),
        );
        genome.n_chr_real = genome.chr_name.len() as u32;
        genome.chr_name_index.clear();
        for (ii, name) in genome.chr_name.iter().enumerate() {
            genome.chr_name_index.insert(name.clone(), ii as u64);
        }

        chr_start1[0].pop();
        genome.chr_start =
            genome_transformgenome_l17_concatenatevectors(&chr_start1[0], &chr_start1[1]);
        genome.chr_length =
            genome_transformgenome_l17_concatenatevectors(&chr_length1[0], &chr_length1[1]);
        gtf.exon_loci =
            genome_transformgenome_l17_concatenatevectors(&exon_loci1[0], &exon_loci1[1])
                .into_iter()
                .map(|exon| [exon[0], exon[1], exon[2], exon[3]])
                .collect();
        let transform_blocks1 = transform_blocks[1].clone();
        genome_transformgenome_l10_appendvector(&mut transform_blocks[0], &transform_blocks1);
        result.transform_blocks_tsv =
            genome_transformgenome_l271_genome_transformblockswrite(&transform_blocks[0]);
    }

    genome.g = g_new;
    genome.n_genome = n_genome1;
    result
}

#[doc = "Original `Genome::transformChrLenStart` at STAR/source/Genome_transformGenome.cpp:171. Args: vcfVariants: map<string,vector<VariantInfo>>, chrStart1: vector<uint64>, chrLength1: vector<uint64>"]
pub fn genome_transformgenome_l171_genome_transformchrlenstart(
    genome: &crate::genome::Genome,
    vcf_variants: &mut std::collections::BTreeMap<
        String,
        Vec<crate::variation::VariantInfo>,
    >,
    chr_start1: &mut Vec<u64>,
    chr_length1: &mut Vec<u64>,
) -> String {
    let mut log_main = String::new();
    *chr_start1 = genome.chr_start.clone();
    *chr_length1 = genome.chr_length.clone();

    for ichr in 0..genome.chr_length.len() {
        if !vcf_variants.contains_key(&genome.chr_name[ichr]) {
            continue;
        }

        let vv = vcf_variants.get_mut(&genome.chr_name[ichr]).unwrap();
        vv.sort_by_key(|vi| vi.pos);

        let mut vv1 = Vec::with_capacity(vv.len());
        let mut g0 = 0u64;
        for v in vv.iter() {
            if v.pos >= g0 {
                vv1.push(v.clone());
            }
            g0 = std::cmp::max(g0, v.pos + v.seq[0].len() as u64);
        }
        log_main.push_str(&format!(
            "{}: filtered out overlapping variants = {}; remaining variants = {}\n",
            genome.chr_name[ichr],
            vv.len() as i64 - vv1.len() as i64,
            vv1.len()
        ));
        *vv = vv1;

        for v in vv.iter() {
            chr_length1[ichr] =
                ((chr_length1[ichr] as i64) + v.seq[1].len() as i64 - v.seq[0].len() as i64) as u64;
        }
        log_main.push_str(&format!(
            "Transformed chr length difference: {} {}\n",
            genome.chr_name[ichr],
            chr_length1[ichr] as i64 - genome.chr_length[ichr] as i64
        ));
    }

    chr_start1[0] = 0;
    for ichr in 0..genome.chr_length.len() {
        chr_start1[ichr + 1] = chr_start1[ichr]
            + ((chr_length1[ichr] + 1) / genome.genome_chr_bin_nbases as u64 + 1)
                * genome.genome_chr_bin_nbases as u64;
        log_main.push_str(&format!(
            "Transformed chr start difference: {} {}\n",
            genome.chr_name[ichr],
            chr_start1[ichr] as i64 - genome.chr_start[ichr] as i64
        ));
    }

    log_main
}

#[doc = "Original `Genome::transformGandBlocks` at STAR/source/Genome_transformGenome.cpp:215. Args: vcfVariants: map<string,vector<VariantInfo>>, chrStart1: vector<uint64>, chrLength1: vector<uint64>, transformBlocks: vector<array<uint64,3>>, Gnew: char"]
pub fn genome_transformgenome_l215_genome_transformgandblocks(
    genome: &crate::genome::Genome,
    vcf_variants: &std::collections::BTreeMap<String, Vec<crate::variation::VariantInfo>>,
    chr_start1: &[u64],
    chr_length1: &[u64],
    transform_blocks: &mut Vec<[u64; 3]>,
    g_new: &mut [u8],
) -> Option<String> {
    let mut debug = String::new();
    for ichr in 0..genome.chr_length.len() {
        if !vcf_variants.contains_key(&genome.chr_name[ichr]) {
            let dst = chr_start1[ichr] as usize;
            let src = genome.chr_start[ichr] as usize;
            let len = genome.chr_length[ichr] as usize;
            g_new[dst..dst + len].copy_from_slice(&genome.g[src..src + len]);
            transform_blocks.push([
                genome.chr_start[ichr],
                genome.chr_length[ichr],
                chr_start1[ichr],
            ]);
            continue;
        }

        let vv = &vcf_variants[&genome.chr_name[ichr]];
        let mut iv = 0usize;
        let mut g1 = chr_start1[ichr];
        let mut g0 = genome.chr_start[ichr];
        transform_blocks.push([g0, 0, g1]);

        while g0 < genome.chr_start[ichr] + genome.chr_length[ichr] {
            if g0 == vv[iv].pos - 1 + genome.chr_start[ichr] {
                let seq = &vv[iv].seq;
                let mut s0 = vec![0u8; seq[0].len()];
                sequencefuns_l131_convertnucleotidestonumbers(
                    seq[0].as_bytes(),
                    &mut s0,
                    seq[0].len() as u64,
                );
                if genome.g[g0 as usize..g0 as usize + seq[0].len()] != s0[..] {
                    debug.push_str(&format!("{} {}\n", g0, seq[0]));
                }

                let mut s1 = vec![0u8; seq[1].len()];
                sequencefuns_l131_convertnucleotidestonumbers(
                    seq[1].as_bytes(),
                    &mut s1,
                    seq[1].len() as u64,
                );
                g_new[g1 as usize..g1 as usize + s1.len()].copy_from_slice(&s1);
                g0 += seq[0].len() as u64;
                g1 += seq[1].len() as u64;

                if vv[iv].len != 0 {
                    let last = transform_blocks.last_mut().unwrap();
                    last[1] = g0 - seq[0].len() as u64
                        + std::cmp::min(seq[0].len(), seq[1].len()) as u64
                        - last[0];
                    transform_blocks.push([g0, 0, g1]);
                }

                if iv < vv.len() - 1 {
                    iv += 1;
                }
            } else {
                g_new[g1 as usize] = genome.g[g0 as usize];
                g0 += 1;
                g1 += 1;
            }
        }

        if transform_blocks.last().unwrap()[1] == 0 {
            let last = transform_blocks.last_mut().unwrap();
            last[1] = g0 - last[0];
        }

        if g1 != chr_start1[ichr] + chr_length1[ichr] {
            debug.push_str(&format!(
                "{} {}\n",
                g1,
                chr_start1[ichr] + chr_length1[ichr]
            ));
        }
    }

    if debug.is_empty() { None } else { Some(debug) }
}

#[doc = "Original `Genome::transformBlocksWrite` at STAR/source/Genome_transformGenome.cpp:271. Args: transformBlocks: vector<array<uint64,3>>"]
pub fn genome_transformgenome_l271_genome_transformblockswrite(
    transform_blocks: &[[u64; 3]],
) -> String {
    let mut conv_stream = format!("{}\t-1\n", transform_blocks.len());
    for tb in transform_blocks {
        conv_stream.push_str(&format!("{}\t{}\t{}\n", tb[2], tb[1], tb[0]));
    }
    conv_stream
}

#[doc = "Original `Genome::transformExonLoci` at STAR/source/Genome_transformGenome.cpp:282. Args: exonLoci: vector<array<uint64,exL>>, transformBlocks: vector<array<uint64,3>>"]
pub fn genome_transformgenome_l282_genome_transformexonloci(
    exon_loci: &mut Vec<[u64; 5]>,
    transform_blocks: &[[u64; 3]],
) -> String {
    let mut exon_loci1 = Vec::new();
    for exon in exon_loci.iter_mut() {
        let exon_s = exon[GENOME_EX_S];
        let exon_e = exon[GENOME_EX_E];

        let mut t_bit = transform_blocks.partition_point(|tb| tb[0] < exon_s);
        t_bit -= 1;
        let mut tb = transform_blocks[t_bit];

        if exon_s < tb[0] + tb[1] {
            exon[GENOME_EX_S] = tb[2] + exon_s - tb[0];
        } else {
            exon[GENOME_EX_S] = transform_blocks[t_bit + 1][2];
        }

        while exon_e > transform_blocks[t_bit][0] + transform_blocks[t_bit][1] {
            t_bit += 1;
        }

        tb = transform_blocks[t_bit];
        if exon_e >= tb[0] {
            exon[GENOME_EX_E] = tb[2] + exon_e - tb[0];
        } else {
            exon[GENOME_EX_E] = transform_blocks[t_bit - 1][2] + transform_blocks[t_bit - 1][1] - 1;
        }

        if exon[GENOME_EX_S] <= exon[GENOME_EX_E] {
            exon_loci1.push(*exon);
        }
    }

    let log_main = format!(
        "Transform exons: removed {}\n",
        exon_loci.len() - exon_loci1.len()
    );
    *exon_loci = exon_loci1;
    log_main
}
