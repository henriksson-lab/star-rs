#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `insertSeqSA` at STAR/source/insertSeqSA.cpp:18. Args: SA: PackedArray, SA1: PackedArray, SAi: PackedArray, G: char, G1: char, nG: uint64, nG1: uint64, nG2: uint64, P: Parameters, mapGen: Genome"]
pub fn insertseqsa_l18_insertseqsa(
    sa: &mut crate::packed_array::PackedArray,
    sa1: &mut crate::packed_array::PackedArray,
    _sai: &mut crate::packed_array::PackedArray,
    _g: &[u8],
    g1: &[u8],
    n_g: u64,
    n_g1: u64,
    n_g2: u64,
    _p: &crate::parameters_chimeric::Parameters,
    map_gen: &crate::genome::Genome,
) -> Result<u64, String> {
    let mut gstrand_bit1 = ((n_g + n_g1) as f64).ln() / 2.0f64.ln();
    gstrand_bit1 = gstrand_bit1.floor() + 1.0;
    let mut gstrand_bit1 = gstrand_bit1 as u64;
    if gstrand_bit1 < 32 {
        gstrand_bit1 = 32;
    }
    if gstrand_bit1 + 1 > sa.word_length {
        return Err("EXITING because of FATAL ERROR: cannot insert sequence on the fly because of strand GstrandBit problem\nSOLUTION: please contact STAR author at https://groups.google.com/forum/#!forum/rna-star\n".to_string());
    }

    let n2bit = 1u64 << (sa.word_length - 1);
    let strand_mask = !n2bit;
    for isa in 0..sa.length {
        let mut ind1 = packedarray_h18_packedarray_index(sa, isa);
        if (ind1 & n2bit) > 0 {
            if (ind1 & strand_mask) >= n_g2 {
                ind1 += n_g1;
                packedarray_l17_packedarray_writepacked(sa, isa, ind1);
            }
        } else if ind1 >= n_g {
            ind1 += n_g1;
            packedarray_l17_packedarray_writepacked(sa, isa, ind1);
        }
    }

    const GENOME_END_FILL_L: usize = 16;
    let n_g1_usize = n_g1 as usize;
    let mut seqq = vec![GENOME_SPACING_CHAR; 4 * n_g1_usize + 3 * GENOME_END_FILL_L];
    let seq0_start = GENOME_END_FILL_L;
    let seq1_start = 2 * GENOME_END_FILL_L + 2 * n_g1_usize;

    seqq[seq0_start..seq0_start + n_g1_usize].copy_from_slice(&g1[..n_g1_usize]);
    for ii in 0..n_g1_usize {
        let value = seqq[seq0_start + ii];
        seqq[seq0_start + 2 * n_g1_usize - 1 - ii] = if value < 4 { 3 - value } else { value };
    }
    let seq0 = seqq[seq0_start..seq0_start + 2 * n_g1_usize].to_vec();
    sequencefuns_l4_complementseqnumbers(
        &seq0,
        &mut seqq[seq1_start..seq1_start + 2 * n_g1_usize],
        (2 * n_g1) as u64,
    );

    let seq0 = &seqq[seq0_start..];
    let seq1 = &seqq[seq1_start..];
    let mut ind_array = vec![0u64; n_g1_usize * 4 + 2];
    for ii in 0..(2 * n_g1_usize) {
        if seq0[ii] > 3 {
            ind_array[ii * 2] = u64::MAX;
        } else {
            let mut l = 0u64;
            ind_array[ii * 2] = suffixarrayfuns_l297_suffixarraysearch1(
                map_gen,
                [seq0, seq1],
                ii as u64,
                10000,
                n_g,
                ii < n_g1_usize,
                0,
                sa.length - 1,
                &mut l,
            );
            ind_array[ii * 2 + 1] = ii as u64;
        }
    }

    let mut n_ind = 0usize;
    for ii in 0..(2 * n_g1_usize) {
        if ind_array[ii * 2] != u64::MAX {
            ind_array[n_ind * 2] = ind_array[ii * 2];
            ind_array[n_ind * 2 + 1] = ind_array[ii * 2 + 1];
            n_ind += 1;
        }
    }

    let sort_len = map_gen.p_ge.g_suffix_length_max as usize / std::mem::size_of::<u64>();
    ind_array[..2 * n_ind].chunks_exact_mut(2).for_each(|_| {});
    let mut pairs: Vec<[u64; 2]> = ind_array[..2 * n_ind]
        .chunks_exact(2)
        .map(|chunk| [chunk[0], chunk[1]])
        .collect();
    pairs.sort_by(|a, b| {
        let cmp = funcompareuintandsuffixesmemcmp_l7_funcompareuintandsuffixesmemcmp(
            a, b, seq0, sort_len,
        );
        cmp.cmp(&0)
    });
    for (ii, pair) in pairs.iter().enumerate() {
        ind_array[2 * ii] = pair[0];
        ind_array[2 * ii + 1] = pair[1];
    }
    ind_array[2 * n_ind] = u64::MAX - 998;
    ind_array[2 * n_ind + 1] = u64::MAX - 998;

    packedarray_l8_packedarray_definebits(sa1, sa.word_length, sa.length + n_ind as u64);
    packedarray_l31_packedarray_allocatearray(sa1);

    let mut isa1 = 0usize;
    let mut isa2 = 0u64;
    for isa in 0..sa.length {
        while isa == ind_array[isa1 * 2] {
            let mut ind1 = ind_array[isa1 * 2 + 1];
            if ind1 < n_g1 {
                ind1 += n_g;
            } else {
                ind1 = (ind1 - n_g1 + n_g2) | n2bit;
            }
            packedarray_l17_packedarray_writepacked(sa1, isa2, ind1);
            isa2 += 1;
            isa1 += 1;
        }
        let ind1 = packedarray_h18_packedarray_index(sa, isa);
        packedarray_l17_packedarray_writepacked(sa1, isa2, ind1);
        isa2 += 1;
    }
    while isa1 < n_ind {
        let mut ind1 = ind_array[isa1 * 2 + 1];
        if ind1 < n_g1 {
            ind1 += n_g;
        } else {
            ind1 = (ind1 - n_g1 + n_g2) | n2bit;
        }
        packedarray_l17_packedarray_writepacked(sa1, isa2, ind1);
        isa2 += 1;
        isa1 += 1;
    }

    Ok(n_ind as u64)
}
