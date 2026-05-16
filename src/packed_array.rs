#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

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

#[doc = "Original `PackedArray::PackedArray` at STAR/source/PackedArray.cpp:3. Args: "]
pub fn packedarray_l3_packedarray_packedarray() -> crate::packed_array::PackedArray {
    crate::packed_array::PackedArray {
        array_allocated: false,
        ..Default::default()
    }
}

#[doc = "Original `PackedArray::defineBits` at STAR/source/PackedArray.cpp:8. Args: Nbits: uint, lengthIn: uint"]
pub fn packedarray_l8_packedarray_definebits(
    packed: &mut crate::packed_array::PackedArray,
    nbits: u64,
    length_in: u64,
) {
    packed.word_length = nbits;
    packed.word_comp_length = (std::mem::size_of::<u64>() as u64) * 8 - packed.word_length;
    packed.bit_rec_mask = (!0u64) >> packed.word_comp_length;
    packed.length = length_in;
    packed.length_byte = if packed.length == 0 {
        0
    } else {
        (packed.length - 1) * packed.word_length / 8 + std::mem::size_of::<u64>() as u64
    };
}

#[doc = "Original inline `PackedArray::operator[]` at STAR/source/PackedArray.h:18. Args: ii: uint"]
#[inline(always)]
pub fn packedarray_h18_packedarray_index(
    packed: &crate::packed_array::PackedArray,
    ii: u64,
) -> u64 {
    let b = ii * packed.word_length;
    let byte = (b / 8) as usize;
    let shift = (b % 8) as u32;
    let a1 = unsafe {
        debug_assert!(byte + std::mem::size_of::<u64>() <= packed.char_array.len());
        std::ptr::read_unaligned(packed.char_array.as_ptr().add(byte) as *const u64)
    };
    ((a1 >> shift) << packed.word_comp_length) >> packed.word_comp_length
}

#[doc = "Original `PackedArray::writePacked` at STAR/source/PackedArray.cpp:17. Args: jj: uint, x: uint"]
pub fn packedarray_l17_packedarray_writepacked(
    packed: &mut crate::packed_array::PackedArray,
    jj: u64,
    x: u64,
) {
    let b = jj * packed.word_length;
    let B = (b / 8) as usize;
    let S = (b % 8) as u32;

    let x = x << S;
    let mut a1_bytes = [0u8; 8];
    a1_bytes.copy_from_slice(&packed.char_array[B..B + 8]);
    let a1 = u64::from_ne_bytes(a1_bytes);
    let a1 = (a1 & !(packed.bit_rec_mask << S)) | x;
    packed.char_array[B..B + 8].copy_from_slice(&a1.to_ne_bytes());
}

#[doc = "Original `PackedArray::pointArray` at STAR/source/PackedArray.cpp:27. Args: pointerCharIn: char"]
pub fn packedarray_l27_packedarray_pointarray(
    packed: &mut crate::packed_array::PackedArray,
    pointer_char_in: &[u8],
) {
    packed.char_array.clear();
    packed.char_array.extend_from_slice(pointer_char_in);
}

#[doc = "Original `PackedArray::allocateArray` at STAR/source/PackedArray.cpp:31. Args: "]
pub fn packedarray_l31_packedarray_allocatearray(
    packed: &mut crate::packed_array::PackedArray,
) {
    packed.char_array = vec![0u8; packed.length_byte as usize];
    let tail_start = (packed.length_byte as usize).saturating_sub(std::mem::size_of::<u64>());
    for ii in tail_start..packed.length_byte as usize {
        packed.char_array[ii] = 0;
    }
    packed.array_allocated = true;
}

#[doc = "Original `PackedArray::deallocateArray` at STAR/source/PackedArray.cpp:37. Args: "]
pub fn packedarray_l37_packedarray_deallocatearray(
    packed: &mut crate::packed_array::PackedArray,
) {
    if packed.array_allocated {
        packed.char_array.clear();
        packed.array_allocated = false;
    }
    packed.char_array.clear();
}
