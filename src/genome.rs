#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `Genome` at STAR/source/Genome.h:13."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genome {
    pub g: Vec<u8>,
    pub sa: Vec<u64>,
    pub sa_packed: PackedArray,
    pub sa_insert: PackedArray,
    pub sa_pass2: Vec<u64>,
    pub sai: Vec<u64>,
    pub sai_packed: PackedArray,
    pub genome_sa_index_start: Vec<u64>,
    pub n_genome: u64,
    pub n_sa: u64,
    pub n_sa_byte: u64,
    pub n_chr_real: u32,
    pub genome_chr_bin_nbases: u32,
    pub chr_bin_n: u32,
    pub gstrand_bit: u32,
    pub gstrand_mask: u64,
    pub sai_mark_absent_mask_c: u64,
    pub sai_mark_nmask: u64,
    pub sai_mark_nmask_c: u64,
    pub sjdb_overhang: u32,
    pub sjdb_length: u32,
    pub sj_gstart: u64,
    pub sj_dstart: Vec<u64>,
    pub sj_astart: Vec<u64>,
    pub sjdb_start: Vec<u64>,
    pub sjdb_end: Vec<u64>,
    pub sjdb_motif: Vec<u8>,
    pub sjdb_shift_left: Vec<u8>,
    pub sjdb_shift_right: Vec<u8>,
    pub sjdb_strand: Vec<u8>,
    pub sjdb_n: u32,
    pub p_ge: ParametersGenome,
    pub chr_start: Vec<u64>,
    pub chr_length: Vec<u64>,
    pub chr_bin: Vec<u32>,
    pub chr_name: Vec<String>,
    pub chr_name_all: Vec<String>,
    pub chr_name_index: std::collections::BTreeMap<String, u64>,
    pub chr_length_all: Vec<u32>,
    pub genome_insert_l: u64,
    pub n_sa_insert: u64,
    pub var: Variation,
    pub genome_out: GenomeOut,
    pub align_sjdb_overhang_min: u32,
    pub align_intron_min: u32,
}

impl Genome {
    /// SAi[idx]. Matches C++ STAR which keeps SAi as a packed array and
    /// unpacks on each access. At production load time `sai_packed` is
    /// populated and `sai` is left empty (saves ~2.86 GB for a human
    /// genomeSAindexNbases=14). Some leaf tests construct `sai` directly
    /// without packing; honor those by falling back to the Vec.
    #[inline]
    pub fn sai_value(&self, idx: u64) -> u64 {
        if self.sai_packed.array_allocated {
            crate::packed_array::packedarray_h18_packedarray_index(&self.sai_packed, idx)
        } else {
            self.sai[idx as usize]
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GenomeOut {
    pub conv_yes: bool,
    pub gaps_are_junctions: bool,
    pub conv_file: String,
    pub conv_blocks: Vec<[u64; 3]>,
    pub n_minus_strand_offset: u64,
}

#[doc = "Original `Genome::Genome` at STAR/source/Genome.cpp:15. Args: P: Parameters, pGe: ParametersGenome"]
pub fn genome_l15_genome_genome(
    p_ge: crate::parameters_genome::ParametersGenome,
) -> crate::genome::Genome {
    let sjdb_overhang = p_ge.sjdb_overhang;
    let sjdb_length = if p_ge.sjdb_overhang == 0 {
        0
    } else {
        p_ge.sjdb_overhang * 2 + 1
    };
    let genome_chr_bin_nbases = 1_u32.checked_shl(p_ge.g_chr_bin_nbits).unwrap_or(0);
    crate::genome::Genome {
        p_ge,
        genome_chr_bin_nbases,
        sjdb_overhang,
        sjdb_length,
        ..Default::default()
    }
}

#[doc = "Original `Genome::freeMemory` at STAR/source/Genome.cpp:33. Args: "]
pub fn genome_l33_genome_freememory(map_gen: &mut crate::genome::Genome) -> bool {
    if map_gen.p_ge.g_load == "NoSharedMemory" {
        map_gen.g.clear();
        map_gen.sa.clear();
        map_gen.sa_pass2.clear();
        map_gen.sai.clear();
        true
    } else {
        false
    }
}

#[doc = "Original `Genome::OpenStream` at STAR/source/Genome.cpp:48. Args: name: string, stream: ifstream, size: uint"]
pub fn genome_l48_genome_openstream(
    map_gen: &crate::genome::Genome,
    name: &str,
    size: u32,
    p: &crate::parameters_chimeric::Parameters,
    raw_time: libc::time_t,
) -> Result<(std::fs::File, u32, String), crate::out_sj::ExitWithErrorResult> {
    use std::io::{Seek, SeekFrom};

    let file_name = format!("{}/{}", map_gen.p_ge.g_dir, name);
    let mut stream = match std::fs::File::open(&file_name) {
        Ok(stream) => stream,
        Err(_) => {
            let err_out = format!(
                "EXITING because of FATAL ERROR: could not open genome file: {}\nSOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permissions\n",
                file_name
            );
            return Err(errorwarning_l8_exitwitherror(
                &err_out,
                true,
                true,
                EXIT_CODE_GENOME_FILES,
                p,
                raw_time,
            ));
        }
    };

    if size > 0 {
        Ok((
            stream,
            size,
            format!("{}: size given as a parameter = {}\n", name, size),
        ))
    } else {
        let size1 = stream.metadata().map(|m| m.len() as i64).unwrap_or(-1);
        if size1 <= 0 {
            let err_out = format!(
                "EXITING because of FATAL ERROR: failed reading from genome file: {}\nSOLUTION: re-generate the genome index\n",
                file_name
            );
            return Err(errorwarning_l8_exitwitherror(
                &err_out, true, true, 1, p, raw_time,
            ));
        }
        let size_out = size1 as u32;
        let _ = stream.seek(SeekFrom::Start(0));
        Ok((
            stream,
            size_out,
            format!(
                "Checking {} sizefile size: {} bytes; state: good=1 eof=0 fail=0 bad=0\n",
                name, size_out
            ),
        ))
    }
}

#[doc = "Original `Genome::HandleSharedMemoryException` at STAR/source/Genome.cpp:83. Args: exc: SharedMemoryException, shmSize: uint64"]
pub fn genome_l83_genome_handlesharedmemoryexception(
    exc: &crate::shared_memory::SharedMemoryException,
    shm_size: u64,
    p: &crate::parameters_chimeric::Parameters,
    raw_time: libc::time_t,
) -> crate::out_sj::ExitWithErrorResult {
    let mut err_out = format!(
        "Shared memory error: {}, errno: {}({})\n",
        exc.error_code, exc.error_detail, exc.error_detail
    );

    let mut exit_code = EXIT_CODE_SHM;
    match exc.error_code {
        SHM_EOPENFAILED => {
            err_out.push_str("EXITING because of FATAL ERROR: problems with shared memory: error from shmget() or shm_open().\n");
            err_out.push_str("SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n");
        }
        SHM_EEXISTS => {
            err_out.push_str(
                "EXITING: fatal error from shmget() trying to allocate shared memory piece.\n",
            );
            err_out.push_str(&format!(
                "Possible cause 1: not enough RAM. Check if you have enough RAM of at least {} bytes\n",
                shm_size + 2_000_000_000
            ));
            err_out.push_str(&format!(
                "Possible cause 2: not enough virtual memory allowed with ulimit. SOLUTION: run ulimit -v {}\n",
                shm_size + 2_000_000_000
            ));
            err_out.push_str("Possible cause 3: allowed shared memory size is not large enough. SOLUTIONS: (i) consult STAR manual on how to increase shared memory allocation; (ii) ask your system administrator to increase shared memory allocation; (iii) run STAR with --genomeLoad NoSharedMemory\n");
        }
        SHM_EFTRUNCATE => {
            err_out.push_str("EXITING: fatal error from ftruncate() error shared memory.\n");
            err_out.push_str(&format!(
                "Possible cause 1: not enough RAM. Check if you have enough RAM of at least {} bytes\n",
                shm_size + 2_000_000_000
            ));
            exit_code = EXIT_CODE_MEMORY_ALLOCATION;
        }
        SHM_EMAPFAILED => {
            err_out.push_str("EXITING because of FATAL ERROR: problems with shared memory: error from shmat() while trying to get address of the shared memory piece.\n");
            err_out.push_str("SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n");
        }
        SHM_ECLOSE => {
            err_out.push_str(
                "EXITING because of FATAL ERROR: could not close the shared memory object.\n",
            );
        }
        SHM_EUNLINK => {
            err_out.push_str("EXITING because of FATAL ERROR: problems with shared memory: error from shmctl() while trying to remove shared memory piece.\n");
            err_out.push_str("SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n");
        }
        _ => {
            err_out.push_str("EXITING because of FATAL ERROR: There was an issue with the shared memory allocation. Try running STAR with --genomeLoad NoSharedMemory to avoid using shared memory.");
        }
    }

    errorwarning_l8_exitwitherror(&err_out, true, true, exit_code, p, raw_time)
}

#[doc = "Original `Genome::chrInfoLoad` at STAR/source/Genome.cpp:139. Args: "]
pub fn genome_l139_genome_chrinfoload(
    map_gen: &mut crate::genome::Genome,
    p: &crate::parameters_chimeric::Parameters,
    raw_time: libc::time_t,
) -> Result<String, crate::out_sj::ExitWithErrorResult> {
    let chr_name_file = format!("{}/chrName.txt", map_gen.p_ge.g_dir);
    let chr_name_contents = match std::fs::read_to_string(&chr_name_file) {
        Ok(contents) => contents,
        Err(_) => {
            let err_out = format!(
                "EXITING because of FATAL error, could not open file {}\nSOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n",
                chr_name_file
            );
            return Err(errorwarning_l8_exitwitherror(
                &err_out,
                true,
                true,
                EXIT_CODE_INPUT_FILES,
                p,
                raw_time,
            ));
        }
    };

    map_gen.chr_name.clear();
    for line in chr_name_contents.lines() {
        if line.is_empty() {
            break;
        }
        map_gen.chr_name.push(line.to_string());
    }
    map_gen.n_chr_real = map_gen.chr_name.len() as u32;
    map_gen.chr_start.resize(map_gen.n_chr_real as usize + 1, 0);
    map_gen.chr_length.resize(map_gen.n_chr_real as usize, 0);

    let chr_length_file = format!("{}/chrLength.txt", map_gen.p_ge.g_dir);
    let chr_length_contents = match std::fs::read_to_string(&chr_length_file) {
        Ok(contents) => contents,
        Err(_) => {
            let err_out = format!(
                "EXITING because of FATAL error, could not open file {}\nSOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n",
                chr_length_file
            );
            return Err(errorwarning_l8_exitwitherror(
                &err_out,
                true,
                true,
                EXIT_CODE_INPUT_FILES,
                p,
                raw_time,
            ));
        }
    };
    for (ii, value) in chr_length_contents
        .split_whitespace()
        .take(map_gen.n_chr_real as usize)
        .enumerate()
    {
        map_gen.chr_length[ii] = value.parse::<u64>().unwrap_or(0);
    }

    let chr_start_file = format!("{}/chrStart.txt", map_gen.p_ge.g_dir);
    let chr_start_contents = match std::fs::read_to_string(&chr_start_file) {
        Ok(contents) => contents,
        Err(_) => {
            let err_out = format!(
                "EXITING because of FATAL error, could not open file {}\nSOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n",
                chr_start_file
            );
            return Err(errorwarning_l8_exitwitherror(
                &err_out,
                true,
                true,
                EXIT_CODE_INPUT_FILES,
                p,
                raw_time,
            ));
        }
    };
    for (ii, value) in chr_start_contents
        .split_whitespace()
        .take(map_gen.n_chr_real as usize + 1)
        .enumerate()
    {
        map_gen.chr_start[ii] = value.parse::<u64>().unwrap_or(0);
    }

    let mut log = format!(
        "Number of real (reference) chromosomes= {}\n",
        map_gen.n_chr_real
    );
    map_gen.chr_name_index.clear();
    for ii in 0..map_gen.n_chr_real as usize {
        log.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            ii + 1,
            map_gen.chr_name[ii],
            map_gen.chr_length[ii],
            map_gen.chr_start[ii]
        ));
        map_gen
            .chr_name_index
            .insert(map_gen.chr_name[ii].clone(), ii as u64);
    }

    map_gen.p_ge.chr_set_mito.clear();
    for cm in &map_gen.p_ge.chr_set_mito_strings {
        let ind1 = map_gen
            .chr_name
            .iter()
            .position(|chr_name| chr_name == cm)
            .unwrap_or(map_gen.chr_name.len()) as u64;
        map_gen.p_ge.chr_set_mito.insert(ind1);
    }

    Ok(log)
}

#[doc = "Original `Genome::chrBinFill` at STAR/source/Genome.cpp:209. Args: "]
pub fn genome_l209_genome_chrbinfill(map_gen: &mut crate::genome::Genome) {
    let n_chr_real = if map_gen.n_chr_real == 0 {
        map_gen.chr_name.len() as u32
    } else {
        map_gen.n_chr_real
    };
    let genome_chr_bin_nbases = if map_gen.genome_chr_bin_nbases == 0 {
        1_u32.checked_shl(map_gen.p_ge.g_chr_bin_nbits).unwrap_or(0)
    } else {
        map_gen.genome_chr_bin_nbases
    };

    map_gen.n_chr_real = n_chr_real;
    map_gen.genome_chr_bin_nbases = genome_chr_bin_nbases;
    map_gen.chr_bin_n =
        (map_gen.chr_start[n_chr_real as usize] / genome_chr_bin_nbases as u64 + 1) as u32;
    map_gen.chr_bin.resize(map_gen.chr_bin_n as usize, 0);

    let mut ichr = 1_usize;
    for ii in 0..map_gen.chr_bin_n as usize {
        if (ii as u64) * genome_chr_bin_nbases as u64 >= map_gen.chr_start[ichr] {
            ichr += 1;
        }
        map_gen.chr_bin[ii] = (ichr - 1) as u32;
    }
}

#[doc = "Original `Genome::genomeSequenceAllocate` at STAR/source/Genome.cpp:219. Args: nGenomeIn: uint64, nG1allocOut: uint64, Gout: char, G1out: char"]
pub fn genome_l219_genome_genomesequenceallocate(
    n_genome_in: u64,
    p: &crate::parameters_chimeric::Parameters,
    raw_time: libc::time_t,
) -> Result<(u64, Vec<u8>, usize), crate::out_sj::ExitWithErrorResult> {
    let n_g1_alloc_out = (n_genome_in + 100) * 2;
    let required = n_g1_alloc_out + n_g1_alloc_out / 3;

    if p.limit_genome_generate_ram < required {
        let err_out = format!(
            "EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM={}is too small for your genome\nSOLUTION: please specify --limitGenomeGenerateRAM not less than {} and make that much RAM available \n",
            p.limit_genome_generate_ram, required
        );
        return Err(errorwarning_l8_exitwitherror(
            &err_out,
            true,
            true,
            EXIT_CODE_INPUT_FILES,
            p,
            raw_time,
        ));
    }

    Ok((
        n_g1_alloc_out,
        vec![GENOME_SPACING_CHAR; n_g1_alloc_out as usize],
        100,
    ))
}
