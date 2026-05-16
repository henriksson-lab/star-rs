#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `createDirectory` at STAR/source/streamFuns.cpp:10. Args: dirPathIn: string, dirPerm: mode_t, dirParameter: string, P: Parameters"]
pub fn streamfuns_l10_createdirectory(
    dir_path_in: &str,
    dir_perm: u32,
    dir_parameter: &str,
) -> Result<String, String> {
    let dir_path = match dir_path_in.rfind('/') {
        Some(pos) => &dir_path_in[..=pos],
        None => "",
    };
    if let Err(err) = streamfuns_create_dir(dir_path, dir_perm) {
        if err.kind() == std::io::ErrorKind::AlreadyExists {
            Ok(format!(
                "{} directory exists and will be overwritten: {}\n",
                dir_parameter, dir_path
            ))
        } else {
            if dir_path.is_empty() {
                return Err(format!(
                    "EXITING because of fatal OUTPUT FILE error: could not create output directory:  for {} {}\n ERROR: {}\nSOLUTION: check the path and permissions.\n",
                    dir_parameter, dir_path_in, err
                ));
            }
            let mut i1 = dir_path[1..]
                .find('/')
                .map(|pos| pos + 1)
                .unwrap_or(dir_path.len());
            while i1 < dir_path.len() {
                let dir_path1 = &dir_path[..i1];
                if let Err(err1) = streamfuns_create_dir(dir_path1, dir_perm) {
                    if err1.kind() != std::io::ErrorKind::AlreadyExists {
                        return Err(format!(
                            "EXITING because of fatal OUTPUT FILE error: could not create output directory: {} for {} {}\n ERROR: {}\nSOLUTION: check the path and permissions.\n",
                            dir_path1, dir_parameter, dir_path_in, err1
                        ));
                    }
                }
                i1 = dir_path[i1 + 1..]
                    .find('/')
                    .map(|pos| pos + i1 + 1)
                    .unwrap_or(dir_path.len());
            }
            Ok(format!(
                "{} directory and its parents created: {}\n",
                dir_parameter, dir_path
            ))
        }
    } else {
        Ok(format!(
            "{} directory created: {}\n",
            dir_parameter, dir_path
        ))
    }
}

#[doc = "Original `fstreamReadBig` at STAR/source/streamFuns.cpp:39. Args: S: std::ifstream, A: char, N: unsigned long long"]
pub fn streamfuns_l39_fstreamreadbig<R: std::io::Read>(
    stream: &mut R,
    a: &mut [u8],
    n: u64,
) -> u64 {
    const FSTREAM_CHUNK_MAX: usize = 2_147_483_647;
    let mut c = 0usize;
    let full_chunks = n / FSTREAM_CHUNK_MAX as u64;
    for _ in 0..full_chunks {
        match stream.read(&mut a[c..c + FSTREAM_CHUNK_MAX]) {
            Ok(nread) => {
                c += nread;
                if nread < FSTREAM_CHUNK_MAX {
                    return c as u64;
                }
            }
            Err(_) => return c as u64,
        }
    }
    let rem = (n % FSTREAM_CHUNK_MAX as u64) as usize;
    if rem > 0 {
        if let Ok(nread) = stream.read(&mut a[c..c + rem]) {
            c += nread;
        }
    } else if n == 0 {
        let _ = stream.read(&mut []);
    }
    c as u64
}

#[doc = "Original `fstreamWriteBig` at STAR/source/streamFuns.cpp:51. Args: S: std::ofstream, A: char, N: unsigned long long, fileName: std::string, errorID: std::string, P: Parameters"]
pub fn streamfuns_l51_fstreamwritebig<W: std::io::Write>(
    stream: &mut W,
    a: &[u8],
    n: u64,
) -> std::io::Result<()> {
    const FSTREAM_CHUNK_MAX: usize = 2_147_483_647;
    let mut c = 0usize;
    let full_chunks = n / FSTREAM_CHUNK_MAX as u64;
    for _ in 0..full_chunks {
        stream.write_all(&a[c..c + FSTREAM_CHUNK_MAX])?;
        c += FSTREAM_CHUNK_MAX;
    }
    let rem = (n % FSTREAM_CHUNK_MAX as u64) as usize;
    if rem > 0 {
        stream.write_all(&a[c..c + rem])?;
    }
    Ok(())
}

#[doc = "Original `ofstrOpen` at STAR/source/streamFuns.cpp:91. Args: fileName: std::string, errorID: std::string, P: Parameters"]
pub fn streamfuns_l91_ofstropen(file_name: &str, error_id: &str) -> Result<std::fs::File, String> {
    std::fs::OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(file_name)
        .map_err(|_| {
            format!(
                "{}: exiting because of *OUTPUT FILE* error: could not create output file {}\nSOLUTION: check that the path exists and you have write permission for this file. Also check \"ulimit -n\" and increase it to allow more open files.\n",
                error_id, file_name
            )
        })
}

#[doc = "Original `fstrOpen` at STAR/source/streamFuns.cpp:102. Args: fileName: std::string, errorID: std::string, P: Parameters, flagDelete: bool"]
pub fn streamfuns_l102_fstropen(
    file_name: &str,
    error_id: &str,
    flag_delete: bool,
) -> Result<std::fs::File, String> {
    let opened = if flag_delete {
        std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(file_name)
    } else {
        std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(file_name)
            .or_else(|_| {
                std::fs::OpenOptions::new()
                    .read(true)
                    .write(true)
                    .create(true)
                    .truncate(true)
                    .open(file_name)
            })
    };
    opened.map_err(|_| {
        format!(
            "{}: exiting because of *OUTPUT FILE* error: could not create input/output file {}\nSolution: check that the path exists and you have write permission for this file\n",
            error_id, file_name
        )
    })
}

#[doc = "Original `ifstrOpen` at STAR/source/streamFuns.cpp:124. Args: fileName: std::string, errorID: std::string, solutionString: std::string, P: Parameters"]
pub fn streamfuns_l124_ifstropen(
    file_name: &str,
    error_id: &str,
    solution_string: &str,
) -> Result<std::fs::File, String> {
    std::fs::File::open(file_name).map_err(|_| {
        let mut err = format!(
            "{}: exiting because of *INPUT FILE* error: could not open input file {}\nSolution: check that the file exists and you have read permission for this file\n",
            error_id, file_name
        );
        if !solution_string.is_empty() {
            err.push_str(&format!("          {}\n", solution_string));
        }
        err
    })
}

#[doc = "Original `ifstrOpenGenomeFile` at STAR/source/streamFuns.cpp:139. Args: fileName: std::string, errorID: std::string, P: Parameters"]
pub fn streamfuns_l139_ifstropengenomefile(
    file_name: &str,
    error_id: &str,
    genome_dir: &str,
) -> Result<std::fs::File, String> {
    streamfuns_l124_ifstropen(
        &format!("{}/{}", genome_dir, file_name),
        error_id,
        "if this file is missing from the genome directory, you will need to *re-generate the genome*",
    )
}

#[doc = "Original `copyFile` at STAR/source/streamFuns.cpp:144. Args: fileIn: string, fileOut: string"]
pub fn streamfuns_l144_copyfile(file_in: &str, file_out: &str) -> std::io::Result<u64> {
    std::fs::copy(file_in, file_out)
}
