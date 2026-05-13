use std::io::{Read, Seek};
use std::path::Path;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

fn file_is_gzip(file: &mut std::fs::File) -> std::io::Result<bool> {
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    file.rewind()?;
    Ok(n == magic.len() && magic == GZIP_MAGIC)
}

pub fn read_bytes_auto_gzip(path: impl AsRef<Path>) -> std::io::Result<Vec<u8>> {
    let mut file = std::fs::File::open(path)?;
    let mut contents = Vec::new();
    if file_is_gzip(&mut file)? {
        let mut decoder = flate2::read::MultiGzDecoder::new(file);
        decoder.read_to_end(&mut contents)?;
    } else {
        let len = file.metadata()?.len() as usize;
        contents.reserve(len);
        file.read_to_end(&mut contents)?;
        if contents.len() != len {
            return Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "file length changed while reading",
            ));
        }
    }
    Ok(contents)
}

pub fn read_to_string_auto_gzip(path: impl AsRef<Path>) -> std::io::Result<String> {
    let mut file = std::fs::File::open(path)?;
    let mut contents = String::new();
    if file_is_gzip(&mut file)? {
        let mut decoder = flate2::read::MultiGzDecoder::new(file);
        decoder.read_to_string(&mut contents)?;
    } else {
        file.read_to_string(&mut contents)?;
    }
    Ok(contents)
}
