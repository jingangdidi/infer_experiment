use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::GzDecoder;

use crate::error::MyError;

/// 读取bed或bed.gz，这样可以返回bed或bed.gz的reader
/// 参考：https://users.rust-lang.org/t/write-to-normal-or-gzip-file-transparently/35561/2
/// 参考：https://github.com/rust-lang/flate2-rs/issues/393
pub fn my_reader(file: &Path) -> Result<Box<dyn BufRead>, MyError> {
    let opened_file = File::open(&file).map_err(|e| MyError::ReadFileError{file: file.to_str().unwrap().to_string(), error: e})?;
    if file.extension().unwrap() == "gz" {
        //Box::new(BufReader::with_capacity(8 * 1024, GzDecoder::new(opened_file)))
        Ok(Box::new(BufReader::new(GzDecoder::new(opened_file))))
    } else {
        //Box::new(BufReader::with_capacity(8 * 1024, opened_file))
        Ok(Box::new(BufReader::new(opened_file)))
    }
}
