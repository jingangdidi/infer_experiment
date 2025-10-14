use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use bio::{
    data_structures::interval_tree::IntervalTree,
    io::{
        bed::Reader,
        gff,
    },
    utils::Interval,
};
use bio_types::strand::Strand;

use crate::{
    error::MyError,
    utils::my_reader,
};

/// 读取指定参考基因bed文件
/// infer_experiment.py使用bx-python的IntervalTree存储bed位置，这里使用rust-bio的IntervalTree
/// https://github.com/rust-bio/rust-bio/issues/459
/// https://docs.rs/bio/latest/bio/data_structures/interval_tree/struct.IntervalTree.html
/// https://docs.rs/bio/2.2.0/bio/io/bed/index.html
pub fn load_bed(ref_bed: &Path) -> Result<HashMap<String, IntervalTree<u64, String>>, MyError> {
    // 读取bed文件
    //let mut bed_reader = Reader::from_file(ref_bed).map_err(|e| MyError::ReadBedError{file: ref_bed.to_str().unwrap().to_string(), error: e.into()})?;
    let mut bed_reader = Reader::new(my_reader(ref_bed)?); // 使用my_reader支持读取bed或bed.gz
    let mut records = bed_reader.records();
    // 存储bed位置和链信息，相同chr存储在一起
    let mut gene_ranges: HashMap<String, IntervalTree<u64, String>> = HashMap::new(); // key: chr, value: IntervalTree
    // 遍历每个record
    let mut chr: String;
    while let Some(Ok(record)) = records.next() {
        chr = record.chrom().to_string();
        if !gene_ranges.contains_key(&chr) {
            gene_ranges.insert(chr.clone(), IntervalTree::new());
        }
        let strand = match record.strand() {
            Some(s) => match s {
                Strand::Forward => "+".to_string(),
                Strand::Reverse => "-".to_string(),
                Strand::Unknown => "*".to_string(),
            },
            None => "*".to_string(),
        };
        let tree = gene_ranges.get_mut(&chr).unwrap();
        tree.insert(Interval::new(record.start()..record.end()).unwrap(), strand);
    }
    Ok(gene_ranges)
}

/// load gtf file, get bed6
pub fn load_gtf(gtf: &Path, feature: &str) -> Result<HashMap<String, IntervalTree<u64, String>>, MyError> {
    let mut gene_ranges: HashMap<String, IntervalTree<u64, String>> = HashMap::new(); // key: chr, value: IntervalTree
    let f = File::open(gtf).map_err(|e| MyError::ReadFileError{file: gtf.to_str().unwrap().to_string(), error: e})?;
    let mut reader = gff::Reader::new(f, gff::GffType::GTF2);
    let mut strand: &str;
    for record in reader.records() {
        let rec = record.map_err(|e| MyError::GtfRecordError{file: gtf.to_str().unwrap().to_string(), error: e.into()})?;
        if rec.feature_type() == feature {
            strand = match rec.strand() {
                Some(Strand::Forward) => "+",
                Some(Strand::Reverse) => "-",
                Some(Strand::Unknown) => "",
                None => "",
            };
            let tree = gene_ranges.get_mut(rec.seqname()).unwrap();
            tree.insert(Interval::new((*rec.start())..(*rec.end())).unwrap(), strand.to_string());
        }
    }
    Ok(gene_ranges)
}
