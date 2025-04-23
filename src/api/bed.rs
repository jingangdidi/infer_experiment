use std::collections::HashMap;
use std::path::Path;

use bio::{
    data_structures::interval_tree::IntervalTree,
    io::bed::Reader,
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
