use std::collections::{HashMap, HashSet};
use std::path::Path;

use rust_htslib::bam::{Read, Reader, Record, Header, HeaderView};

use crate::{
    bed::load_bed,
    error::MyError,
};

/// 开始分析
pub fn run_infer(bam_file: &Path, ref_bed: &Path, sample_size: usize, q_cut: u8) -> Result<(), MyError> {
    // 读取bam文件
    //let mut bam_reader = rust_htslib::bam::IndexedReader::from_path(bam_file).map_err(|e| MyError::ReadBamError{file: bam_file.to_str().unwrap().to_string(), error: e})?; // 这样读取会导致下面循环读取record时返回None退出循环
    let mut bam_reader = Reader::from_path(bam_file).map_err(|e| MyError::ReadBamError{file: bam_file.to_str().unwrap().to_string(), error: e})?;
    // 需要根据tid从header中获取chr，参考：https://github.com/rust-bio/rust-htslib/issues/288
    let header = Header::from_template(&bam_reader.header());
    let head_view = HeaderView::from_header(&header);
    // 读取指定参考基因bed文件
    let gene_ranges = load_bed(ref_bed)?; // HashMap<String, IntervalTree>
    // 记录record数量，达到指定的sample_size数量则停止
    let mut count: usize = 0;
    let mut p_strandness: HashMap<String, f64> = HashMap::new(); // key: read_id(1/2) + map_strand(+/-) + strand_from_gene(“:”拼接的基因strand), value: count
    let mut s_strandness: HashMap<String, f64> = HashMap::new(); // key: map_strand(+/-) + strand_from_gene(“:”拼接的基因strand), value: count
    let mut read_id_map_strand_gene_strand: String;
    let mut p_strandness_sum: f64 = 0.0;
    let mut s_strandness_sum: f64 = 0.0;
    let mut read_start: u64;
    let mut read_end: u64;
    // 先声明一个record，后面使用bam的read方法往里面写入新record，避免每次重新声明record，更高效
    let mut record = Record::new();
    // 遍历bam每条read
    while let Some(result) = bam_reader.read(&mut record) {
        // 达到指定的sample_size数量则停止
        if count >= sample_size {
            break
        }
        // 这里result是`Result<()>`解析下没报错就行
        result.map_err(|e| MyError::ReadBamRecordError{file: bam_file.to_str().unwrap().to_string(), error: e})?;
        // skip low quanlity
        if record.is_quality_check_failed() {
            continue
        }
        // skip duplicate read
        if record.is_duplicate() {
            continue
        }
        // skip non primary hit
        if record.is_secondary() {
            continue
        }
        // skip unmap read
        if record.is_unmapped() {
            continue
        }
        // 舍弃质量分数低于自定义阈值的record
        if record.mapq() < q_cut {
            continue
        }
        // 根据tid从header中获取当前record的chr
        let chrom = String::from_utf8(head_view.tid2name(record.tid() as u32).to_vec()).unwrap();
        // 统计
        //println!("{}", count);
        if record.is_paired() {
            //read_id_map_strand_gene_strand = if is_in_flag(record.flags(), SamFlag::FIRST_IN_PAIR) { // First in pair, 64
            read_id_map_strand_gene_strand = if record.is_first_in_template() { // First in pair, 64
                "1".to_string()
            //} else if is_in_flag(record.flags(), SamFlag::SECOND_IN_PAIR) { // Second in pair, 128
            } else if record.is_last_in_template() { // Second in pair, 128
                "2".to_string()
            } else {
                "0".to_string() // 不应该运行到这里
            };
            if record.is_reverse() {
                read_id_map_strand_gene_strand += "-";
            } else {
                read_id_map_strand_gene_strand += "+";
            }
            read_start = record.pos() as u64;
            read_end = read_start + record.seq_len() as u64; // 这里加上原始read长度，比如150
            if let Some(tree) = gene_ranges.get(&chrom) {
                let tmp = tree.find(read_start..read_end);
                if tmp.clone().count() == 0 {
                    continue
                }
                let tmp_str = tmp.map(|s| s.data().clone()).collect::<HashSet<_>>().into_iter().collect::<Vec<String>>().join(":");
                read_id_map_strand_gene_strand += &tmp_str;
                match p_strandness.get_mut(&read_id_map_strand_gene_strand) {
                    Some(s) => *s += 1.0,
                    None => {
                        p_strandness.insert(read_id_map_strand_gene_strand, 1.0);
                    },
                }
                p_strandness_sum += 1.0;
                count += 1;
            }
        } else {
            read_id_map_strand_gene_strand = if record.is_reverse() {
                "-".to_string()
            } else {
                "+".to_string()
            };
            read_start = record.pos() as u64;
            read_end = read_start + record.seq_len() as u64; // 这里加上原始read长度，比如150
            if let Some(tree) = gene_ranges.get(&chrom) {
                let tmp = tree.find(read_start..read_end);
                if tmp.clone().count() == 0 {
                    continue
                }
                let tmp_str = tmp.map(|s| s.data().clone()).collect::<HashSet<_>>().into_iter().collect::<Vec<String>>().join(":");
                read_id_map_strand_gene_strand += &tmp_str;
                match s_strandness.get_mut(&read_id_map_strand_gene_strand) {
                    Some(s) => *s += 1.0,
                    None => {
                        s_strandness.insert(read_id_map_strand_gene_strand, 1.0);
                    },
                }
                s_strandness_sum += 1.0;
                count += 1;
            }
        }
    }

    //println!("p_keys: {:?}", p_strandness.keys());
    //println!("s_keys: {:?}", s_strandness.keys());
    // 最后统计，并打印结果
    println!("Total {} usable reads were sampled", count);
    let mut spec1: f64 = 0.0;
    let mut spec1_each: [f64;4] = [0.0, 0.0, 0.0, 0.0];
    let mut spec2: f64 = 0.0;
    let mut spec2_each: [f64;4] = [0.0, 0.0, 0.0, 0.0];
    let other: f64;
    if p_strandness.len() > 0 && s_strandness.len() == 0 {
        //spec1 = (p_strandness.get("1++").unwrap() + p_strandness.get("1--").unwrap() + p_strandness.get("2+-").unwrap() + p_strandness.get("2-+").unwrap()) / p_strandness_sum;
        if let Some(v) = p_strandness.get("1++") {
            spec1 += v;
            spec1_each[0] = v / p_strandness_sum;
        }
        if let Some(v) = p_strandness.get("1--") {
            spec1 += v;
            spec1_each[1] = v / p_strandness_sum;
        }
        if let Some(v) = p_strandness.get("2+-") {
            spec1 += v;
            spec1_each[2] = v / p_strandness_sum;
        }
        if let Some(v) = p_strandness.get("2-+") {
            spec1 += v;
            spec1_each[3] = v / p_strandness_sum;
        }
        spec1 /= p_strandness_sum;
        //spec2 = (p_strandness.get("1+-").unwrap() + p_strandness.get("1-+").unwrap() + p_strandness.get("2++").unwrap() + p_strandness.get("2--").unwrap()) / p_strandness_sum;
        if let Some(v) = p_strandness.get("1+-") {
            spec2 += v;
            spec2_each[0] = v / p_strandness_sum;
        }
        if let Some(v) = p_strandness.get("1-+") {
            spec2 += v;
            spec2_each[1] = v / p_strandness_sum;
        }
        if let Some(v) = p_strandness.get("2++") {
            spec2 += v;
            spec2_each[2] = v / p_strandness_sum;
        }
        if let Some(v) = p_strandness.get("2--") {
            spec2 += v;
            spec2_each[3] = v / p_strandness_sum;
        }
        spec2 /= p_strandness_sum;
        other = if 1.0 - spec1 - spec2 < 0.0 {
            0.0
        } else {
            1.0 - spec1 - spec2
        };
        println!("This is PairEnd Data");
        println!("Fraction of reads failed to determine: {:.4}", other);
        println!("Fraction of reads explained by \"1++,1--,2+-,2-+\": {:.4} ({:.4}, {:.4}, {:.4}, {:.4})", spec1, spec1_each[0], spec1_each[1], spec1_each[2], spec1_each[3]);
        println!("Fraction of reads explained by \"1+-,1-+,2++,2--\": {:.4} ({:.4}, {:.4}, {:.4}, {:.4})", spec2, spec2_each[0], spec2_each[1], spec2_each[2], spec2_each[3]);
    } else if s_strandness.len() > 0 && p_strandness.len() == 0 {
        //spec1 = (s_strandness.get("++").unwrap() + s_strandness.get("--").unwrap()) / s_strandness_sum;
        if let Some(v) = p_strandness.get("++") {
            spec1 += v;
            spec1_each[0] = v / s_strandness_sum;
        }
        if let Some(v) = p_strandness.get("--") {
            spec1 += v;
            spec1_each[1] = v / s_strandness_sum;
        }
        spec1 /= s_strandness_sum;
        //spec2 = (s_strandness.get("+-").unwrap() + s_strandness.get("-+").unwrap()) / s_strandness_sum;
        if let Some(v) = p_strandness.get("+-") {
            spec2 += v;
            spec2_each[0] = v / s_strandness_sum;
        }
        if let Some(v) = p_strandness.get("-+") {
            spec2 += v;
            spec2_each[1] = v / s_strandness_sum;
        }
        spec2 /= s_strandness_sum;
        other = if 1.0 - spec1 - spec2 < 0.0 {
            0.0
        } else {
            1.0 - spec1 - spec2
        };
        println!("This is SingleEnd Data");
        println!("Fraction of reads failed to determine: {:.4}", other);
        println!("Fraction of reads explained by \"++,--\": {:.4} ({:.4}, {:.4})", spec1, spec1_each[0], spec1_each[1]);
        println!("Fraction of reads explained by \"+-,-+\": {:.4} ({:.4}, {:.4})", spec2, spec2_each[0], spec2_each[1]);
    } else {
        println!("Unknown Data type");
    }
    Ok(())
}

/*
// 通过flag是否含有64或128来判断1或2：https://github.com/rust-bio/rust-htslib/pull/440
// 检查flag是否含有指定flag值，这段代码还没加到包中，因此写到这里进行调用
// 代码参考自：https://github.com/rust-bio/rust-htslib/pull/440/commits/bfb2db10aa66beb494eb5993608282b1d83e9d98
struct SamFlag;

impl SamFlag {
    //const PAIRED:          u16 = 1;
    //const PROPERLY_PAIRED: u16 = 2;
    //const READ_UNMAPPED:   u16 = 4;
    //const MATE_UNMAPPED:   u16 = 8;
    //const READ_RERVERSE:   u16 = 16;
    //const MATE_REVERSE:    u16 = 32;
    const FIRST_IN_PAIR:   u16 = 64;
    const SECOND_IN_PAIR:  u16 = 128;
    //const NOT_PRIMARY_ALN: u16 = 256;
    //const FAIL_QC:         u16 = 512;
    //const DUPLICATE:       u16 = 1024;
    //const SUPPLEMENTARY:   u16 = 2048;
}

/// This function uses bitwise operations to test if flags are set
/// # Arguments
/// * `flag`: u16 - The record flag you want to test
/// * `in_`: u16 - The flags you want to check if they are set (use 0 for no test)
/// # Usage:
/// let test if a flag is both paired and first in pair
/// ```
/// use rust_htslib::bam::flags::{Flag, is_in_flag};
/// let read_flag = 65;
/// assert_eq!(is_in_flag(read_flag, Flag::PAIRED + Flag::FIRST_IN_PAIR), true);
/// ```
fn is_in_flag(flag: u16, in_: u16) -> bool {
    if (in_ & flag) != in_ {
        return false;
    }
    true
}
*/
