use std::path::PathBuf;

use argh::FromArgs;

/// error: 定义的错误类型，用于错误传递
use crate::{
    error::MyError,
};

#[derive(FromArgs)]
/// infer experiment
struct Paras {
    /// input alignment file in SAM or BAM format
    #[argh(option, short = 'i')]
    input_file: String,

    /// reference gene model in bed fomat
    #[argh(option, short = 'r')]
    refgene: String,

    /// number of reads sampled from SAM/BAM file. default=200000
    #[argh(option, short = 's')]
    sample_size: Option<usize>,

    /// minimum mapping quality (phred scaled) for an alignment to be considered as \"uniquely mapped\". default=30
    #[argh(option, short = 'q')]
    mapq: Option<u8>,
}

/// 存储解析后的命令行参数
///#[derive(Debug, Default)]
pub struct ParsedParas {
    pub input_file:  PathBuf, // bam比对文件
    pub refgene:     PathBuf, // 相应物种基因bed文件
    pub sample_size: usize,   // 对bam前几个符合筛选条件的record进行统计，默认1000
    pub mapq:        u8,      // mapq阈值，默认30
}

/// 解析参数
pub fn parse_para() -> Result<ParsedParas, MyError> {
    let para: Paras = argh::from_env();
    let out: ParsedParas = ParsedParas{
        input_file: {
            let tmp_bam = PathBuf::from(&para.input_file);
            if !(tmp_bam.exists() && tmp_bam.is_file()) {
                return Err(MyError::FileNotExistError{file: para.input_file})
            }
            tmp_bam
        },
        refgene: {
            let tmp_bed = PathBuf::from(&para.refgene);
            if !(tmp_bed.exists() && tmp_bed.is_file()) {
                return Err(MyError::FileNotExistError{file: para.refgene})
            }
            tmp_bed
        },
        sample_size: match para.sample_size {
            Some(n) => {
                if n < 1000 {
                    println!("Warn: Sample Size ({}) too small to give a accurate estimation", n);
                }
                n
            },
            None => 200000,
        },
        mapq: match para.mapq {
            Some(m) => m,
            None => 30,
        },
    };
    Ok(out)
}
