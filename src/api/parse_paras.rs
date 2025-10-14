use std::path::PathBuf;

use argh::FromArgs;

/// error: 定义的错误类型，用于错误传递
use crate::{
    error::MyError,
};

#[derive(FromArgs)]
#[argh(help_triggers("-h", "--help"))] // https://github.com/google/argh/pull/106
/// infer experiment
struct Paras {
    /// input alignment file in SAM or BAM format
    #[argh(option, short = 'i')]
    input_file: String,

    /// reference gene model in bed fomat
    #[argh(option, short = 'r')]
    refgene: Option<String>,

    /// reference gtf file
    #[argh(option, short = 'g')]
    gtf: Option<String>,

    /// gtf feature, default: gene
    #[argh(option, short = 'f')]
    feature: Option<String>,

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
    pub input_file:  PathBuf,         // bam比对文件
    pub refgene:     Option<PathBuf>, // 相应物种基因bed文件
    pub gtf:         Option<PathBuf>, // gtf file
    pub feature:     String,          // gtf feature
    pub sample_size: usize,           // 对bam前几个符合筛选条件的record进行统计，默认1000
    pub mapq:        u8,              // mapq阈值，默认30
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
        refgene: match para.refgene {
            Some(bed) => {
                let tmp_bed = PathBuf::from(&bed);
                if !(tmp_bed.exists() && tmp_bed.is_file()) {
                    return Err(MyError::FileNotExistError{file: bed})
                }
                // -f only valid for -g
                if para.feature.is_some() {
                    println!("Warning - -f only valid for -g");
                }
                Some(tmp_bed)
            },
            None => None,
        },
        gtf: match para.gtf {
            Some(g) => {
                let tmp_gtf = PathBuf::from(&g);
                if !(tmp_gtf.exists() && tmp_gtf.is_file()) {
                    return Err(MyError::FileNotExistError{file: g})
                }
                Some(tmp_gtf)
            },
            None => None,
        },
        feature: match para.feature {
            Some(f) => f,
            None => "gene".to_string(),
        },
        sample_size: match para.sample_size {
            Some(n) => {
                if n < 1000 {
                    println!("Warning - Sample Size ({}) too small to give a accurate estimation", n);
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
    // must specify -r or -g, but cannot use -r and -g simultaneously
    match (out.refgene.is_some(), out.gtf.is_some()) {
        (true, true) => return Err(MyError::ParaError{para: "couldn't specify -r and -g simultaneously".to_string()}),
        (true, false) => (),
        (false, true) => (),
        (false, false) => return Err(MyError::ParaError{para: "you must specify -r or -g".to_string()}),
    }

    Ok(out)
}
