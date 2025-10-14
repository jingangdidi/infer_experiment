use infer_experiment::{
    parse_paras::parse_para,
    error::MyError,
    infer::run_infer,
};

fn main() {
    if let Err(e) = run() {
        println!("{}", e); // 这里不要用`{:?}`，会打印结构体而不是打印指定的错误信息
    }
}

fn run() -> Result<(), MyError> {
    // 解析参数
    let paras = parse_para()?;

    // 开始统计
    run_infer(&paras.input_file, paras.refgene, paras.gtf, &paras.feature, paras.sample_size, paras.mapq)
}
