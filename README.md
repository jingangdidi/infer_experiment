# infer_experiment
RSeQC infer_experiment written in rust, compiles into a single executable.

# download releases
[download release](https://github.com/jingangdidi/infer_experiment/releases)

# usage
```
Usage: infer_experiment -i <input-file> [-r <refgene>] [-g <gtf>] [-f <feature>] [-s <sample-size>] [-q <mapq>]

infer experiment

Options:
  -i, --input-file  input alignment file in SAM or BAM format
  -r, --refgene     reference gene model in bed fomat
  -g, --gtf         reference gtf file
  -f, --feature     gtf feature, default: gene
  -s, --sample-size number of reads sampled from SAM/BAM file. default=200000
  -q, --mapq        minimum mapping quality (phred scaled) for an alignment to be considered as "uniquely mapped". default=30
  -h, --help        display usage information
```

# example
1. use `-i <bam>` and `-r <bed>`:
```
./infer_experiment -i test.bam -r hg38_GENCODE_V42_Basic.bed
```
output result:
```
Total 200000 usable reads were sampled
This is PairEnd Data
Fraction of reads failed to determine: 0.0769
Fraction of reads explained by "1++,1--,2+-,2-+": 0.8897 (0.2380, 0.2076, 0.2069, 0.2371)
Fraction of reads explained by "1+-,1-+,2++,2--": 0.0334 (0.0074, 0.0091, 0.0092, 0.0077)
```
2. use `-i <bam>` and `-g <gtf>`:
```
./infer_experiment -i test.bam -g hg38.gtf
```
output result:
```
Total 200000 usable reads were sampled
This is PairEnd Data
Fraction of reads failed to determine: 0.0769
Fraction of reads explained by "1++,1--,2+-,2-+": 0.8897 (0.2380, 0.2076, 0.2069, 0.2371)
Fraction of reads explained by "1+-,1-+,2++,2--": 0.0334 (0.0074, 0.0091, 0.0092, 0.0077)
```

# Building from source
```
git clone https://github.com/jingangdidi/infer_experiment.git
cd infer_experiment
cargo build --release
```

# Related tools
[RSeQC](https://rseqc.sourceforge.net/#infer-experiment-py)
