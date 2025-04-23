# infer_experiment
RSeQC infer_experiment written in rust, compiles into a single executable.

# download releases
[download release](https://github.com/jingangdidi/infer_experiment/releases/tag/v0.1.0)

# usage
```
Usage: infer_experiment -i <input-file> -r <refgene> [-n <sample-size>] [-q <mapq>]

infer experiment

Options:
  -i, --input-file  input alignment file in SAM or BAM format
  -r, --refgene     reference gene model in bed or bed.gz fomat
  -s, --sample-size number of reads sampled from SAM/BAM file. default=200000
  -q, --mapq        minimum mapping quality (phred scaled) for an alignment to be considered as "uniquely mapped". default=30
  --help, help      display usage information
```

# example
```
./infer_experiment -i test.bam -r hg38_GENCODE_V42_Basic.bed
```

# result
```
Total 200000 usable reads were sampled
This is PairEnd Data
Fraction of reads failed to determine: 0.0769
Fraction of reads explained by "1++,1--,2+-,2-+": 0.8897 (0.2380, 0.2076, 0.2069, 0.2371)
Fraction of reads explained by "1+-,1-+,2++,2--": 0.0334 (0.0074, 0.0091, 0.0092, 0.0077)
```

# Related tools
[RSeQC](https://rseqc.sourceforge.net/#infer-experiment-py)
