

# Blastn-Python

| File         | Desc               |
| ------------ | ------------------ |
| Module.py    | Smith-Waterman算法 |
| blastn.ipynb | BLASTN算法实现     |

## 人类基因组 数据集下载

https://blog.csdn.net/qq_53947118/article/details/122571107

```sh
blastn -query query.fasta -db GRCh37 -out result.txt -outfmt 6 
blastn -query query.fasta -db GRCh38 -out result.txt -outfmt 6 

cd data/
wget https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
mv human_g1k_v37.fasta GRCh37.fasta
makeblastdb -in GRCh37.fasta -dbtype nucl -out GRCh37


wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta GRCh38.fasta
makeblastdb -in GRCh38.fasta -dbtype nucl -out GRCh38
```

# 微调参数对齐结果


0,13881,
160600871,160614799,

13881,33970,
160625899,160646057,

33970,56952,
160623879,160647005,

56952,61475,
160630365,160634928,

61475,67948,
160640473,160647022,

67948,84499,
160630382,160647052,

84499,90532,
160635962,160642017,

90533,122386,
160614698,160646656,

122387,126786,
160618927,160623340,

126788,145063,
160634437,160652815