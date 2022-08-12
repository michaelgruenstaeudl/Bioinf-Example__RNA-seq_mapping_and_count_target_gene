Bioinformatics Example: Mapping RNA-Seq data and inferring count for gene of interest
=====================================================================================

### STEP 1. Select and download RNA-Seq sequence data
This example will be performed on two medium-sized paired-end RNA-Seq sequence data sets of human origin. Hence, I select NCBI SRA runs [SRR10677729](https://www.ncbi.nlm.nih.gov/sra/?term=SRR10677729) and [SRR10677731](https://www.ncbi.nlm.nih.gov/sra/?term=SRR10677731) as example input data. This data was selected at random and is used here only for demonstrative purposes.

#####  Import the paired-end RNA-Seq samples from NCBI SRA using the [NCBI SRA-Toolkit](https://github.com/ncbi/sra-tools)
```
module load SRA-Toolkit
```

Download data of sample 1
```
wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR10677729/SRR10677729.1
fastq-dump --split-files --gzip --skip-technical --read-filter pass --origfmt --readids --clip SRR10677729.1 &
mv SRR10677729.1_pass_1.fastq.gz sample1_R1.fastq.gz
mv SRR10677729.1_pass_2.fastq.gz sample1_R2.fastq.gz
```

Download data of sample 2
```
wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR10677731/SRR10677731.1
fastq-dump --split-files --gzip --skip-technical --read-filter pass --origfmt --readids --clip SRR10677731.1 &
mv SRR10677731.1_pass_1.fastq.gz sample2_R1.fastq.gz
mv SRR10677731.1_pass_2.fastq.gz sample2_R2.fastq.gz
```
