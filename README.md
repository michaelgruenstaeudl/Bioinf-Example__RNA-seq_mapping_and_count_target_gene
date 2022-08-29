Bioinformatics Example: Mapping RNA-Seq data and inferring count for gene of interest
=====================================================================================

### STEP 1. Select and download RNA-Seq sequence data
This example will be performed on two medium-sized paired-end RNA-Seq sequence data sets of human origin. Hence, I select NCBI SRA runs [SRR10677729](https://www.ncbi.nlm.nih.gov/sra/?term=SRR10677729) and [SRR10677731](https://www.ncbi.nlm.nih.gov/sra/?term=SRR10677731) as example input data. This data was selected at random and is used here only for demonstrative purposes.

##### Import the paired-end RNA-Seq samples from NCBI SRA using the [NCBI SRA-Toolkit](https://github.com/ncbi/sra-tools)

Environmental modules needed:
```
module load SRA-Toolkit
```

###### Sample 1
```
wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR10677729/SRR10677729.1
fastq-dump --split-files --gzip --skip-technical --read-filter pass --origfmt --readids --clip SRR10677729.1 &
mv SRR10677729.1_pass_1.fastq.gz sample1_R1.fastq.gz
mv SRR10677729.1_pass_2.fastq.gz sample1_R2.fastq.gz
```

###### Sample 2
```
wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR10677731/SRR10677731.1
fastq-dump --split-files --gzip --skip-technical --read-filter pass --origfmt --readids --clip SRR10677731.1 &
mv SRR10677731.1_pass_1.fastq.gz sample2_R1.fastq.gz
mv SRR10677731.1_pass_2.fastq.gz sample2_R2.fastq.gz
```

### STEP 2. Perform Quality control step with FastQC and MultiQC

Note: FastQC cannot handle paired-end data; QC must be run on R1 and R2 separately

###### Sample 1
```
module load FastQC
mkdir -p fastqc_output
fastqc -t 2 -f fastq -o fastqc_output *.fastq.gz 1>fastqc_runtime.log 2>&1 &
```

###### Sample 2
```
module load MultiQC
mkdir -p multiqc_output
multiqc -o multiqc_output ./fastqc_output 1>multiqc_runtime.log 2>&1 &
```

### STEP 3. Perform trimming and adapter removal step with cutadapt

Environmental modules needed:
```
module load cutadapt
```

###### Checking if FASTQ files properly paired
```
cutadapt -o /dev/null -p /dev/null -j 4 sample1_R1.fastq.gz sample1_R2.fastq.gz 1>checkPairing_sample1_runtime.log 2>&1 &
cutadapt -o /dev/null -p /dev/null -j 4 sample2_R1.fastq.gz sample2_R2.fastq.gz 1>checkPairing_sample2_runtime.log 2>&1 &
```
Note: According to the NCBI SRA records, the Truseq standard mRNA Sample Prep Kit (Illumina) was employed for library prep. The correct adapteers, thus, are: R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA; R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT; other adapater seqs: R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC; R2: TGATCGTCGGACTGTAGAACTCTGAACGTGTAGA

###### Sample 1
```
cutadapt -m 22 -O 4 -j 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o sample1_R1_trimmed.fastq.gz sample1_R1.fastq.gz 1>cutadapt_sample1_R1_runtime.log 2>&1 &
cutadapt -m 22 -O 4 -j 4 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o sample1_R2_trimmed.fastq.gz sample1_R2.fastq.gz 1>cutadapt_sample1_R2_runtime.log 2>&1 &
```

###### Sample 2
```
cutadapt -m 22 -O 4 -j 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o sample2_R1_trimmed.fastq.gz sample2_R1.fastq.gz 1>cutadapt_sample2_R1_runtime.log 2>&1 &
cutadapt -m 22 -O 4 -j 4 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o sample2_R2_trimmed.fastq.gz sample2_R2.fastq.gz 1>cutadapt_sample2_R2_runtime.log 2>&1 &
```

### STEP 4. Perform mapping with RNA STAR

Environmental modules needed:
```
module load STAR
```

```
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.p13.genome.fa.gz
gunzip -k gencode.v40.annotation.gtf.gz
gunzip -k GRCh38.p13.genome.fa.gz
```

##### Build index for chromosome 17

```
mkdir -p GRCh38_chr17_index
STAR --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir GRCh38_chr17_index \
  --genomeFastaFiles GRCh38.p13.genome_chr17.fa \
  --sjdbGTFfile gencode.v40.annotation_chr17.gtf \
  --sjdbOverhang 149 \
  --outFileNamePrefix chr17 \
  1>STAR_buildingIndex_runtime.log 2>&1 &
```

Explanations for the above command: 
*  --runThreadN: number of threads
*  --runMode: genomeGenerate mode
*  --genomeDir: /path/to/store/genome_indices
*  --genomeFastaFiles: /path/to/FASTA_file
*  --sjdbGTFfile: /path/to/GTF_file
*  --sjdbOverhang: readlength -1

##### Perform mapping itself using the software RNA STAR
```
STAR --genomeDir GRCh38_chr17_index \
  --readFilesIn sample2_R1_trimmed.fastq.gz sample2_R2_trimmed.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outFileNamePrefix GRCh38_chr17_mapping/sample2_ \
  1>STAR_mapping_sample2_runtime.log 2>&1 &
```


##### Optional: Removing duplicates
Note: The need for duplicate removal was indicated by the results of MultiQC (see above)

Environmental modules needed:
```
module load picard
```

```
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  I=input.bam \
  O=marked_duplicates.bam \
  M=marked_dup_metrics.txt \
  REMOVE_DUPLICATE=true
```
