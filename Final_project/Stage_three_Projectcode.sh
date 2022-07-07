#!/usr/bin/bash

# downloading the datasets
mkdir reads
cd reads

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR987/006/SRR9873306/SRR9873306_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR987/006/SRR9873306/SRR9873306_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR133/095/SRR13342195/SRR13342195_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR133/095/SRR13342195/SRR13342195_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/021/SRR14160921/SRR14160921_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/021/SRR14160921/SRR14160921_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/094/SRR17859194/SRR17859194_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/094/SRR17859194/SRR17859194_2.fastq.gz
wget
wget

cd ..

cat > list.txt
SRR9873306
SRR13342195
SRR14160921
SRR17859194

CTRL + D

# quality checking of the dataset

mkdir -p fastqc_results 

for sample in `cat list.txt`
do
  fastqc reads/${sample}*.fastq.gz -o fastqc_results
done

multiqc fastqc_results -o multiqc_results

# Trimming the reads

mkdir -p trimmed_reads

for sample in `cat list.txt`
do
      fastp
-i "reads/${sample}_1.fastq.gz"
-I "reads/${sample}_2.fastq.gz"
-o "trimmed_reads/${sample}.trimmed_1.fastq.gz"
-O "trimmed_reads/${sample}.trimmed_2.fastq.gz"
-h "trimmed_reads/${sample}.fastp.html"
done

#assembling the reads

For sample in `cat list.txt`
do
spades.py -1 "reads/${sample}_1.fastq.gz" -2 "reads/${sample}_2.fastq.gz" --careful --cov-cutoff auto -o spades_report
done
