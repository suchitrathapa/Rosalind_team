#!/usr/bin/bash

# logging into the server
ssh einstein@7.tcp.eu.ngrok.io -p 18467

einstein

# creating a folder of my name and entering it
cd Rosalind
mkdir Suchitra
cd Suchitra

# downloading the data into a newly created folder name raw_data

mkdir -p raw_data 
cd raw_data
	
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz	
cd ..
# downloading the reference sequence
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#unzipping the reference sequence

gunzip hg19.chr5_12_17.fa.gz

# create a txt file and input the sequence name for later use

cat > list.txt

SLGFSK-N_231335
SLGFSK-T_231336

ctrl + D

# quality checking using fastQC 

mkdir -p Fastqc_Reports 

for sample in `cat list.txt`
do
	fastqc raw_data/${sample}*.fastq.gz -o Fastqc_Reports
done

multiqc Fastqc_Reports -o Multiqc_Reports

# Trimming the reads

mkdir -p trimmed_reads
cd trimmed_reads
mkdir trimmed_fastqc_results
cd ..

for sample in `cat list.txt`
do
       trimmomatic PE -threads 8 raw_data/${sample}_r1_chr5_12_17.fastq.gz raw_data/${sample}_r2_chr5_12_17.fastq.gz \
               trimmed_reads/${sample}_r1_paired.fq.gz trimmed_reads/${sample}_r1_unpaired.fq.gz \
               trimmed_reads/${sample}_r2_paired.fq.gz trimmed_reads/${sample}_r2_unpaired.fq.gz \
               ILLUMINACLIP:$HOME/sample_tut/TruSeq3-PE.fa:2:30:10:8:keepBothReads \
               LEADING:3 TRAILING:10 MINLEN:25
done

fastqc  trimmed_reads/${sample}_r1_paired.fq.gz  trimmed_reads/${sample}_r2_paired.fq.gz \
                 -o trimmed_reads/Fastqc_trimresults
                 
multiqc  trimmed_reads/Fastqc_trimresults  -o trimmed_reads/Multiqc_trimresults

#Index reference file	
bwa index hg19.chr5_12_17.fa 

# Mapping or Aligning the reads 

mkdir Mapping

bwa mem -R '@RG\tID:231335\tSM:Normal' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-N_231335_r1_paired.fq.gz \
       trimmed_reads/SLGFSK-N_231335_r2_paired.fq.gz > Mapping/SLGFSK-N_231335.sam

bwa mem -R '@RG\tID:231336\tSM:Tumor' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-T_231336_r1_paired.fq.gz \
        trimmed_reads/SLGFSK-T_231336_r2_paired.fq.gz > Mapping/SLGFSK-T_231336.sam
        
        
        
        
        
        
        