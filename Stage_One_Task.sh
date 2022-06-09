#!/usr/bin/bash

#Bash code for counting the sequence
mkdir task2
cd task2
wget https://raw.githubusercontent.com/HackBio-Internship/wale-home-tasks/main/DNA.fa

# every sequence data has 4 line so whatever the number comes by total number of lines it should be divided by 4 to get the number of sequence
cat DNA.fa | wc –l

# or each sequence line starts with “>” sign. So counting the lines with “>” will list the number of line s that contains the sequence
cat DNA.fa | grep '>' | wc -l

# for printing the total number of each nucleotide in the sequence A T G C
grep -o 'A' DNA.fa | wc -l && grep -o 'T' DNA.fa | wc -l && grep -o 'G' DNA.fa | wc -l && grep -o 'C' DNA.fa | wc -l

#to install miniconda in your terminal 
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh

chmod +x Miniconda3-py38_4.12.0-Linux-x86_64.sh

./Miniconda3-py38_4.12.0-Linux-x86_64.sh

conda activate base

# code for installing fastqc , fastp and spades one by one
conda install -c bioconda fastqc
conda install -c bioconda fastp
conda install -c bioconda spades

# uploading forward and reverse read dataset from github link
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R1.fastq.gz?raw=true -O Alsen_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R2.fastq.gz?raw=true -O Alsen_R2.fastq.gz

# creating an output folder to send analysis output into it
mkdir output

# running fastqc on reads that was downloaded earlier and sending the report file into “output” folder
fastqc Alsen_R1.fastq.gz -o output/
fastqc Alsen_R2.fastq.gz -o output/

# running fastp on reads that was downloaded earlier and saving the output as new fastp file 
fastp -i Alsen_R1.fastq.gz -o Alsen_R1.fastp.gz
fastp -i Alsen_R2.fastq.gz -o Alsen_R2.fastp.gz

# running spades on reads that was downloaded earlier and sending the report file into “spades_report” folder 
spades.py -1 Alsen_R1.fastq.gz -2 Alsen_R2.fastq.gz --careful --cov-cutoff auto -o spades_report
