#!/usr/bin/bash
# Bash code for printing name
NAME="Suchitra"
SURNAME="Thapa"
echo "$NAME $SURNAME"
echo "$NAME "
echo "$SURNAME"
#script for story 1
mkdir suchitra
mkdir biocomputing && cd biocomputing
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna 
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
mv wildtype.fna ../suchitra/
rm wildtype.gbk.1
cd suchitra
grep tatatata wildtype.fna
grep -n tatatata wildtype.fna > mutanttype.fna
clear && history
ls && ../biocomputing
cd


#script for story 2
#Answer 1
sudo apt install figlet toilet
figlet Suchitra
figlet -c Suchitra Thapa

#Answer 3
mkdir compare
cd compare
wget https://www.bioinformatics.babraham.ac.uk/training/Introduction%20to%20Unix/>
ls
gunzip unix_intro_data.tar.gz
ls
tar -xvf unix_intro_data.tar
ls
cd seqmonk_genomes/
ls
cd Saccharomyces\ cerevisiae/
ls
EF4/
cd EF4/
ls
grep rRNA Mito.dat
cp Mito.dat ~/compare
nano Mito.dat
mv Mito.dat Mitochondrion.txt
ls
cd
#Answer 4
cd compare
ls
cd FastQ_Data
ls
zcat lane8_DD_P4_TTAGGC_L008_R1.fastq.gz | wc –l
zcat -n *.fastq.gz | wc -l > totalcount
cd
exit
