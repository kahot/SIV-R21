#!/bin/bash

#This script for the data from Mac251 R21 project

mkdir Output/Run_4_15_Animal_3516

#0 adapter trimming
bbduk.sh in1=/Volumes/Kaho_Data/SIV_DATA/mac251_No_PID_Run_1-118926115/FASTQ_Generation_2019-02-25_18_22_24Z-163551411/15_L001-ds.6e407d39e55e43d9830d6d6114b79deb/15_S15_L001_R1_001.fastq.gz in2=/Volumes/Kaho_Data/SIV_DATA/mac251_No_PID_Run_1-118926115/FASTQ_Generation_2019-02-25_18_22_24Z-163551411/15_L001-ds.6e407d39e55e43d9830d6d6114b79deb/15_S15_L001_R2_001.fastq.gz  out=Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=Output/Run_4_15_Animal_3516/stats30_0.txt

#1 Trim reads at both ends at Score<30
bbduk.sh in=Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_adp.trimmed.fastq out=Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_trimmed.q30.fastq qtrim=rl trimq=30 stats=Output/Run_4_15_Animal_3516/stats30_1.txt

#2. Kmer filtering
bbduk.sh in=Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_trimmed.q30.fastq out=Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_unmatched.q30.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=Output/Run_4_15_Animal_3516/stats30_2.txt

#3.Remove reads with ave score <30 (Phred score =20 99% accuracy 1% chance of mistake)
bbduk.sh in=Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_unmatched.q30.fq out=Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_clean.q30.fq maq=30 stats=Output/Run_4_15_Animal_3516/stats30_3.txt

#4. Align the file using bwa to the reference 
#bwa index -p SIV -a is Data/SIVMM239_ENV.fasta

bwa mem -t 4 -k 15 -a SIV Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_clean.q30.fq  > Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_BWAmapped.sam

#5. convert sam to bam
samtools view -S -b Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_BWAmapped.sam > Output/bam/Run_4_15_Animal_3516_BWAmapped.bam


rm -r Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_BWAmapped.sam
rm -r Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_unmatched.q30.fq
rm -r Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_trimmed.q30.fastq
rm -r Output/Run_4_15_Animal_3516/Run_4_15_Animal_3516_adp.trimmed.fastq

