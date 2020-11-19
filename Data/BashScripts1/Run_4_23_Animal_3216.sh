#!/bin/bash

#This script for the data from Mac251 R21 project

mkdir Output/Run_4_23_Animal_3216

#0 adapter trimming
bbduk.sh in1=/Volumes/Kaho_Data/SIV_DATA/mac251_No_PID_Run_1-118926115/FASTQ_Generation_2019-02-25_18_22_24Z-163551411/23_L001-ds.5b9c411bac2544d2ae5b6957f54c6dd1/23_S23_L001_R1_001.fastq.gz in2=/Volumes/Kaho_Data/SIV_DATA/mac251_No_PID_Run_1-118926115/FASTQ_Generation_2019-02-25_18_22_24Z-163551411/23_L001-ds.5b9c411bac2544d2ae5b6957f54c6dd1/23_S23_L001_R2_001.fastq.gz  out=Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=Output/Run_4_23_Animal_3216/stats30_0.txt

#1 Trim reads at both ends at Score<30
bbduk.sh in=Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_adp.trimmed.fastq out=Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_trimmed.q30.fastq qtrim=rl trimq=30 stats=Output/Run_4_23_Animal_3216/stats30_1.txt

#2. Kmer filtering
bbduk.sh in=Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_trimmed.q30.fastq out=Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_unmatched.q30.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=Output/Run_4_23_Animal_3216/stats30_2.txt

#3.Remove reads with ave score <30 (Phred score =20 99% accuracy 1% chance of mistake)
bbduk.sh in=Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_unmatched.q30.fq out=Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_clean.q30.fq maq=30 stats=Output/Run_4_23_Animal_3216/stats30_3.txt

#4. Align the file using bwa to the reference 
#bwa index -p SIV -a is Data/SIVMM239_ENV.fasta

bwa mem -t 4 -k 15 -a SIV Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_clean.q30.fq  > Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_BWAmapped.sam

#5. convert sam to bam
samtools view -S -b Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_BWAmapped.sam > Output/bam/Run_4_23_Animal_3216_BWAmapped.bam


rm -r Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_BWAmapped.sam
rm -r Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_unmatched.q30.fq
rm -r Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_trimmed.q30.fastq
rm -r Output/Run_4_23_Animal_3216/Run_4_23_Animal_3216_adp.trimmed.fastq

