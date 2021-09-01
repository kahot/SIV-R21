#!/bin/bash

#This script for the data from Mac251 R21 project

mkdir Output/SAMPLE

#0 adapter trimming
bbduk.sh in1=FASTQFile1 in2=FASTQFile2  out=Output/SAMPLE/SAMPLE_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=Output/SAMPLE/stats30_0.txt

#1 Trim reads at both ends at Score<30
bbduk.sh in=Output/SAMPLE/SAMPLE_adp.trimmed.fastq out=Output/SAMPLE/SAMPLE_trimmed.q30.fastq qtrim=rl trimq=30 stats=Output/SAMPLE/stats30_1.txt

#2. Kmer filtering
bbduk.sh in=Output/SAMPLE/SAMPLE_trimmed.q30.fastq out=Output/SAMPLE/SAMPLE_unmatched.q30.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=Output/SAMPLE/stats30_2.txt

#3.Remove reads with ave score <30 (Phred score =20 99% accuracy 1% chance of mistake)
bbduk.sh in=Output/SAMPLE/SAMPLE_unmatched.q30.fq out=Output/SAMPLE/SAMPLE_clean.q30.fq maq=30 stats=Output/SAMPLE/stats30_3.txt

#4. Align the file using bwa to the reference 
#bwa index -p SIV -a is Data/SIVMM239_ENV.fasta

bwa mem -t 4 -k 15 -a SIV Output/SAMPLE/SAMPLE_clean.q30.fq  > Output/SAMPLE/SAMPLE_BWAmapped.sam

#5. convert sam to bam
samtools view -S -b Output/SAMPLE/SAMPLE_BWAmapped.sam > Output/bam/SAMPLE_BWAmapped.bam


rm -r Output/SAMPLE/SAMPLE_BWAmapped.sam
rm -r Output/SAMPLE/SAMPLE_unmatched.q30.fq
rm -r Output/SAMPLE/SAMPLE_trimmed.q30.fastq
rm -r Output/SAMPLE/SAMPLE_adp.trimmed.fastq

