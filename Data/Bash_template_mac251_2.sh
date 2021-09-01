#!/bin/bash

#This script for the data from Mac251 R21 project


#1. map to its own consensus
bwa index -p SAMPLE -a is Output/Consensus/SAMPLE_Consensus.fasta

bwa mem -t 4 -k 15 -a SAMPLE Output/SAMPLE/SAMPLE_clean.q30.fq  > Output/SAMPLE/SAMPLE_BWAmapped2.sam

#5. convert sam to bam
samtools view -S -b Output/SAMPLE/SAMPLE_BWAmapped2.sam > Output/SAMPLE/SAMPLE_BWAmapped2.bam


#6. sort the bam file
samtools sort  Output/SAMPLE/SAMPLE_BWAmapped2.bam -o  Output/bam2/SAMPLE_BWA_sort.bam

#7. index the bam file
samtools index  Output/bam2/SAMPLE_BWA_sort.bam  Output/bam2/SAMPLE_BWA_sort.bam.bai

rm -r Output/SAMPLE/SAMPLE_BWAmapped2.bam
rm -r Output/SAMPLE/SAMPLE_BWAmapped2.sam
