#!/bin/bash

#This script for the data from Mac251 R21 project


#1. map to its own consensus
bwa index -p Run_5_24_Animal_4016 -a is Output/Consensus/Run_5_24_Animal_4016_Consensus.fasta

bwa mem -t 4 -k 15 -a Run_5_24_Animal_4016 Output/Run_5_24_Animal_4016/Run_5_24_Animal_4016_clean.q30.fq  > Output/Run_5_24_Animal_4016/Run_5_24_Animal_4016_BWAmapped2.sam

#5. convert sam to bam
samtools view -S -b Output/Run_5_24_Animal_4016/Run_5_24_Animal_4016_BWAmapped2.sam > Output/Run_5_24_Animal_4016/Run_5_24_Animal_4016_BWAmapped2.bam


#6. sort the bam file
samtools sort  Output/Run_5_24_Animal_4016/Run_5_24_Animal_4016_BWAmapped2.bam -o  Output/bam2/Run_5_24_Animal_4016_BWA_sort.bam

#7. index the bam file
samtools index  Output/bam2/Run_5_24_Animal_4016_BWA_sort.bam  Output/bam2/Run_5_24_Animal_4016_BWA_sort.bam.bai

rm -r Output/Run_5_24_Animal_4016/Run_5_24_Animal_4016_BWAmapped2.bam
rm -r Output/Run_5_24_Animal_4016/Run_5_24_Animal_4016_BWAmapped2.sam
