#!/bin/bash

#This script for the data from Mac251 R21 project


#1. map to its own consensus
bwa index -p Run_4_10_Animal_20615 -a is Output/Consensus/Run_4_10_Animal_20615_Consensus.fasta

bwa mem -t 4 -k 15 -a Run_4_10_Animal_20615 Output/Run_4_10_Animal_20615/Run_4_10_Animal_20615_clean.q30.fq  > Output/Run_4_10_Animal_20615/Run_4_10_Animal_20615_BWAmapped2.sam

#5. convert sam to bam
samtools view -S -b Output/Run_4_10_Animal_20615/Run_4_10_Animal_20615_BWAmapped2.sam > Output/Run_4_10_Animal_20615/Run_4_10_Animal_20615_BWAmapped2.bam


#6. sort the bam file
samtools sort  Output/Run_4_10_Animal_20615/Run_4_10_Animal_20615_BWAmapped2.bam -o  Output/bam2/Run_4_10_Animal_20615_BWA_sort.bam

#7. index the bam file
samtools index  Output/bam2/Run_4_10_Animal_20615_BWA_sort.bam  Output/bam2/Run_4_10_Animal_20615_BWA_sort.bam.bai

rm -r Output/Run_4_10_Animal_20615/Run_4_10_Animal_20615_BWAmapped2.bam
rm -r Output/Run_4_10_Animal_20615/Run_4_10_Animal_20615_BWAmapped2.sam
