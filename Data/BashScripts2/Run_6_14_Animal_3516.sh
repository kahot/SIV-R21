#!/bin/bash

#This script for the data from Mac251 R21 project


#1. map to its own consensus
bwa index -p Run_6_14_Animal_3516 -a is Output/Consensus/Run_6_14_Animal_3516_Consensus.fasta

bwa mem -t 4 -k 15 -a Run_6_14_Animal_3516 Output/Run_6_14_Animal_3516/Run_6_14_Animal_3516_clean.q30.fq  > Output/Run_6_14_Animal_3516/Run_6_14_Animal_3516_BWAmapped2.sam

#5. convert sam to bam
samtools view -S -b Output/Run_6_14_Animal_3516/Run_6_14_Animal_3516_BWAmapped2.sam > Output/Run_6_14_Animal_3516/Run_6_14_Animal_3516_BWAmapped2.bam


#6. sort the bam file
samtools sort  Output/Run_6_14_Animal_3516/Run_6_14_Animal_3516_BWAmapped2.bam -o  Output/bam2/Run_6_14_Animal_3516_BWA_sort.bam

#7. index the bam file
samtools index  Output/bam2/Run_6_14_Animal_3516_BWA_sort.bam  Output/bam2/Run_6_14_Animal_3516_BWA_sort.bam.bai

rm -r Output/Run_6_14_Animal_3516/Run_6_14_Animal_3516_BWAmapped2.bam
rm -r Output/Run_6_14_Animal_3516/Run_6_14_Animal_3516_BWAmapped2.sam
