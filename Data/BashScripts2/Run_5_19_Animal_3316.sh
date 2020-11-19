#!/bin/bash

#This script for the data from Mac251 R21 project


#1. map to its own consensus
bwa index -p Run_5_19_Animal_3316 -a is Output/Consensus/Run_5_19_Animal_3316_Consensus.fasta

bwa mem -t 4 -k 15 -a Run_5_19_Animal_3316 Output/Run_5_19_Animal_3316/Run_5_19_Animal_3316_clean.q30.fq  > Output/Run_5_19_Animal_3316/Run_5_19_Animal_3316_BWAmapped2.sam

#5. convert sam to bam
samtools view -S -b Output/Run_5_19_Animal_3316/Run_5_19_Animal_3316_BWAmapped2.sam > Output/Run_5_19_Animal_3316/Run_5_19_Animal_3316_BWAmapped2.bam


#6. sort the bam file
samtools sort  Output/Run_5_19_Animal_3316/Run_5_19_Animal_3316_BWAmapped2.bam -o  Output/bam2/Run_5_19_Animal_3316_BWA_sort.bam

#7. index the bam file
samtools index  Output/bam2/Run_5_19_Animal_3316_BWA_sort.bam  Output/bam2/Run_5_19_Animal_3316_BWA_sort.bam.bai

rm -r Output/Run_5_19_Animal_3316/Run_5_19_Animal_3316_BWAmapped2.bam
rm -r Output/Run_5_19_Animal_3316/Run_5_19_Animal_3316_BWAmapped2.sam
