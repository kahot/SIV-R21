#!/bin/bash

#This script for the data from Mac251 R21 project


#1. map to its own consensus
bwa index -p Run_5_21_Animal_3816 -a is Output/Consensus/Run_5_21_Animal_3816_Consensus.fasta

bwa mem -t 4 -k 15 -a Run_5_21_Animal_3816 Output/Run_5_21_Animal_3816/Run_5_21_Animal_3816_clean.q30.fq  > Output/Run_5_21_Animal_3816/Run_5_21_Animal_3816_BWAmapped2.sam

#5. convert sam to bam
samtools view -S -b Output/Run_5_21_Animal_3816/Run_5_21_Animal_3816_BWAmapped2.sam > Output/Run_5_21_Animal_3816/Run_5_21_Animal_3816_BWAmapped2.bam


#6. sort the bam file
samtools sort  Output/Run_5_21_Animal_3816/Run_5_21_Animal_3816_BWAmapped2.bam -o  Output/bam2/Run_5_21_Animal_3816_BWA_sort.bam

#7. index the bam file
samtools index  Output/bam2/Run_5_21_Animal_3816_BWA_sort.bam  Output/bam2/Run_5_21_Animal_3816_BWA_sort.bam.bai

rm -r Output/Run_5_21_Animal_3816/Run_5_21_Animal_3816_BWAmapped2.bam
rm -r Output/Run_5_21_Animal_3816/Run_5_21_Animal_3816_BWAmapped2.sam