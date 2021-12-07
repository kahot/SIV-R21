# SIV-R21

## Analysis of SIV R21 SIV/Mtb Co-infection Study


* Data directory contains sample information. SampleSheet_Mac251.csv has the most comprehensive sample information.
* Output directory contains all output files, except the bam files (due to the memory limit)
* Rscripts directory contains all scripts for the analysis. 


## Analysis steps

### 1. Create bash files to process FASTQ files with CreateBashFiles.R 

* Output bash files are stored in Data/BashScripts1/
* Run all bash files to trim/filter/map to create bam files
* After creating a consensus sequence for each file using Geneious (or another tool of your choice), run the files in Data/BashScripts2/ 
	
```bash
#Run all bash scripts in a directory
for f in Data/BashScripts1/*.sh; do
	bash "$f" -H   || break 
done
```
* Need bam and bai (index) files for the next step  
    
### 2. After mapping is done, start the analysis in R by following the file numbers

* Scripts 1 to 3 create frequency tables for each bam file.
* Script 4 analyze the stock and control files 
* Script 5 assess transmitted/founder variants 
* Script 6 assess SIV diversity across tissues/time/cohorts and compare with RNA copy nunbers and CD4/CD8 frequency
* Script 7 analyze the indels 
* Script 8 assess the known immune escape mutations 
* Script 9 assess an interesting mutation at AA112
* Script 10 analyze the effects of genetic drift between tissues
* Script 11 conducts glmm to understand the effects of different factors on SIV diversity
    
