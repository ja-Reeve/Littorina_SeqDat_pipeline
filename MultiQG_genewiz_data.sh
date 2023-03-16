#!/bin/bash

### Generating read quality reports for WGS  Littorina saxatilis samples from
### arcoss the North Atlantic. Reads were sequenced at Azenta's GeneWiz 
### facility in Leipzig at 10x target coverage.  Reports are generated for each 
### sample using MultiQC.
### James Reeve - Göteborgs Universitet
### 28/10/2022

### Job parameters
#SBATCH -A snic2022-5-266
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 20:00:00
#SBATCH -J BS_MultiQC_GeneWiz
#SBATCH --output logfiles/BS_multiQC_GeneWiz.log 

### Required software
module load bioinfo-tools
module load FastQC/0.11.9
module load MultiQC/1.10

### Command
DIR=/proj/snic2020-6-155/James

fastqc --dir $TMPDIR -t 20 $DIR/fastq/Broad_scale/GeneWiz/*.fastq.gz --outdir $DIR/fastqc/Broad_scale/GeneWiz
multiqc $DIR/fastqc/Broad_scale/GeneWiz -o $DIR/MultiQC
