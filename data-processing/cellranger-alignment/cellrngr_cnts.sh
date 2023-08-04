#!/usr/bin/env bash

pthread=$1
echo $pthead

export PATH=/projects/$USER/software/cellranger-6.1.2:$PATH ### Ensure path is correct

#for the command below, add your sample name to the "id", "fastqs", and "sample" options And ensure all paths are correct!
#if you get an error that fastqs aren't in directory make sure the path is correct and the sample name matches
cellranger count --id=run_count_<sample name> \
                 --fastqs=/scratch/summit/$USER/project_01/01_input/<sample name> \
                 --sample=<sample name> \
                 --transcriptome=/projects/$USER/references/K9_ref_genome \
                 --localcores=$pthread \
                 --expect-cells=5000 ### MODIFY as needed
