#!/usr/bin/env bash

#this command will be populated from the SLURM directives and the .sbatch script
#make sure the fastqs are in subdirectories within ../01_input and everything shoudl go smoothly

cellranger count --id=${1} \
                 --fastqs=../01_input/${1}/ \
                 --sample=${2} \
                 --transcriptome=/projects/$USER/references/K9_ref_genome \
                 --expect-cells=5000 ### MODIFY as needed
