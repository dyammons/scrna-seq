#!/usr/bin/env bash

#SBATCH --job-name=kb-cnt

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=04:00:80

#SBATCH --partition=amilan
#SBATCH --qos=normal

#SBATCH --mail-type=END
#SBATCH --mail-user=dyammons@colostate.edu
#SBATCH --output=kb-cnt-%j.log

##### Notes #####
#This sbatch file is designed to run on sample though kb count using kb-python.
#Pending successful run (started 11/06/2023 @ 6:40pm) this script will be coverted into an array script


##### Call bash script #####

#ensure the node is clear of pre-loaded software & load conda env with kb-python
module purge
module load anaconda
conda activate dsci512

#excute kb count
kb count -i /projects/$USER/references/canine/canFam3.1.kb-python.idx -g /projects/dyammons@colostate.edu/references/canine/transcripts_to_genes.txt -x 10xv3 \
--h5ad -t $SLURM_NTASKS --tmp /scratch/alpine/dyammons@colostate.edu/tmp/tmp/ \
-o ../03_output/2023-11-06_kb-count_output/duod_norm_1/ \
/scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/01_input/duod_norm_1/duod_norm_1_CKDL220032323-1A_HMLG2DSX5_S6_L002_R1_001.fastq.gz /scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/01_input/duod_norm_1/duod_norm_1_CKDL220032323-1A_HMLG2DSX5_S6_L002_R2_001.fastq.gz
