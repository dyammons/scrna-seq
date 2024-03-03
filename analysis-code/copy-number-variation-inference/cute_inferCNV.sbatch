#!/usr/bin/env bash
#SBATCH --job-name=infCNV
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:60:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mail-type=END
#SBATCH --mail-user=dyammons@colostate.edu
#SBATCH --output=%x-%A_%a.log
#SBATCH --array=0-7 #set this to 0-(# of samples - 1), so the example is for 8 samples -- if you are only running 1 sample, then you can set it to 0-0

#load conda env
source /curc/sw/anaconda3/latest
conda activate r_env

#run the Rscript of interest
Rscript seurat_inferCNVindv.R -t $SLURM_NTASKS -i $SLURM_ARRAY_TASK_ID
