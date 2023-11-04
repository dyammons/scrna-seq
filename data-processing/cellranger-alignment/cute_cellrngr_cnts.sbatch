#!/usr/bin/env bash

#SBATCH --job-name=cellrngr_cnt
#SBATCH --ntasks=24       # modify this number to reflect how many cores you want to use (up to 64)
#SBATCH --nodes=1         # this script is designed to run on one node
#SBATCH --time=06:00:00   # set time; default = 4 hours

#SBATCH --partition=amilan  # modify this to reflect which queue you want to use. Either 'shas' or 'shas-testing'
#SBATCH --qos=normal      # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'

#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=dyammons@colostate.edu ### change to your email ###

#SBATCH --output=cellrngr_cnt_%j.log  #modify as desired - will output a log file where the "%j" inserts the job ID number
#SBATCH --array=0-7

##### Load cellranger #####
module purge
module load cellranger


##### Load in sample names to run #####
samples=$(ls -lh ../01_input/ | grep "^d" | awk '{print $9}')
declare -a StringArray=($samples)

#if you have extra dirs or only want to run select samples, then store the sample names in the StringArray variable
#declare -a StringArray=("sample1" "sample2" "sample3")

##### Call bash script #####
bash rngr_cnts.sh ${StringArray[${SLURM_ARRAY_TASK_ID}]} $(ls ../01_input/${StringArray[${SLURM_ARRAY_TASK_ID}]}/ | grep "fastq.gz" | head -n1 | awk -F "_S" '{print $1}')
