#!/usr/bin/env bash

#SBATCH --job-name=kal_cnt
#SBATCH --ntasks=8       # modify this number to reflect how many cores you want to use (up to 64)
#SBATCH --nodes=1         # this script is designed to run on one node
#SBATCH --time=00:40:00   # set time; default = 4 hours

#SBATCH --partition=amilan  # modify this to reflect which queue you want to use. Either 'shas' or 'shas-testing'
#SBATCH --qos=normal      # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'

#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=dyammons@colostate.edu ### change to your email ###

#SBATCH --output=kal_cnt_%A_%a.log  #modify as desired - will output a log file where the "%A" inserts the job ID number and the %a

#SBATCH --array=0-7 #set this to 0-(# of samples - 1), so the example is for 8 samples -- if you are only running 1 sample, then you can set it to 0-0 (although you may want to consider not using an array)



####### BEGIN USER PREFERENCES #######

#pwd to kallisto index
pwd_idx="/projects/$USER/references/canine/canFam3.1.kb-python.idx"

#pwd to transcript to gene file
pwd_tr2g="/projects/$USER/references/canine/transcripts_to_genes.txt"

#chemistry used in the assay
chem="10xv3"

#set the tmp dir -- it should not already exist
tmp_dir="/scratch/alpine/$USER/tmp/kallisto/"

#set the output directory
DATE="11-07-2023"
out_dir="../03_output/$DATE-kb-count_output/"

#set sample names in a sting array
#samples=$(ls -lh ../01_input/ | grep "^d" | awk '{print $9}')
#declare -a StringArray=($samples)

#if you have extra dirs or only want to run select samples, then store the sample names in the StringArray variable
declare -a StringArray=("duod_norm_2" "duod_norm_3")

#Ensure the node is clear of pre-loaded software & load conda env with kb-python
module purge
module load anaconda
conda activate kb-tools #change the env to match the env with the required software installed

####### END USER PREFERENCES #######



####### BEGIN CODE #######

### Extract sample variables
sampleName=$(ls ../01_input/${StringArray[${SLURM_ARRAY_TASK_ID}]}/ | grep "fastq.gz" | head -n1 | awk -F "_S" '{print $1}')
R1=$(ls ../01_input/${StringArray[${SLURM_ARRAY_TASK_ID}]}/ | grep "_R1_")
R2=$(ls ../01_input/${StringArray[${SLURM_ARRAY_TASK_ID}]}/ | grep "_R2_")

### Update user
echo -e ">>> INPUT: This script will process the sample files into the following names: "
echo -e "\tSAMPLE1\tSAMPLE2\tNAMES"
echo -e "\t$R1\t$R2\t$sampleName\n"

### Excute kb count
cmd1="kb count -i $pwd_idx \
-g $pwd_tr2g \
-x $chem \
-t $SLURM_NTASKS \
--tmp $tmp_dir \
-o $out_dir/${StringArray[${SLURM_ARRAY_TASK_ID}]}/ \
--h5ad \
../01_input/${StringArray[${SLURM_ARRAY_TASK_ID}]}/$R1 ../01_input/${StringArray[${SLURM_ARRAY_TASK_ID}]}/$R2"

echo $cmd1
echo -e "\t$ $cmd1"
time eval $cmd1

####### END CODE #######
