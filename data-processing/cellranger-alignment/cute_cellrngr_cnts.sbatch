#!/usr/bin/env bash

#SBATCH --job-name=cellrngr_cnt
#SBATCH --nodes=1         # this script is designed to run on one node
#SBATCH --time=00:05:00

#SBATCH --partition=atesting  # modify this to reflect which queue you want to use. Either 'amilan' or 'atesting'
#SBATCH --qos=testing      # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'

#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=dyammons@colostate.edu ### change to your email ###

#SBATCH --output=cellrngr_cnt_%A_%a.log  #modify as desired - will output a log file where the "%A" inserts the job ID number and the %a

#SBATCH --array=0-7 #set this to 0-(# of samples - 1), so the example is for 8 samples -- if you are only running 1 sample, then you can set it to 0-0



##### USER SETTINGS #####

#specify path to reference
ref_pwd=/projects/dyammons@colostate.edu/references/canine/canine_ref_genome_cellranger_7_1_0_ROS_Cfam_1_0_110_base


#specify sample names to run
samples=$(ls -lh ../01_input/ | grep "^d" | awk '{print $9}')
declare -a StringArray=($samples)

#if you have extra dirs or only want to run select samples, then store the sample names in the StringArray variable
#declare -a StringArray=("sample1" "sample2" "sample3")

##### END SETTINGS #####


##### SET ENV #####

### Load cellranger
module purge
module load cellranger/7.1.0
cellranger --version


##### BEGIN CODE #####

### Excute cellranger count
sampleName=$(ls ../01_input/${StringArray[${SLURM_ARRAY_TASK_ID}]}/ | grep "fastq.gz" | head -n1 | awk -F "_S" '{print $1}')

cmd1="cellranger count --id=${StringArray[${SLURM_ARRAY_TASK_ID}]} \
                       --fastqs=../01_input/${StringArray[${SLURM_ARRAY_TASK_ID}]}/ \
                       --sample=${sampleName} \
                       --transcriptome=${ref_pwd} \
                       --expect-cells=5000 \
                       --dry" ### some add --include-introns=false \
echo ${cmd1}
time eval $cmd1

##### END CODE #####
