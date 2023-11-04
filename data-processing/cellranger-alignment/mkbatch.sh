#!/usr/bin/env bash

########################################################################
#  Function: automate job script creation for multiple scRNA samples!  #
#                                                                      #
#  Useage: bash mkbatch.sh                                             #
#                                                                      #
#  Assumptions:                                                        #
#    - cellranger is to be loaded with module load cellranger          #
#    - run in a directory in which input relative path is ../01_input  #
#                                                                      #
#  Requirments: update the user preferences below                      #
#                                                                      #
#  Created: October 2023                                               #
#  Updated: November 22, 2023 - by DA                                  #
########################################################################


##### set user preferences #####
numNode=1
numTasks=24
email="dyammons@colostate.edu"
runTime="06:00:00"
refpwd="/projects/$USER/references/canine/k9_ref_genome"
partition="amilan"

##### make the files #####
#retrieve samples based on what is present in the 01_input folder
dirs=$(ls -lh ../01_input/ | grep "^d" | awk '{print $9}')
declare -a StringArray=($dirs)

> jobList.txt
for val in "${StringArray[@]}"; do

	#create sbatch file
	> cute_cnts_$val.sbatch
	echo "#!/usr/bin/env bash" >> cute_cnts_$val.sbatch
	echo "" >> cute_cnts_$val.sbatch
	echo "#SBATCH --job-name=cnt_$val" >> cute_cnts_$val.sbatch
 
	echo "" >> cute_cnts_$val.sbatch
	echo "#SBATCH --nodes=$numNode" >> cute_cnts_$val.sbatch                       
	echo "#SBATCH --ntasks=$numTasks" >> cute_cnts_$val.sbatch
	echo "#SBATCH --time=$runTime" >> cute_cnts_$val.sbatch
 
	echo "" >> cute_cnts_$val.sbatch
	echo "#SBATCH --partition=$partition" >> cute_cnts_$val.sbatch
	echo "#SBATCH --qos=normal" >> cute_cnts_$val.sbatch
 
	echo "" >> cute_cnts_$val.sbatch
	echo "#SBATCH --mail-type=END" >> cute_cnts_$val.sbatch
	echo "#SBATCH --mail-user=$email" >> cute_cnts_$val.sbatch
	echo "#SBATCH --output=cellRngr_$val-%j.log" >> cute_cnts_$val.sbatch
 
	echo "" >> cute_cnts_$val.sbatch
	echo "module purge" >> cute_cnts_$val.sbatch
	echo "module load cellranger" >> cute_cnts_$val.sbatch
 
	echo "" >> cute_cnts_$val.sbatch
	echo "##### Call bash script #####" >> cute_cnts_$val.sbatch
	echo "" >> cute_cnts_$val.sbatch
	echo "bash cnts_$val.sh" >> cute_cnts_$val.sbatch
	echo "" >> cute_cnts_$val.sbatch
	
	#create cellranger counts script
	> cnts_$val.sh
	echo "#!/usr/bin/env bash" >> cnts_$val.sh
	echo "" >> cnts_$val.sh
	
 
	echo "" >> cnts_$val.sh
	echo "##### Call cellranger count #####" >> cnts_$val.sh
	echo "cellranger count --id=run_count_$val \\" >> cnts_$val.sh
	echo "--fastqs=../01_input/$val/ \\" >> cnts_$val.sh
	samName=$(ls ../01_input/$val/ | grep "fastq.gz" | head -n1 | awk -F "_S" '{print $1}') #this needs a sed cmd probs
	echo "--sample=$samName \\" >> cnts_$val.sh
	echo "--transcriptome=$refpwd \\" >> cnts_$val.sh
	echo "--expect-cells=5000" >> cnts_$val.sh
	
	#create list of cmds    
	echo "sbatch cute_cnts_$val.sbatch" >> jobList.txt
done 
