#!/usr/bin/env bash

#SBATCH --job-name=cellranger_mkref

#SBATCH --nodes=1         # this script is designed to run on one node
#SBATCH --ntasks=8
#SBATCH --time=01:00:00   # set time; default = 4 hours; 0 = no restriction

#SBATCH --partition=amilan  # modify this to reflect which queue you want to use. Either 'shas' or 'shas-testing'
#SBATCH --qos=normal      # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'

#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=dyammons@colostate.edu ### change to your email address ###

#SBATCH --output=cellrngr_mkref_%j.log  # this will capture all output in a logfile with %j as the job #

######### Instructions ###########

#make reference genome
bash mkref.sh
