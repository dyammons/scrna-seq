#!/usr/bin/env bash

### Script to move output count matricies from scratch to pl
### Input: absolute path to folder containing all of the cellranger for outputs
### Output: copies required files into subsdiretories named based on sample
###
### Usage: mv-matrix.sh path/to/CellRanger/outputs
###
### Known bugs: if there are extra directories in $1 then the script will likely fail
###
### Written: unknown - DA
### Upated: 01/31/2024


### Begin code


echo "INFO: Retrieving paths to Cell Ranger output files."
dirs=$( ls -l $1 | grep "^d" | awk '{print $9}' )
declare -a StringArray=($dirs)

echo "INFO: The following samples will be transerfered: ${StringArray[@]} to $PWD"

#read in the list and copy the data required to run ESAT
for val in "${StringArray[@]}"; do
  folder="./$val/"
  mkdir $folder
  
  filez="${1}/$val/outs/filtered_feature_bc_matrix/*"
  cp $filez $folder

done

echo "INFO: End script. ${#StringArray[@]} samples were copied."


### End code
