#!/usr/bin/env bash

###CMD to run:
# bash mkCB.sh meta.csv

### MODIFY ###
#set to parent dir
export CBDATAROOT='/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/can_duod_update_2023_07_25'

#name of output dir
outDir="../can_duod_n3n4"

#command that will be used to generate the browser -- remove -r if only working with 1 dataset
cmd1="cbBuild -r -i ./cellbrowser.conf -o ${outDir}"

#CB settings:
name="duod_n3n4"
tags="10x"
organism="Canis lupus familiaris"
project="Canine duodenum atlas"
dotSize=3

title="A single-cell RNA sequencing atlas of canine duodenum in health and CIE"
submitter="Dylan Ammons"
version=1
submission_date="2023-07-25"

### END MODIFY ###

#generate config files:
while read line
do

    subName=( $(echo $line | cut -f1 -d',' --output-delimiter=' ') )
    defaultClus=( $(echo $line | cut -f2 -d',' --output-delimiter=' ') )
    pertyName=( $(echo $line | cut -f3 -d',' --output-delimiter=' ') )    

    echo 'shortLabel= '${pertyName} > $subName/cellbrowser.conf

    echo 'enumFields = ["'${defaultClus}'"]' >> $subName/cellbrowser.conf ### modify this to make universal

    echo 'coords=[{"file":"umap.coords.tsv", "shortLabel":"UMAP"}, ]' >> $subName/cellbrowser.conf


    echo 'clusterField = "'${defaultClus}'"' >> $subName/cellbrowser.conf
    echo 'labelField = "'${defaultClus}'"' >> $subName/cellbrowser.conf


    echo 'organisms = ["'${organism}'"]' >> $subName/cellbrowser.conf
    echo 'projects = ["'${project}'"]' >> $subName/cellbrowser.conf
    echo "radius = ${dotSize}" >> $subName/cellbrowser.conf

    echo 'quickGenesFile = "quickGenes.csv"' >> $subName/cellbrowser.conf

    #immutables:
    echo 'exprMatrix="exprMatrix.tsv.gz"' >> $subName/cellbrowser.conf
    echo 'geneIdType="raw"' >> $subName/cellbrowser.conf
    echo 'meta="meta.tsv"' >> $subName/cellbrowser.conf
    echo 'markers=[{"file":"markers.tsv", "shortLabel":"Cluster-specific markers"}]' >> $subName/cellbrowser.conf
    echo 'unit = "log of read count/UMI"' >> $subName/cellbrowser.conf
    echo 'matrixType="auto"' >> $subName/cellbrowser.conf
    
done < $1

#make main conf file
echo 'name = "'${name}'"' > cellbrowser.conf
echo 'shortLabel = "'${project}'"' >> cellbrowser.conf
echo 'tags = ["'${tag}'"]' >> cellbrowser.conf

#make main desc.conf
echo 'title = "'${title}'"' > desc.conf
echo 'submitter = "'${submitter}'"' >> desc.conf
echo 'version = "'${version}'"' >> desc.conf
echo 'submission_date = "'${submission_date}'"' >> desc.conf

#run cb
echo -e "\t$ ${cmd1}"
time eval $cmd1

cp $outDir/index.html $outDir/index.aspx
