#!/usr/bin/env bash

###CMD to run:
# bash mkCB.sh meta.csv

#set to parent dir
export CBDATAROOT='/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/cb_input'
outDir="../cb_output"
cmd1="cbBuild -r -i ./cellbrowser.conf -o ${outDir}"

#generate configs:

### setting:
name="canine-duodenum-CIE-healthy"
tags="10x"
organism="Canis lupus familiaris"
project="Canine duodenum scRNA-seq"
geo_series="GSE254005"

title="Single cell transcriptomic analysis of the canine duodenum in chronic inflammatory enteropathy and health "
submitter="Dylan Ammons"
version=1
submission_date="2024-05-13"

# abstract of paper or dataset summary text
showAbstract=TRUE
showMethods=TRUE

### code
#create md5sum hash for all data files
find -type f -name *tsv* -not -name "*ipynb*" -not -name "*-checkpoint*" -exec md5sum '{}' \; > md5sum.txt

while read line
do

    subName=$(echo "$line" | cut -f1 -d, --output-delimiter=' ')
    defaultClus=$(echo "$line" | cut -f2 -d, --output-delimiter=' ')
    enumClus=$(echo "$line" | cut -f3 -d, --output-delimiter=' ')
    pertyName=$(echo "$line" | cut -f4 -d, --output-delimiter=' ') 

    echo 'shortLabel = "'${pertyName}'"' > $subName/cellbrowser.conf

    echo 'enumFields = ["'${defaultClus}'"]' >> $subName/cellbrowser.conf ### modify this to make universal

    echo 'coords=[{"file":"umap.coords.tsv", "shortLabel":"UMAP"}, ]' >> $subName/cellbrowser.conf


    echo 'clusterField = "'${defaultClus}'"' >> $subName/cellbrowser.conf
    echo 'labelField = "'${defaultClus}'"' >> $subName/cellbrowser.conf


    echo 'organisms = ["'${organism}'"]' >> $subName/cellbrowser.conf
    echo 'projects = ["'${project}'"]' >> $subName/cellbrowser.conf

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
echo 'institution= "Colorado State University"' >> desc.conf
echo 'geo_series = "'${geo_series}'"' >> desc.conf
echo 'version = "'${version}'"' >> desc.conf
echo 'submission_date = "'${submission_date}'"' >> desc.conf

if [ $showAbstract == TRUE ]
then
    echo 'abstract = "abstract.html"' >> desc.conf
fi

if [ $showMethods == TRUE ]
then
    echo 'methods = "methods.html"' >> desc.conf
fi



#run cb
echo -e "\t$ ${cmd1}"
time eval $cmd1

cp $outDir/index.html $outDir/index.aspx
