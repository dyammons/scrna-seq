#!/usr/bin/env bash

##########################################################################################
#  Function: this short script will collect data from Novogene's server,                 #
#            check md5sum to ensure integrity, and put each sample in its own directory  #
#                                                                                        #
#  Useage: bash getNovoData.sh url.txt                                                   #
#                                                                                        #
#  Input: the only input is url.txt                                                      #
#         - this file should be a list of urls with each url on a new line.              #
#                                                                                        #
#  Created: October 2023                                                                 #
#  Updated: August 4, 2023 - by DA                                                       #
##########################################################################################


#read in the urls from Novogene and convert urls to lists
url_list=()
while read line
do
    echo $line | cut -d "/" -f5 >> samples.tmp
    if echo $line | grep -q "MD5.txt"
    then
        echo $line >> getMd5.tmp
    else
        url_list+=" \"$line\""
    fi
done < $1


#make output dir(s) for each sample being retrieved
sample_list=$(sort -u samples.tmp)
mkdir $sample_list

#get and move md5sum files to appropriate dir
#this was a little tricky b/c all of the checksum files have the same name (MDS.txt), to this code gets that sorted
md5s=$(sort -u getMd5.tmp)

for mSum in $md5s
do
    wget -q $mSum
    getName=$(echo $mSum | cut -d "/" -f5)
    cp MD5.txt ./$getName/
    sed -i "s/  /  .\/$getName\//g" MD5.txt
    cat MD5.txt >> md5Chk.tmp
    rm MD5.txt
done

#use parallelization to get the raw data - uses 8 workers
echo $url_list | xargs -n 1 -P 8 wget -q -nc

#move the sample(s) to their cooresponding dir(s)
for sample in $sample_list
do
   mv $sample* ./$sample/ 2>/dev/null
done

#check the md5sums for each file and record results in "md5_verify.log"
md5sum -c md5Chk.tmp > md5_verify.log

#clean up any .tmp files generated
rm *.tmp
