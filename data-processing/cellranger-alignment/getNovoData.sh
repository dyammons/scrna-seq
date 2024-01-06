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


echo "INFO: Processing url list"
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

echo "INFO: Concatinating MD5.txt files"
md5s=$(sort -u getMd5.tmp)
for mSum in $md5s
do
    wget -q $mSum
    cat MD5.txt >> md5Chk.tmp
    rm MD5.txt
done

echo "INFO: Downloading data to ${PWD}"
#use parallelization to get the raw data - uses 8 workers
echo $url_list | xargs -n 1 -P 8 wget -q -nc

echo "INFO: Ensuring file integrity"
md5sum -c md5Chk.tmp > md5_verify.log

echo "INFO: Data transfer complete. Cleaning up directory."
cp md5Chk.tmp MD5.txt
rm *.tmp

echo "INFO: Script complete. Please review the md5_verify.log file to ensure files transfered properly."
