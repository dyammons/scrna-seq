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

#get and cat md5sum files for each sample
md5s=$(sort -u getMd5.tmp)
for mSum in $md5s
do
    wget -q $mSum
    cat MD5.txt >> md5Chk.tmp
    rm MD5.txt
done

#use parallelization to get the raw data - uses 8 workers
echo $url_list | xargs -n 1 -P 8 wget -q -nc

#check the md5sums for each file and record results in "md5_verify.log"
md5sum -c md5Chk.tmp > md5_verify.log

#clean up any .tmp files generated
cp md5Chk.tmp MD5.txt
rm *.tmp
