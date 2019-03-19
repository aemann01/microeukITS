#!/bin/bash
# Script for demultiplexing ITS data, gather some simple summary statistics on the raw reads
# Use: ./demultiplex.sh
###################
#Demultiplexing
###################
#after logging onto deep thought
cd /home/user/mann
bcl2fastq -R /media/Synology_Data/Seq_Run_Data/LMAMR51_Mar_2019/ -o fastq --sample-sheet sampleSheet.csv --use-bases-mask Y*,I*,I*,Y*

#get file name and raw count
cd fastq
ls *gz > file
ls *gz | while read line; do zcat $line | wc -l ; done >> raw
awk '{print $1/4}' raw > reads
paste file reads > readcount.txt
rm file raw reads
#get summary stats of raw reads
average=$(awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' readcount.txt)
median=$(awk '{print $2}' readcount.txt | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else pri
nt a[x-1]; }')
max=$(awk '{print $2}' readcount.txt | sort -n | tail -n1)
min=$(awk '{print $2}' readcount.txt | sort -n | head -n1)
echo "Average:" $average >> readcount.stats.txt
echo "Median:" $median >> readcount.stats.txt
echo "Minimum read count:" $min >> readcount.stats.txt
echo "Maximum read count:" $max >> readcount.stats.txt

