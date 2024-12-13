#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
time_stamp=$(date +"%Y-%m-%d %H:%M:%S")

sample_id="ERR1019038"
read_1_file_name="ERR1019038_1.fastq.gz"
read_2_file_name="ERR1019038_2.fastq.gz"

csv_file="fastq-list-$sample_id.csv"
touch $csv_file

printf "[$time_stamp]: "
printf "Writing to FASTQ list file... \n"
printf "RGID,RGSM,RGLB,Lane,Read1File,Read2File\n" >> $csv_file
printf "$sample_id,$sample_id,RGLB,1,$read_1_file_name,$read_2_file_name\n" >> $csv_file