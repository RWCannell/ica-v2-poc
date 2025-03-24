#!/bin/bash

file_path="/dataA/1000G/SRR1295540_1.fastq.gz"

compressed_file_size=$(gzip -t $file_path)
echo $compressed_file_size