#!/bin/bash

file_path="Documents/test_data/fastq/1_control_trnL_2019_minq7.fastq"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"

uploaded_file_response=$(icav2 projectdata upload $file_path --project-id $project_id)

file_id=$(echo $uploaded_file_response | jq -r '.id')
echo $file_id