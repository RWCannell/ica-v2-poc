#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
pipeline_code="basic_pipeline"
pipeline_id="bfecca03-6443-45bd-b313-e4f555cd0748"
user_reference="regan_test_analysis_03"
storage_size="Small"
file_name="NZ_GG704940.fa"

projectanalyses_list_response=$(icav2 projectanalyses list --project-id $project_id)
analysis_id=$(echo $projectanalyses_response | jq -r ".items[] | select(.userReference == \"$user_reference\").id")

analysis_response=$(icav2 projectanalyses input $analysis_id)

analysis_data_id=$(echo $analysis_response | jq -r ".items[].analysisData[] | select(.name == \"$file_name\").dataId")
analysis_data_code=$(echo $analysis_response | jq -r ".items[].code")

file_ref="$analysis_data_code:$analysis_data_id"

icav2 projectpipelines start nextflow $pipeline_id \
--user-reference $user_reference \
--project-id $project_id \
--storage-size $storage_size \
--input $file_ref
