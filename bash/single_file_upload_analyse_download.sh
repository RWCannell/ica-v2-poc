#!/bin/bash

file_path="$HOME/Documents/ica_data_uploads/fasta/Yokenella_regensburgei_ATCC_43003_uid65133/NZ_AGCL00000000.scaffold.fa"
file_name="NZ_JH417859.fa"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
pipeline_code="basic_pipeline"
pipeline_id="bfecca03-6443-45bd-b313-e4f555cd0748"
user_reference="regan_test_analysis_03"
storage_size="Small"
local_download_path="${HOME}/Documents/ica_data_downloads/"

time_stamp=$(date +"%Y-%m-%d_%H:%M:%S")

upload_file_path="$file_path/$file_name"

# We need the fileId and analysisId to start Nextflow analysis
file_id=""
analysis_id=""

printf "[$time_stamp]: "
printf "Uploading '$file_name'... \n"

upload_file_response=$(icav2 projectdata upload $upload_file_path --project-id $project_id)

if [[ $? != 0 ]]; then
    printf "Failed to upload file '$file_name'. \n"
    exit 1
else
    temp_file_name="$file_name.txt"
    touch $temp_file_name
    printf "$upload_file_response\n" > $temp_file_name
    # id of file starts with 'fil.'
    file_id=$(grep -i "\"id\": \"fil\." $temp_file_name | grep -o 'fil[^\"]*')

    printf "[$time_stamp]: "
    printf "Uploaded file with ID '$file_id' \n"

    printf "Deleting file '$temp_file_name' \n"
    rm $temp_file_name
fi

printf "Fetching list of project analyses... \n"
projectanalyses_list_response=$(icav2 projectanalyses list --project-id $project_id)

if [[ $? != 0 ]]; then
    printf "Failed to fetch projectanalyses list. \n"
    exit 1
else
    analysis_id=$(echo $projectanalyses_list_response | jq -r ".items[] | select(.userReference == \"$user_reference\").id")
    printf "Fetching project analysis with id '$analysis_id'... \n"
    analysis_response=$(icav2 projectanalyses input $analysis_id)
    if [[ $? != 0 ]]; then
        printf "Failed to fetch projectanalysis with id '$analysis_id'. \n"
        exit 1
    else
        analysis_data_code=$(echo $analysis_response | jq -r ".items[].code")

        printf "Constructing 'analysisCode:fileId' from analysis response and file details... \n"
        file_ref="$analysis_data_code:$file_id"

        printf "Starting Nextflow analysis... \n"
        icav2 projectpipelines start nextflow $pipeline_id \
            --user-reference $user_reference \
            --project-id $project_id \
            --storage-size $storage_size \
            --input $file_ref
    fi
fi

printf "$time_stamp: "
printf "Fetching analysis output of analysis with id '$analysis_id'... \n"
analysis_output=$(icav2 projectanalyses output $analysis_id --project-id $project_id)

if [[ $? != 0 ]]; then
    printf "Failed to fetch analysis output. \n"
else
    # analysis_output_data=$(echo $analysis_output | jq -r ".items[].data[]")
    printf "Fetching analysis output folder id... \n"
    analysis_output_folder_id=$(echo $analysis_output | jq -r ".items[].data[].dataId")

    printf "Downloading analysis output folder with id '$analysis_output_folder_id' to '$local_download_path'... \n"
    analysis_output_download_response=$(icav2 projectdata download $analysis_output_folder_id $local_download_path)
fi