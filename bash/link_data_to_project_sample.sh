#!/bin/bash

folder_path="$HOME/Documents/ica_data_uploads/fasta"
folder_name="Citrobacter_freundii_4_7_47CFAA_uid46379"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
sample_id="4bef7ead-314c-4651-8bb6-da0656dc8332"
sample_name="Citrobacter_freundii_sample"
sample_description="a sample containing Citrobacter_freundii_4_7_47CFAA fasta files"
time_stamp=$(date +"%Y-%m-%d_%H:%M:%S")

upload_folder_path="$folder_path/$folder_name"

printf "[$time_stamp]: "
printf "Fetching list of samples... \n"
project_sample_list_response=$(icav2 projectsamples list --project-id $project_id)
if [[ $? != 0 ]]; then
    printf "[$time_stamp]: "
    printf "Failed to fetch list of project samples. \n"
    exit 1
else
    project_sample_id=$(echo $project_sample_list_response | jq -r ".items[].sample | select(.name == \"$sample_name\").id")
    project_sample_status=$(echo $project_sample_list_response | jq -r ".items[].sample | select(.name == \"$sample_name\").id")
    printf "[$time_stamp]: "
    printf "Found project sample '$sample_name' with id '$project_sample_id'. \n"
fi

printf "Fetching path of folder '$folder_name'... \n"
projectdata_list_response=$(icav2 projectdata list --project-id $project_id)

if [[ $? != 0 ]]; then
    printf "[$time_stamp]: "
    printf "Failed to fetch list of project data. \n"
    exit 1
else
    folder_path=$(echo $projectdata_list_response | jq -r ".items[].details | select(.name == \"$folder_name\").path")
    printf "[$time_stamp]: "
    printf "Found path '$folder_path' of folder '$folder_name'. \n"

    printf "[$time_stamp]: " 
    printf "Fetching details of folder '$folder_name'... \n"
    folder_details_response=$(icav2 projectdata get $folder_path --project-id $project_id)
    if [[ $? != 0 ]]; then
        printf "[$time_stamp]: "
        printf "Failed to fetch details of folder '$folder_name'. \n"
        exit 1
    else
        folder_id=$(echo $folder_details_response | jq -r ".id")
        printf "[$time_stamp]: " 
        printf "Fetched id '$folder_id' of folder '$folder_name'... \n"

        printf "[$time_stamp]: "
        printf "Linking sample '$sample_name' with folder '$folder_name'... \n"

        project_sample_link_response=$(icav2 projectsamples link $project_sample_id \
            --project-id $project_id \
            --data-id $folder_id)
        if [[ $? != 0 ]]; then
            printf "[$time_stamp]: "
            printf "Failed to link project sample '$sample_name' with folder '$folder_name'. \n"
            exit 1
        else
            printf "[$time_stamp]: "
            printf "Successfully linked project sample '$sample_name' with folder '$folder_name'. \n"
        fi
    fi
fi
