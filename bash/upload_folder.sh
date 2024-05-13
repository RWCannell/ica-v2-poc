#!/bin/bash

folder_path="$HOME/Documents/ica_data_uploads/fasta"
folder_name="Yokenella_regensburgei_ATCC_43003_uid65133"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"

upload_folder_path="$folder_path/$folder_name"

printf "Uploading '$folder_name'... \n"

upload_folder_response=$(icav2 projectdata upload $upload_folder_path --project-id $project_id)

if [[ $? != 0 ]]; then
    printf "Failed to upload folder '$folder_name'. \n"
    exit 1
else
    temp_file_name="$folder_name.txt"
    touch $temp_file_name
    printf "$upload_folder_response\n" > $temp_file_name
    # id of folder starts with 'fol.'
    folder_id=$(grep -i "\"id\": \"fol\." $temp_file_name | grep -o 'fol[^\"]*')
    printf "Uploaded folder with ID '$folder_id' \n"
fi

rm $temp_file_name