#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
sample_name="Citrobacter_30_2_sample"
sample_description="a sample containing Citrobacter_30_2 fasta files"

printf "Creating project sample '$sample_name'... \n"

projectsample_create_response=$(icav2 projectsamples create $sample_name \
    --project-id $project_id \
    --description $sample_description \
    --user-tag "regan_citrobacter_sample_01" \
    --technical-tag "citrobacter_fasta_sample")


if ! [ $? != 0 ]; then
    printf "Failed to create project sample '$sample_name'. \n"
    exit 1
else
    projectsample_id=$(echo $projectsample_create_response | jq -r ".sample.id")
    printf "Created project sample with id '$projectsample_id'. \n"
fi
