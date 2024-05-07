#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
analysis_user_reference="regan_test_analysis_01"

projectanalyses_response=$(icav2 projectanalyses list --project-id $project_id)

analysis_id=$(echo $projectanalyses_response | jq -r ".items[] | select(.userReference == \"$analysis_user_reference\").id")
analysis_status=$(echo $projectanalyses_response | jq -r ".items[] | select(.userReference == \"$analysis_user_reference\").status")

echo $projectanalyses_response
echo $analysis_id

