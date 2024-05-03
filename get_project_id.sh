#!/bin/bash

project_name="SGDP"

projects=$(icav2 projects list)
project_id=$(echo $projects | jq -r ".items[] | select(.name == \"$project_name\").id")
echo $project_id
