PROJECTS=$(icav2 projects list)
PROJECT_ID=$(echo $PROJECTS | jq -r '.items[] | select(.name == "SGDP").id')
echo $PROJECT_ID
