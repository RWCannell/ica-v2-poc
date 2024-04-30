import requests

# The base API endpoint
projects_base_url = "https://ica.illumina.com/ica/rest/api/projects/"

# The project id
project_id = "049307d6-85dd-4cdc-b88d-a740e4e9e550"

# The analysis id
analysis_id = "387b5178-732f-4706-9c41-e67a0cd00dc6"

# The headers to be sent with the request
request_headers = { "accept: application/vnd.illumina.v4+json",
                    "X-API-Key: XXXXXXXXXXXXXXXX"
                  }

# A GET request to the API for a project by id
async def get_project_by_id(project_id, headers):
    project_request_url = project_id
    response = requests.get(project_request_url, headers)
    response_json = response.json()
    return response_json

# A GET request to the API for an analysis by id
async def get_analysis_by_id(project_id, analysis_id, headers):
    analysis_request_url = project_id + "/analyses" + analysis_id
    response = requests.get(analysis_request_url, headers)
    response_json = response.json()
    return response_json

# A GET request to the API for the status of analysis
async def get_analysis_status(project_id, analysis_id, headers):
    analysis_by_id = get_analysis_by_id(project_id, analysis_id, headers)
    status_of_analysis = analysis_by_id['status']
    if status_of_analysis == 'SUCCEEDED':
        print("Analysis is complete and has SUCCEEDED")
    return status_of_analysis

# Print results
project_by_id = get_project_by_id("049307d6-85dd-4cdc-b88d-a740e4e9e550", request_headers)
analysis_by_id = get_analysis_by_id("049307d6-85dd-4cdc-b88d-a740e4e9e550", "387b5178-732f-4706-9c41-e67a0cd00dc6", request_headers)

print(project_by_id)
print(analysis_by_id)
