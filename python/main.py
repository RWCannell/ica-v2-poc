import requests

# The base API endpoint
base_url = "https://ica.illumina.com/ica/rest/api/"

# The project id
project_id = "049307d6-85dd-4cdc-b88d-a740e4e9e550"

# The analysis id
analysis_id = "387b5178-732f-4706-9c41-e67a0cd00dc6"

# The request url
request_url = base_url + "projects/" + project_id + "/analyses" + analysis_id

# The headers to be sent with the request
request_headers = { "accept: application/vnd.illumina.v4+json",
                    "X-API-Key: XXXXXXXXXXXXXXXX"
                  }

# A GET request to the API
response = requests.get(request_url)

# Print the response
response_json = response.json()
print(response_json)
