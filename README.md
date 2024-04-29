# ica-v2-poc   
## Introduction   
This is a simple proof of concept (POC) for automating certain processes that use the _Illumina Connected Analytics (ICA)_ API. The main processes that we wish to automate are:   

- uploading files for analysis
- run Nextflow pipelines for analysis
- trigger download of output file(s)
- delete output file after download succeeds   

We can use a combination of both the [API](https://ica.illumina.com/ica/api/swagger/index.html#/) and the [CLI](https://help.ica.illumina.com/command-line-interface/cli-indexcommands).   

Before we can begin, we need to have an existing project or create a new project. For the rest of this README, we'll be referring to the existing project, **SGDP**.   

## Authentication
Authentication is required in order to use the API or the CLI. After logging in to the UI, an API key needs to be created. Instructions for generating an API key can be found over [here](https://help.ica.illumina.com/account-management/am-iam#api-keys).   

There are two ways to authenticate in order to make use of the API:
1. API key + JWT for the entire API, except for the POST `/tokens` endpoint.
2. API key + Basic Authentication (username/email and password) for the POST `/tokens` endpoint.    

When using the CLI, authentication takes place when running the command:
```bash
icav2 config set
```
There will be prompts. The defaults can be used by simply pressing `Enter` or `Return`. When the API key is prompted, provide the value that has been generated in the UI. 
```bash
icav2 config set
Creating /Users/regan/.icav2/config.yaml
Initialize configuration settings [default]
server-url [ica.illumina.com]: 
x-api-key : myAPIKey
output-format (allowed values table,yaml,json defaults to table) : 
colormode (allowed values none,dark,light defaults to none) :
```
The `/Users/regan/.icav2/config.yaml` file can be modified if the default settings are wished to be changed.

## Project and Project Data   
A project can be created in the UI. After a project is created, the project object can be obtained by using the following `curl` command with the API:
```bash
curl -X 'GET' \
  'https://ica.illumina.com/ica/rest/api/projects?includeHiddenProjects=false' \
  -H 'accept: application/vnd.illumina.v3+json' \
  -H 'X-API-Key: V21shKiHQvX24fF7vhVs0lUIWPZI5g'
```
The response body will contain a JSON object with a lot of details about the project. (The object is too big to be displayed here).   

Data for the project can be created by using the following `curl` command:
```bash
 curl -X 'POST' \
 'https://ica.illumina.com/ica/rest/api/projects/{projectId}/data' \
 -H 'accept: application/vnd.illumina.v3+json' \
 -H 'X-API-Key: XXXXXXXXXXXXXXXX' \
 -H 'Content-Type: application/vnd.illumina.v3+json' \
 -d '{
 "name": "tempFile.txt",
 "folderId": "fol.579eda846f1b4f6e2d1e08db91408069",
 "dataType": "FILE"
 }'
```
The example above generated a partial file called `tempFile.txt` within our project, and is situated inside a folder with the folder ID `fol.579eda846f1b4f6e2d1e08db91408069`. The project, file, or folder IDs can be accessed either by logging into the ICA web UI or by using the ICA V2 CLI.   

To get details about the data of the project, the `id` of the project needs to be sent as a path parameter in the following request:
```bash
curl -X 'GET' \
  'https://ica.illumina.com/ica/rest/api/projects/{projectId}/data?filePathMatchMode=STARTS_WITH_CASE_INSENSITIVE' \
  -H 'accept: application/vnd.illumina.v3+json' \
  -H 'X-API-Key: V21shKiHQvX24fF7vhVs0lUIWPZI5g'
```

## Uploading Files to Illumina Connected Analytics   
We can use the UI to upload a file to a project.   
![Uploading File using UI](public/assets/images/upload_file_using_ui.png "Uploading File using UI")

The CLI can also be used to upload a file with the following command:
```bash
icav2 projectdata upload /path/to/file/or/folder
```
To upload using the API, we can use `curl` and the correct REST API endpoint:
```bash

``` 

## Run Nextflow Pipelines   


## Trigger Download of Output File(s)   


## Delete Output File   

