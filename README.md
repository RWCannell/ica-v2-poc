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
Creating /Users/regancannell/.icav2/config.yaml
Initialize configuration settings [default]
server-url [ica.illumina.com]: 
x-api-key : myAPIKey
output-format (allowed values table,yaml,json defaults to table) : 
colormode (allowed values none,dark,light defaults to none) :
```
The `/Users/regancannell/.icav2/config.yaml` file can be modified if the default settings are wished to be changed.

## Project and Project Data   
A project can be created in the UI. After a project is created, the project object can be obtained by using the following `curl` command with the API:
```bash
curl -X 'GET' \
  'https://ica.illumina.com/ica/rest/api/projects?includeHiddenProjects=false' \
  -H 'accept: application/vnd.illumina.v3+json' \
  -H 'X-API-Key: XXXXXXXXXXXXXXXX'
```
The response body will contain a JSON object with a lot of details about the project. (The object is too big to be displayed here).   

## Uploading Files to Illumina Connected Analytics   
We can use the UI to upload a file to a project.   
![Uploading File using UI](public/assets/images/upload_file_using_ui.png "Uploading File using UI")

The CLI can also be used to upload a file. We can get a list of projects with the following CLI command:
```bash
icav2 projects list
```
From the returned object (which is JSON by default), the `projectId` should be noted, since it will be required in subsequent commands/requests. To set the project context, run the following command:
```bash
icav2 projects enter <projectId>
```
If the command runs successfully, then a response like
```bash
"Context switched to project <projectId>"
```
should be received. Now the terminal session will be in the specified project's context. It's possible to now upload a file to the project with:
```bash
icav2 projectdata upload <localFileFolder>
```
If a remote path for uploading is not specified, the data will be uploaded to the top level of the project's storage folder.   

The API can also be used for uploading data. Data for the project can be created by using the following `curl` command:
```bash
 curl -X 'POST' \
 'https://ica.illumina.com/ica/rest/api/projects/{projectId}/data' \
 -H 'accept: application/vnd.illumina.v3+json' \
 -H 'X-API-Key: XXXXXXXXXXXXXXXX' \
 -H 'Content-Type: application/vnd.illumina.v3+json' \
 -d '{
 "name": "gencode.v45.lncRNA_transcripts.fa",
 "folderId": "fol.579eda846f1b4f6e2d1e08db91408069",
 "dataType": "FILE"
 }'
```
The example above generated a partial file called `gencode.v45.lncRNA_transcripts.fa` within our project, and is situated inside a folder with the folder ID `fol.579eda846f1b4f6e2d1e08db91408069`. The project, file, or folder IDs can be accessed either by logging into the ICA web UI or by using the ICA V2 CLI.   

To get details about the data of the project, the `id` of the project needs to be sent as a path parameter in the following request:
```bash
curl -X 'GET' \
  'https://ica.illumina.com/ica/rest/api/projects/{projectId}/data?filePathMatchMode=STARTS_WITH_CASE_INSENSITIVE' \
  -H 'accept: application/vnd.illumina.v3+json' \
  -H 'X-API-Key: XXXXXXXXXXXXXXXX'
```

## Downloading Files from Illumina Connected Analytics   
Files can be downloaded from the ICA storage using the CLI with the following CLI command:
```bash
icav2 projectdata list # take note of the dataId
icav2 projectdata download <dataId>
```

## Run Nextflow Pipelines   
A Nextflow pipeline can be created in the web UI. A tutorial on creating a Nextflow pipeline and running an analysis through the web UI can be found over [here](https://help.ica.illumina.com/tutorials/nextflow). The analysis from the example in the tutorial takes about 30 minutes to complete. When completed, there is output data in the ICA storage.

## Trigger Download of Output File(s)   
The output data files can be downloaded from the ICA storage using the web UI. A download can even be scheduled through the web UI.   
![Schedule Download through Web UI](public/assets/images/download_scheduled_with_ui.png "Schedule Download through Web UI")

## Delete Output File   

