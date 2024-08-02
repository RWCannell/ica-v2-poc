#!/usr/bin/env nextflow
nextflow.enable.dsl=2
projectId = params.projectId
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval
analysisStatusCheckLimit = params.analysisStatusCheckLimit
read1FileId = params.read1FileId
read2FileId = params.read2FileId
analysisId = params.analysisId
readsFileUploadPath = params.readsFileUploadPath
referenceFileId = params.referenceFileId
readsPairFilesUploadPath = params.readsPairFilesUploadPath
referenceFileUploadPath = params.referenceFileUploadPath
localDownloadPath = params.localDownloadPath

process downloadAnalysisOutput {
    debug true
    
    input:
    val(analysisId)
    val(localDownloadPath)

    output:
    stdout

    script:
    """
    #!/bin/bash
    printf "[\${time_stamp}]: "
    printf "Fetching analysis output response...\n"
    analysis_output_response=\$(icav2 projectanalyses output ${analysisId})
    analysis_output_folder_id=\$(echo \${analysis_output_response} | jq -r ".items[].data[].dataId")
    printf "Analysis output folder ID is '\${analysis_output_folder_id}'\n"

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Downloading analysis output folder with ID '\${analysis_output_folder_id}' to '${localDownloadPath}'...\n"

    icav2 projectdata download \${analysis_output_folder_id} ${localDownloadPath}
    """
}

process deleteData {
    debug true

    input:
    val(read1FileId)
    val(read2FileId)
    val(analysisId)

    output:
    stdout

    script:
    """
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded read 1 file with ID '${read1FileId}'...\n"
    icav2 projectdata delete ${read1FileId}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded read 2 file with ID '${read2FileId}'...\n"
    icav2 projectdata delete ${read2FileId}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    analysis_output_response=\$(icav2 projectanalyses output ${analysisId})
    analysis_output_folder_id=\$(echo \${analysis_output_response} | jq -r ".items[].data[].dataId")
    printf "[\${timeStamp}]: Deleting analysis output folder with ID '\${analysis_output_folder_id}'...\n"
    icav2 projectdata delete \${analysis_output_folder_id}

    printf "Uploaded file and analysis output folder successfully deleted.\n"
    """
}

workflow {
    downloadAnalysisOutput(params.analysisId, params.localDownloadPath)
    deleteData(params.read1FileId, params.read2FileId, params.analysisId)
}
