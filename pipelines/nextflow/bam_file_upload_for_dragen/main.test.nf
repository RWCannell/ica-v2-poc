#!/usr/bin/env nextflow
nextflow.enable.dsl=2
projectId = params.projectId
bamAnalysisDataCode = params.bamAnalysisDataCode
referenceAnalysisDataCode = params.referenceAnalysisDataCode
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval
bamFilesUploadPath = params.bamFilesUploadPath
bamFilePairsUploadPath = params.bamFilePairsUploadPath
referenceFileUploadPath = params.referenceFileUploadPath
referenceFileIcaPath = params.referenceFileIcaPath
referenceFileId = params.referenceFileId
localDownloadPath = params.localDownloadPath
icaUploadPath = params.icaUploadPath

process checkAnalysisStatus {
    debug true
    
    input:
    path(analysisResponse)
    val(analysisStatusCheckInterval)

    output:
    stdout

    script:
    """
    #!/bin/bash

    analysisStatusCheckCount=0
    analysisStatusCheckLimit=10
    analysisStatus="REQUESTED"

    analysisId=\$(cat ${analysisResponse} | jq -r ".id")
    analysisRef=\$(cat ${analysisResponse} | jq -r ".reference")
    
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\${timeStamp}]: Checking status of analysis with id '\${analysisId}' every ${analysisStatusCheckInterval} seconds, until status is 'SUCCEEDED'..."
    while true;
    do
        ((\${StatusCheckCount}+=1))
        updatedAnalysisResponse=\$(icav2 projectanalyses get \${analysisId})

        echo "Checking status of analysis with reference '\${analysisRef}'..."
        analysisStatus=\$(echo \${updatedAnalysisResponse} | jq -r ".status")

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        echo "[\${timeStamp}]: Current status of analysis is '\${analysisStatus}'..."

        if [[ \${analysisStatus} == "SUCCEEDED" ]]; then
            echo "Analysis SUCCEEDED"
            echo "Fetching analysis output response..."
            analysisOutputResponse=\$(icav2 projectanalyses output \$analysisId)
            analysisOutputFolderId=\$(echo \${analysisOutputResponse} | jq -r ".items[].data[].dataId")
            echo "Analysis output folder ID is '\${analysisOutputFolderId}'"

            touch analysisOutputFolderId.txt
            echo "\${analysisOutputFolderId}" > analysisOutputFolderId.txt
            break;

        elif [[ \${analysisStatus} == "FAILED" ]]; then
            echo "Analysis FAILED \n"
            break;

        elif [[ \${analysisStatus} == "FAILED_FINAL" ]]; then
            echo "Analysis FAILED_FINAL"
            break;

        elif [[ \${analysisStatus} == "ABORTED" ]]; then
            echo "Analysis ABORTED"
            break;

        elif [[ \${analysisStatusCheckCount} -gt \${analysisStatusCheckLimit} ]]; then
            echo "Analysis status has been checked more than \${analysisStatusCheckLimit} times. Stopping..."
            break;

        else
            echo "Analysis still in progress..."
        fi

        sleep ${analysisStatusCheckInterval};
    done
    """
}

workflow {
  bamFilePairsChannel = Channel.fromFilePairs(params.bamFilePairsUploadPath, checkIfExists:true) { 
    file -> file.name.replaceAll(/.bam|.bai$/,'') 
  }
  analysisResponsePath = Channel.fromPath("analysis_response_example.txt", checkIfExists:true)
  checkAnalysisStatus(analysisResponsePath, analysisStatusCheckInterval)
}

