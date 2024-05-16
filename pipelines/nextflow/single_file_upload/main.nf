#!/usr/bin/env nextflow
nextflow.enable.dsl=2
filePath = Channel.fromPath("ica_data_uploads/fasta/Citrobacter_30_2_uid32453/NZ_ACDJ00000000.scaffold.fa/NZ_GG657370.fa", checkIfExists: true)
projectId = params.projectId
analysisDataCode = params.analysisDataCode
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
fileStatusCheckInterval = params.fileStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval

process uploadFile {
    debug true
    input:
    path(filePath)
    val(projectId)

    output:
    path "${fileName}.txt", emit: fileUploadResponse

    script:
    fileName = filePath.baseName
    """
    #!/bin/bash

    echo "Uploading file '${fileName}' to project with id '${projectId}'..."

    touch ${fileName}.txt

    fileUploadResponse=\$(icav2 projectdata upload ${filePath} --project-id ${projectId})

    echo "Successfully uploaded file '${fileName}' to project with id '${projectId}'."

    echo "\${fileUploadResponse}" > ${fileName}.txt
    """
}

process checkFileUploadStatus {
    debug true
    
    input:
    path(fileUploadResponse)
    val(analysisDataCode)
    val(fileStatusCheckInterval)

    output:
    path "fileReference.txt", emit: fileRef

    script:
    fileReference = ""
    """
    #!/bin/bash
    
    fileId=\$(cat ${fileUploadResponse} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    fileStatusCheckCount=0
    fileStatusCheckLimit=10

    while true;
    do
        echo "Checking status of uploaded file..."
        ((fileStatusCheckCount+=1))
        fileResponse=\$(icav2 projectdata get \${fileId} --project-id ${projectId})
        fileStatus=\$(echo \$fileResponse | jq -r ".details.status")
        if [[ "\${fileStatus}" == "AVAILABLE" ]]; then
            echo "File is AVAILABLE"

            fileReference="${analysisDataCode}:\${fileId}"

            echo "File Reference: "

            echo \${fileReference}

            touch fileReference.txt

            echo "\${fileReference}" > fileReference.txt
            break;
        elif
            [[ fileStatusCheckLimit > 10 ]]
            echo "File status has been checked more than \${fileStatusCheckLimit} times. Stopping..."
        else
            printf "File '\${fileId}' is still not AVAILABLE... \n"
        fi
        sleep \$fileStatusCheckInterval;
    done
    """
}

process startAnalysis {
    debug true
    
    input:
    path(fileRef)

    output:
    val analysisResponse, emit: analysisResponse

    script:
    analysisResponse = ""

    """
    #!/bin/bash
    echo "Starting Nextflow analysis..."

    fileReference=\$(cat ${fileRef})

    echo "File Ref: \${fileReference}"

    analysisResponse=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference ${userReference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${fileReference})

    """
}

process checkAnalysisStatus {
    debug true
    
    input:
    val(analysisResponse)
    val(analysisStatusCheckInterval)

    output:
    val analysisOutputFolderId, emit: analysisOutputFolderId

    script:
    analysisOutputFolderId = ""
    """
    #!/bin/bash

    analysisStatusCheckCount=0
    analysisStatusCheckLimit=10
    analysisStatus="REQUESTED"

    analysisId=\$(echo ${analysisResponse} | jq -r ".id")
    analysisRef=\$(echo ${analysisResponse} | jq -r ".reference")
    
    echo "Checking status of analysis with id '\${analysisId}' every ${analysisStatusCheckInterval} seconds, until status is 'SUCCEEDED'..."
    while true;
    do
        updatedAnalysisResponse=\$(icav2 projectanalyses get \${analysisId})
        analysisStatus=$(echo \${updatedAnalysisResponse} | jq -r ".items[] | select(.reference == "\${analysisRef}").status")

        if [ \${analysisStatus} == "SUCCEEDED" ]; then
            echo "Analysis SUCCEEDED"
            echo "Fetching analysis output response..."
            analysisOutputResponse=\$(icav2 projectanalyses output \$analysisId)
            analysisOutputFolderId=$(echo \${analysisOutputResponse} | jq -r ".items[].data[].dataId")
            echo "Analysis output folder ID is '\${analysisOutputFolderId}'"
            break;

        elif [ \${analysisStatus} == "FAILED" ]; then
            echo "Analysis FAILED \n"
            break;

        elif [ \${analysisStatus} == "FAILED_FINAL" ]; then
            echo "Analysis FAILED_FINAL"
            break;

        elif [ \${analysisStatus} == "ABORTED" ]; then
            echo "Analysis ABORTED"
            break;

        elif [[ fileStatusCheckLimit > 10 ]]
            echo "Analysis status has been checked more than \${analysisStatusCheckLimit} times. Stopping..."

        else
            echo "Analysis still in progress..."
        fi

        sleep ${analysisStatusCheckInterval};
    done
    """
}

process downloadAnalysisOutput {
    debug true
    
    input:
    val(analysisOutputFolderId)
    val(localDownloadPath)

    output:
    stdout

    script:
    """
    #!/bin/bash

    echo "Downloading analysis output folder with ID '${analysisOutputFolderId}' to '${localDownloadPath}'..."

    icav2 projectdata download ${analysisOutputFolderId} ${local_download_path}
    """
}

workflow {
    uploadFile(filePath, params.projectId)
    
    checkFileUploadStatus(uploadFile.out.fileUploadResponse.view(), params.analysisDataCode, params.fileStatusCheckInterval)

    startAnalysis(checkFileUploadStatus.out.fileRef.view())

    checkAnalysisStatus(startAnalysis.out.analysisResponse.view(), params.analysisStatusCheckInterval)

    downloadAnalysisOutput(checkAnalysisStatus.out.analysisOutputFolderId, params.localDownloadPath)
}
