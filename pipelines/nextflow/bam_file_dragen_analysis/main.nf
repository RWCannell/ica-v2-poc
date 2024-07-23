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
bamFileId = params.bamFileId
referenceFileId = params.referenceFileId
localDownloadPath = params.localDownloadPath
icaUploadPath = params.icaUploadPath

process startDragenAnalysis {
    debug true
    
    input:
    path(dataFile)

    output:
    path "analysisResponse.txt", emit: analysisResponse

    script:
    analysisResponse = ""

    """
    #!/bin/bash
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Extracting bam file and reference file id from data file, respectively...\n"

    bam_analysis_code

    printf "[\${timeStamp}]: Starting Nextflow analysis...\n"

    fileReference=\$(cat ${fileRef})

    echo "File Ref: \${fileReference}"

    analysisResponse=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference ${userReference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${bam_analysis_code} \
        --input \${reference_analysis_code} \
        --parameters enable-variant-caller:true \
        --parameters RGID:Illumina_RGID \
        --parameters RGSM:\${sample_id} \
        --parameters output-directory:\${output_directory} \
        --parameters output-file-prefix:\${sample_id} 

    echo "\${analysisResponse}" > analysisResponse.txt
    """

    scp regan@ui.core.wits.ac.za:~/dataG/SGDVP/ERR1019034 /Users/regancannell/Documents/ica_data_uploads/fastq
}

process checkAnalysisStatus {
    debug true
    
    input:
    path(analysisResponse)
    val(analysisStatusCheckInterval)

    output:
    path "analysisOutputFolderId.txt", emit: analysisOutputFolderId

    script:
    analysisOutputFolderId = ""
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
  dataFileChannel = Channel.fromFilePairs("data.txt", checkIfExists:true)

  startDragenAnalysis(dataFileChannel)
}
