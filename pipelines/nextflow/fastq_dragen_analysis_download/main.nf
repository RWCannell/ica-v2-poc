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

process checkAnalysisStatus {
    debug true
    
    input:
    val(analysisId)
    val(analysisStatusCheckInterval)

    output:
    path "data.txt", emit: dataFile

    script:
    def dataFile = "data.txt" 
    """
    #!/bin/bash

    touch ${dataFile}
    printf "read1:${read1FileId}\n" >> ${dataFile}
    printf "read2:${read2FileId}\n" >> ${dataFile}
    printf "analysisId:${analysisId}\n" >> ${dataFile}

    analysis_status_check_count=0
    analysis_status="REQUESTED"

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Checking status of analysis with id '${analysisId}' every ${analysisStatusCheckInterval} seconds, until status is 'SUCCEEDED'...\n"
    while true;
    do
        ((analysis_status_check_count +=1 ))
        updated_analysis_response=\$(icav2 projectanalyses get ${analysisId})
        analysis_status=\$(echo \${updated_analysis_response} | jq -r ".status")

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Current status of analysis is '\${analysis_status}'...\n"

        if [[ \${analysis_status} == "SUCCEEDED" ]]; then
            printf "Analysis SUCCEEDED\n"
            printf "analysisStatus:SUCCEEDED\n" >> ${dataFile}
            break;

        elif [[ \${analysis_status} == "FAILED" ]]; then
            printf "Analysis FAILED \n"
            printf "analysisStatus:FAILED\n" >> ${dataFile}
            break;

        elif [[ \${analysis_status} == "FAILED_FINAL" ]]; then
            printf "Analysis FAILED_FINAL\n"
            printf "analysisStatus:FAILED_FINAL\n" >> ${dataFile}
            break;

        elif [[ \${analysis_status} == "ABORTED" ]]; then
            printf "Analysis ABORTED\n"
            printf "analysisStatus:ABORTED\n" >> ${dataFile}
            break;

        elif [[ \${analysis_status_check_count} -gt ${analysisStatusCheckLimit} ]]; then
            printf "Analysis status has been checked more than ${analysisStatusCheckLimit} times. Stopping...\n"
            printf "analysisStatus:TIMEOUT\n" >> ${dataFile}
            break;

        else
            printf "Analysis still in progress...\n"
        fi

        sleep ${analysisStatusCheckInterval};
    done
    """
}

process downloadAnalysisOutput {
    debug true
    
    input:
    path(dataFile)
    val(localDownloadPath)

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash
    analysis_id=\$(cat ${dataFile} | grep -o 'analysisId:.*' | cut -f2- -d:)

    printf "[\${time_stamp}]: "
    printf "Fetching analysis output response...\n"
    analysis_output_response=\$(icav2 projectanalyses output \${analysis_id})
    analysis_output_folder_id=\$(echo \${analysis_output_response} | jq -r ".items[].data[].dataId")
    printf "Analysis output folder ID is '\${analysis_output_folder_id}'\n"
    printf "Writing id of analysis output folder to existing data file...\n"
    printf "outputFolderId:\${analysis_output_folder_id}\n" >> ${dataFile}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Downloading analysis output folder with ID '\${analysis_output_folder_id}' to '${localDownloadPath}'...\n"

    icav2 projectdata download \${analysis_output_folder_id} ${localDownloadPath}

    printf "downloadComplete:true\n" >> ${dataFile}
    """
}

process deleteData {
    debug true

    input:
    path(dataFile)

    output:
    stdout

    script:
    """
    read_1_file_id=\$(cat ${dataFile} | grep -o 'read1:.*' | cut -f2- -d:)
    read_2_file_id=\$(cat ${dataFile} | grep -o 'read2:.*' | cut -f2- -d:)
    analysis_output_folder_id=\$(cat ${dataFile} | grep -o 'outputFolderId:.*' | cut -f2- -d:)
    download_complete=\$(cat ${dataFile} | grep -o 'downloadComplete:.*' | cut -f2- -d:)

    if [ "\${download_complete}" = true ]; then  
        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting uploaded read 1 file with ID '\${read_1_file_id}'...\n"
        icav2 projectdata delete \${read_1_file_id}
        
        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting uploaded read 2 file with ID '\${read_2_file_id}'...\n"
        icav2 projectdata delete \${read_2_file_id}

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting analysis output folder with ID '\${analysis_output_folder_id}'...\n"
        icav2 projectdata delete \${analysis_output_folder_id}
    fi
    printf "Uploaded file and analysis output folder successfully deleted.\n"
    """
}

workflow {
    checkAnalysisStatus(params.analysisId, params.analysisStatusCheckInterval)
    downloadAnalysisOutput(checkAnalysisStatus.out.dataFile, params.localDownloadPath)
    deleteData(downloadAnalysisOutput.out.dataFile)
}
