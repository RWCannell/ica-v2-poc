#!/usr/bin/env nextflow
nextflow.enable.dsl=2
projectId = params.projectId
analysisDataCode = params.analysisDataCode
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
fileStatusCheckInterval = params.fileStatusCheckInterval
fileStatusCheckLimit = params.fileStatusCheckLimit
analysisStatusCheckInterval = params.analysisStatusCheckInterval
analysisStatusCheckLimit = params.analysisStatusCheckLimit
localUploadPath = params.localUploadPath
localDownloadPath = params.localDownloadPath

process uploadFile {
    debug true
    input:
    path(filePath)
    val(projectId)

    output:
    path "data.txt", emit: dataFile

    script:
    sampleId = filePath.baseName
    dataFile="data.txt"
    """
    #!/bin/bash
    time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
    file_id=""
    ica_upload_path="/fasta/${sampleId}/"

    upload_response_file="file_upload_response.txt"
    touch \${upload_response_file}

    printf "[\${time_stamp}]: "
    printf "Uploading file '${sampleId}'... \n"
    file_upload_response=\$(icav2 projectdata upload ${filePath} \${ica_upload_path} --project-id ${projectId})
    echo "\${file_upload_response}" > \${upload_response_file}

    # id of file starts with 'fil.'
    file_id=\$(cat \${upload_response_file} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    touch ${dataFile}

    printf "[\${time_stamp}]: "
    printf "Writing file data to existing data file...\n"

    printf "sampleId:${sampleId}\n" >> ${dataFile}
    printf "in:\${file_id}\n" >> ${dataFile}
    """
}

process startAnalysis {
    debug true
    
    input:
    path(dataFile)

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash
    file_status_check_count=0
    file_status="PARTIAL"
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")

    sample_id=\$(cat ${dataFile} | grep -o 'sampleId:.*' | cut -f2- -d:)
    file_id=\$(cat ${dataFile} | grep -o 'in:.*' | cut -f2- -d:)
    file_analysis_code=\$(cat ${dataFile} | grep -E "in")
    while true;
    do
        ((file_status_check_count +=1 ))
        file_data_response=\$(icav2 projectdata get \${file_id})
        file_status=\$(echo \${file_data_response} | jq -r ".details.status")

        if [[ \${file_status} == "AVAILABLE" ]]; then
            printf "File is AVAILABLE\n"
            break;
        elif [[ \${file_status} == "PARTIAL" ]]; then
            printf "File upload is PARTIAL\n"
        elif [[ \${file_status_check_count} -gt ${fileStatusCheckLimit} ]]; then
            printf "File status has been checked more than ${fileStatusCheckLimit} times. Stopping...\n"
            printf "analysisStatus:TIMEOUT\n" >> ${dataFile}
            break;
        else
            printf "File is not AVAILABLE\n"
        fi
        sleep ${fileStatusCheckInterval};
    done

    echo "[\${timeStamp}]: Starting Nextflow analysis..."
    analysis_response=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference ${userReference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${file_analysis_code})
    
    analysis_response_file="analysis_response.txt"
    touch \${analysis_response_file}
    echo "\${analysis_response}" > \${analysis_response_file}

    analysis_id=\$(cat \${analysis_response_file} | jq -r ".id")
    analysis_ref=\$(cat \${analysis_response_file} | jq -r ".reference")

    printf "[\${time_stamp}]: "
    printf "Writing id of analysis '\${analysis_ref}' to existing data file...\n"
    printf "analysisId:\${analysis_id}\n" >> ${dataFile}

    printf "Writing reference of analysis '\${analysis_ref}' to existing data file...\n"
    printf "analysisRef:\${analysis_ref}\n" >> ${dataFile}
    """
}

process checkAnalysisStatus {
    debug true
    
    input:
    path(dataFile)
    val(analysisStatusCheckInterval)

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash

    analysis_status_check_count=0
    analysis_status="REQUESTED"

    analysis_id=\$(cat ${dataFile} | grep -o 'analysisId:.*' | cut -f2- -d:)
    analysis_ref=\$(cat ${dataFile} | grep -o 'analysisRef:.*' | cut -f2- -d:)
    
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Checking status of analysis with id '\${analysis_id}' every ${analysisStatusCheckInterval} seconds, until status is 'SUCCEEDED'...\n"
    while true;
    do
        ((analysis_status_check_count +=1 ))
        updated_analysis_response=\$(icav2 projectanalyses get \${analysis_id})

        printf "Checking status of analysis with reference '\${analysis_ref}'...\n"
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
    sample_id=\$(cat ${dataFile} | grep -o 'sampleId:.*' | cut -f2- -d:)
    file_id=\$(cat ${dataFile} | grep -o 'in:.*' | cut -f2- -d:)
    analysis_output_folder_id=\$(cat ${dataFile} | grep -o 'outputFolderId:.*' | cut -f2- -d:)

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded file with ID '\${file_id}'...\n"
    icav2 projectdata delete \${file_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting analysis output folder with ID '\${analysis_output_folder_id}'...\n"
    icav2 projectdata delete \${analysis_output_folder_id}

    printf "Uploaded file and analysis output folder successfully deleted.\n"
    """

}

workflow {
    filePath = Channel.fromPath(params.localUploadPath, checkIfExists: true)

    uploadFile(filePath, params.projectId)
    
    startAnalysis(uploadFile.out.dataFile)

    checkAnalysisStatus(startAnalysis.out.dataFile, params.analysisStatusCheckInterval)

    downloadAnalysisOutput(checkAnalysisStatus.out.dataFile, params.localDownloadPath)

    deleteData(downloadAnalysisOutput.out.dataFile)
}
