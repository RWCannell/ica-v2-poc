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
read1AnalysisDataCode = params.read1AnalysisDataCode
read2AnalysisDataCode = params.read2AnalysisDataCode
fastqListDataCode = params.fastqListDataCode
sampleId = params.sampleId
read1FileId = params.read1FileId
read2FileId = params.read2FileId
fastqListFileId = params.fastqListFileId
outputFolderId = params.outputFolderId
analysisId = params.analysisId
referenceFileId = params.referenceFileId
readsPairFilesUploadPath = params.readsPairFilesUploadPath
referenceFileUploadPath = params.referenceFileUploadPath
localDownloadPath = params.localDownloadPath

process createDataFile {
    debug true

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Creating data file...\n"
    
    data_file="data.txt"
    touch \${data_file}

    printf "[\${time_stamp}]: "
    printf "Writing file data to existing data file...\n"

    printf "sampleId:${sampleId}\n" >> \${data_file}
    printf "${read1AnalysisDataCode}:${read1FileId}\n" >> \${data_file}
    printf "${read2AnalysisDataCode}:${read2FileId}\n" >> \${data_file}
    printf "${fastqListDataCode}:${fastqListFileId}\n" >> \${data_file}
    """
}

process downloadAnalysisOutput {
    debug true
    
    input:
    path(dataFile)

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash
    printf "[\${time_stamp}]: "
    printf "Writing id of analysis output folder to existing data file...\n"
    printf "outputFolderId:${outputFolderId}\n" >> ${dataFile}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Downloading analysis output folder with ID '${outputFolderId}' to '${localDownloadPath}'...\n"

    icav2 projectdata download ${outputFolderId} ${localDownloadPath}

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
    sample_id=\$(cat ${dataFile} | grep -o 'sampleId:.*' | cut -f2- -d:)
    read_1_file_id=\$(cat ${dataFile} | grep -o 'read1:.*' | cut -f2- -d:)
    read_2_file_id=\$(cat ${dataFile} | grep -o 'read2:.*' | cut -f2- -d:)
    fastq_list_file_id=\$(cat ${dataFile} | grep -o 'fastq_list:.*' | cut -f2- -d:)
    analysis_output_folder_id=\$(cat ${dataFile} | grep -o 'outputFolderId:.*' | cut -f2- -d:)

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded read 1 file with ID '\${read_1_file_id}'...\n"
    icav2 projectdata delete \${read_1_file_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded read 2 file with ID '\${read_2_file_id}'...\n"
    icav2 projectdata delete \${read_2_file_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded CSV file with ID '\${fastq_list_file_id}'...\n"
    icav2 projectdata delete \${fastq_list_file_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting analysis output folder with ID '\${analysis_output_folder_id}'...\n"
    icav2 projectdata delete \${analysis_output_folder_id}

    printf "Uploaded file and analysis output folder successfully deleted.\n"
    """
}

workflow {
    createDataFile()
    downloadAnalysisOutput(createDataFile.out.dataFile)
    deleteData(downloadAnalysisOutput.out.dataFile)
}
