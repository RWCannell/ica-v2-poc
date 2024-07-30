#!/usr/bin/env nextflow
nextflow.enable.dsl=2
projectId = params.projectId
read1AnalysisDataCode = params.read1AnalysisDataCode
read2AnalysisDataCode = params.read2AnalysisDataCode
referenceAnalysisDataCode = params.referenceAnalysisDataCode
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
hashTableConfigFile = params.hashTableConfigFile
referenceDirectory = params.referenceDirectory
intermediateResultsDirectory = params.intermediateResultsDirectory
fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval
analysisStatusCheckLimit = params.analysisStatusCheckLimit
readsFileUploadPath = params.readsFileUploadPath
referenceFileId = params.referenceFileId
readsPairFilesUploadPath = params.readsPairFilesUploadPath
referenceFileUploadPath = params.referenceFileUploadPath
localDownloadPath = params.localDownloadPath

process uploadFastqFilePairs {
    debug true
    maxForks 2
    input:
    tuple val(sampleId), file(reads)
    val(projectId)

    output:
    path "data.txt", emit: dataFile

    script:
    def (read_1_file, read_2_file) = reads
    """
    #!/bin/bash
    time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
    read_1_file_id=""
    read_2_file_id=""
    ica_upload_path="/fastq/$sampleId/"

    printf "sampleId: $sampleId \n"
    printf "reads: $reads \n"
    printf "read_1_file: $read_1_file \n"
    printf "read_2_file: $read_2_file \n"

    read_1_file_response="read_1_file_response.txt"
    read_2_file_response="read_2_file_response.txt"

    touch \${read_1_file_response}
    touch \${read_2_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading read 1 file '${read_1_file}'... \n"
    read_1_upload_response=\$(icav2 projectdata upload ${read_1_file} \${ica_upload_path} --project-id ${projectId})
    echo "\${read_1_upload_response}" > \${read_1_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading read 2 file '${read_2_file}'... \n"
    read_2_upload_response=\$(icav2 projectdata upload ${read_2_file} \${ica_upload_path} --project-id ${projectId})
    echo "\${read_2_upload_response}" > \${read_2_file_response}

    # id of file starts with 'fil.'
    read_1_file_id=\$(cat \${read_1_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
    read_2_file_id=\$(cat \${read_2_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    read_1_uploaded_file_data_response=\$(icav2 projectdata get \${read_1_file_id})
    if [[ \$? != 0 ]]; then
        printf "Failed to fetch data about file with id '\${read_1_file_id}'. \n"
        exit 1
    else
        read_1_uploaded_file_path=\$(echo \${read_1_uploaded_file_data_response} | jq -r ".details.path")
        printf "Path of uploaded file is '\${read_1_uploaded_file_path}'. \n"
    fi

    read_2_uploaded_file_data_response=\$(icav2 projectdata get \${read_2_file_id})
    if [[ \$? != 0 ]]; then
        printf "Failed to fetch data about file with id '\${read_2_file_id}'. \n"
        exit 1
    else
        read_2_uploaded_file_path=\$(echo \${read_2_uploaded_file_data_response} | jq -r ".details.path")
        printf "Path of uploaded file is '\${read_2_uploaded_file_path}'. \n"
    fi

    data_file="data.txt"

    if ! [ -f \${data_file} ]; then
        echo "Data file does not exist. Creating one..."
        touch \${data_file}
    fi

    printf "[\${time_stamp}]: "
    printf "Writing file data to existing data file...\n"

    printf "sampleId:${sampleId}\n" >> \${data_file}
    printf "${read1AnalysisDataCode}:\${read_1_file_id}\n" >> \${data_file}
    printf "${read2AnalysisDataCode}:\${read_2_file_id}\n" >> \${data_file}
    """
}

process getReferenceFile {
  debug true

  input:
  path(dataFile)

  output:
  path "data.txt", emit: dataFile

  script:
  def reference_file_name = ""
  """
  #!/bin/bash
  time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")

  get_reference_file_response_file="get_reference_file_response.txt"
  touch \${get_reference_file_response_file}

  printf "[\${time_stamp}]: "
  printf "Getting reference file data JSON response...\n"
  get_reference_file_response=\$(icav2 projectdata get ${referenceFileId} --project-id ${projectId})
  reference_file_name=\$(echo \${get_reference_file_response} | jq -r ".details.name")
  echo "\${get_reference_file_response}" > \${get_reference_file_response_file}

  printf "[\${time_stamp}]: "
  printf "Writing reference file '\${reference_file_name}' id to existing data file...\n"

  printf "${referenceAnalysisDataCode}:${referenceFileId}\n" >> ${dataFile}
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

    sample_id=\$(cat ${dataFile} | grep -o 'sampleId:.*' | cut -f2- -d:)
    read_1_file_id=\$(cat ${dataFile} | grep -o 'read1:.*' | cut -f2- -d:)
    read_2_file_id=\$(cat ${dataFile} | grep -o 'read2:.*' | cut -f2- -d:)

    read_1_analysis_code=\$(cat ${dataFile} | grep -E "read1")
    read_2_analysis_code=\$(cat ${dataFile} | grep -E "read2")
    reference_analysis_code=\$(cat ${dataFile} | grep -E "ref_tar")

    output_directory="/output/\${sample_id}/"

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Starting Nextflow analysis...\n"

    analysisResponse=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference ${userReference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${reference_analysis_code} \
        --input fastqs:"\${read_1_file_id},\${read_2_file_id}" \
        --parameters enable_map_align:true \
        --parameters enable_map_align_output:false \
        --parameters output_format:BAM \
        --parameters enable_variant_caller:true \
        --parameters vc_emit_ref_confidence:BP_RESOLUTION \
        --parameters enable_cnv:false \
        --parameters enable_sv:false \
        --parameters repeat_genotype_enable:false \
        --parameters enable_hla:false \
        --parameters enable_variant_annotation:false \
        --parameters output_file_prefix:"\${sample_id}")

    analysis_response_file="analysis_response.txt"
    touch \${analysis_response_file}
    echo "\${analysis_response}" > \${analysis_response_file}

    analysis_id=\$(cat \${analysis_response} | jq -r ".id")
    analysis_ref=\$(cat \${analysis_response} | jq -r ".reference")

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
            break;

        elif [[ \${analysis_status} == "FAILED" ]]; then
            printf "Analysis FAILED \n"
            break;

        elif [[ \${analysis_status} == "FAILED_FINAL" ]]; then
            printf "Analysis FAILED_FINAL\n"
            break;

        elif [[ \${analysis_status} == "ABORTED" ]]; then
            printf "Analysis ABORTED\n"
            break;

        elif [[ \${analysis_status_check_count} -gt ${analysisStatusCheckLimit} ]]; then
            printf "Analysis status has been checked more than ${analysisStatusCheckLimit} times. Stopping...\n"
            break;

        else
            printf "Analysis still in progress...\n"
        fi

        printf "Fetching analysis output response...\n"
        analysis_output_response=\$(icav2 projectanalyses output \${analysis_id})
        analysis_output_folder_id=\$(echo \${analysis_output_response} | jq -r ".items[].data[].dataId")
        printf "Analysis output folder ID is '\${analysis_output_folder_id}'\n"

        printf "[\${time_stamp}]: "
        printf "Writing id of analysis output folder to existing data file...\n"
        printf "outputFolderId:\${analysis_output_folder_id}\n" >> ${dataFile}

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
    stdout

    script:
    """
    #!/bin/bash

    output_folder_id=\$(cat ${dataFile} | grep -o 'outputFolderId:.*' | cut -f2- -d:)

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Downloading analysis output folder with ID '\${output_folder_id}' to '${localDownloadPath}'...\n"

    icav2 projectdata download \${output_folder_id} ${localDownloadPath}
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
    output_folder_id=\$(cat ${dataFile} | grep -o 'outputFolderId:.*' | cut -f2- -d:)

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded read 1 file with ID '\${read_1_file_id}'...\n"
    icav2 projectdata delete \${read_1_file_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded read 2 file with ID '\${read_2_file_id}'...\n"
    icav2 projectdata delete \${read_2_file_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting analysis output folder with ID '\${output_folder_id}'...\n"
    icav2 projectdata delete \${output_folder_id}

    printf "Uploaded file and analysis output folder successfully deleted.\n"
    """
}

workflow {
    fastqFilePairs = Channel.fromFilePairs(readsPairFilesUploadPath, checkIfExists:true)
    referenceFilePath = Channel.fromPath(params.referenceFileUploadPath, checkIfExists: true)

    uploadFastqFilePairs(fastqFilePairs, params.projectId)
    getReferenceFile(uploadFastqFilePairs.out.dataFile)
    startAnalysis(getReferenceFile.out.dataFile)
    checkAnalysisStatus(startAnalysis.out.dataFile, params.analysisStatusCheckInterval)
    downloadAnalysisOutput(checkAnalysisStatus.out.dataFile, params.localDownloadPath)
    deleteData(checkAnalysisStatus.out.dataFile)
}
