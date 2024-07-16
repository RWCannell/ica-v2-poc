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
fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval
readsFileUploadPath = params.readsFileUploadPath
readsPairFilesUploadPath = params.readsPairFilesUploadPath
referenceFileUploadPath = params.referenceFileUploadPath
localDownloadPath = params.localDownloadPath
icaUploadPath = params.icaUploadPath

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

process uploadReferenceFile {
  debug true
  maxForks 3

  input:
  path(dataFile)
  path(referenceFileUploadPath)

  output:
  path "data.txt", emit: dataFile

  script:
  def reference_file = referenceFileUploadPath.baseName
  """
  #!/bin/bash
  time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
  reference_file_id=""
  reference_file_ica_path="/fastq/${reference_file}/"

  get_reference_file_response_file="get_reference_file_response.txt"
  reference_file_upload_response_file="reference_file_upload_response.txt"

  touch \${reference_file_response}

  printf "[\${time_stamp}]: "
  printf "Checking if reference file already exists in ICA...\n"
  get_reference_file_response=\$(icav2 projectdata get \${reference_file_ica_path} --project-id ${projectId})
  echo "\${get_reference_file_response}" > \${get_reference_file_response_file}

  if grep -iq "No data found for path" \${get_reference_file_response_file};then
      printf "[\${time_stamp}]: "
      printf "Reference file not found in ICA. Uploading reference file '${reference_file}'... \n"
      reference_file_upload_response=\$(icav2 projectdata upload ${referenceFileUploadPath} \${reference_file_ica_path} --project-id ${projectId})
      echo "\${reference_file_upload_response}" > \${reference_file_upload_response_file}

      # id of file starts with 'fil.'
      printf "[\${time_stamp}]: "
      printf "Extracting file_id of reference file '${reference_file}' from upload response... \n"
      reference_file_id=\$(cat \${reference_file_upload_response_file} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
  else
      # id of file starts with 'fil.'
      printf "[\${time_stamp}]: "
      printf "Extracting file_id of reference file '${reference_file}' from get response... \n"
      reference_file_id=\$(cat \${get_reference_file_response_file} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
  fi

  printf "[\${time_stamp}]: "
  printf "Writing file data to existing data file...\n"

  printf "${referenceAnalysisDataCode}:\${reference_file_id}\n" >> ${dataFile}
  """
}

process startAnalysis {
    debug true
    
    input:
    path(dataFile)

    output:
    path "data.txt", emit: dataFile
    path "analysisResponse.txt", emit: analysisResponse

    script:
    analysisResponse = ""

    """
    #!/bin/bash

    sample_id=\$(cat ${dataFile} | grep -E "sampleId")
    read1_analysis_code=\$(cat ${dataFile} | grep -E "read1")
    read2_analysis_code=\$(cat ${dataFile} | grep -E "read2")
    reference_analysis_code=\$(cat ${dataFile} | grep -E "ref_tar")
    output_directory="/output/\${sample_id}/"

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\${timeStamp}]: Starting Nextflow analysis...\n"

    analysisResponse=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference ${userReference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${read1_analysis_code} \
        --input \${read2_analysis_code} \
        --input \${reference_analysis_code} \
        --parameters enable-variant-caller:true \
        --parameters RGID:Illumina_RGID \
        --parameters RGSM:\${sample_id} \
        --parameters output-directory:\${output_directory} \
        --parameters output-file-prefix:\${sample_id} 

    touch "analysisResponse.txt"
    echo "\${analysisResponse}" > analysisResponse.txt
    """
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
    fastqFilePairs = Channel.fromFilePairs(readsPairFilesUploadPath, checkIfExists:true)
    referenceFilePath = Channel.fromPath(params.referenceFileUploadPath, checkIfExists: true)

    uploadFastqFilePairs(fastqFilePairs, params.projectId)
    uploadReferenceFile(uploadFastqFilePairs.out.dataFile, referenceFilePath)
    startAnalysis(uploadReferenceFile.out.dataFile)
    // checkAnalysisStatus(startAnalysis.out.analysisResponse, params.analysisStatusCheckInterval)
    // downloadAnalysisOutput(checkAnalysisStatus.out.analysisOutputFolderId, params.localDownloadPath)
    // deleteData(uploadFile.out.fileUploadResponse.view(), downloadAnalysisOutput.out.outputFolderId)
}
