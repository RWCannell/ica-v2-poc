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

process uploadReferenceFile {
  debug true

  input:
  path(dataFilePath)
  val(referenceFileId)
  path(referenceFileUploadPath)

  output:
  path "data.txt", emit: dataFile

  script:
  def reference_file = referenceFileUploadPath.baseName
  """
  #!/bin/bash
  time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
  reference_file_id=${referenceFileId}
  reference_file_ica_path=${referenceFileIcaPath}

  get_reference_file_response_file="get_reference_file_response.txt"
  reference_file_upload_response_file="reference_file_upload_response.txt"

  touch \${reference_file_response}

  printf "[\${time_stamp}]: "
  printf "Checking if reference file id has been set in params.json...\n"
  if [ -z "${referenceFileId}" ]; then 
    printf "[\${time_stamp}]: "
    printf "Reference file id has not been set in params.json. Checking if reference file exists in ICA path...\n"

    get_reference_file_response=\$(icav2 projectdata get ${referenceFileIcaPath} --project-id ${projectId})
    echo "\${get_reference_file_response}" > \${get_reference_file_response_file}
  else
    printf "[\${time_stamp}]: "
    printf "Reference file id has been set in params.json. Getting reference file data JSON response...\n"
    get_reference_file_response=\$(icav2 projectdata get ${referenceFileId} --project-id ${projectId})
    echo "\${get_reference_file_response}" > \${get_reference_file_response_file}
  fi

  if grep -iq "No data found" \${get_reference_file_response_file}; then
    printf "[\${time_stamp}]: "
    printf "Reference file not found in ICA. Uploading reference file '${reference_file}'... \n"
    reference_file_upload_response=\$(icav2 projectdata upload ${referenceFileUploadPath} ${referenceFileIcaPath} --project-id ${projectId})
    echo "\${reference_file_upload_response}" > \${reference_file_upload_response_file}

    printf "[\${time_stamp}]: "
    printf "Extracting file_id of reference file '${reference_file}' from upload response... \n"
    reference_file_id=\$(cat \${reference_file_upload_response_file} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
  else
      printf "[\${time_stamp}]: "
      printf "Extracting file_id of reference file '${reference_file}' from get response... \n"
      reference_file_id=\$(cat \${get_reference_file_response_file} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
  fi

  printf "[\${time_stamp}]: "
  printf "Writing file data to existing data file...\n"

  printf "${referenceAnalysisDataCode}:\${reference_file_id}\n" >> ${dataFilePath}
  """
}

workflow {
  bamFilePairsChannel = Channel.fromFilePairs(params.bamFilePairsUploadPath, checkIfExists:true) { 
    file -> file.name.replaceAll(/.bam|.bai$/,'') 
  }
  dataFilePath = Channel.fromPath("data.txt", checkIfExists:true)
  uploadReferenceFile(dataFilePath, referenceFileId, referenceFileIcaPath)
}

