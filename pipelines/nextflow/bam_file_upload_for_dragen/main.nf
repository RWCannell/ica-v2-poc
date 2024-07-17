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

process uploadBamFiles {
  debug true
  maxForks 3
  tag "$sampleId"

  input:
  tuple val(sampleId), file(bam)

  output:
  path "manifest.tsv", emit: manifestFile

  script:
  def bam_file = sampleId + ".bam"
  def bai_file = sampleId + ".bam.bai"
  """
  #!/bin/bash
  time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
  bam_file_id=""
  bai_file_id=""
  ica_upload_path="/bam/$sampleId/"

  printf "sampleId: $sampleId \n"
  printf "bam_file: $bam_file \n"
  printf "bai_file: $bai_file \n"

  bam_file_response="bam_file_response.txt"
  bai_file_response="bai_file_response.txt"

  touch \${bam_file_response}
  touch \${bai_file_response}

  printf "[\${time_stamp}]: "
  printf "Uploading .bam file '${bam_file}'... \n"
  bam_file_upload_response=\$(icav2 projectdata upload ${bam_file} \${ica_upload_path} --project-id ${projectId})
  echo "\${bam_file_upload_response}" > \${bam_file_response}

  printf "[\${time_stamp}]: "
  printf "Uploading .bam.bai file '${bai_file}'... \n"
  bai_file_upload_response=\$(icav2 projectdata upload ${bai_file} \${ica_upload_path} --project-id ${projectId})
  echo "\${bai_file_upload_response}" > \${bai_file_response}

  # id of file starts with 'fil.'
  bam_file_id=\$(cat \${bam_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
  bai_file_id=\$(cat \${bai_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

  bam_file_uploaded_file_data_response=\$(icav2 projectdata get \${bam_file_id})
  if [[ \$? != 0 ]]; then
      printf "Failed to fetch data about file with id '\${bam_file_id}'. \n"
      exit 1
  else
      bam_file_uploaded_file_path=\$(echo \${bam_file_uploaded_file_data_response} | jq -r ".details.path")
      printf "Path of uploaded .bam file is '\${bam_file_uploaded_file_path}'. \n"
  fi

  bai_file_uploaded_file_data_response=\$(icav2 projectdata get \${bai_file_id})
  if [[ \$? != 0 ]]; then
      printf "Failed to fetch data about file with id '\${bai_file_id}'. \n"
      exit 1
  else
      bai_file_uploaded_file_path=\$(echo \${bai_file_uploaded_file_data_response} | jq -r ".details.path")
      printf "Path of uploaded .bai file is '\${bai_file_uploaded_file_path}'. \n"
  fi

  data_file="data.txt"

  if ! [ -f \${data_file} ]; then
      echo "Data file does not exist. Creating one..."
      touch \${data_file}
  fi

  printf "[\${time_stamp}]: "
  printf "Writing file data to existing data file...\n"

  printf "sampleId:${sampleId}\n" >> \${data_file}
  printf "${bamAnalysisDataCode}:\${bam_file_id}\n" >> \${data_file}
  """
}

process uploadReferenceFile {
  debug true

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
  reference_file_id=${referenceFileId}
  reference_file_ica_path=${referenceFileIcaPath}

  get_reference_file_response_file="get_reference_file_response.txt"
  reference_file_upload_response_file="reference_file_upload_response.txt"

  touch \${reference_file_response}

  printf "[\${time_stamp}]: "
  printf "Checking if reference file id has been set in params.json...\n"
  if [ -z "${referenceFileId}" ]; then 
    printf "[\${time_stamp}]: "
    printf "Reference file id has been set in params.json. Getting reference file data JSON response...\n"
    get_reference_file_response=\$(icav2 projectdata get ${referenceFileId} --project-id ${projectId})
    echo "\${get_reference_file_response}" > \${get_reference_file_response_file}
  else 
    printf "[\${time_stamp}]: "
    printf "Reference file id has not been set in params.json. Checking if reference file exists in ICA path...\n"

    get_reference_file_response=\$(icav2 projectdata get ${referenceFileUploadPath} --project-id ${projectId})
    echo "\${get_reference_file_response}" > \${get_reference_file_response_file}
  fi

  if grep -iq "No data found for path" \${get_reference_file_response_file}; then
    printf "[\${time_stamp}]: "
    printf "Reference file not found in ICA. Uploading reference file '${reference_file}'... \n"
    reference_file_upload_response=\$(icav2 projectdata upload ${referenceFileUploadPath} \${reference_file_ica_path} --project-id ${projectId})
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

  printf "${referenceAnalysisDataCode}:\${reference_file_id}\n" >> ${dataFile}
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
  bamFilePairsChannel = Channel.fromFilePairs(params.bamFilePairsUploadPath, checkIfExists:true) { 
    file -> file.name.replaceAll(/.bam|.bai$/,'') 
  }
  uploadBamFiles(bamFilePairsChannel)
  uploadReferenceFile(uploadBamFiles.out.manifestFile)
}

// nextflow.enable.dsl=2
// projectId = params.projectId
// bamAnalysisDataCode = params.bamAnalysisDataCode
// referenceAnalysisDataCode = params.referenceAnalysisDataCode
// pipelineId = params.pipelineId
// pipelineCode = params.pipelineCode
// userReference = params.userReference
// storageSize = params.storageSize
// fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
// analysisStatusCheckInterval = params.analysisStatusCheckInterval
// bamFilesUploadPath = params.bamFilesUploadPath
// bamFilePairsUploadPath = params.bamFilePairsUploadPath
// referenceFileUploadPath = params.referenceFileUploadPath
// referenceFileIcaPath = params.referenceFileIcaPath
// localDownloadPath = params.localDownloadPath
// icaUploadPath = params.icaUploadPath

// process writeLine1 {
//   debug true

//   input:
//   val(projectId)

//   output:
//   path "project.txt", emit: projectFile

//   script:
//   def project_text = "Project id1 is " + projectId
//   """
//   #!/bin/bash
//   time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
//   project_file="project.txt"

//   printf "Printing project text to project file...\n"
//   printf "${project_text}\n" >> \${project_file}
//   """
// }

// process writeLine2 {
//   debug true

//   input:
//   path(projectFile)

//   output:
//   path "project.txt", emit: projectFile

//   script:
//   def project_text = "Project id2 is " + projectId
//   """
//   #!/bin/bash
//   time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")

//   printf "Printing line 2 of project text to project file...\n"
//   printf "${project_text}\n" >> ${projectFile}
//   """
// }
// workflow {
//   writeLine1(params.projectId)
    
//   writeLine2(writeLine1.out.projectFile)
// }
