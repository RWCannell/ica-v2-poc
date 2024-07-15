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

  manifest_file="manifest.tsv"
  manifest_tab=\$(printf "\t")

  if ! [ -f \${manifest_file} ]; then
      echo "Manifest file does not exist. Creating one..."
      touch \${manifest_file}

      manifest_header_1="file_name"
      manifest_header_2="file_analysis_code"

      printf "[\${time_stamp}]: "
      printf "Writing file data to manifest...\n"

      printf "\${manifest_header_1} \${manifest_tab} \${manifest_header_2}\n" >> \${manifest_file}
      printf "\${bam_file_id} \${manifest_tab} ${bamAnalysisDataCode}:\${bam_file_id}\n" >> \${manifest_file}
  else
      printf "[\${time_stamp}]: "
      printf "Writing file data to existing manifest...\n"

      printf "\${bam_file_id} \${manifest_tab} ${bamAnalysisDataCode}:\${bam_file_id}\n" >> \${manifest_file}
  fi
  """
}

process uploadReferenceFile {
  debug true
  maxForks 3

  input:
  path(manifestFile)
  path(referenceFileUploadPath)
  path(referenceFileIcaPath)

  output:
  path "manifest.tsv", emit: manifestFile

  script:
  def reference_file = referenceFileUploadPath.baseName
  """
  #!/bin/bash
  time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
  reference_file_id=""
  reference_file_ica_path="/bam/$reference_file/"

  get_reference_file_response_file="get_reference_file_response.txt"
  reference_file_upload_response_file="reference_file_upload_response.txt"

  touch \${reference_file_response}

  printf "[\${time_stamp}]: "
  printf "Checking if reference file already exists in ICA...\n"
  get_reference_file_response=\$(icav2 projectdata get ${referenceFileIcaPath} --project-id ${projectId})
  echo "\${get_reference_file_response}" > \${get_reference_file_response_file}

  if grep -iq "No data found for path" \${get_reference_file_response_file};then
      printf "[\${time_stamp}]: "
      printf "Reference file not found in ICA. Uploading reference file '${reference_file}'... \n"
      reference_file_upload_response=\$(icav2 projectdata upload ${reference_file} ${referenceFileIcaPath} --project-id ${projectId})
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

  manifest_tab=\$(printf "\t")

  printf "[\${time_stamp}]: "
  printf "Writing file data to existing manifest...\n"

  printf "\${reference_file_id} \${manifest_tab} ${referenceAnalysisDataCode}:\${reference_file_id}\n" >> ${manifestFile}
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
