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
bamFilesUploadPath = params.bamFilesUploadPath
bamFilePairsUploadPath = params.bamFilePairsUploadPath
referenceFileUploadPath = params.referenceFileUploadPath
localDownloadPath = params.localDownloadPath
icaUploadPath = params.icaUploadPath

process uploadBamFiles {
  debug true
  tag "$sampleId"

  input:
  tuple val(sampleId), file(bam)

  script:
  """
  echo "sampleId: ${sampleId}"
  echo "bam: ${sampleId}.bam"
  echo "bam.bai: ${sampleId}.bam.bai"
  """
}

workflow {
  bamFilePairsChannel = Channel.fromFilePairs(params.bamFilePairsUploadPath, checkIfExists:true) { 
    file -> file.name.replaceAll(/.bam|.bai$/,'') 
  }
  uploadBamFiles(bamFilePairsChannel)
}
