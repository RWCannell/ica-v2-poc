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

process startAnalysis {
    debug true
    
    input:
    path(dataFile)

    output:
    path "output.txt"

    script:

    """
    #!/bin/bash

    sample_id=\$(cat ${dataFile} | grep -E "sampleId")
    read1_analysis_code=\$(cat ${dataFile} | grep -E "read1")
    read2_analysis_code=\$(cat ${dataFile} | grep -E "read2")
    reference_analysis_code=\$(cat ${dataFile} | grep -E "ref_tar")
    output_directory="/output/\${sample_id}/"

    output="output.txt"
    touch \${output}

    printf "read1_analysis_code - \${read1_analysis_code}\n"
    printf "read2_analysis_code - \${read2_analysis_code}\n"
    printf "reference_analysis_code - \${reference_analysis_code}\n"
    printf "output_directory - \${output_directory}\n"

    printf "read1_analysis_code - \${read1_analysis_code}\n" >> \${output}
    printf "read2_analysis_code - \${read2_analysis_code}\n" >> \${output}
    printf "reference_analysis_code - \${reference_analysis_code}\n" >> \${output}
    printf "output_directory - \${output_directory}\n" >> \${output}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "Starting Nextflow analysis..."
    """
}

workflow {
    // fastqFilePairs = Channel.fromFilePairs(readsPairFilesUploadPath, checkIfExists:true)
    // referenceFilePath = Channel.fromPath(params.referenceFileUploadPath, checkIfExists: true)

    // uploadFastqFilePairs(fastqFilePairs, params.projectId)
    // uploadReferenceFile(uploadFastqFilePairs.out.dataFile, referenceFilePath)
    data = Channel.fromPath("/Users/regancannell/Documents/RWCannell/ica-v2-poc/pipelines/nextflow/fastq_pairwise_upload_for_dragen/data_file.txt", checkIfExists: true)
    startAnalysis(data)
}
