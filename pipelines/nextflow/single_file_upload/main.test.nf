#!/usr/bin/env nextflow
nextflow.enable.dsl=2
filePath = Channel.fromPath("ica_data_uploads/mixed.txt", checkIfExists: true)
// filePath = Channel.fromPath("ica_data_uploads/fasta/Citrobacter_30_2_uid32453/NZ_ACDJ00000000.scaffold.fa/NZ_GG657383.fa", checkIfExists: true)
projectId = params.projectId
analysisDataCode = params.analysisDataCode
pipelineId = params.pipelineId
userReference = params.userReference
storageSize = params.storageSize

process uploadFile {
    debug true
    input:
    path(filePath)
    val(projectId)

    output:
    path "${fileName}.txt", emit: fileUploadResponse

    script:
    fileName = filePath.baseName
    """
    #!/bin/bash

    fileUploadResponse=\$(cat ${filePath})

    echo "\${fileUploadResponse}" > ${fileName}.txt
    """
}

process constructFileReference {
    debug true
    
    input:
    path(fileUploadResponse)
    val(analysisDataCode)

    output:
    path "fileReference.txt", emit: fileRef

    script:
    fileReference = ""
    """
    #!/bin/bash
    
    fileId=\$(cat ${fileUploadResponse} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    fileReference="${analysisDataCode}:\${fileId}"

    echo "File Reference: "

    echo \${fileReference}

    touch fileReference.txt

    echo "\${fileReference}" > fileReference.txt
    """
}

process startAnalysis {
    debug true
    
    input:
    val(fileRef)

    // output:
    // val(fileReference), emit: fileRef

    script:
    """
    #!/bin/bash
    echo "Starting Nextflow analysis..."

    fileRef=\$(cat ${fileRef})

    echo "File Ref: \${fileRef}"

    icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference ${userReference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${fileRef}
    """
}

workflow {
    uploadFile(filePath, params.projectId)
    
    constructFileReference(uploadFile.out.fileUploadResponse.view(), params.analysisDataCode)

    startAnalysis(constructFileReference.out.fileRef.view())
}
