#!/usr/bin/env nextflow
nextflow.enable.dsl=2
filePath = Channel.fromPath("ica_data_uploads/fasta/Citrobacter_freundii_4_7_47CFAA_uid46379/NZ_ADLG00000000.scaffold.fa/NZ_JH414877.fa", checkIfExists: true)
projectId = params.projectId
analysisDataCode = params.analysisDataCode

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

    echo "Uploading file '${fileName}' to project with id '${projectId}'..."

    touch ${fileName}.txt

    fileUploadResponse=\$(icav2 projectdata upload ${filePath} --project-id ${projectId})

    echo "Successfully uploaded file '${fileName}' to project with id '${projectId}'."

    echo "\${fileUploadResponse}" > ${fileName}.txt
    """
}

process constructFileReference {
    debug true
    
    input:
    path(fileUploadResponse)
    val(analysisDataCode)

    output:
    val(fileReference), emit: fileRef

    script:
    fileReference = ""
    """
    #!/bin/bash
    
    fileId=\$(cat ${fileUploadResponse} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    fileReference="${analysisDataCode}:\${fileId}"

    echo "File Reference: "

    echo \${fileReference}
    """
}

workflow {
    uploadFile(filePath, params.projectId)
    
    constructFileReference(uploadFile.out.fileUploadResponse.view(), params.analysisDataCode)
}
