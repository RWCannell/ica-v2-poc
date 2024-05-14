#!/usr/bin/env nextflow
nextflow.enable.dsl=2

filePath = Channel.fromPath("ica_data_uploads/fasta/Yokenella_regensburgei_ATCC_43003_uid65133/NZ_AGCL00000000.scaffold.fa/NZ_JH417861.fa", checkIfExists: true)
projectId = params.projectId

process uploadFile {
    debug true
    input:
        path(filePath)
        val(projectId)
    output:
        val(fileUploadResponse)
        // stdout emit: fileId

    script:
    fileName = filePath.baseName
    fileUploadResponse = ''
    """
    #!/bin/bash

    echo "Uploading file '${fileName}' to project with id '${projectId}'..."

    fileUploadResponse=\$(icav2 projectdata upload ${filePath} --project-id ${projectId})

    echo "Successfully uploaded file '${fileName}' to project with id '${projectId}'."
    """
}

workflow {
    uploadFile(filePath, params.projectId)
}