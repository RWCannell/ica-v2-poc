#!/usr/bin/env nextflow
nextflow.enable.dsl=2
filePath = Channel.fromPath("ica_data_uploads/fasta/Citrobacter_30_2_uid32453/NZ_ACDJ00000000.scaffold.fa/NZ_GG657368.fa", checkIfExists: true)
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
    path "fileReference.txt", emit: fileRef

    script:
    fileReference = ""
    fileIsReady = false
    """
    #!/bin/bash
    
    fileId=\$(cat ${fileUploadResponse} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    fileResponse=\$(icav2 projectdata get \${fileId} --project-id ${projectId})

    intervalInSeconds=30
    statusCheckCount=0
    statusCheckLimit=10

    while true;
    do
        echo "Checking status of uploaded file..."
        ((statusCheckCount+=1))
        fileStatus=\$(echo \$fileResponse | jq -r ".details.status")
        if [[ "\${fileStatus}" == "AVAILABLE" ]]; then
            echo "File is AVAILABLE"

            fileReference="${analysisDataCode}:\${fileId}"

            echo "File Reference: "

            echo \${fileReference}

            touch fileReference.txt

            echo "\${fileReference}" > fileReference.txt
            break;
        else
            printf "File '\${fileId}' is still not AVAILABLE... \n"
        fi
        sleep \$intervalInSeconds;
    done
    """
}

process startAnalysis {
    debug true
    
    input:
    path(fileRef)

    script:
    """
    #!/bin/bash
    echo "Starting Nextflow analysis..."

    fileReference=\$(cat ${fileRef})

    echo "File Ref: \${fileReference}"

    icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference ${userReference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${fileReference}
    """
}

workflow {
    uploadFile(filePath, params.projectId)
    
    constructFileReference(uploadFile.out.fileUploadResponse.view(), params.analysisDataCode)

    startAnalysis(constructFileReference.out.fileRef.view())
}
