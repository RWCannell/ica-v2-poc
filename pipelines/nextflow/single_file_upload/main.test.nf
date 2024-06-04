#!/usr/bin/env nextflow
nextflow.enable.dsl=2
nextflow.preview.recursion=true
// filePath = Channel.fromPath("ica_data_uploads/fasta/Citrobacter_30_2_uid32453/NZ_ACDJ00000000.scaffold.fa/NZ_GG657383.fa", checkIfExists: true)
projectId = params.projectId
analysisDataCode = params.analysisDataCode
pipelineId = params.pipelineId
userReference = params.userReference
storageSize = params.storageSize

process getAnalysisId {
    debug true

    input:
    path(filePath)

    output:
    stdout

    script:
    """
    #!/bin/bash

    analysisResponse=\$(cat ${filePath})

    analysisId=\$(echo \${analysisResponse} | jq -r ".id")
    echo "\${analysisId}"
    """
}

process getStatus {
    debug true
    
    input:
    val analysisId

    output:
    stdout

    script:
    """
    #!/bin/bash

    echo "analysisId: ${analysisId}"  

    updatedAnalysisResponse=\$(icav2 projectanalyses get ${analysisId})
    
    analysisStatus=\$(echo \${updatedAnalysisResponse} | jq -r ".status")
    
    echo "\${analysisStatus}"
    """
}

process readStatus {
    debug true
    
    input:
    val(finalStatus)

    script:
    """
    #!/bin/bash

    echo "Final Status: ${finalStatus}"
    """
}

workflow {
    filePath = Channel.fromPath("pipeline_start_response.txt", checkIfExists: true)
    getAnalysisId(filePath)
    // getStatus.recurse(getAnalysisId.out.analysisId.view().collect()).times(5)
    getStatus(getAnalysisId.out.view().first())
    // .recurse(getAnalysisId.out.analysisId.view().first())
    // .times(5)
    // .until { it -> it.size() > 100 }
    readStatus(getStatus.out.view())
    // getStatus
    // .out
    // .finalStatus
    // .view()
}
