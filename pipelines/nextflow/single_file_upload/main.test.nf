#!/usr/bin/env nextflow
nextflow.enable.dsl=2
nextflow.preview.recursion=true
filePath = Channel.fromPath("pipeline_start_response.txt", checkIfExists: true)
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
    val analysisId, emit: analysisId

    script:
    analysisId = ""
    """
    #!/bin/bash

    analysisResponse=\$(cat ${filePath})

    echo "analysisResponse: \${analysisResponse}" 

    analysisId=\$(echo \${analysisResponse} | jq -r ".id")
    echo "analysisId: \$analysisId"
    """
}

process readStatus {
    debug true
    
    input:
    val(analysisId)

    output:
    val(finalStatus), emit: finalStatus

    script:
    finalStatus = "REQUESTED"
    """
    #!/bin/bash

    echo "analysisId: ${analysisId}"  

    updatedAnalysisResponse=\$(icav2 projectanalyses get ${analysisId})
    
    analysisStatus=\$(echo \${updatedAnalysisResponse} | jq -r ".status")
    
    finalStatus=\${analysisStatus}

    echo "Final Status: \${analysisStatus}"
    """
}

workflow {
    // ch_collected = filePath.collect()
    getAnalysisId(filePath)
    // readStatus.recurse(getAnalysisId.out.analysisId.view().collect()).times(5)
    readStatus(getAnalysisId.out.analysisId.view())
    // .recurse(getAnalysisId.out.analysisId.view().first())
    // .times(5)
    // .until { it -> it.size() > 100 }

    // readStatus
    // .out
    // .view(it)
}
