#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
pipeline_id="858708e9-6cdd-42ae-ad97-8b86ee07cdc5"
user_reference="dragen_germline_whole_genome"
storage_size="Medium"
local_download_path="${HOME}/Documents/ica_data_downloads/"
interval_in_seconds=120
sample_id="SRR622459"
output_directory="/"
read1_file_id="fil.68c7bd2bd2dd423d1eeb08dd49010d1c"
read2_file_id="fil.0a7f991a58154b69f17c08dd4847c564"
fastq_list_file_id="fil.050ab2d7e51a4716ebe308dd499c30b4"
read1_file_name="SRR622459_1.fastq.gz"
read2_file_name="SRR622459_2.fastq.gz"
read1_analysis_code="read1:fil.68c7bd2bd2dd423d1eeb08dd49010d1c"
read2_analysis_code="read2:fil.0a7f991a58154b69f17c08dd4847c564"
reference_analysis_code="ref_tar:fil.2e3fd8d802ee4963da2208dc484ea8f0"
hash_table_config_file="/scratch/reference/hash_table.cfg"

analysis_response=$(icav2 projectpipelines start nextflow $pipeline_id \
--user-reference $user_reference \
--project-id $project_id \
--storage-size $storage_size \
--input $reference_analysis_code \
--input fastqs:$read1_file_id,$read2_file_id \
--input fastq_list:$fastq_list_file_id \
--parameters enable_map_align:true \
--parameters enable_map_align_output:true \
--parameters output_format:CRAM \
--parameters enable_duplicate_marking:true \
--parameters enable_variant_caller:true \
--parameters vc_emit_ref_confidence:GVCF \
--parameters enable_cnv:true \
--parameters enable_sv:true \
--parameters repeat_genotype_enable:false \
--parameters enable_hla:false \
--parameters enable_variant_annotation:false \
--parameters output_file_prefix:"$sample_id")

echo $analysis_response


