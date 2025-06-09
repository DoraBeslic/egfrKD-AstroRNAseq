#!/usr/bin/env nextflow

process FEATURE_COUNTS {
    tag "$sampleId"

    publishDir "s3://nextflow-data-bucket-1146/results/feature_counts", mode: 'copy'
    
    input:
    tuple val(sampleId), path(bam)
    path annotation

    output:
    tuple val(sampleId), path("${sampleId}.counts.txt"), emit: counts

    script:
    """
    # Create the publishDir if it doesn't exist
    mkdir -p results/feature_counts

    featureCounts -T ${task.cpus} \
        -a ${annotation} \
        -o ${sampleId}.counts.txt \
        $bam
    """
}
