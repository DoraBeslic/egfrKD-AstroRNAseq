#!/usr/bin/env nextflow

process FASTP {

    publishDir "s3://nextflow-data-bucket-1146/results/trimming", mode: 'copy'
    tag "$sampleId"

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path("${sampleId}_trimmed.fq.gz"), emit: trimmed_reads
    tuple val(sampleId), path("${sampleId}_fastp.html"), emit: html
    tuple val(sampleId), path("${sampleId}_fastp.json"), emit: json

    script:
    """
    # Create the publishDir if it doesn't exist
    mkdir -p results/trimming
    
    fastp \
        -i ${reads} \
        -o ${sampleId}_trimmed.fq.gz \
        --html ${sampleId}_fastp.html \
        --json ${sampleId}_fastp.json \
        --thread ${task.cpus}
    """
}
