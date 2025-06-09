#!/usr/bin/env nextflow

process FASTQC {

    publishDir "s3://nextflow-data-bucket-1146/results/fastqc", mode: 'copy'

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path("${sampleId}_fastqc.zip"), emit: zip
    tuple val(sampleId), path("${sampleId}_fastqc.html"), emit: html

    script:
    """
    # Create the publishDir if it doesn't exist
    mkdir -p results/fastqc

    fastqc $reads
    """
}
