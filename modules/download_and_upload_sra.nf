#!/usr/bin/env nextflow

process DOWNLOAD_AND_UPLOAD_SRA {
    tag "$srr_id"

    publishDir 's3://nextflow-data-bucket-1146/data/fastq/', mode: 'copy'

    input:
    val srr_id

    output:
    tuple val(srr_id), path("${srr_id}.fastq.gz"), emit: fastq

    script:
    """
    prefetch ${srr_id}
    fasterq-dump ${srr_id} -O .
    gzip ${srr_id}.fastq
    """
}
