#!/usr/bin/env nextflow

process STAR_ALIGN {

    publishDir "s3://nextflow-data-bucket-1146/results/star_align", mode: 'copy'
    tag "$sampleId"

    input:
    tuple val(sampleId), path(reads)
    path genomeDir

    output:
    tuple val(sampleId), path("${sampleId}_Aligned.sortedByCoord.out.bam"), emit: aligned_bam
    tuple val(sampleId), path("${sampleId}_Log.out"), emit: log

    script:
    """
    STAR \
         --runThreadN ${task.cpus} \
         --genomeDir ${genomeDir} \
         --readFilesIn ${reads} \
         --readFilesCommand gunzip -c \
         --outFileNamePrefix ${sampleId}_ \
         --outSAMtype BAM SortedByCoordinate
    """
}
