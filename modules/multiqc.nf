#!/usr/bin/env nextflow

process MULTIQC {

    publishDir 's3://nextflow-data-bucket-1146/results/multiqc', mode: 'copy'

    input:
    tuple val(sampleId), path(qc_files)

    output:
    path("${sampleId}_multiqc.html"), emit: report
    path("${sampleId}_multiqc_data"), emit: data

    script:
    """
    mkdir -p results/multiqc

    # Copy QC files into a temporary directory to avoid conflicts
    mkdir -p work_qc
    cp ${qc_files.join(' ')} work_qc/

    # Run MultiQC on the collected files
    multiqc work_qc -n "${sampleId}_multiqc.html" --outdir .
    """
}
