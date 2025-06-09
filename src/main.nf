#!/usr/bin/env nextflow

// Load modules
modules_base = '/Users/isidorabeslic/Workspace/egfrKD-AstroRNAseq/modules'
include { DOWNLOAD_AND_UPLOAD_SRA } from file("${modules_base}/download_and_upload_sra.nf")
include { FASTP }                   from file("${modules_base}/modules/fastp.nf")
include { FASTQC }                  from file("${modules_base}/modules/fastqc.nf")
include { STAR_ALIGN }              from file("${modules_base}/modules/star_align.nf")
include { FEATURE_COUNTS }          from file("${modules_base}/modules/featurecounts.nf")
include { MULTIQC }                 from file("${modules_base}/modules/multiqc.nf")


workflow {

    // Read SRRs
    srr_ch = Channel.fromPath(params.srr_file)
        .flatMap { file -> file.readLines().findAll { it.trim() } }

    // Download SRA to S3
    DOWNLOAD_AND_UPLOAD_SRA(srr_ch)

    // Get fastq files
    fastq_ch = DOWNLOAD_AND_UPLOAD_SRA.out.fastq

    genome_dir_ch = Channel.value(params.genome_dir)
    gtf_ch        = Channel.value(params.gtf)

    // Run core modules
    FASTQC(fastq_ch)
    FASTP(fastq_ch)
    STAR_ALIGN(FASTP.out.trimmed_reads, genome_dir_ch)
    FEATURE_COUNTS(STAR_ALIGN.out.aligned_bam, gtf_ch)

    // MultiQC by sample
    all_qc_by_sample = FASTQC.out.html
        .mix(FASTQC.out.zip)
        .mix(FASTP.out.json)
        .mix(STAR_ALIGN.out.log)
        .groupTuple()
    MULTIQC(all_qc_by_sample)

}

workflow.onComplete {
    if (workflow.success) {
        log.info "Workflow finished successfully. Deleting workDir from S3..."

        def delete_cmd = "aws s3 rm s3://nextflow-data-bucket-1146/nextflow-workdir/ --recursive"
        log.info "Running command: ${delete_cmd}"
        def result = delete_cmd.execute().text
        log.info "S3 cleanup result:\n${result}"
    } else {
        log.warn "Workflow failed — workDir not deleted."
    }
}
