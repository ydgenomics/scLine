process FASTQC {
    tag "$sample_id"
    conda 'fastqc=0.12.1'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    """
    fastqc -t ${task.cpus} $r1 $r2
    """
}