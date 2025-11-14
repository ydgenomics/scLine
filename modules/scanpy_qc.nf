process SCANPY_QC {
    tag "$sample_id"
    conda 'scanpy=1.9.3 pandas=2.0'
    publishDir "${params.outdir}/scanpy_qc", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_qc.h5ad"), emit: h5ad
    path "${sample_id}_qc_report.html",               emit: html

    script:
    """
    scanpy_qc.py \\
        --sample_id $sample_id \\
        --r1 $r1 \\
        --r2 $r2 \\
        --min_genes ${params.min_genes} \\
        --min_cells ${params.min_cells} \\
        --max_mt_percent ${params.max_mt_percent} \\
        --outdir ./
    """
}