process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    conda 'multiqc=1.19'      // 与 fastqc 同频道，兼容性好
    container 'ewels/multiqc:1.19'  // docker/singularity 备用

    input:
    path('*_fastqc.zip') qc_zips   // 通配任意 FastQC zip，也可以传 file('*') 全收

    output:
    path('multiqc_report.html') , emit: html
    path('multiqc_data')        , emit: data_dir

    script:
    """
    multiqc . -n multiqc_report.html
    """
}