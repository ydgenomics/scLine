params.metaneighbor_sh = "${projectDir}/bin/04_metaneighbor/metaneighbor.sh"
params.metaNeighbor_r = "${projectDir}/bin/04_metaneighbor/metaNeighbor.R"
params.sanky_plot_py = "${projectDir}/bin/04_metaneighbor/sanky_plot.py"

process MetaNeighbor {
    tag "$batch"
    publishDir "$outdir/metaneighbor", mode: 'copy', overwrite: true
    
    input:
    path inputFile
    val prefix
    val batch_key
    val cluster_key
    val seq
    val slimit
    val outdir
    
    output:
    path("*.csv"), emit: csv
    path("*.pdf"), emit: pdf
    path("*.html"), emit: html
    
    script:
    """
    sh ${params.metaneighbor_sh} "$inputFile" "$prefix" "$batch_key" \
    "$cluster_key" "$seq" "$slimit" \
    "${params.metaNeighbor_r}" "${params.sanky_plot_py}" >> METANEIGHBOR.log 2>&1
    """
}