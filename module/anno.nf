// anno

params.anno_sh = "${projectDir}/bin/02_anno/anno.sh"
params.anno_singler_r = "${projectDir}/bin/02_anno/anno_singler.R"
params.plot_r = "${projectDir}/bin/02_anno/plot.R"
params.anno_sctype_r = "${projectDir}/bin/02_anno/anno_sctype.R"

process ANNOR {
    tag "batch"
    
    publishDir "$outdir/anno", mode: 'copy', overwrite: true

    input:
    path input_query_rds, stageAs: 'input/*'
    path rdsORcsv, stageAs: 'input/*'
    val cluster_key
    val umap_name
    val ref_cluster_key
    val outdir

    output:
    path("*.pdf"), emit: pdf
    path("*.csv"), emit: csv
    path("*anno.csv"), emit: anno
    path("*.rds"), emit: rds
    
    script:
    """
    sh ${params.anno_sh} "$input_query_rds" "$cluster_key" "$rdsORcsv" \
    "$umap_name" "$ref_cluster_key" "${params.anno_singler_r}" \
    "${params.plot_r}" "${params.anno_sctype_r}" >> ANNOR.log 2>&1
    """
}
