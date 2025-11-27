params.dea_sh = "${projectDir}/bin/05_dea/dea.sh"
params.findmarkers_r = "${projectDir}/bin/05_dea/findmarkers.R"
params.single_volcano_r = "${projectDir}/bin/05_dea/single_volcano.R"

process FindMarkers {
    tag "$batch"
    publishDir "$outdir/dea", mode: 'copy', overwrite: true
    
    input:
    path inputFile
    val cluster_key
    val ident_1
    val ident_2
    val outdir
    
    output:
    path("*.csv"), emit: csv
    
    script:
    """
    sh ${params.dea_sh} $inputFile "$ident_1" "$ident_2" "$cluster_key" \
    ${params.findmarkers_r} ${params.single_volcano_r} >> DEA.log 2>&1
    """
}