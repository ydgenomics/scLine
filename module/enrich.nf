params.build_orgdb_r = "${projectDir}/bin/06_enrich/build_orgdb.R"
params.ko_json="${projectDir}/bin/06_enrich/ko00001.json"
params.go_obo="${projectDir}/bin/06_enrich/go_obo_result.csv"

process makeOrgPackage {
    tag "$batch"
    publishDir "$outdir/enrich", mode: 'copy', overwrite: true
    
    input:
    path emapper_xlsx, stageAs: 'input/*'
    val genus
    val species
    val outdir
    
    output:
    path("org.*.eg.db"), emit: orgdb
    path("kegg_info.RData"), emit: kegg_info_RData
    
    script:
    """
    Rscript ${params.build_orgdb_r} --emapper_xlsx "$emapper_xlsx" \
    --ko_json "${params.ko_json}" --go_obo "${params.go_obo}" \
    --taxid "1111" --genus "$genus" --species "$species" >> makeOrgPackage.log 2>&1
    """
}

params.enrich_r = "${projectDir}/bin/06_enrich/enrich.R"

process clusterProfiler {
    tag "$batch"
    publishDir "$outdir/enrich", mode: 'copy', overwrite: true
    
    input:
    path db, stageAs: 'input/*'
    path gene_csv, stageAs: 'input/*'
    path kegg_info_RData, stageAs: 'input/*'
    val genus
    val species
    val minp
    val outdir
    
    output:
    path("*.pdf"), emit: pdf
    path("*.txt"), emit: txt
    
    script:
    """
    Rscript ${params.enrich_r} --gene_csv "$gene_csv" \
    --kegg_info_RData "$kegg_info_RData" --db "$db" \
    --minp $minp --genus $genus --species $species >> clusterProfiler.log 2>&1
    """
}