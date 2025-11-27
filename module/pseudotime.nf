params.cytotrace_py = "${projectDir}/bin/07_pseudotime/CytoTRACE.py"

process CYTOTRACE {
    tag "$batch"
    publishDir "$outdir/cytotrace", mode: 'copy', overwrite: true
    
    input:
    path inputFile, stageAs: 'input/*'
    val batch_key
    val cluster_key
    val cluster_value
    val outdir
    
    output:
    path("*.h5ad"), emit: h5ad
    path("figures/*.pdf"), emit: pdf1
    path("*.pdf"), emit: pdf2
    
    script:
    """
    python ${params.cytotrace_py} --input_h5ad "$inputFile" \
    --batch_key "$batch_key" --cluster_key "$cluster_key" \
    --cluster_value "$cluster_value" >> CYTOTRACE.log 2>&1
    """
}