params.convert_sh = "${projectDir}/bin/convert/convert.sh"
params.convert_rdsAh5ad2_r = "${projectDir}/bin/convert/convert_rdsAh5ad2.R"
params.deal_layers_ydgenomics_py = "${projectDir}/bin/convert/deal_layers_ydgenomics.py"
params.convert_rdsAh5ad_r = "${projectDir}/bin/convert/convert_rdsAh5ad.R"

process CONVERT_0 {
    tag "batch"
    
    publishDir "$outdir", mode: 'copy'

    input:
    path inputFile, stageAs: 'input/*'
    val assay
    val python_env
    val outdir

    output:
    path("*r*"), emit: result
    
    script:
    """
    source /opt/software/miniconda3/bin/activate
    conda activate convert
    sh ${params.convert_sh} $inputFile "multi_convert" "$assay" ${params.convert_rdsAh5ad2_r} \
    ${params.deal_layers_ydgenomics_py} ${params.convert_rdsAh5ad_r} "$python_env" >> CONVERT_0.log 2>&1
    """
}

process CONVERT_1 {
    tag "batch"
    
    // publishDir "${params.outdir}/convert", mode: 'copy'

    input:
    path inputFile, stageAs: 'input/*'

    output:
    path("*r*"), emit: result
    
    script:
    """
    source /opt/software/miniconda3/bin/activate
    conda activate convert
    sh ${params.convert_sh} $inputFile "multi_convert" "RNA" ${params.convert_rdsAh5ad2_r} \
    ${params.deal_layers_ydgenomics_py} ${params.convert_rdsAh5ad_r} ${params.python_env} >> CONVERT_1.log 2>&1
    """
}

process CONVERT_2 {
    tag "batch"
    
    // publishDir "${params.outdir}/convert", mode: 'copy'

    input:
    path inputFile, stageAs: 'input/*'

    output:
    path("*r*"), emit: result
    
    script:
    """
    source /opt/software/miniconda3/bin/activate
    conda activate convert
    sh ${params.convert_sh} $inputFile "multi_convert" "SCT" ${params.convert_rdsAh5ad2_r} \
    ${params.deal_layers_ydgenomics_py} ${params.convert_rdsAh5ad_r} ${params.python_env} >> CONVERT_2.log 2>&1
    """
}