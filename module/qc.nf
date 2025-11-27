#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.soupx_sh = "${projectDir}/bin/01_qc/soupx.sh"
params.run_soupx_r = "${projectDir}/bin/01_qc/run_soupx.R"

process SOUPX {
    tag "batch"
    
    // publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
     path raw, stageAs: 'raw/*'
     path filter, stageAs: 'filter/*'
     val sample
     val tfidfMin // int(default: 1)

    output:
    path("$sample"), emit: matrix
    path("*.txt"), emit: txt
    
    script:
    """
    sh ${params.soupx_sh} "$raw" "$filter" "$sample" 100 $tfidfMin "${params.run_soupx_r}" >> SOUPX.log 2>&1
    """
}



params.scrublet_sh = "${projectDir}/bin/01_qc/scrublet.sh"
params.run_scrublet_py = "${projectDir}/bin/01_qc/run_scrublet.py"

process SCRUBLET {
    tag "batch"
    
    publishDir "$outdir", mode: 'copy', overwrite: true

    input:
    path filterDirs,  stageAs: 'filter/*'
    path rhoDirs,  stageAs: 'rho/*'
    val sampleList
    val biosampleList
    val rlst
    val outdir

    output:
    path("qc"), emit: result
    path("qc/*/*.h5ad"), emit: h5ad
    path("qc/*/markers.csv/*.csv"), emit: csv
    
    script:
    """
    sh ${params.scrublet_sh} --biosample_value "${biosampleList.join(' ')}" --group_key 'sample' \
    --filter_list "$filterDirs" --splice_list "$filterDirs" \
    --unsplice_list "$filterDirs" --sample_list "${sampleList.join(' ')}" \
    --input_mingenes 100 --input_mincells 3 --mitogenes_csv "${sampleList[0]}" \
    --mito_threshold 5 --n_hvg 3000 --rlst "$rlst" \
    --doublet_threshold 2 --rhotxt_list "$rhoDirs" \
    --run_scrublet_py "${params.run_scrublet_py}" >> SCRUBLET.log 2>&1
    """
}