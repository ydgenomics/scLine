// ----- concat -----

params.concat_py = "${projectDir}/bin/03_integrate/concat.py"

process CONCAT {
    tag "batch"
    
    // publishDir "${params.outdir}/integrate", mode: 'copy'

    input:
    path rawDirs,     stageAs: 'raw/*'
    path annoDirs,     stageAs: 'csv/*'
    val biosampleList
    val prefix

    output:
    path("*.h5ad"), emit: h5ad
    
    script:
    """
    python ${params.concat_py} "$rawDirs" "$annoDirs" "${biosampleList.join(' ')}" "$prefix" "biosample" >> CONCAT.log 2>&1
    """
}

params.concat2_py = "${projectDir}/bin/03_integrate/concat2.py"

process CONCAT2 {
    tag "batch"
    
    // publishDir "${params.outdir}/integrate", mode: 'copy'

    input:
    path rawDirs,     stageAs: 'raw/*'
    path annoDirs,     stageAs: 'csv/*'
    val biosampleList
    val prefix

    output:
    path("*.h5ad"), emit: h5ad
    
    script:
    """
    python ${params.concat2_py} "$rawDirs" "$annoDirs" "${biosampleList.join(' ')}" "$prefix" "biosample" >> CONCAT2.log 2>&1
    """
}


// ----- scvi -----

params.scVI_integration_py = "${projectDir}/bin/03_integrate/scVI_integration.py"

process SCVI {
    tag "batch"
    
    // publishDir "$outdir/integrate", mode: 'copy'

    input:
    path input_file,     stageAs: 'input/*'
    val prefix
    val resolutionSet
    val outdir

    output:
    path("*.h5ad"), emit: h5ad
    path("*.pdf"), emit: pdf
    
    script:
    """
    python ${params.scVI_integration_py} "$input_file" \
    "${prefix}_scVI_integrated.h5ad" "${prefix}_scVI_integrated_UMAP.pdf" \
    --batch_key "biosample" --sample_key "sample" --cluster_key "celltype" --resolution_set $resolutionSet >> SCVI.log 2>&1
    """
}


params.harmony_integration_py = "${projectDir}/bin/03_integrate/harmony_integration.py"

process HARMONY {
    tag "batch"
    
    // publishDir "$outdir/integrate", mode: 'copy'

    input:
    path input_file,     stageAs: 'input/*'
    val prefix
    val resolutionSet
    val outdir
    
    output:
    path("*.h5ad"), emit: h5ad
    path("*.pdf"), emit: pdf
    
    script:
    """
    python ${params.harmony_integration_py} "$input_file" \
    "${prefix}_harmony_integrated.h5ad" "${prefix}_harmony_integrated_UMAP.pdf" \
    --batch_key "biosample" --sample_key "sample" --cluster_key "celltype" --resolution_set $resolutionSet >> HARMONY.log 2>&1
    """
}

params.unintegration_py = "${projectDir}/bin/03_integrate/unintegration.py"

process UNINTEGRATION {
    tag "batch"
    
    // publishDir "$outdir/integrate", mode: 'copy'

    input:
    path input_file,     stageAs: 'input/*'
    val prefix
    val resolutionSet
    val outdir

    output:
    path("*.h5ad"), emit: h5ad
    path("*.pdf"), emit: pdf
    
    script:
    """
    python ${params.unintegration_py} "$input_file" \
    "${prefix}_unintegration_integrated.h5ad" "${prefix}_unintegration_integrated_UMAP.pdf" \
    --batch_key "biosample" --sample_key "sample" --cluster_key "celltype" --resolution_set $resolutionSet >> UNINTEGRATION.log 2>&1
    """
}

params.CCA_integration_r = "${projectDir}/bin/03_integrate/SCTransform.CCA_integration.R"

process SCT_CCA {
    tag "batch"
    
    // publishDir "$outdir/integrate", mode: 'copy'

    input:
    path input_file,     stageAs: 'input/*'
    val prefix
    val resolutionSet
    val outdir

    output:
    path("*.rds"), emit: rds
    path("*.pdf"), emit: pdf
    
    script:
    """
    Rscript ${params.CCA_integration_r} --input_rds "$input_file" \
    --out_rds "${prefix}_SCTransform.CCA_integrated.rds" \
    --out_UMAP "${prefix}_SCTransform.CCA_integrated_UMAP.pdf" --batch_key "biosample" \
    --sample_key "sample" --cluster_key "celltype" --resolution_set $resolutionSet >> SCT_CCA.log 2>&1
    """
}

params.harmony_integration_r = "${projectDir}/bin/03_integrate/SCTransform.CCA_integration.R"

process SCT_HARMONY {
    tag "batch"
    
    // publishDir "$outdir/integrate", mode: 'copy'

    input:
    path input_file,     stageAs: 'input/*'
    val prefix
    val resolutionSet
    val outdir

    output:
    path("*.rds"), emit: rds
    path("*.pdf"), emit: pdf
    
    script:
    """
    # sct.harmony
    Rscript ${params.harmony_integration_r} --input_rds "$input_file" \
    --out_rds "${prefix}_SCTransform.harmony_integrated.rds" \
    --out_UMAP "${prefix}_SCTransform.harmony_integrated_UMAP.pdf" \
    --batch_key "biosample" --sample_key "sample" --cluster_key "celltype" --resolution_set $resolutionSet >> SCT_HARMONY.log 2>&1
    """
}


params.rliger_r = "${projectDir}/bin/03_integrate/rliger.INMF_integration.R"

process RLIGER {
    tag "batch"
    
    // publishDir "$outdir/integrate", mode: 'copy'

    input:
    path input_file,     stageAs: 'input/*'
    val prefix
    val resolutionSet
    val outdir

    output:
    path("*.rds"), emit: rds
    path("*.pdf"), emit: pdf
    
    script:
    """
    Rscript ${params.rliger_r} --input_rds "$input_file" --out_rds "${prefix}_rliger.INMF_integrated.rds" \
    --out_UMAP "${prefix}_rliger.INMF_integrated_UMAP.pdf" --batch_key "biosample" \
    --sample_key "sample" --cluster_key "celltype" --resolution_set $resolutionSet >> RLIGER.log 2>&1
    """
}


params.dealplus_sh = "${projectDir}/bin/03_integrate/dealplus.sh"
params.dealplus_r = "${projectDir}/bin/03_integrate/dealplus.R"
params.dealplus_py = "${projectDir}/bin/03_integrate/dealplus.py"

process DEALPLUS {
    tag "batch"
    
    publishDir "$outdir/integrate", mode: 'copy'

    input:
    path input_file,     stageAs: 'input/*'
    val other1_key
    val other2_key
    val outdir

    output:
    path("*_dealplus*"), emit: result
    path("*.pdf"), emit: pdf
    
    script:
    """
    sh ${params.dealplus_sh} $input_file $other1_key $other2_key \
    ${params.dealplus_r} ${params.dealplus_py} >> RLIGER.log 2>&1
    """
}

params.methods_file = "harmony rliger.INMF SCTransform.CCA SCTransform.harmony" // 无意义
params.pcas_file = "X_pca_harmony X_inmfnorm X_pca X_harmony" // 无意义
params.batch_key = "biosample"
params.n_jobs = 8 // match with the number of cpu in nextflow.config
params.scIB_py="${projectDir}/bin/03_integrate/scIB.py"

process SCIB {
    tag "batch"
    
    publishDir "$outdir/integrate/scib", mode: 'copy', overwrite: true

    input:
    path unint_h5ad, stageAs: 'input/*'
    path int_list, stageAs: 'input/*'
    val label_key
    val prefix
    val outdir

    output:
    path("*.h5ad"), emit: h5ad
    path("*_integrated*"), emit: best_h5ad
    path("*.csv"), emit: csv
    path("*.pdf"), emit: pdf
    
    script:
    """
    python ${params.scIB_py} \
    --unintegrated_h5ad "$unint_h5ad" \
    --integrated_file "$int_list" \
    --methods_file "${params.methods_file}" \
    --pcas_file "${params.pcas_file}" \
    --deals_file "N N N N N N" \
    --tests_file "true true true true true true true true true true" \
    --batch_key "biosample" \
    --label_key "$label_key" \
    --n_jobs "${params.n_jobs}" \
    --prefix "$prefix" >> SCIB.log 2>&1
    """
}

/*
workflow {
    CONCAT(ch_raws, ch_csvs, ch_biosamples)
    CONCAT.out.view()
    // SCVI(CONCAT.out.h5ad)
    HARMONY(CONCAT.out.h5ad)
    UNINTEGRATION(CONCAT.out.h5ad) 
    CONVERT(CONCAT.out.h5ad)
    SCT_CCA(CONVERT.out.result)
    SCT_HARMONY(CONVERT.out.result)
    RLIGER(CONVERT.out.result)
    
    ch_convert_input = RLIGER.out.rds
        .mix(SCT_CCA.out.rds)
        .mix(SCT_HARMONY.out.rds)
        .view { "CONVERT_INT input: $it" }
    CONVERT_INT(ch_convert_input)
    ch_final_input = HARMONY.out.h5ad
        .mix(CONVERT_INT.out.result.collect().flatten())
        .collect()
        .view { "Final input list (${it.size()} files): $it" }
    SCIB(UNINTEGRATION.out.h5ad, ch_final_input)
}
*/