#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = "${projectDir}/output"

// params.rawPaths = '/data/work/scline/example/SRR31330969_raw_sampled|/data/work/scline/example/SRR31330970_raw_sampled'
// params.filterPaths = '/data/work/scline/example/SRR31330969_filter_sampled|/data/work/scline/example/SRR31330970_filter_sampled'
// params.biosampleValues = 'group1|group2'
// params.sampleValues = 'SRR31330969|SRR31330970'
// params.rdsORcsv = "/data/work/scline/example/markergene.csv"

params.rawPaths = '/data/work/scline/input/sample1/sample1_raw|/data/work/scline/input/sample2/sample2_raw'
params.filterPaths = '/data/work/scline/input/sample1/sample1_filter|/data/work/scline/input/sample2/sample2_filter'
params.biosampleValues = 'group1|group2'
params.sampleValues = 'sample1|sample2'
params.rdsORcsv = "/data/work/scline/input/markergene.csv"

def rawList = params.rawPaths.split('\\|')
def filterList = params.filterPaths.split('\\|')
def biosampleList = params.biosampleValues.split('\\|')
def sampleList = params.sampleValues.split('\\|')

ch_raws     = Channel.value(rawList.collect{ p -> file(p) })   // List[Path]
ch_filters  = Channel.value(filterList.collect{ p -> file(p) })
ch_samples  = Channel.value(sampleList)                       // List<String> 即可
ch_biosamples  = Channel.value(biosampleList)                       // List<String> 即可

// qc

params.soupx_sh = "${projectDir}/bin/01_qc/soupx.sh"
params.run_soupx_r = "${projectDir}/bin/01_qc/run_soupx.R"

params.scrublet_sh = "${projectDir}/bin/01_qc/scrublet.sh"
params.run_scrublet_py = "${projectDir}/bin/01_qc/run_scrublet.py"


process QC {
    tag "batch"
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path rawDirs,     stageAs: 'raw/*'      // 整包软链到 raw/
    path filterDirs,  stageAs: 'filter/*'   // 整包软链到 filter/
    val sampleList
    val biosampleList

    output:
    path("qc"), emit: result
    path("qc/*/*.h5ad"), emit: h5ad
    
    script:
    """
    source /opt/software/miniconda3/bin/activate
    conda activate scline
    
    sh ${params.soupx_sh} "$rawDirs" "$filterDirs" "${sampleList.join(' ')}" 100 1 "${params.run_soupx_r}" >> QC.log 2>&1
    
    sh ${params.scrublet_sh} --biosample_value "${biosampleList.join(' ')}" --group_key 'sample' \
    --filter_list "${sampleList.join(' ')}" --splice_list "${sampleList.join(' ')}" \
    --unsplice_list "${sampleList.join(' ')}" --sample_list "${sampleList.join(' ')}" \
    --input_mingenes 100 --input_mincells 3 --mitogenes_csv "${sampleList[0]}" \
    --mito_threshold 5 --n_hvg 3000 --rlst '0.2,0.5,0.8,1.0' \
    --doublet_threshold 2 --rhotxt_list "${sampleList[0]}" \
    --run_scrublet_py "${params.run_scrublet_py}" >> QC.log 2>&1
    """
}

// convert

params.convert_sh = "${projectDir}/bin/convert/convert.sh"
params.convert_rdsAh5ad2_r = "${projectDir}/bin/convert/convert_rdsAh5ad2.R"
params.deal_layers_ydgenomics_py = "${projectDir}/bin/convert/deal_layers_ydgenomics.py"
params.convert_rdsAh5ad_r = "${projectDir}/bin/convert/convert_rdsAh5ad.R"
params.python_env = "/opt/software/miniconda3/envs/scline/bin/python"

process CONVERT {
    tag "batch"
    
    // publishDir "${params.outdir}/convert", mode: 'copy', overwrite: true

    input:
    path inputFile, stageAs: 'input/*'

    output:
    path("*r*"), emit: result
    
    script:
    """
    source /opt/software/miniconda3/bin/activate
    conda activate convert
    sh ${params.convert_sh} $inputFile "multi_convert" "RNA" ${params.convert_rdsAh5ad2_r} \
    ${params.deal_layers_ydgenomics_py} ${params.convert_rdsAh5ad_r} ${params.python_env}  >> CONVERT.log 2>&1
    """
}

// anno

params.anno_sh = "${projectDir}/bin/02_anno/anno.sh"
params.cluster_key = "leiden_res_0.20"
params.umap_name = "Xumap_"
params.ref_cluster_key = "sctype"
params.anno_singler_r = "${projectDir}/bin/02_anno/anno_singler.R"
params.plot_r = "${projectDir}/bin/02_anno/plot.R"
params.anno_sctype_r = "${projectDir}/bin/02_anno/anno_sctype.R"

process ANNO {
    tag "batch"
    
    publishDir "${params.outdir}/anno", mode: 'copy', overwrite: true

    input:
    tuple path(input_query_rds, stageAs: 'input/*'), path(rdsORcsv, stageAs: 'input/*')

    output:
    path("*.pdf"), emit: pdf
    path("*.csv"), emit: csv
    path("*anno.csv"), emit: anno
    path("*.rds"), emit: rds
    
    script:
    """
    source /opt/software/miniconda3/bin/activate
    conda activate scline
    sh ${params.anno_sh} "$input_query_rds" "${params.cluster_key}" "$rdsORcsv" \
    "${params.umap_name}" "${params.ref_cluster_key}" "${params.anno_singler_r}" \
    "${params.plot_r}" "${params.anno_sctype_r}" >> ANNO.log 2>&1
    """
}

workflow {
    ch_raws.view { "RAW: $it" }
    QC(ch_raws,ch_filters,ch_samples,ch_biosamples)
    QC.out.h5ad.view { "QC h5ad output: $it" }
    
    /*
    QC.out.h5ad
        .flatten()
        .view { "inputFile: $it" }
        .set { ch_single_files }
    */
    
    CONVERT(QC.out.h5ad.flatten())
    CONVERT.out.result.view { "CONVERT result: $it" }
    ch_reference = Channel.fromPath(params.rdsORcsv)
    CONVERT.out.result
        .combine(ch_reference)
        .view { "query: $it[0], ref: $it[1]" }
        .set { ch_anno_input }
    ANNO(ch_anno_input)
}