#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = "${projectDir}/output"

params.rawPaths = '/data/work/scline/input/sample1/sample1_raw|/data/work/scline/input/sample2/sample2_raw'
params.filterPaths = '/data/work/scline/input/sample1/sample1_filter|/data/work/scline/input/sample2/sample2_filter'
params.biosampleValues = 'group1|group2'
params.sampleValues = 'sample1|sample2'


println "rawPaths: ${params.rawPaths}"
println "filterPaths: ${params.filterPaths}"
println "sampleValues: ${params.sampleValues}"


def rawList = params.rawPaths.split('\\|')
def filterList = params.filterPaths.split('\\|')
def biosampleList = params.biosampleValues.split('\\|')
def sampleList = params.sampleValues.split('\\|')

ch_raws     = Channel.value(rawList.collect{ p -> file(p) })   // List[Path]
ch_filters  = Channel.value(filterList.collect{ p -> file(p) })
ch_samples  = Channel.value(sampleList)                       // List<String> 即可
ch_biosamples  = Channel.value(biosampleList)                       // List<String> 即可

params.soupx_sh = "${projectDir}/bin/01_qc/soupx.sh"
params.run_soupx_r = "${projectDir}/bin/01_qc/run_soupx.R"

params.scrublet_sh = "${projectDir}/bin/01_qc/scrublet.sh"
params.run_scrublet_py = "${projectDir}/bin/01_qc/run_scrublet.py"

process QC {
    tag "batch"
    
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path rawDirs,     stageAs: 'raw/*'      // 整包软链到 raw/
    path filterDirs,  stageAs: 'filter/*'   // 整包软链到 filter/
    val sampleList
    val biosampleList

    output:
    path("result"), emit: result
    path("result/*/*.h5ad"), emit: h5ad
    
    script:
    """
    source /opt/software/miniconda3/bin/activate
    conda activate scline
    
    sh ${params.soupx_sh} "${rawDirs.join('|')}" "${filterDirs.join('|')}" "${sampleList.join('|')}" 100 1 "${params.run_soupx_r}" >> qc.log 2>&1
    
    bash ${params.scrublet_sh} --biosample_value "${biosampleList.join('|')}" --group_key 'sample' \
    --filter_list "${sampleList.join('|')}" --splice_list "${sampleList.join('|')}" \
    --unsplice_list "${sampleList.join('|')}" --sample_list "${sampleList.join('|')}" \
    --input_mingenes 100 --input_mincells 3 --mitogenes_csv "${sampleList[0]}" \
    --mito_threshold 5 --n_hvg 3000 --rlst '0.2,0.5,0.8,1.0' \
    --doublet_threshold 2 --rhotxt_list "${sampleList[0]}" \
    --run_scrublet_py "${params.run_scrublet_py}" >> qc.log 2>&1
    """
}

workflow {
    // 把所有整包通道一起发（单任务）
    ch_raws.view { "RAW: $it" }
    QC(
        ch_raws,
        ch_filters,
        ch_samples,
        ch_biosamples
    )
}