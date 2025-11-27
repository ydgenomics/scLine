#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = "${projectDir}/output"

params.outdir = "${projectDir}/output"

params.rawPaths = '/data/work/scline/input/sample1/sample1_raw|/data/work/scline/input/sample2/sample2_raw'
params.filterPaths = '/data/work/scline/input/sample1/sample1_filter|/data/work/scline/input/sample2/sample2_filter'
params.biosampleValues = 'group1|group1'
params.sampleValues = 'raw|filter'


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

params.soupxR = "${projectDir}/bin/test.sh"

process SOUPX_BATCH {
    tag "batch"
    
    publishDir "${params.outdir}/test", mode: 'copy'

    input:
    path rawDirs,     stageAs: 'raw/*'      // 整包软链到 raw/
    path filterDirs,  stageAs: 'filter/*'   // 整包软链到 filter/
    val sampleList

    output:
    path "test.txt"
    
    script:
    """
    biosamples=(${sampleList.join(' ')})

    # 2. 复制所有样本目录
    mkdir -p result
    cp -r "${sampleList[@]}" result/
    echo "${rawDirs.join('|')}" > test.txt
    # sh "${params.soupxR}" "${rawDirs.join('|')}"
    """
}

workflow {
    // 把所有整包通道一起发（单任务）
    ch_raws.view { "RAW: $it" }
    SOUPX_BATCH(
        ch_raws,
        ch_filters,
        ch_samples
    )
}