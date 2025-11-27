#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = "${projectDir}/output"

params.rawPaths = '/data/work/scline/input/sample1/raw|/data/work/scline/input/sample2/raw'
params.filterPaths = '/data/work/scline/input/sample1/filter|/data/work/scline/input/sample2/filter'
params.splicePaths = '/data/work/scline/input/sample1/splice|/data/work/scline/input/sample2/splice'
params.unsplicePaths = '/data/work/scline/input/sample1/unsplice|/data/work/scline/input/sample2/unsplice'
params.biosampleValues = 'group1|group1'
params.sampleValues = 'sample1|sample2'


println "rawPaths: ${params.rawPaths}"
println "filterPaths: ${params.filterPaths}"
println "sampleValues: ${params.sampleValues}"


def rawList = params.rawPaths.split('\\|')
def filterList = params.filterPaths.split('\\|')
def spliceList = params.splicePaths.split('\\|')
def unspliceList = params.unsplicePaths.split('\\|')
def biosampleList = params.biosampleValues.split('\\|')
def sampleList = params.sampleValues.split('\\|')


ch_input = Channel.from(
    [biosampleList, sampleList, rawList, filterList, spliceList, unspliceList].transpose()
).map { biosample, sample, raw, filter, splice, unsplice ->
    tuple(biosample,
          sample,
          file(raw, checkIfExists: true),
          file(filter, checkIfExists: true),
          file(splice, checkIfExists: true),
          file(unsplice, checkIfExists: true)
          )
}

params.soupxR = "${projectDir}/bin/01_qc/run_soupx.R"

process QC {
    tag "$sample"
    
    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

    input:
    tuple val(biosample), val(sample), path(raw), path(filter), path(splice), path(unsplice)
    
    output:
    path("${sample}"), emit: matrix
    path("${sample}_*.txt"), emit: txt
    

    script:
    """
    source /opt/software/miniconda3/bin/activate
    
    conda activate scline
    
    Rscript "${params.soupxR}" --raw_path $raw --filter_path $filter --sample_name $sample --input_mingenes 100 --tfidfMin 1
    """
}

process ANNO {
    tag "$sample"
    
    publishDir "${params.outdir}/scrublet", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(raw), path(filter)
    
    output:
    tuple val(sample), path("${sample}")

    script:
    """
    
    """
}

workflow {
    println "输入样本数量: ${sampleList.size()}"
    println "输出目录: ${params.outdir}"
    QC(ch_input)
    QC.out.view()
}