#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// 第1段：参数定义
params.txt = 'data/*.txt'
params.outdir = 'results'
println "参数定义正常"

// 第2段：通道定义
ch_txt = Channel.fromPath(params.txt)
                .map { file -> tuple(file.simpleName, file) }
println "通道定义正常"

// 第3段：第一个进程
process TEST {
    tag "$sample"
    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(text)
    
    output:
    tuple val(sample), path("${sample}_copy.txt")
    
    script:
    """
    cp "$text" "${sample}_copy.txt"
    """
}
println "TEST进程定义正常"

// 第4段：第二个进程
process ADD {
    tag "$sample"
    publishDir "${params.outdir}/add", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(orig)
    
    output:
    tuple val(sample), path("${sample}_withHello.txt")
    
    script:
    """
    cp "$orig" "${sample}_withHello.txt"
    echo "hello world" >> "${sample}_withHello.txt"
    """
}
println "ADD进程定义正常"

// 第5段：工作流
workflow {
    TEST(ch_txt)
    ADD(TEST.out)
}
