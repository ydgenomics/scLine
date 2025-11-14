#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FASTQC    } from './modules/fastqc.nf'
include { MULTIQC   } from './modules/multiqc.nf'
include { SCANPY_QC } from './modules/scanpy_qc.nf'

workflow {

    // 读取样本表
    ch_input = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.r1), file(row.r2)) }

    // 1. FastQC 原始数据质控
    FASTQC(ch_input)

    // 2. Scanpy 细胞级质控
    SCANPY_QC(ch_input)

    // 3. 汇总报告
    MULTIQC(FASTQC.out.collect())
}