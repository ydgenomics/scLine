// 读取并解析CSV文件
ch_input_csv = Channel.value(file(params.input, checkIfExists: true))

// 解析CSV文件的函数
def parseInputCSV() {
    def csvFile = file(params.input)
    def lines = csvFile.readLines()
    
    if (lines.size() < 2) {
        throw new IllegalArgumentException("CSV must include 4 columns and 1 row")
    }
    
    def header = lines[0].split(',')
    def rawPaths = []
    def filterPaths = []
    def biosampleValues = []
    def sampleValues = []
    
    // 检查必需的列
    def requiredColumns = ['rawPaths', 'filterPaths', 'biosampleValues', 'sampleValues']
    requiredColumns.each { col ->
        if (!header.contains(col)) {
            throw new IllegalArgumentException("CSV lacks this column: ${col}")
        }
    }
    
    // 解析数据行
    for (int i = 1; i < lines.size(); i++) {
        def line = lines[i].trim()
        if (line.isEmpty()) continue
        
        def values = line.split(',')
        if (values.size() < header.size()) {
            throw new IllegalArgumentException("The number of columns in the ${i+1}$th row does not match the header.")
        }
        
        // 创建映射
        def rowMap = [:]
        header.eachWithIndex { col, idx ->
            rowMap[col] = values[idx]?.trim()
        }
        
        // 收集数据
        rawPaths.add(rowMap.rawPaths)
        filterPaths.add(rowMap.filterPaths)
        biosampleValues.add(rowMap.biosampleValues)
        sampleValues.add(rowMap.sampleValues)
    }
    
    return [
        rawPaths: rawPaths,
        filterPaths: filterPaths,
        biosampleValues: biosampleValues,
        sampleValues: sampleValues
    ]
}

// 检查并处理biosample分组的辅助函数
def checkBiosampleGroups() {
    def parsedData = parseInputCSV()
    def biosampleGroups = [:]
    def sampleToBiosample = [:]
    
    // 构建biosample到sample的映射
    parsedData.biosampleValues.eachWithIndex { biosample, idx ->
        def sample = parsedData.sampleValues[idx]
        if (!biosampleGroups.containsKey(biosample)) {
            biosampleGroups[biosample] = []
        }
        biosampleGroups[biosample].add(sample)
        sampleToBiosample[sample] = biosample
    }
    
    println "=== biosample分组信息 ==="
    biosampleGroups.each { biosample, samples ->
        println "${biosample}: ${samples}"
    }
    println "======================"
    
    return [
        groups: biosampleGroups,
        mapping: sampleToBiosample
    ]
}


workflow {
    // 解析输入CSV
    def parsedData = parseInputCSV()
    
    // 设置全局参数
    params.rawPaths = parsedData.rawPaths
    params.filterPaths = parsedData.filterPaths
    params.biosampleValues = parsedData.biosampleValues
    params.sampleValues = parsedData.sampleValues
    
    // 使用辅助函数检查分组
    def biosampleInfo = checkBiosampleGroups()
    println "Group info: ${biosampleInfo}"
    QC(params.rawPaths, params.filterPaths, params.biosampleValues, params.sampleValues, params.rlst, params.runsoupx, params.tfidfMin, params.outdir)
    CONVERT(QC.out.results.flatten(), "RNA", params.python_env, params.outdir)
    METANEIGHBOR(CONVERT.out.results.collect(), params.prefix, "biosample", "leiden_res_0.50", params.biosampleValues, params.slimit, params.outdir)
    CONCAT_2(QC.out.results.collect(), METANEIGHBOR.out.results.collect(), params.biosampleValues, params.prefix)
    INTEGRATE(CONCAT_2.out.results, params.prefix, params.resolutionSet, params.runsct, "biosample", "anno", params.outdir)
}

// ----- qc -----
include { SOUPX } from './module/qc'
include { SCRUBLET } from './module/qc'

workflow QC {
    take: 
        rawList
        filterList
        biosampleList
        sampleList
        rlst
        runsoupx
        tfidfMin
        outdir
    
    main:
        def flatSampleList = sampleList.flatten()
        ch_samples1 = Channel.fromList(flatSampleList)
        ch_raws1 = Channel.fromList(rawList.collect{ p -> file(p) })
        ch_filters1 = Channel.fromList(filterList.collect{ p -> file(p) })
        
        ch_samples2 = Channel.value(sampleList)
        ch_biosamples2 = Channel.value(biosampleList)
        ch_filters2 = Channel.value(filterList.collect{ p -> file(p) })
        
        if (runsoupx == 'yes') {
            SOUPX(ch_raws1, ch_filters1, ch_samples1, tfidfMin)
            SCRUBLET(SOUPX.out.matrix.collect(), SOUPX.out.txt.collect(), ch_samples2, ch_biosamples2, rlst, outdir)
        } else {
            SCRUBLET(ch_filters2, ch_filters2, ch_samples2, ch_biosamples2, rlst, outdir)
        }
    
    emit:
        results = SCRUBLET.out.h5ad
        csv = SCRUBLET.out.csv
}

// ----- convert -----

include { CONVERT_0 } from './module/convert'

workflow CONVERT {
    take: 
        inputFile
        assay
        python_env
        outdir
    main:
        CONVERT_0(inputFile, assay, python_env, outdir)
    emit:
        results = CONVERT_0.out.result
}

workflow CONVERT2 {
    take: 
        inputFile
        assay
        python_env
        outdir
    main:
        CONVERT_0(inputFile, assay, python_env, outdir)
    emit:
        results = CONVERT_0.out.result
}

// ----- integrate -----
include { CONCAT2 } from './module/integrate'

workflow CONCAT_2 {
    take: 
        ch_h5ad
        ch_csv
        biosampleValues
        prefix
    main:
        // --- concat ---
        ch_biosamples = Channel.value(biosampleValues)
        CONCAT2(ch_h5ad, ch_csv, ch_biosamples, prefix)
    emit:
        results = CONCAT2.out.h5ad
}

include { SCVI } from './module/integrate'
include { HARMONY } from './module/integrate'
include { UNINTEGRATION } from './module/integrate'
include { CONVERT_1 } from './module/convert'
include { SCT_CCA } from './module/integrate'
include { SCT_HARMONY } from './module/integrate'
include { RLIGER } from './module/integrate'
include { DEALPLUS } from './module/integrate'
include { SCIB } from './module/integrate'
include { CONVERT_1 as CONVERT_1_int } from './module/convert'
include { CONVERT_2 as CONVERT_2_int } from './module/convert'

workflow INTEGRATE {
    take: 
        inputFile
        prefix
        resolutionSet
        runsct
        other1_key
        other2_key
        outdir
    main:
        // --- scvi ---
        SCVI(inputFile, prefix, resolutionSet, outdir)
        // --- harmony ---
        HARMONY(inputFile, prefix, resolutionSet, outdir)
        // --- unintegration ---
        UNINTEGRATION(inputFile, prefix, resolutionSet, outdir)
        CONVERT_1(inputFile)
        // --- rliger.inmf ---
        RLIGER(CONVERT_1.out.result, prefix, resolutionSet, outdir)
        CONVERT_1_int(RLIGER.out.rds)
        
        if (runsct == 'yes') {
            // --- sct.transform ---
            SCT_CCA(CONVERT_1.out.result, prefix, resolutionSet, outdir)
            // --- sct.harmony ---
            SCT_HARMONY(CONVERT_1.out.result, prefix, resolutionSet, outdir)
            ch_dealplus = HARMONY.out.h5ad
                .mix(SCVI.out.h5ad)
                .mix(UNINTEGRATION.out.h5ad)
                .mix(RLIGER.out.rds)
                .mix(SCT_CCA.out.rds)
                .mix(SCT_HARMONY.out.rds)
            ch_convert2_int = SCT_CCA.out.rds
                .mix(SCT_HARMONY.out.rds)
            CONVERT_2_int(ch_convert2_int.flatten())
            ch_scib_input = HARMONY.out.h5ad
                .mix(SCVI.out.h5ad)
                .mix(CONVERT_1_int.out.result)
                .mix(CONVERT_2_int.out.result.flatten())
                .collect()
        } else {
            ch_dealplus = HARMONY.out.h5ad
                .mix(SCVI.out.h5ad)
                .mix(UNINTEGRATION.out.h5ad)
                .mix(RLIGER.out.rds)
            ch_scib_input = HARMONY.out.h5ad
                .mix(SCVI.out.h5ad)
                .mix(CONVERT_1_int.out.result)
                .collect()
        }
        
        DEALPLUS(ch_dealplus.flatten(), other1_key, other2_key, outdir)
        SCIB(UNINTEGRATION.out.h5ad, ch_scib_input, other2_key, prefix, outdir)
    emit:
        best_h5ad = SCIB.out.best_h5ad
        rds = RLIGER.out.rds
}

// ----- metaneighbor -----

include { MetaNeighbor } from './module/metaneighbor'

workflow METANEIGHBOR {
    take: 
        inputFile
        prefix
        batch_key
        cluster_key
        seq
        slimit
        outdir
    main:
        MetaNeighbor(inputFile, prefix, batch_key, cluster_key, seq, slimit, outdir)
    emit:
        results = MetaNeighbor.out.anno
}

// ----- dea -----

include { FindMarkers } from './module/dea'

workflow DEA {
    take: 
        inputFile
        cluster_key
        ident_1
        ident_2
        outdir
    main:
        FindMarkers(inputFile, cluster_key, ident_1, ident_2, outdir)
    emit:
        results = FindMarkers.out.csv
}


// ----- pseudotime -----
include { CYTOTRACE } from './module/pseudotime'

workflow PSEUDOTIME {
    take: 
        inputFile
        batch_key
        cluster_key
        cluster_value
        outdir
    main:
        CYTOTRACE(inputFile, batch_key, cluster_key, cluster_value, outdir)
    emit:
        results = CYTOTRACE.out.h5ad
}