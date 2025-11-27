#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 参数验证函数
def validateParams() {
    def errors = []
    
    // 检查模式
    if (!params.mode) {
        errors << "必须指定 --mode 参数，可选值: qc, convert, anno, integrate, metaneighbor, dea, enrich, pseudotime, all"
        return errors
    }
    
    if (!(params.mode in ['qc', 'convert', 'anno', 'integrate', 'metaneighbor', 'dea', 'enrich', 'pseudotime', 'all'])) {
        errors << "无效的模式: ${params.mode}，可选值: qc, convert, anno, integrate, metaneighbor, dea, enrich, pseudotime, all"
        return errors
    }
    
    // 模式特定的验证
    switch(params.mode) {
        case 'qc':
            if (!params.rawPaths) errors << "QC模式必须指定 --rawPaths 参数"
            if (!params.filterPaths) errors << "QC模式必须指定 --filterPaths 参数"
            if (!params.biosampleValues) errors << "QC模式必须指定 --biosampleValues 参数"
            if (!params.sampleValues) errors << "QC模式必须指定 --sampleValues 参数"
            break
            
        case 'convert':
            if (!params.inputFile) errors << "Convert模式必须指定 --inputFile 参数"
            if (!params.assay) errors << "Convert模式必须指定 --assay 参数"
            if (!params.python_env) errors << "Convert模式必须指定 --python_env 参数"
            break
            
        case 'anno':
            if (!params.input_query_rds) errors << "Annotation模式必须指定 --input_query_rds 参数"
            if (!params.rdsORcsv) errors << "Annotation模式必须指定 --rdsORcsv 参数"
            if (!params.cluster_key) errors << "Annotation模式必须指定 --cluster_key 参数"
            break
            
        case 'integrate':
            if (!params.inputH5ad) errors << "Integration模式必须指定 --inputH5ad 参数"
            if (!params.biosampleValues) errors << "Integration模式必须指定 --biosampleValues 参数"
            if (!params.prefix) errors << "Integration模式必须指定 --prefix 参数"
            break
            
        case 'metaneighbor':
            if (!params.inputFile) errors << "MetaNeighbor模式必须指定 --inputFile 参数"
            if (!params.prefix) errors << "MetaNeighbor模式必须指定 --prefix 参数"
            if (!params.batch_key) errors << "MetaNeighbor模式必须指定 --batch_key 参数"
            if (!params.cluster_key) errors << "MetaNeighbor模式必须指定 --cluster_key 参数"
            if (!params.seq) errors << "MetaNeighbor模式必须指定 --seq 参数"
            break
            
        case 'dea':
            if (!params.inputFile) errors << "DEA模式必须指定 --inputFile 参数"
            if (!params.cluster_key) errors << "DEA模式必须指定 --cluster_key 参数"
            if (params.ident_1 == null) errors << "DEA模式必须指定 --ident_1 参数"
            if (params.ident_2 == null) errors << "DEA模式必须指定 --ident_2 参数"
            break
            
        case 'enrich':
            if (!params.emapper_xlsx) errors << "Enrichment模式必须指定 --emapper_xlsx 参数"
            if (!params.gene_csv) errors << "Enrichment模式必须指定 --gene_csv 参数"
            break
            
        case 'pseudotime':
            if (!params.inputFile) errors << "Pseudotime模式必须指定 --inputFile 参数"
            if (!params.batch_key) errors << "Pseudotime模式必须指定 --batch_key 参数"
            if (!params.cluster_key) errors << "Pseudotime模式必须指定 --cluster_key 参数"
            if (!params.cluster_value) errors << "Pseudotime模式必须指定 --cluster_value 参数"
            break
        
        case 'all':
            if (!params.rawPaths) errors << "All模式必须指定 --rawPaths 参数"
            if (!params.filterPaths) errors << "All模式必须指定 --filterPaths 参数"
            if (!params.biosampleValues) errors << "All模式必须指定 --biosampleValues 参数"
            if (!params.sampleValues) errors << "All模式必须指定 --sampleValues 参数"
            break
    }
    
    return errors
}

// 显示帮助信息并退出
if (params.help) {
    switch(params.mode) {
        case 'qc': println qcHelpMessage(); break
        case 'convert': println convertHelpMessage(); break
        case 'anno': println annoHelpMessage(); break
        case 'integrate': println integrateHelpMessage(); break
        case 'metaneighbor': println metaneighborHelpMessage(); break
        case 'dea': println deaHelpMessage(); break
        case 'enrich': println enrichHelpMessage(); break
        case 'pseudotime': println pseudotimeHelpMessage(); break
        case 'all': println allHelpMessage(); break
        default: println helpMessage()
    }
    System.exit(0)
}

def helpMessage() {
    """
    ================================================================================
    scline: single-cell analysis pipeline
    ================================================================================
    
    Usage:
    nextflow run main.nf --mode <MODE> [OPTIONS]
    
    Available Modes:
    - qc          : Quality control and preprocessing
    - convert     : Data format conversion
    - anno        : Cell type annotation
    - integrate   : Data integration and batch correction
    - metaneighbor: MetaNeighbor analysis
    - dea         : Differential expression analysis
    - enrich      : Functional enrichment analysis
    - pseudotime  : Pseudotime trajectory analysis
    
    For mode-specific help:
    nextflow run main.nf --mode <MODE> --help
    
    ================================================================================
    """
}

def qcHelpMessage() {
    """
    ================================================================================
    QC Mode - Quality Control and Preprocessing
    ================================================================================
    
    Usage:
    nextflow run main.nf --mode qc [OPTIONS]
    
    Required Parameters:
    --rawPaths str                Pipe-separated list of raw data paths
    --filterPaths str             Pipe-separated list of filtered data paths
    --biosampleValues str         Pipe-separated biosample group names
    --sampleValues str            Pipe-separated sample names
    
    Optional Parameters:
    --tfidfMin int/float          Minimum value of tfidf to accept for a marker gene (SoupX)
                                  Default: 1 (If SoupX is not working properly, try decreasing it.)
    --rlst str                    Pipe-separated list of rlst values for SoupX
                                  Default: "0.2|0.5|0.8|1.0"
    --runsoupx str                Run SoupX contamination removal? 'yes' or 'no'
                                  Default: 'no'
    --outdir str                  Output directory path
                                  Default: './results'
    --qc_cpu str                  CPU resource [SOUPX|SCRUBLET]
                                  Default: '2|1'                      
    --qc_mem str                  Memory resource (GB) [SOUPX|SCRUBLET]
                                  Default: '4|2'
                                  
    Example:
    nextflow run main.nf --mode qc \\
        --rawPaths '/data/work/scline/input/sample1/sample1_raw|/data/work/scline/input/sample2/sample2_raw' \\
        --filterPaths '/data/work/scline/input/sample1/sample1_filter|/data/work/scline/input/sample2/sample2_filter' \\
        --biosampleValues 'group1|group2' \\
        --sampleValues 'sample1|sample2' \\
        --runsoupx 'yes' \\
        --outdir './results/qc'
    
    ================================================================================
    """
}

def convertHelpMessage() {
    """
    ================================================================================
    Convert Mode - Data Format Conversion
    ================================================================================
    
    Usage:
    nextflow run main.nf --mode convert [OPTIONS]
    
    Required Parameters:
    --inputFile str               Input file path for conversion
    --assay str                   Assay type will be converted
    --python_env str              Python environment to use
    
    Optional Parameters:
    --outdir str                  Output directory path
                                  Default: './results'
    --convert_cpu int             CPU resource
                                  Default: 1                    
    --convert_mem int             Memory resource (GB)
                                  Default: 2
                                  
    Example:
    nextflow run main.nf --mode convert \\
        --inputFile '/data/work/scline/results/qc/qc/group1/group1.h5ad' \\
        --assay 'RNA' \\
        --python_env '/opt/software/miniconda3/envs/scline/bin/python' \\
        --outdir './results/convert'
    
    ================================================================================
    """
}

def annoHelpMessage() {
    """
    ================================================================================
    Annotation Mode - Cell Type Annotation
    ================================================================================
    
    Usage:
    nextflow run main.nf --mode anno [OPTIONS]
    
    Required Parameters:
    --input_query_rds str         Input query .rds file (unannotated)
    --rdsORcsv str                Annotation reference file for scType(.csv) or singleR(.rds)
    --cluster_key str             Cluster key in the data for annotation
    --umap_name str               UMAP name for plotting
                                  Default: 'Xumap_' (must exist in .obsm)
    --ref_cluster_key str         Reference cluster key in rdsORcsv(.rds)
                                  Default: 'celltype' (when run singleR is required!)
    
    Optional Parameters:
    --outdir str                  Output directory path
                                  Default: './results'
    --anno_cpu int                CPU resource
                                  Default: 1                    
    --anno_mem int                Memory resource (GB)
                                  Default: 2
    
    Example:
    nextflow run main.nf --mode anno \\
        --input_query_rds '/data/work/scline/results/convert/group1.hr.rds' \\
        --rdsORcsv '/data/work/scline/input/markergene.csv' \\
        --cluster_key 'leiden_res_0.50' \\
        --umap_name 'Xumap_' \\
        --outdir './results/anno'
    
    ================================================================================
    """
}

def integrateHelpMessage() {
    """
    ================================================================================
    Integration Mode - Data Integration and Batch Correction
    ================================================================================
    
    Usage:
    nextflow run main.nf --mode integrate [OPTIONS]
    
    Required Parameters:
    --inputH5ad str               Pipe-separated list of input H5AD files
    --annoCsv str                 Annotation CSV file (two columns: 'barcode' & 'anno')
                                  (if not exist, could use the same input of '--inputH5ad')
    --biosampleValues str         Pipe-separated biosample names
    --prefix str                  Output prefix
    --other1_key str              Additional metadata key 1 (must exist in .obs.column)
                                  Default: 'biosample'
    --other2_key str              Additional metadata key 2 (must exist in .obs.column and best with annotation info)
                                  Default: 'leiden_res_0.50'
    --python_env str              Python environment to use                              
    
    Optional Parameters:
    --resolutionSet float         Resolution parameter for clustering after integration
                                  Default: 0.5
    --runsct str                  Run SCTransform(SCT.CCA and SCT.harmony)? 'yes' or 'no'
                                  Default: 'no' (Running it will spend many resource)
    --outdir str                  Output directory path
                                  Default: './results'
    --inte_cpu str                CPU resource [CONCAT|SCVI|HARMONY|UNINTEGRATION|CONVERT_1|RLIGER|CONVERT_1_int|SCT_CCA|SCT_HARMONY|CONVERT_2_int|DEALPLUS|SCIB]
                                  Default: '1|8|2|2|1|2|1|4|4|1|1|8'                    
    --inte_mem str                Memory resource (GB) [CONCAT|SCVI|HARMONY|UNINTEGRATION|CONVERT_1|RLIGER|CONVERT_1_int|SCT_CCA|SCT_HARMONY|CONVERT_2_int|DEALPLUS|SCIB]
                                  Default: '2|4|4|4|4|4|4|8|8|4|4|8'
                                  
    Example:
    nextflow run main.nf --mode integrate \\
        --inputH5ad '/data/work/scline/results/qc/qc/group1/group1.h5ad|/data/work/scline/results/qc/qc/group2/group2.h5ad' \\
        --annoCsv '/data/work/scline/results/qc/qc/group1/group1.h5ad|/data/work/scline/results/qc/qc/group2/group2.h5ad' \\
        --biosampleValues 'group1|group2' \\
        --prefix 'cotton' \\
        --other1_key 'biosample' \\
        --other2_key 'leiden_res_0.50' \\
        --python_env '/opt/software/miniconda3/envs/scline/bin/python' \\
        --runsct 'yes' \\
        --outdir './results/integrate'
    
    ================================================================================
    """
}

def metaneighborHelpMessage() {
    """
    ================================================================================
    MetaNeighbor Mode - Cell Type Similarity Analysis
    ================================================================================
    
    Usage:
    nextflow run main.nf --mode metaneighbor [OPTIONS]
    
    Required Parameters:
    --inputFile str               Input integrated RDS file
    --prefix str                  Output prefix
    --batch_key str               Batch key in metadata
    --cluster_key str             Cluster key in metadata
    --seq str                     Pipe-separated sequence identifiers (values of batch_key combined with '|')
    
    Optional Parameters:
    --slimit float                Similarity limit threshold
                                  Default: 0.8
    --outdir str                  Output directory path
                                  Default: './results'
    --meta_cpu int                CPU resource
                                  Default: 1                    
    --meta_mem int                Memory resource (GB)
                                  Default: 2
                                  
    Example:
    nextflow run main.nf --mode metaneighbor \\
        --inputFile '/data/work/scline/results/integrate/cotton_rliger.INMF_integrated_dealplus.rds' \\
        --prefix 'cotton' \\
        --batch_key 'biosample' \\
        --cluster_key 'leiden_res_0.50' \\
        --seq 'group1|group2' \\
        --outdir './results/metaneighbor'
    
    ================================================================================
    """
}

def deaHelpMessage() {
    """
    ================================================================================
    DEA Mode - Differential Expression Analysis
    ================================================================================
    
    Usage:
    nextflow run main.nf --mode dea [OPTIONS]
    
    Required Parameters:
    --inputFile str               Input RDS file
    --cluster_key str             Cluster key in metadata
    --ident_1 str                 Identity 1 for comparison
    --ident_2 str                 Identity 2 for comparison
    
    Optional Parameters:
    --outdir str                  Output directory path
                                  Default: './results'
    --dea_cpu int                 CPU resource
                                  Default: 1                    
    --dea_mem int                 Memory resource (GB)
                                  Default: 2
                                  
    Example:
    nextflow run main.nf --mode dea \\
        --inputFile '/data/work/scline/results/integrate/cotton_rliger.INMF_integrated_dealplus.rds' \\
        --cluster_key 'leiden_res_0.50' --ident_1 '0' --ident_2 '1' \\
        --outdir './results/dea'
    
    ================================================================================
    """
}

def enrichHelpMessage() {
    """
    ================================================================================
    ENRICH Mode - Functional Enrichment Analysis
    ================================================================================
    
    Usage:
    nextflow run main.nf --mode enrich [OPTIONS]
    
    Required Parameters:
    --emapper_xlsx str            EggNOG mapper annotations file
    --gene_csv str                Gene markers CSV file included these columns (cluster/gene/) 
    
    Optional Parameters:
    --genus str                   Genus name
                                  Default: "Genus"
    --species str                 Species name
                                  Default: "species"
    --minp float                  Minimum p-value threshold
                                  Default: 0.05
    --outdir str                  Output directory path
                                  Default: './results'
    --enrich_cpu int              CPU resource
                                  Default: 1                    
    --enrich_mem int              Memory resource (GB)
                                  Default: 2
                                  
    Example:
    nextflow run main.nf --mode enrich \\
        --emapper_xlsx '/data/work/scline/bin/06_enrich/out.emapper.annotations.xlsx' \\
        --gene_csv '/data/work/scline/bin/06_enrich/leiden_res_0.50.markers.csv' \\
        --outdir './results/enrich'
    
    ================================================================================
    """
}

def pseudotimeHelpMessage() {
    """
    ================================================================================
    Pseudotime Mode - Trajectory Analysis
    ================================================================================
    
    Usage:
    nextflow run main.nf --mode pseudotime [OPTIONS]
    
    Required Parameters:
    --inputFile str               Input H5AD file
    --batch_key str               Batch key in metadata
    --cluster_key str             Cluster key in metadata
    --cluster_value str           Pipe-separated cluster values for pseudotime (values of cluster_key combined with '|')
    
    Optional Parameters:
    --outdir str                  Output directory path
                                  Default: './results'
    --time_cpu int                CPU resource
                                  Default: 1                    
    --time_mem int                Memory resource (GB)
                                  Default: 4
                                  
    Example:
    nextflow run main.nf --mode pseudotime \\
        --inputFile '/data/work/scline/results/integrate/cotton_harmony_integrated_dealplus.h5ad' \\
        --batch_key 'biosample' \\
        --cluster_key 'celltype' \\
        --cluster_value '0|1' \\
        --outdir './results/pseudotime'
    
    ================================================================================
    """
}

def allHelpMessage() {
    """
    ================================================================================
    ALL Mode - A pipeline analysis for existing parts
    ================================================================================
    
    Usage:
    nextflow run main.nf --mode all [OPTIONS]
    
    Required Parameters:
    --rawPaths str                Pipe-separated list of raw data paths
    --filterPaths str             Pipe-separated list of filtered data paths
    --biosampleValues str         Pipe-separated biosample group names
    --sampleValues str            Pipe-separated sample names
    --python_env str              Python environment to use
    --rdsORcsv str                Annotation reference file for scType(.csv) or singleR(.rds)
    --cluster_key str             Cluster key in the data for annotation
    --ref_cluster_key str         Reference cluster key in rdsORcsv(.rds)
                                  Default: 'celltype' (when run singleR is required!)
    --prefix str                  Output prefix                              
    --other1_key str              Additional metadata key 1 (must exist in .obs.column)
                                  Default: 'biosample'
    --other2_key str              Additional metadata key 2 (must exist in .obs.column and best with annotation info)
                                  Default: 'leiden_res_0.50'
    --ident_1 str                 Identity 1 for comparison
    --ident_2 str                 Identity 2 for comparison
    --cluster_value str           Pipe-separated cluster values for pseudotime (values of cluster_key combined with '|')
    
    Optional Parameters:
    --tfidfMin int/float          Minimum value of tfidf to accept for a marker gene (SoupX)
                                  Default: 1 (If SoupX is not working properly, try decreasing it.)
    --rlst str                    Pipe-separated list of rlst values for SoupX
                                  Default: "0.2|0.5|0.8|1.0"
    --runsoupx str                Run SoupX contamination removal? 'yes' or 'no'
                                  Default: 'no'
    --resolutionSet float         Resolution parameter for clustering after integration
                                  Default: 0.5
    --runsct str                  Run SCTransform(SCT.CCA and SCT.harmony)? 'yes' or 'no'
                                  Default: 'no' (Running it will spend many resource)
    --slimit float                Similarity limit threshold
                                  Default: 0.8
    --emapper_xlsx str            EggNOG mapper annotations file
    --genus str                   Genus name
                                  Default: "Genus"
    --species str                 Species name
                                  Default: "species"
    --minp float                  Minimum p-value threshold
    --qc_cpu str                  CPU resource [SOUPX|SCRUBLET]
                                  Default: '2|1'                      
    --qc_mem str                  Memory resource (GB) [SOUPX|SCRUBLET]
                                  Default: '4|2'
    --convert_cpu int             CPU resource
                                  Default: 1                    
    --convert_mem int             Memory resource (GB)
                                  Default: 2
    --anno_cpu int                CPU resource
                                  Default: 1                    
    --anno_mem int                Memory resource (GB)
                                  Default: 2
    --inte_cpu str                CPU resource [CONCAT|SCVI|HARMONY|UNINTEGRATION|CONVERT_1|RLIGER|CONVERT_1_int|SCT_CCA|SCT_HARMONY|CONVERT_2_int|DEALPLUS|SCIB]
                                  Default: '1|8|2|2|1|2|1|4|4|1|1|8'                    
    --inte_mem str                Memory resource (GB) [CONCAT|SCVI|HARMONY|UNINTEGRATION|CONVERT_1|RLIGER|CONVERT_1_int|SCT_CCA|SCT_HARMONY|CONVERT_2_int|DEALPLUS|SCIB]
                                  Default: '2|4|4|4|4|4|4|8|8|4|4|8'
    --meta_cpu int                CPU resource
                                  Default: 1                    
    --meta_mem int                Memory resource (GB)
                                  Default: 2
    --dea_cpu int                 CPU resource
                                  Default: 1                    
    --dea_mem int                 Memory resource (GB)
                                  Default: 2                              
    --enrich_cpu int              CPU resource
                                  Default: 1                    
    --enrich_mem int              Memory resource (GB)
                                  Default: 2
    --time_cpu int                CPU resource
                                  Default: 1                    
    --time_mem int                Memory resource (GB)
                                  Default: 4
    --outdir str                  Output directory path
                                  Default: './results'
    
    Example:
    nextflow run main.nf --mode all \\
        --rawPaths '/data/work/scline/input/sample1/sample1_raw|/data/work/scline/input/sample2/sample2_raw' \\
        --filterPaths '/data/work/scline/input/sample1/sample1_filter|/data/work/scline/input/sample2/sample2_filter' \\
        --biosampleValues 'group1|group2' \\
        --sampleValues 'sample1|sample2' \\
        --python_env '/opt/software/miniconda3/envs/scline/bin/python' \\
        --rdsORcsv '/data/work/scline/input/markergene.csv' \\
        --cluster_key 'leiden_res_0.50' \\
        --prefix 'cotton' \\
        --other1_key 'biosample' \\
        --other2_key 'anno' \\
        --ident_1 'celltype2' --ident_2 'celltype3' \\
        --cluster_value 'celltype2|celltype3' \\
        --runsoupx 'yes' \\
        --runsct 'yes' \\
        --outdir 'results2'
    
    ================================================================================
    """
}


// 参数验证（在help检查之后）
def validationErrors = validateParams()
if (validationErrors) {
    println "参数验证失败:"
    validationErrors.each { error ->
        println " - $error"
    }
    println ""
    
    // 显示对应的帮助信息
    switch(params.mode) {
        case 'qc': println qcHelpMessage(); break
        case 'convert': println convertHelpMessage(); break
        case 'anno': println annoHelpMessage(); break
        case 'integrate': println integrateHelpMessage(); break
        case 'metaneighbor': println metaneighborHelpMessage(); break
        case 'dea': println deaHelpMessage(); break
        case 'enrich': println enrichHelpMessage(); break
        case 'pseudotime': println pseudotimeHelpMessage(); break
        case 'all': println allHelpMessage(); break
        default: println helpMessage()
    }
    System.exit(1)
}

// 在 workflow 外部定义辅助函数
def checkBiosampleGroups() {
    def biosamples = params.biosampleValues.split('\\|') as List
    def uniqueBiosamples = biosamples.unique()  // 现在 biosamples 是 List，有 unique() 方法
    return [
        biosamples: biosamples,
        uniqueBiosamples: uniqueBiosamples,
        isSingleGroup: uniqueBiosamples.size() == 1,
        groupCount: uniqueBiosamples.size()
    ]
}

workflow {
    // 使用辅助函数
    def biosampleInfo = checkBiosampleGroups()
    println "分组信息: ${biosampleInfo}"
    
    if (params.mode == 'qc') {
        QC(params.rawPaths, params.filterPaths, params.biosampleValues, params.sampleValues, params.rlst, params.runsoupx, params.tfidfMin, params.outdir)
    }
    else if (params.mode == 'convert') {
        CONVERT(params.inputFile, params.assay, params.python_env, params.outdir)
    }
    else if (params.mode == 'anno') {
        ANNO(params.input_query_rds, params.rdsORcsv, params.cluster_key, params.umap_name, params.ref_cluster_key, params.outdir)
    }
    else if (params.mode == 'integrate') {
        CONCAT_1(params.inputH5ad, params.final_annoCsv, params.biosampleValues, params.prefix)
        INTEGRATE(CONCAT_1.out.results, params.prefix, params.resolutionSet, params.runsct, params.other1_key, params.other2_key, params.outdir)
    }
    else if (params.mode == 'metaneighbor') {
        METANEIGHBOR(params.inputFile, params.prefix, params.batch_key, params.cluster_key, params.seq, params.slimit, params.outdir)
    }
    else if (params.mode == 'dea') {
        DEA(params.inputFile, params.cluster_key, params.ident_1, params.ident_2, params.outdir)
    }
    else if (params.mode == 'enrich') {
        ENRICH(params.emapper_xlsx, params.gene_csv, params.genus, params.species, params.minp, params.outdir)
    }
    else if (params.mode == 'pseudotime') {
        PSEUDOTIME(params.inputFile, params.batch_key, params.cluster_key, params.cluster_value, params.outdir)
    }
    else if (params.mode == 'all') {
        if (biosampleInfo.isSingleGroup) {
            println "单组样本: ${biosampleInfo.uniqueBiosamples[0]}"
            QC(params.rawPaths, params.filterPaths, params.biosampleValues, params.sampleValues, params.rlst, params.runsoupx, params.tfidfMin, params.outdir)
            CONVERT(QC.out.results.flatten(), "RNA", params.python_env, params.outdir)
            ANNO(CONVERT.out.results.flatten(), params.rdsORcsv, params.cluster_key, "Xumap_", params.ref_cluster_key, params.outdir)
            CONCAT_2(QC.out.results.collect(), ANNO.out.results.collect(), params.biosampleValues, params.prefix)
            CONVERT2(CONCAT_2.out.results, "RNA", params.python_env, params.outdir)
            DEA(CONVERT2.out.results, params.other2_key, params.ident_1, params.ident_2, params.outdir)
            if (!params.emapper_xlsx) {
                println "Don't run ENRICH"
            } else {
                ch_enrich_input = DEA.out.results.mix(QC.out.csv.flatten()).collect()
                ENRICH(params.emapper_xlsx, ch_enrich_input, params.genus, params.species, params.minp, params.outdir)
            }
            PSEUDOTIME(CONCAT_2.out.results, "sample", "anno", params.cluster_value, params.outdir)
        } else {
            println "多组样本: ${biosampleInfo.uniqueBiosamples}"
            QC(params.rawPaths, params.filterPaths, params.biosampleValues, params.sampleValues, params.rlst, params.runsoupx, params.tfidfMin, params.outdir)
            CONVERT(QC.out.results.flatten(), "RNA", params.python_env, params.outdir)
            ANNO(CONVERT.out.results.flatten(), params.rdsORcsv, params.cluster_key, "Xumap_", params.ref_cluster_key, params.outdir)
            CONCAT_2(QC.out.results.collect(), ANNO.out.results.collect(), params.biosampleValues, params.prefix)
            INTEGRATE(CONCAT_2.out.results, params.prefix, params.resolutionSet, params.runsct, "biosample", "anno", params.outdir)
            METANEIGHBOR(INTEGRATE.out.rds, params.prefix, "biosample", params.cluster_key, params.biosampleValues, params.slimit, params.outdir)
            DEA(INTEGRATE.out.rds, params.other2_key, params.ident_1, params.ident_2, params.outdir)
            if (!params.emapper_xlsx) {
                println "Don't run ENRICH"
            } else {
                ch_enrich_input = DEA.out.results.mix(QC.out.csv.flatten()).collect()
                ENRICH(params.emapper_xlsx, ch_enrich_input, params.genus, params.species, params.minp, params.outdir)
                // ENRICH(params.emapper_xlsx, QC.out.csv.flatten(), params.genus, params.species, params.minp, params.outdir)
                // ENRICH(params.emapper_xlsx, DEA.out.results, params.genus, params.species, params.minp, params.outdir)
            }
            PSEUDOTIME(INTEGRATE.out.best_h5ad, "biosample", "anno", params.cluster_value, params.outdir)
        }
    }
    else if (params.mode) {
        println "ERROR: Unknown mode: ${params.mode}"
        println helpMessage()
        System.exit(1)
    }
    else {
        println helpMessage()
        System.exit(0)
    }
}


// ----- qc -----

include { SOUPX } from './module/qc'
include { SCRUBLET } from './module/qc'

workflow QC {
    take: 
        rawPaths
        filterPaths
        biosampleValues
        sampleValues
        rlst
        runsoupx
        tfidfMin
        outdir
    
    main:
        def rawList = rawPaths.split('\\|')
        def filterList = filterPaths.split('\\|')
        def biosampleList = biosampleValues.split('\\|')
        def sampleList = sampleValues.split('\\|')

        // 添加调试信息
        println "=== QC工作流调试信息 ==="
        println "sampleList: ${sampleList}"
        println "sampleList class: ${sampleList.getClass()}"
        println "sampleList size: ${sampleList.size()}"
        
        // 确保 sampleList 是平面列表
        def flatSampleList = sampleList.flatten()
        println "flatSampleList: ${flatSampleList}"
        println "======================"

        // 使用确保是平面列表的版本
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
// ----- anno -----

include { ANNOR } from './module/anno'

workflow ANNO {
    take: 
        input_query_rds
        rdsORcsv
        cluster_key
        umap_name
        ref_cluster_key
        outdir
    main:
        ANNOR(input_query_rds, rdsORcsv, cluster_key, umap_name, ref_cluster_key, outdir)
    emit:
        results = ANNOR.out.anno
}

// ----- integrate -----

include { CONCAT } from './module/integrate'
include { CONCAT2 } from './module/integrate'

workflow CONCAT_1 {
    take: 
        inputH5ad
        annoCsv
        biosampleValues
        prefix
    main:
        // --- concat ---
        def inputH5ad = params.inputH5ad.split('\\|')
        def annoCsv = params.annoCsv.split('\\|')
        def biosampleValues = params.biosampleValues.split('\\|')
        ch_h5ad = Channel.value(inputH5ad.collect{ p -> file(p) })
        ch_csv = Channel.value(annoCsv.collect{ p -> file(p) })
        ch_biosamples = Channel.value(biosampleValues)
        CONCAT(ch_h5ad, ch_csv, ch_biosamples, prefix)
    emit:
        results = CONCAT.out.h5ad
}

workflow CONCAT_2 {
    take: 
        ch_h5ad
        ch_csv
        biosampleValues
        prefix
    main:
        // --- concat ---
        def biosampleValues = params.biosampleValues.split('\\|')
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

/*
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
*/

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
                .mix(UNINTEGRATION.out.h5ad)
                .mix(RLIGER.out.rds)
                .mix(SCT_CCA.out.rds)
                .mix(SCT_HARMONY.out.rds)
            ch_convert2_int = SCT_CCA.out.rds
                .mix(SCT_HARMONY.out.rds)
            CONVERT_2_int(ch_convert2_int.flatten())
            ch_scib_input = HARMONY.out.h5ad
                .mix(CONVERT_1_int.out.result)
                .mix(CONVERT_2_int.out.result.flatten())
                .collect()
        } else {
            ch_dealplus = HARMONY.out.h5ad
                .mix(UNINTEGRATION.out.h5ad)
                .mix(RLIGER.out.rds)
            ch_scib_input = HARMONY.out.h5ad
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
        results = MetaNeighbor.out.csv
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

// ----- enrich -----

include { makeOrgPackage } from './module/enrich'
include { clusterProfiler } from './module/enrich'

workflow ENRICH {
    take: 
        emapper_xlsx
        gene_csv
        genus
        species
        minp
        outdir
    main:
        makeOrgPackage(emapper_xlsx, genus, species, outdir)
        clusterProfiler(makeOrgPackage.out.orgdb, gene_csv, makeOrgPackage.out.kegg_info_RData, genus, species, minp, outdir)
    emit:
        results = makeOrgPackage.out.orgdb
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