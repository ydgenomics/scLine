# nextflow tutorial
nextflow至少需要一个*.nf文件作为运行文件，代码`nextflow run main.nf`。main.nf主要有两部分组成(process和workflow)，process对应子任务，每个子任务可以包含以下部分(input, output, script)，其中至少包括script这部分；workflow定义子任务的运行关系和输入。对于流程外部的参数可以通过`params.varaible_name`进行定义，process使用参数的方法有两种，一是直接在input定义后做传参，另一种是直接用`${params.variable_name}`的方式直接调用，不用在workflow调用子process时显示输入。
逐个元素 → Channel.fromList(...)
整体列表 → Channel.value(...)

`nextflow.config` 是 Nextflow 的核心配置文件，它定义了工作流的全局行为、参数、资源配置等。下面详细讲解它的作用：

## 1. 基本结构和位置

### 文件位置
- 默认在工作流目录下
- Nextflow 会自动加载
- 可以有多个（项目级、用户级、系统级）

### 基本结构
```nextflow
// nextflow.config

// 1. 参数定义
params {
    outdir = 'results'
    input = 'data/*.fastq'
    threads = 4
}

// 2. 进程配置
process {
    cpus = 2
    memory = '4 GB'
    time = '1 hour'
}

// 3. 执行配置
executor {
    name = 'local'
    queueSize = 100
}

// 4. 其他配置
docker {
    enabled = true
}
```

## 2. 参数定义 (params)

### 默认参数
```nextflow
params {
    // 基本参数
    input = null
    outdir = './results'
    
    // 流程参数
    minReads = 1000
    qualityThreshold = 20
    
    // 文件路径
    reference = '/data/genomes/hg38.fa'
    annotation = '/data/annotations/genes.gtf'
    
    // 布尔参数
    skipQC = false
    saveIntermediate = true
}
```

### 参数验证
```nextflow
// 参数验证和默认值
if (!params.input) {
    error "必须提供 --input 参数"
}

if (!params.outdir) {
    params.outdir = "results_${new Date().format('yyyyMMdd_HHmmss')}"
}
```

## 3. 进程配置 (process)

### 全局进程配置
```nextflow
process {
    // 资源默认值
    cpus = 1
    memory = '2 GB'
    time = '1 h'
    
    // 执行器配置
    executor = 'local'
    queue = 'default'
    
    // 错误处理
    errorStrategy = 'retry'
    maxRetries = 3
    maxErrors = '-1'
    
    // 调试
    echo = true
    shell = ['/bin/bash', '-euo', 'pipefail']
}
```

### 进程特定配置
```nextflow
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = '4 GB'
        time = '30 min'
    }
    
    withName: 'BWA_MEM' {
        cpus = 8
        memory = '16 GB'
        time = '2 hours'
        container = 'biocontainers/bwa:0.7.17'
    }
    
    withName: 'SALMON' {
        cpus = 12
        memory = '32 GB'
        time = '4 hours'
        errorStrategy = 'ignore'
    }
    
    // 使用正则表达式匹配多个进程
    withName: '.*QC' {
        memory = '8 GB'
    }
}
```

## 4. 执行器配置 (executor)

### 本地执行
```nextflow
executor {
    name = 'local'
    queueSize = 100
    pollInterval = '5 sec'
}
```

### SLURM 集群
```nextflow
executor {
    name = 'slurm'
    queueSize = 100
    pollInterval = '1 min'
    
    // SLURM 特定配置
    queue = 'normal'
    account = 'my_project'
}
```

### AWS Batch
```nextflow
executor {
    name = 'awsbatch'
    queueSize = 100
    awsRegion = 'us-east-1'
    jobRole = 'arn:aws:iam::123456789:role/NextflowRole'
}
```

## 5. 环境配置

### Conda 配置
```nextflow
conda {
    cacheDir = '/scratch/conda_cache'
    useMamba = true
    createTimeout = '30 min'
}

// 为特定进程指定 conda 环境
process {
    withName: 'SOUPX' {
        conda = '/opt/software/miniconda3/envs/scline'
    }
    
    withName: 'SEURAT' {
        conda = 'environment.yml'
    }
}
```

### Docker 配置
```nextflow
docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
    temp = 'auto'
}

process {
    withName: 'CELLRANGER' {
        container = 'nfcore/cellranger:7.1.0'
    }
}
```

### 环境变量
```nextflow
env {
    PYTHONPATH = '/opt/software/lib/python'
    R_LIBS_USER = '/opt/software/R/library'
    TMPDIR = '/scratch/tmp'
}
```

## 6. 日志和报告配置

### 详细的日志配置
```nextflow
// 时间线
timeline {
    enabled = true
    file = "${params.outdir}/logs/timeline.html"
}

// 执行跟踪
trace {
    enabled = true
    file = "${params.outdir}/logs/trace.txt"
    fields = 'task_id,hash,process,status,exit,module,container,cpus,time,disk,memory,attempt,submit,start,complete,duration,realtime,percent_cpu,rss,vmem,peak_rss,peak_vmem,rchar,wchar'
}

// 报告
report {
    enabled = true
    file = "${params.outdir}/logs/report.html"
}

// DAG 图
dag {
    enabled = true
    file = "${params.outdir}/logs/pipeline_dag.html"
}

// 邮件通知
notification {
    enabled = true
    to = 'user@example.com'
    from = 'nextflow@server.com'
}
```

## 7. 资源配置

### 内存和 CPU
```nextflow
process {
    memory = { 
        2.GB * task.attempt 
    }
    cpus = { 
        Math.min( task.attempt * 2, 16 ) 
    }
    time = { 
        task.attempt == 1 ? '1h' : '2h' 
    }
}
```

### 动态资源分配
```nextflow
process {
    withName: 'ALIGNMENT' {
        memory = { 
            // 根据输入文件大小动态分配内存
            def size = task.getInputSize()
            if (size < 1.GB) '8.GB'
            else if (size < 10.GB) '16.GB' 
            else '32.GB'
        }
    }
}
```

## 8. 配置文件继承和覆盖

### 多环境配置
```nextflow
// nextflow.config - 基础配置
params {
    outdir = './results'
    input = null
}

profiles {
    // 开发环境
    dev {
        process.executor = 'local'
        process.memory = '4.GB'
        params.input = 'test_data/*.fastq'
    }
    
    // 生产环境 - SLURM
    production {
        process.executor = 'slurm'
        process.memory = '32.GB'
        process.time = '24.h'
        process.queue = 'highmem'
        params.outdir = '/project/results'
    }
    
    // AWS 环境
    aws {
        process.executor = 'awsbatch'
        process.container = 'myorg/mypipeline:latest'
        aws.region = 'us-east-1'
    }
}
```

### 使用配置文件
```bash
# 使用开发配置
nextflow run pipeline.nf -profile dev

# 使用生产配置
nextflow run pipeline.nf -profile production

# 组合多个配置
nextflow run pipeline.nf -profile production,docker
```

## 9. 完整的实际示例

```nextflow
// nextflow.config

// 参数定义
params {
    // 输入输出
    reads = null
    outdir = './results'
    
    // 分析参数
    min_genes = 100
    min_cells = 3
    nfeatures = 2000
    
    // 参考数据
    reference = '/data/references/hg38'
    annotation = '/data/annotations/gencode.v38.gtf'
    
    // 流程控制
    skip_qc = false
    save_intermediate = true
}

// 参数验证
if (!params.reads) {
    error "必须提供 --reads 参数，包含输入文件路径"
}

// 进程配置
process {
    // 全局默认值
    cpus = 1
    memory = '4.GB'
    time = '1.h'
    errorStrategy = 'retry'
    maxRetries = 3
    
    // 进程特定配置
    withName: 'FASTQC' {
        cpus = 2
        memory = '8.GB'
        time = '30.min'
    }
    
    withName: 'CELL_RANGER' {
        cpus = 16
        memory = '64.GB'
        time = '24.h'
        container = 'nfcore/cellranger:7.1.0'
    }
    
    withName: 'SEURAT_PROCESSING' {
        cpus = 8
        memory = '32.GB'
        time = '4.h'
        conda = '/opt/conda/envs/seurat'
    }
}

// 执行器配置
executor {
    name = 'slurm'
    queueSize = 100
    pollInterval = '30 sec'
}

// 日志配置
timeline {
    enabled = true
    file = "${params.outdir}/logs/timeline.html"
}

report {
    enabled = true
    file = "${params.outdir}/logs/report.html"
}

trace {
    enabled = true
    file = "${params.outdir}/logs/trace.txt"
}

// 环境配置
env {
    TMPDIR = '/scratch/tmp'
    R_LIBS_USER = '/opt/software/R/library'
}

// 配置文件
profiles {
    local {
        process.executor = 'local'
        process.memory = '8.GB'
    }
    
    slurm {
        process.executor = 'slurm'
        process.queue = 'normal'
        process.clusterOptions = '--account=my_project'
    }
}
```

## 10. 配置文件的优先级

Nextflow 按以下顺序加载配置：
1. 系统默认配置
2. 用户主目录的 `~/.nextflow/config`
3. 项目目录的 `nextflow.config`
4. 命令行指定的配置文件 (`-c file`)
5. 命令行参数 (`--param value`)

**总结**：`nextflow.config` 是 Nextflow 工作流的核心，它提供了：
- 参数管理和验证
- 资源分配和优化
- 执行环境配置
- 多环境支持
- 监控和日志记录
- 错误处理和重试策略