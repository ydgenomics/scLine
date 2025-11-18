# --- input ---
rawPaths= $1 # 'file1|file2'
filterPaths= $2
sampleValues= $3
inputMinGenes= $4 # 100
tfidfMin= $5 # 1
# --- run ---
IFS='|' read -ra rawPaths <<< "$rawPaths"   # 拆成数组 paths
IFS='|' read -ra filterPaths <<< "$filterPaths"   # 拆成数组 paths
IFS='|' read -ra sampleValues <<< "$sampleValues"   # 拆成数组 paths
for idx in ""${!sampleValues[@]}""; do
    rawPath=${rawPaths[$idx]}
    filterPath=${filterPaths[$idx]}
    sampleValue=${sampleValues[$idx]}
    echo "rawPath=$rawPath  filterPath=$filterPath  sampleValue=$sampleValue"
    /opt/conda/bin/Rscript /bin/01_qc/run_soupx.R \
    --raw_path $rawPath --filter_path $filterPath \
    --sample_name $sampleValue --input_mingenes $inputMinGenes --tfidfMin $tfidfMin
done