#!/bin/bash

# --- input ---
rawPaths=$1 # 'file1|file2'
filterPaths=$2
sampleValues=$3
inputMinGenes=$4 # 100
tfidfMin=$5 # 1
run_soupx_r=$6
# --- run ---
IFS=' ' read -ra raws <<< "$rawPaths"   # 拆成数组 paths
IFS=' ' read -ra filters <<< "$filterPaths"   # 拆成数组 paths
IFS=' ' read -ra samples <<< "$sampleValues"   # 拆成数组 paths

echo "sample counts: ${#samples[@]}" 

for idx in ""${!samples[@]}""; do
    raw=${raws[$idx]}
    filter=${filters[$idx]}
    sample=${samples[$idx]}
    echo "rawPath=$raw  filterPath=$filter  sampleValue=$sample"
    Rscript $run_soupx_r \
    --raw_path $raw --filter_path $filter \
    --sample_name $sample --input_mingenes $inputMinGenes --tfidfMin $tfidfMin
done