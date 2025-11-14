#!/usr/bin/env bash
TEMP=$(getopt -o b:g:f:s:u:l:e:c:m:t:n:r: --long biosample_value:,group_key:,filter_list:,splice_list:,unsplice_list:,sample_list:,input_mingenes:,input_mincells:,mitogenes_csv:,mito_threshold:,n_hvg:,rlst:,doublet_threshold:,rhotxt_list: -n "$0" -- "$@")
eval set -- "$TEMP"

while true; do
  case "$1" in
    -b|--biosample_value) biosample_value=$2; shift 2 ;;
    -g|--group_key) group_key=$2; shift 2 ;;
    -f|--filter_list) filter_list=$2; shift 2 ;;
    -s|--splice_list) splice_list=$2; shift 2 ;;
    -u|--unsplice_list) unsplice_list=$2; shift 2 ;;
    -l|--sample_list) sample_list=$2; shift 2 ;;
    -e|--input_mingenes) input_mingenes=$2; shift 2 ;;
    -c|--input_mincells) input_mincells=$2; shift 2 ;;
    -m|--mitogenes_csv) mitogenes_csv=$2; shift 2 ;;
    -t|--mito_threshold) mito_threshold=$2; shift 2 ;;
    -n|--n_hvg) n_hvg=$2; shift 2 ;;
    -r|--rlst) rlst=$2; shift 2 ;;
    -d|--doublet_threshold) doublet_threshold=$2; shift 2 ;;
    -o|--rhotxt_list) rhotxt_list=$2; shift 2 ;;
    --) shift; break ;;
    *) echo "Internal error!"; exit 1 ;;
  esac
done

# echo "biosample_value=$biosample_value  group_key=$group_key"

# 1. 拆成数组
IFS=, read -ra bv <<<"$biosample_value"
IFS=, read -ra sl  <<<"$sample_list"
IFS=, read -ra fl  <<<"$filter_list"
IFS=, read -ra sp  <<<"$splice_list"
IFS=, read -ra us  <<<"$unsplice_list"
IFS=, read -ra rh  <<<"$rhotxt_list"

# 检查长度是否一致
if (( ${#bv[@]} != ${#sl[@]} || ${#sl[@]} != ${#fl[@]} )); then
  echo "ERROR: The length of inputs of biosample_value, sample_value, and FilterMatrix is inconsistent！" >&2
  echo "  biosample_value (${#bv[@]}): ${bv[*]}" >&2
  echo "  sample_value/sample_list (${#sl[@]}): ${sl[*]}" >&2
  echo "  FilterMatrix/filter_list (${#fl[@]}): ${fl[*]}" >&2
  exit 1
fi

# 2. 按 biosample_value 分组（关联数组的 value ）
declare -A grp_sl grp_fl grp_sp grp_us grp_rh
for i in "${!bv[@]}"; do
    grp_sl["${bv[i]}"]+="${sl[i]},"
    grp_fl["${bv[i]}"]+="${fl[i]},"
    grp_sp["${bv[i]}"]+="${sp[i]},"
    grp_us["${bv[i]}"]+="${us[i]},"
    grp_rh["${bv[i]}"]+="${rh[i]},"
done
  
# 3. 对每组分别运算
for k in "${!grp_sl[@]}"; do
    # 去掉末尾多余的逗号
    sample_list=${grp_sl[$k]%,}
    filter_list=${grp_fl[$k]%,}
    splice_list=${grp_sp[$k]%,}
    # splice_list=${splice_list:-$filter_list} # If it is empty, instead of flter matrix
    unsplice_list=${grp_us[$k]%,}
    # unsplice_list=${unsplice_list:-$filter_list}
    rhotxt_list=${grp_rh[$k]%,}
    # rhotxt_list=${rhotxt_list:-temp.txt}

    echo "===== biosample: $k ====="
    echo "sample_list=$sample_list"
    echo "filter_list=$filter_list"
    echo "splice_list=$splice_list"
    echo "unsplice_list=$unsplice_list"
    echo "rhotxt_list=$rhotxt_list"
    
    mkdir $k; cd $k
    /opt/conda/bin/python /WDL/Dataget/v1.2.3/run_scrublet.py \
    --biosample_value $k --group_key $group_key --filter_list $filter_list --splice_list $splice_list --unsplice_list $unsplice_list \
    --sample_list $sample_list --input_mingenes $input_mingenes --mitogenes_csv $mitogenes_csv --mito_threshold $mito_threshold \
    --input_mincells $input_mincells --n_hvg $n_hvg --rlst $rlst --doublet_threshold $doublet_threshold
    mkdir temp
    for c in $(echo $rhotxt_list | tr ',' ' '); do
        cp $c ./temp
    done
    cat temp/*.txt >> summary.txt
    rm -rf temp
    cd ..
done