# 01_qc
echo "01_qc"

echo "Run SoupX"
sh ./01_qc/soupx.sh $rawPaths $filterPaths $sampleValues $inputMinGenes $tfidfMin

echo "Run scrublet"
sh ./01_qc/scrublet.sh \
--biosample_value $biosample_value --group_key "sample" --filter_list $filter_list --splice_list "$splice_list" \
--unsplice_list "$unsplice_list" --sample_list $sample_list --input_mingenes 100 --mitogenes_csv $mitogenes_csv --mito_threshold 5 \
--input_mincells 3 --n_hvg 3000 --rlst "0.2,0.5,1.0,1.5" --doublet_threshold $doublet_threshold --rhotxt_list "$rhotxt_list"

# 02_anno
echo "02_anno"

echo "Run scType"


echo "Run singleR"


echo "Run Anno"


# 03_enrich



# 04_integration



# 05_metaneighbor


# 06_dea


# 07_pseudotime

# convert

# align


