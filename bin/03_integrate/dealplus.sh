i=$1
other1_key=$2
other2_key=$3
dealplus_r=$4
dealplus_py=$5

outfile=$(basename "$i")
type="${outfile##*.}"
if [ $type == "rds" ]
then
name1=${i%.rds}
name2=${name1##*/}
outrds="${name2}_dealplus.rds"
outpdf="${name2}_plus.pdf"
Rscript $dealplus_r \
--input_rds $i --out_rds $outrds --out_UMAP $outpdf \
--other1_key $other1_key --other2_key $other2_key
elif [ $type == "h5ad" ]
then
name1=${i%.h5ad}
name2=${name1##*/}
outh5ad="${name2}_dealplus.h5ad"
outpdf="${name2}_plus.pdf"
python $dealplus_py $i $outh5ad $outpdf \
--other1_key $other1_key --other2_key $other2_key
fi