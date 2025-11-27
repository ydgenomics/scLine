# emapper_xlsx=$1
# gene_csv=$2
# ko_json=$3
# go_obo=$4
# minp=$5
# genus=$6
# species=$7
# build_orgdb_r=$8
# enrich_r=$9

emapper_xlsx="/data/work/scline/input/Galaxy6.xlsx"
ko_json="/data/work/scline/bin/06_enrich/ko00001.json"
go_obo="/data/work/scline/bin/06_enrich/go_obo_result.csv"
genus="genus"
species="species"
build_orgdb_r="/data/work/scline/bin/06_enrich/build_orgdb.R"


# pwd=$(pwd)
# emapper_xlsx=("${emapper_xlsx[@]/#/$pwd/}")

# genus="$(tr '[:lower:]' '[:upper:]' <<< "${genus:0:1}")${genus:1}" # Ensure that the first letter is capitalized
# echo "$genus"
Rscript $build_orgdb_r \
--emapper_xlsx $emapper_xlsx --ko_json $ko_json --go_obo $go_obo \
--taxid "1111" --genus $genus --species $species
# mv kegg_info.RData ../
# cd ..


genus="genus"
species="species"
kegg_info_RData="kegg_info.RData"
db="/data/work/org.Gspecies.eg.db"
gene_csv="/data/work/scline/results2/qc/group1/markers.csv/leiden_res_0.50.markers.csv"
minp=0.05
enrich_r="/data/work/scline/bin/06_enrich/enrich.R"

Rscript $enrich_r \
--gene_csv $gene_csv --kegg_info_RData $kegg_info_RData --db $db \
--minp $minp --genus $genus --species $species

# input_folder="${species}_enrich"
# folder_name=$(basename "$input_folder")
# tar -czvf "~{prefix}".tar.gz -C "$(dirname "$input_folder")" "$folder_name"