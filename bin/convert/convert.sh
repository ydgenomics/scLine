#!/bin/bash
# Date: 250819
# single_convert似乎有问题，不建议用

input_file=$1
tools=$2
layers=$3

ext="${input_file##*.}"
echo "input file extension is: $ext"
if [ "$tools" == "single_convert" ]; then
    echo "Converting in single_convert..."
    /software/conda/Anaconda/bin/Rscript /WDL/Convert/v1.0.1/convert_rdsAh5ad.R --input_file $input_file
    cp "$input_file" ./
elif [ "$tools" == "multi_convert" ]; then
    if [ "$ext" == "rds" ]; then
        echo "Converting rds to h5ad..."
        /software/conda/Anaconda/bin/Rscript /WDL/Convert/v1.0.1/convert_rdsAh5ad2.R --input_file $input_file --layers $layers
        python /WDL/Convert/v1.0.1/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
        echo "Copying rds file..."
        cp "$input_file" ./
        rm saved_layers.txt saved_paths.txt
    elif [ "$ext" == "h5ad" ]; then
        echo "Converting h5ad to rds..."
        python /WDL/Convert/v1.0.1/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
        /software/conda/Anaconda/bin/Rscript /WDL/Convert/v1.0.1/convert_rdsAh5ad2.R --input_file $input_file --layers $layers
        echo "Copying h5ad file..."
        cp "$input_file" ./
        rm saved_layers.txt saved_paths.txt
    else
        echo "Error: Unsupported file extension '$ext'. Only 'rds' and 'h5ad' are supported."
    fi
else
    echo "Error: Unsupported tools "$tools". Only 'single_convert' and 'multi_convert' are supported."
fi