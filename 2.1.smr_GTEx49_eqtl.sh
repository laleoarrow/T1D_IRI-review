#!/bin/zsh
# This zsh script runs the smr for 49 eqtl as exposure and 1 gwas summary as the outcome.

# Define eqtl data directory (parent folder)
eqtl_parent_dir="/Users/leoarrow/project/ref/eqtl/GTEx/GTEx_V8_cis_eqtl_summary"

# Define smr output directory
output_dir="/Users/leoarrow/project/iridocyclitis/output/smr/eqtl/GTEx49"

# Define GWAS summary data directories
t1d_gwas="/Users/leoarrow/project/iridocyclitis/data/smr/t1d1.ma"
iri_gwas="/Users/leoarrow/project/iridocyclitis/data/smr/iri3.ma"

# Create output directory
mkdir -p "$output_dir"

# run_smr
run_smr() {
    eqtl_dir="$1"
    gwas="$2"
    output_dir="$3"
    eqtl_name=$(basename "$eqtl_dir")
    eqtl_file="$eqtl_dir/${eqtl_name}"
    
    if [[ $(basename "$gwas") == "t1d1.ma" ]]; then
       output_prefix="${eqtl_name}@t1d1"
       else
       output_prefix="${eqtl_name}@iri3"
    fi
    
    output_file="$output_dir/${output_prefix}"
    echo "Running SMR for $eqtl_name and $(basename "$gwas")..."
    smr --bfile /Users/leoarrow/project/ref/1kg.v3/EUR \
        --gwas-summary "$gwas" \
        --beqtl-summary "$eqtl_file" \
        --out "$output_file" \
        --thread-num 2
    echo "SMR for $eqtl_name and $(basename "$gwas") completed."
}
export -f run_smr

echo " >>>>>>>>>>>>>>>>>> Going for t1d >>>>>>>>>>>>>>>>>> "
find "$eqtl_parent_dir" -mindepth 1 -maxdepth 1 -type d | parallel run_smr {} "$t1d_gwas" $output_dir

echo " >>>>>>>>>>>>>>>>>> Going for iri >>>>>>>>>>>>>>>>>> "
find "$eqtl_parent_dir" -mindepth 1 -maxdepth 1 -type d | parallel run_smr {} "$iri_gwas" $output_dir

# run example
# nohup bash ~/project/iridocyclitis/code/2.1.smr_GTEx49_eqtl.sh > parallel_running_records.out 2>&1 &
# jobs
# kill -9 %1
# ps aux | grep smr
# pkill -f run_smr

# single core version---decrapated
# Run smr for type1d GWAS
# echo "Processing type1d GWAS..."
# for eqtl_dir in "$eqtl_parent_dir"/*; do
#     if [ -d "$eqtl_dir" ]; then  # Check if it's a directory
#         eqtl_name=$(basename "$eqtl_dir")
#         eqtl_file="$eqtl_dir/${eqtl_name}"  # Assuming the file name is the same as the folder name without extension
#         output_prefix="${eqtl_name}@t1d1"
#         output_file="$output_dir/${output_prefix}.eqtl"
#         echo "Running SMR for $eqtl_name and type1d..."
#         smr --bfile /Users/leoarrow/project/ref/1kg.v3/EUR \
#             --gwas-summary "$t1d_dir" \
#             --beqtl-summary "$eqtl_file" \
#             --out "$output_file" \
#             --thread-num 10
#         echo "SMR for $eqtl_name and t1d1 completed."
#     fi
# done

# Run smr for iri_meta GWAS
# echo "Processing iri_meta GWAS..."
# for eqtl_dir in "$eqtl_parent_dir"/*; do
#     if [ -d "$eqtl_dir" ]; then  # Check if it's a directory
#         eqtl_name=$(basename "$eqtl_dir")
#         eqtl_file="$eqtl_dir/${eqtl_name}"  # Assuming the file name is the same as the folder name without extension
#         output_prefix="${eqtl_name}@iri3"
#         output_file="$output_dir/${output_prefix}.eqtl"
#         echo "Running SMR for $eqtl_name and iri_meta..."
#         smr --bfile /Users/leoarrow/project/ref/1kg.v3/EUR \
#             --gwas-summary "$iri_dir" \
#             --beqtl-summary "$eqtl_file" \
#             --out "$output_file" \
#             --thread-num 10
#         echo "SMR for $eqtl_name and iri3 completed."
#     fi
# done