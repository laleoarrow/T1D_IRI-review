# !/bin/zsh
# This zsh script runs the smr for mQTL as exposure and GWAS summary as the outcome.

# Define mqtl data directory (parent folder)
mqtl_parent_dir="/Users/leoarrow/project/ref/mqtl/yanglab/Hatton_EUR" # for Hatton big meta
# mqtl_parent_dir="/Users/leoarrow/project/ref/mqtl/yanglab/LBC_BSGS_meta" # for LBC_BSGS_meta

# Define smr output directory
output_dir="/Users/leoarrow/project/iridocyclitis/output/smr/mqtl/Hatton" # for Hatton big meta
# output_dir="/Users/leoarrow/project/iridocyclitis/output/smr/mqtl/LBC" # for LBC_BSGS_meta

# Define GWAS summary data directories
t1d_gwas="/Users/leoarrow/project/iridocyclitis/data/smr/t1d1.ma"
iri_gwas="/Users/leoarrow/project/iridocyclitis/data/smr/iri3.ma"

# Create output directory
mkdir -p "$output_dir"

# Function to run SMR analysis
run_smr_analysis() {
    local gwas_summary=$1
    local gwas_name=$(basename "$gwas_summary" _h19.ma)
    
    for chr_num in {1..22}; do
        mqtl_path=$mqtl_parent_dir/EUR_chr${chr_num} # for Hatton big meta
        # mqtl_path=$mqtl_parent_dir/bl_mqtl_chr${chr_num} # for LBC_BSGS_meta
        tissue_name=$(basename "$mqtl_path")
        echo "############# Processing tissue: $tissue_name #############"
        output_prefix="${tissue_name}_mqtl@${gwas_name}"
        output_path="$output_dir/${output_prefix}"
        echo "Running SMR for $tissue_name and $gwas_name..."
        smr --bfile /Users/leoarrow/project/ref/1kg.v3/EUR \
            --gwas-summary "$gwas_summary" \
            --beqtl-summary "$mqtl_path" \
            --out "$output_path" \
            --thread-num 10
    done
}

# Run smr for type1d GWAS
# echo "Processing type1d GWAS..."
# run_smr_analysis "$t1d_gwas"

# Run smr for iri_meta GWAS
echo "Processing iri_meta GWAS..."
run_smr_analysis "$iri_gwas"