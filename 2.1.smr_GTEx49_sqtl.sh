# !/bin/zsh
# This zsh script runs the smr for sQTL as exposure and GWAS summary as the outcome.
# NOTE: The sqtl parent fold must have only the folders of the tissues with sub-chr smr files.

# Define sQTL data directory (parent folder)
sqtl_parent_dir="/Users/leoarrow/project/ref/sqtl/smr_sqtl_full_data"

# Define smr output directory
output_dir="/Users/leoarrow/project/iridocyclitis/output/smr/sqtl/GTEx49"

# Define GWAS summary data directories
t1d_gwas="/Users/leoarrow/project/iridocyclitis/data/smr/t1d1.ma"
iri_gwas="/Users/leoarrow/project/iridocyclitis/data/smr/iri3.ma"

# Create output directory
mkdir -p "$output_dir"

# Function to run SMR analysis for a specific tissue and chromosome
run_smr() {
    local gwas_summary="$1"
    local tissue_dir="$2"
    local chr_num="$3"
    local output_dir="$4"
    local gwas_name=$(basename "$gwas_summary" .ma)
    local tissue_name=$(basename "$tissue_dir")
    local sqtl_path="$tissue_dir/$tissue_name.query.chr${chr_num}"
    local output_prefix="${tissue_name}_chr${chr_num}@${gwas_name}"
    local output_file="$output_dir/${output_prefix}"

    echo " >>> Running SMR for $tissue_name chr<$chr_num> and $gwas_name..."
    smr --bfile /Users/leoarrow/project/ref/1kg.v3/EUR \
        --gwas-summary "$gwas_summary" \
        --beqtl-summary "$sqtl_path" \
        --out "$output_file" \
        --thread-num 10
    echo " >>> SMR for $tissue_name chr<$chr_num> and $gwas_name completed."
}
export -f run_smr

# Run smr for type1d GWAS
echo ">>>>>>>>>>>>>>>>> Processing t1d1 GWAS... >>>>>>>>>>>>>>>>>"
find "$sqtl_parent_dir" -mindepth 1 -maxdepth 1 -type d | while read tissue_dir; do
    for chr_num in {1..22}; do
        echo "$tissue_dir $chr_num"
    done
done | parallel --colsep ' ' run_smr "$t1d_gwas" {1} {2} "$output_dir"

# Run smr for iri_meta GWAS
echo ">>>>>>>>>>>>>>>>> Processing iri3 GWAS... >>>>>>>>>>>>>>>>>"
find "$sqtl_parent_dir" -mindepth 1 -maxdepth 1 -type d | while read tissue_dir; do
    for chr_num in {1..22}; do
        echo "$tissue_dir $chr_num"
    done
done | parallel --colsep ' ' run_smr "$iri_gwas" {1} {2} $output_dir


'''
nohup bash ~/project/iridocyclitis/code/2.1.smr_GTEx49_sqtl.sh > parallel_running_records_sqtl.out 2>&1 &
jobs
kill -9 %1
ps aux | grep smr
pkill -f run_smr
'''

# single core version---decrapated
'''
# Function to run SMR analysis
run_smr_analysis() {
    local gwas_summary=$1 #第一个参数赋值给局部变量 gwas_summary
    local gwas_name=$(basename "$gwas_summary" _h19.ma)

    for tissue_dir in "$sqtl_parent_dir"/*; do
        tissue_name=$(basename "$tissue_dir")
        echo "############# Processing tissue: $tissue_name #############"
        for chr_num in {1..22}; do
            echo "- Processing chr_num..."
            sqtl_path=$tissue_dir/$tissue_name.query.chr${chr_num}
            output_prefix="${tissue_name}_chr${chr_num}@${gwas_name}"
            output_file="$output_dir/${output_prefix}"
            echo "Running SMR for $tissue_name chr<$chr_num> and $gwas_name..."
            smr --bfile /Users/leoarrow/project/ref/1kg.v3/EUR \
                --gwas-summary "$gwas_summary" \
                --beqtl-summary "$sqtl_path" \
                --out "$output_file" \
                --thread-num 10
        done
    done
}

# Run smr for type1d GWAS
echo ">>>>>>>>>>>>>>>>> Processing t1d1 GWAS... >>>>>>>>>>>>>>>>>"
run_smr_analysis "$t1d_gwas"

# Run smr for iri_meta GWAS
echo ">>>>>>>>>>>>>>>>> Processing iri3 GWAS... >>>>>>>>>>>>>>>>>"
run_smr_analysis "$iri_gwas"
'''

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo ">>>>>>>>>>>>>>>>>!!!!!! All tasked finished !!!!!!>>>>>>>>>>>>>>>>>>"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"