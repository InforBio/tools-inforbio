#!/bin/bash
#


# ============================================================================
# ChIP-seq Analysis Pipeline: MACS2 Peak Calling with IDR For Self-Pseudo Replicates
# ============================================================================

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Start timing
start_time=$(date +%s)

# ============================================================================
# Load Required Modules
# ============================================================================
module load samtools/1.9
module load macs2/2.2.7.1
module load idr/2.0.4.2


# ============================================================================
# Validate Input Arguments
# ============================================================================
if [ $# -ne 2 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <sample_name> <replicate_number>"
    echo "Example: $0 H1 1"
    exit 1
fi

sample_name="$1"
rep_num="$2"

# ============================================================================
# Define Variables
# ============================================================================
input_wt="${sample_name}_WT${rep_num}"
input_ko="${sample_name}_KO${rep_num}"

proj_dir="/shared/projects/chipseq_topo6"
data_dir="$proj_dir/data"
out_dir="$proj_dir/outputs"
tmp_dir="$out_dir/tmp"
log_dir="$proj_dir/logs"

# ============================================================================
# Create Output Directories
# ============================================================================
mkdir -p "$tmp_dir"
mkdir -p "$out_dir/idr"
mkdir -p "$log_dir"

echo "========================================================================"
echo "ChIP-seq Analysis Pipeline"
echo "========================================================================"
echo "Sample: $sample_name"
echo "Replicate: $rep_num"
echo "WT file: ${input_wt}.bam"
echo "KO file: ${input_ko}.bam"
echo "Start time: $(date)"
echo "========================================================================"

# ============================================================================
# Check Input Files Exist
# ============================================================================
echo ""
echo "===> [${sample_name}] Checking input BAM files <==="

for bam_file in "${input_wt}.bam" "${input_ko}.bam"; do
    if [ ! -f "${data_dir}/${bam_file}" ]; then
        echo "ERROR: ${data_dir}/${bam_file} not found!"
        exit 1
    fi
    echo "[OK] Found: ${bam_file}"
done

# ============================================================================
# Create Pseudo-Replicates - WT
# ============================================================================
echo ""
echo "===> [${sample_name}] Creating pseudo-replicates for WT <==="
echo $(date)

samtools collate -o "${tmp_dir}/${input_wt}.collate.bam" "${data_dir}/${input_wt}.bam"
samtools view -H "${tmp_dir}/${input_wt}.collate.bam" > "${tmp_dir}/${input_wt}.collate.header"

nlines=$(samtools view "${tmp_dir}/${input_wt}.collate.bam" | wc -l)
nlines=$(( (nlines + 1) / 2 ))

echo "Splitting ${input_wt} into two pseudo-replicates (${nlines} reads each)..."

samtools view "${tmp_dir}/${input_wt}.collate.bam" | \
    split -d -l ${nlines} - "${tmp_dir}/${input_wt}.collate."

cat "${tmp_dir}/${input_wt}.collate.header" "${tmp_dir}/${input_wt}.collate.00" | \
    samtools view -bS - > "${out_dir}/${input_wt}.collate.00.bam"

cat "${tmp_dir}/${input_wt}.collate.header" "${tmp_dir}/${input_wt}.collate.01" | \
    samtools view -bS - > "${out_dir}/${input_wt}.collate.01.bam"
rm -f \
    "${tmp_dir}/${input_wt}.collate.header" \
    "${tmp_dir}/${input_wt}.collate.00" \
    "${tmp_dir}/${input_wt}.collate.01" \
    "${tmp_dir}/${input_wt}.collate.bam"

echo "[OK] Created pseudo-replicates for ${input_wt}"

# ============================================================================
# Create Pseudo-Replicates - KO
# ============================================================================
echo ""
echo "===> [${sample_name}] Creating pseudo-replicates for KO <==="
echo $(date)

samtools collate -o "${tmp_dir}/${input_ko}.collate.bam" "${data_dir}/${input_ko}.bam"
samtools view -H "${tmp_dir}/${input_ko}.collate.bam" > "${tmp_dir}/${input_ko}.collate.header"

nlines=$(samtools view "${tmp_dir}/${input_ko}.collate.bam" | wc -l)
nlines=$(( (nlines + 1) / 2 ))

echo "Splitting ${input_ko} into two pseudo-replicates (${nlines} reads each)..."

samtools view "${tmp_dir}/${input_ko}.collate.bam" | \
    split -d -l ${nlines} - "${tmp_dir}/${input_ko}.collate."

cat "${tmp_dir}/${input_ko}.collate.header" "${tmp_dir}/${input_ko}.collate.00" | \
    samtools view -bS - > "${out_dir}/${input_ko}.collate.00.bam"

cat "${tmp_dir}/${input_ko}.collate.header" "${tmp_dir}/${input_ko}.collate.01" | \
    samtools view -bS - > "${out_dir}/${input_ko}.collate.01.bam"
rm -f \
    "${tmp_dir}/${input_ko}.collate.header" \
    "${tmp_dir}/${input_ko}.collate.00" \
    "${tmp_dir}/${input_ko}.collate.01" \
    "${tmp_dir}/${input_ko}.collate.bam"


echo "[OK] Created pseudo-replicates for KO"

# ============================================================================
# MACS2 Peak Calling - Pseudo-replicate 00
# ============================================================================
echo ""
echo "===> [${sample_name}] MACS2 peak calling on pseudo-replicate 00 <==="
echo $(date)

macs2 callpeak \
    -t "${out_dir}/${input_wt}.collate.00.bam" \
    -c "${out_dir}/${input_ko}.collate.00.bam" \
    --format BAMPE \
    --gsize 1.87e9 \
    --outdir "${out_dir}/" \
    --name "${input_wt}.selfR00" \
    --bw 300 \
    --mfold 2 50 \
    --pvalue 0.01 \
    --tempdir "${tmp_dir}"

rm -f ${out_dir}/${input_wt}.collate.00.bam ${out_dir}/${input_ko}.collate.00.bam

echo "[OK] MACS2 completed for pseudo-replicate 00"

# ============================================================================
# MACS2 Peak Calling - Pseudo-replicate 01
# ============================================================================
echo ""
echo "===> [${sample_name}] MACS2 peak calling on pseudo-replicate 01 <==="
echo $(date)

macs2 callpeak \
    -t "${out_dir}/${input_wt}.collate.01.bam" \
    -c "${out_dir}/${input_ko}.collate.01.bam" \
    --format BAMPE \
    --gsize 1.87e9 \
    --outdir "${out_dir}/" \
    --name "${input_wt}.selfR01" \
    --bw 300 \
    --mfold 2 50 \
    --pvalue 0.01 \
    --tempdir "${tmp_dir}"

rm -f ${out_dir}/${input_wt}.collate.01.bam ${out_dir}/${input_ko}.collate.01.bam

echo "[OK] MACS2 completed for pseudo-replicate 01"

# ============================================================================
# IDR Analysis
# ============================================================================
echo ""
echo "===> [${sample_name}] Running IDR analysis <==="
echo $(date)

# Sort peaks by significance (column 8: -log10(p-value))
sort -k8,8nr "${out_dir}/${input_wt}.selfR00_peaks.narrowPeak" > \
    "${tmp_dir}/${input_wt}.selfR00.sorted_By_Column_8.narrowPeak"

sort -k8,8nr "${out_dir}/${input_wt}.selfR01_peaks.narrowPeak" > \
    "${tmp_dir}/${input_wt}.selfR01.sorted_By_Column_8.narrowPeak"

# Run IDR
idr --samples \
    "${tmp_dir}/${input_wt}.selfR00.sorted_By_Column_8.narrowPeak" \
    "${tmp_dir}/${input_wt}.selfR01.sorted_By_Column_8.narrowPeak" \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "${out_dir}/idr/${input_wt}_sp_idr.txt" \
    --plot \
    --log-output-file "${out_dir}/idr/${input_wt}_sp_idr.log"

echo "[OK] IDR analysis completed"

# ============================================================================
# Cleanup
# ============================================================================
# echo ""
# echo "===> [${sample_name}] Cleaning up temporary files <==="
# rm -f "${tmp_dir}/${input_wt}."*
# rm -f "${tmp_dir}/${input_ko}."*
# echo "[OK] Cleanup completed (temporary files retained for debugging)"

# ============================================================================
# Summary and Timing
# ============================================================================
echo ""
echo "========================================================================"
echo "Pipeline Completed Successfully"
echo "========================================================================"

end_time=$(date +%s)
elapsed=$((end_time - start_time))
hours=$((elapsed / 3600))
minutes=$(((elapsed % 3600) / 60))
seconds=$((elapsed % 60))

echo "Sample: $sample_name"
echo "End time: $(date)"
printf "Total processing time: %02d:%02d:%02d (HH:MM:SS)\n" $hours $minutes $seconds
echo "========================================================================"
