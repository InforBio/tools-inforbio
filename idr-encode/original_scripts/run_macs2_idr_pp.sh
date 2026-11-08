#!/bin/bash
#


# ============================================================================
# ChIP-seq Analysis Pipeline: MACS2 Peak Calling with IDR For Pooled-Pseudo Replicates
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
if [ $# -ne 1 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <sample_name>"
    echo "Example: $0 H1"
    exit 1
fi

sample_name="$1"

# ============================================================================
# Define Variables
# ============================================================================
input_wt1="${sample_name}_WT1"
input_wt2="${sample_name}_WT2"
input_ko1="${sample_name}_KO1"
input_ko2="${sample_name}_KO2"

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
echo "WT files: ${input_wt1}.bam, ${input_wt2}.bam"
echo "KO file: ${input_ko1}.bam, ${input_ko2}.bam"
echo "Start time: $(date)"
echo "========================================================================"

# ============================================================================
# Check Input Files Exist
# ============================================================================
echo ""
echo "===> [${sample_name}] Checking input BAM files <==="

for bam_file in "${input_wt1}.bam" "${input_wt2}.bam" "${input_ko1}.bam" "${input_ko2}.bam"; do
    if [ ! -f "${data_dir}/${bam_file}" ]; then
        echo "ERROR: ${data_dir}/${bam_file} not found!"
        exit 1
    fi
    echo "[OK] Found: ${bam_file}"
done

# ============================================================================
# Create Pooled-pseudo Replicates - WT
# ============================================================================
echo ""
echo "===> [${sample_name}] Creating pooled-pseudo (PP) replicates for WT <==="
echo $(date)

# Merge WT BAMs
samtools merge -u "${tmp_dir}/${sample_name}_WT_merged.bam" \
    "${data_dir}/${input_wt1}.bam" "${data_dir}/${input_wt2}.bam"
samtools view -H "${tmp_dir}/${sample_name}_WT_merged.bam" > ${tmp_dir}/${sample_name}_WT.header


#Split merged WT BAMs
nlines=$(samtools view ${tmp_dir}/${sample_name}_WT_merged.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view ${tmp_dir}/${sample_name}_WT_merged.bam | \
    shuf - | split -d -l ${nlines} - "${tmp_dir}/${sample_name}_WT" # This will shuffle the lines in the file and split it into two SAM files
cat ${tmp_dir}/${sample_name}_WT.header ${tmp_dir}/${sample_name}_WT00 | \
    samtools view -bS - > ${out_dir}/${sample_name}_WT00.bam
cat ${tmp_dir}/${sample_name}_WT.header ${tmp_dir}/${sample_name}_WT01 | \
    samtools view -bS - > ${out_dir}/${sample_name}_WT01.bam
rm -f \
    ${tmp_dir}/${sample_name}_WT.header \
    ${tmp_dir}/${sample_name}_WT01 ${tmp_dir}/${sample_name}_WT00 \
    ${tmp_dir}/${sample_name}_WT_merged.bam

echo "[OK] Created pseudo-replicates for ${input_wt1} and ${input_wt2}"


# ============================================================================
# Create Pseudo-Replicates - KO
# ============================================================================

echo ""
echo "===> [${sample_name}] Creating pseudo-replicates for KO <==="
echo $(date)

#Merge KO BAMS
samtools merge -u "${tmp_dir}/${sample_name}_KO_merged.bam" \
    "${data_dir}/${input_ko1}.bam" "${data_dir}/${input_ko2}.bam"
samtools view -H "${tmp_dir}/${sample_name}_KO_merged.bam" > ${tmp_dir}/${sample_name}_KO.header

#Split merged KO BAMs
nlines=$(samtools view ${tmp_dir}/${sample_name}_KO_merged.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view ${tmp_dir}/${sample_name}_KO_merged.bam | \
    shuf - | split -d -l ${nlines} - "${tmp_dir}/${sample_name}_KO" # This will shuffle the lines in the file and split it into two SAM files
cat ${tmp_dir}/${sample_name}_KO.header ${tmp_dir}/${sample_name}_KO00 | \
    samtools view -bS - > ${out_dir}/${sample_name}_KO00.bam
cat ${tmp_dir}/${sample_name}_KO.header ${tmp_dir}/${sample_name}_KO01 | \
    samtools view -bS - > ${out_dir}/${sample_name}_KO01.bam
rm -f \
    ${tmp_dir}/${sample_name}_KO.header \
    ${tmp_dir}/${sample_name}_KO01 ${tmp_dir}/${sample_name}_KO00 \
    ${tmp_dir}/${sample_name}_KO_merged.bam

echo "[OK] Created pseudo-replicates for for ${input_ko1} and ${input_ko2}"


# ============================================================================
# MACS2 Peak Calling - Pooled-Pseudo Replicate 00
# ============================================================================
echo ""
echo "===> [${sample_name}] MACS2 peak calling on pooled-pseudo replicate 00 <==="
echo $(date)

macs2 callpeak \
    -t "${out_dir}/${sample_name}_WT00.bam" \
    -c "${out_dir}/${sample_name}_KO00.bam" \
    --format BAMPE \
    --gsize 1.87e9 \
    --outdir "${out_dir}/" \
    --name "${sample_name}.pooledR00" \
    --bw 300 \
    --mfold 2 50 \
    --pvalue 0.01 \
    --tempdir "${tmp_dir}"

rm -f ${out_dir}/${sample_name}_WT00.bam ${out_dir}/${sample_name}_KO00.bam

echo "[OK] MACS2 completed for pooled-pseudo replicate 00"

# ============================================================================
# MACS2 Peak Calling - Pseudo-replicate 01
# ============================================================================
echo ""
echo "===> [${sample_name}] MACS2 peak calling on pseudo-replicate 01 <==="
echo $(date)

macs2 callpeak \
    -t "${out_dir}/${sample_name}_WT01.bam" \
    -c "${out_dir}/${sample_name}_KO01.bam" \
    --format BAMPE \
    --gsize 1.87e9 \
    --outdir "${out_dir}/" \
    --name "${sample_name}.pooledR01" \
    --bw 300 \
    --mfold 2 50 \
    --pvalue 0.01 \
    --tempdir "${tmp_dir}"

rm -f ${out_dir}/${sample_name}_WT01.bam ${out_dir}/${sample_name}_KO01.bam

echo "[OK] MACS2 completed for pooled-pseudo replicate 01"

# ============================================================================
# IDR Analysis
# ============================================================================
echo ""
echo "===> [${sample_name}] Running IDR analysis <==="
echo $(date)

# Sort peaks by significance (column 8: -log10(p-value))
sort -k8,8nr "${out_dir}/${sample_name}.pooledR00_peaks.narrowPeak" > \
    "${tmp_dir}/${sample_name}.pooledR00.sorted_By_Column_8.narrowPeak"

sort -k8,8nr "${out_dir}/${sample_name}.pooledR01_peaks.narrowPeak" > \
    "${tmp_dir}/${sample_name}.pooledR01.sorted_By_Column_8.narrowPeak"

# Run IDR
idr --samples \
    "${tmp_dir}/${sample_name}.pooledR00.sorted_By_Column_8.narrowPeak" \
    "${tmp_dir}/${sample_name}.pooledR01.sorted_By_Column_8.narrowPeak" \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "${out_dir}/idr/${sample_name}_pp_idr.txt" \
    --plot \
    --log-output-file "${out_dir}/idr/${sample_name}_pp_idr.log"

echo "[OK] IDR analysis completed"

# ============================================================================
# Cleanup
# ============================================================================
# echo ""
# echo "===> [${sample_name}] Cleaning up temporary files <==="
# rm -f "${tmp_dir}/${sample_name}.pooled"*
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
