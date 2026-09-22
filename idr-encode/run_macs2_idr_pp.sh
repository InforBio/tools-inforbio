#!/usr/bin/env bash
#
# ChIP-seq Analysis Pipeline: MACS2 Peak Calling with IDR For Pooled-Pseudo Replicates

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_macs2_idr_pp.sh --sample NAME --output-dir DIR
      --treatment-1 BAM --control-1 BAM --treatment-2 BAM --control-2 BAM
      [MACS2 and IDR options]

MACS2 options (original script defaults):
  --format FORMAT             BAMPE
  --genome-size SIZE          1.87e9
  --bandwidth BP              300
  --mfold-low N               2
  --mfold-high N              50
  --pvalue P                  0.01 (mutually exclusive with --qvalue)
  --qvalue Q                  Use an FDR cutoff instead of a p-value cutoff
  --macs2-extra-arg ARG       Append one literal MACS2 argument; repeat as needed

Other options:
  --rank METHOD               p.value, q.value, or signal.value [p.value]
  --threads N                 samtools threads [1]
  --seed N                    pooled record-shuffle seed [0]
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_value() {
    [[ $# -ge 2 ]] || die "option $1 requires a value"
}

split_bam() {
    local input_bam="$1"
    local output_0="$2"
    local output_1="$3"
    local label="$4"
    local split_seed="$5"
    local header="$tmp_dir/${label}.header.sam"
    local body_prefix="$tmp_dir/${label}."
    local body_0="${body_prefix}00"
    local body_1="${body_prefix}01"
    local record_count
    local records_per_half

    samtools view -H "$input_bam" > "$header"
    record_count=$(samtools view "$input_bam" | wc -l)
    [[ "$record_count" -ge 2 ]] || die "cannot split $input_bam: fewer than two alignment records"
    records_per_half=$(( (record_count + 1) / 2 ))

    # shuffle individual SAM alignment records, then split the stream into two
    # equal-sized chunks. Random keys make that strategy repeatable for a seed.
    samtools view "$input_bam" | \
        awk -v seed="$split_seed" 'BEGIN { srand(seed) } { printf "%.17f\t%012d\t%s\n", rand(), NR, $0 }' | \
        LC_ALL=C sort -T "$tmp_dir" -k1,1n -k2,2n | \
        cut -f3- | \
        split -d -l "$records_per_half" - "$body_prefix"

    [[ -s "$body_0" && -s "$body_1" ]] || die "cannot split $input_bam into two nonempty pseudoreplicates"
    { cat "$header" "$body_0"; } | samtools view -@ "$threads" -bS - > "$output_0"
    { cat "$header" "$body_1"; } | samtools view -@ "$threads" -bS - > "$output_1"
}

sample_name=""
output_dir=""
treatment_1=""
treatment_2=""
control_1=""
control_2=""
macs2_format="BAMPE"
genome_size="1.87e9"
bandwidth="300"
mfold_low="2"
mfold_high="50"
cutoff_type="pvalue"
cutoff_value="0.01"
idr_rank="p.value"
threads="${GALAXY_SLOTS:-1}"
seed="0"
macs2_extra_args=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) require_value "$@"; sample_name="$2"; shift 2 ;;
        --output-dir) require_value "$@"; output_dir="$2"; shift 2 ;;
        --treatment-1) require_value "$@"; treatment_1="$2"; shift 2 ;;
        --treatment-2) require_value "$@"; treatment_2="$2"; shift 2 ;;
        --control-1) require_value "$@"; control_1="$2"; shift 2 ;;
        --control-2) require_value "$@"; control_2="$2"; shift 2 ;;
        --format) require_value "$@"; macs2_format="$2"; shift 2 ;;
        --genome-size) require_value "$@"; genome_size="$2"; shift 2 ;;
        --bandwidth) require_value "$@"; bandwidth="$2"; shift 2 ;;
        --mfold-low) require_value "$@"; mfold_low="$2"; shift 2 ;;
        --mfold-high) require_value "$@"; mfold_high="$2"; shift 2 ;;
        --pvalue) require_value "$@"; cutoff_type="pvalue"; cutoff_value="$2"; shift 2 ;;
        --qvalue) require_value "$@"; cutoff_type="qvalue"; cutoff_value="$2"; shift 2 ;;
        --rank) require_value "$@"; idr_rank="$2"; shift 2 ;;
        --threads) require_value "$@"; threads="$2"; shift 2 ;;
        --seed) require_value "$@"; seed="$2"; shift 2 ;;
        --macs2-extra-arg) require_value "$@"; macs2_extra_args+=("$2"); shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done

[[ -n "$sample_name" ]] || die "--sample is required"
[[ "$sample_name" =~ ^[A-Za-z0-9._-]+$ ]] || die "--sample may contain only letters, numbers, dots, underscores, and hyphens"
[[ -n "$output_dir" ]] || die "--output-dir is required"
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
[[ "$seed" =~ ^-?[0-9]+$ ]] || die "--seed must be an integer"

case "$idr_rank" in
    signal.value) rank_column=7 ;;
    p.value) rank_column=8 ;;
    q.value) rank_column=9 ;;
    *) die "--rank must be p.value, q.value, or signal.value" ;;
esac

for value in "$treatment_1" "$control_1" "$treatment_2" "$control_2"; do
    [[ -n "$value" ]] || die "all four BAM options are required"
    [[ -f "$value" ]] || die "BAM file not found: $value"
done
for program in samtools macs2 idr awk sort cut split wc; do
    command -v "$program" >/dev/null 2>&1 || die "$program was not found on PATH"
done

start_time=$(date +%s)
mkdir -p "$output_dir/idr" "$output_dir/tmp"
tmp_dir=$(mktemp -d "$output_dir/tmp/pooled-pseudoreplicates.XXXXXX")
trap 'rm -rf -- "$tmp_dir"' EXIT

echo "========================================================================"
echo "Pooled-pseudoreplicate IDR branch"
echo "Sample: $sample_name"
echo "Start time: $(date)"
echo "========================================================================"

echo "===> [$sample_name] Pooling treatment and control BAMs"
pooled_treatment="$tmp_dir/${sample_name}.treatment.pooled.bam"
pooled_control="$tmp_dir/${sample_name}.control.pooled.bam"
samtools merge -@ "$threads" -f -u "$pooled_treatment" "$treatment_1" "$treatment_2"
samtools merge -@ "$threads" -f -u "$pooled_control" "$control_1" "$control_2"

treatment_00="$tmp_dir/${sample_name}.treatment.pooled.00.bam"
treatment_01="$tmp_dir/${sample_name}.treatment.pooled.01.bam"
control_00="$tmp_dir/${sample_name}.control.pooled.00.bam"
control_01="$tmp_dir/${sample_name}.control.pooled.01.bam"
split_bam "$pooled_treatment" "$treatment_00" "$treatment_01" "pooled-treatment" "$seed"
split_bam "$pooled_control" "$control_00" "$control_01" "pooled-control" "$((seed + 1))"

macs2_common=(
    --format "$macs2_format"
    --gsize "$genome_size"
    --bw "$bandwidth"
    --mfold "$mfold_low" "$mfold_high"
    "--$cutoff_type" "$cutoff_value"
    --tempdir "$tmp_dir"
)
if [[ ${#macs2_extra_args[@]} -gt 0 ]]; then
    macs2_common+=("${macs2_extra_args[@]}")
fi

echo "===> [$sample_name] MACS2 peak calling on pooled pseudoreplicate 00"
macs2 callpeak \
    -t "$treatment_00" \
    -c "$control_00" \
    --outdir "$output_dir" \
    --name "${sample_name}.pooledR00" \
    "${macs2_common[@]}"

echo "===> [$sample_name] MACS2 peak calling on pooled pseudoreplicate 01"
macs2 callpeak \
    -t "$treatment_01" \
    -c "$control_01" \
    --outdir "$output_dir" \
    --name "${sample_name}.pooledR01" \
    "${macs2_common[@]}"

sorted_peak_00="$tmp_dir/${sample_name}.pooledR00.sorted.narrowPeak"
sorted_peak_01="$tmp_dir/${sample_name}.pooledR01.sorted.narrowPeak"
LC_ALL=C sort -k"$rank_column","$rank_column"nr "$output_dir/${sample_name}.pooledR00_peaks.narrowPeak" > "$sorted_peak_00"
LC_ALL=C sort -k"$rank_column","$rank_column"nr "$output_dir/${sample_name}.pooledR01_peaks.narrowPeak" > "$sorted_peak_01"

echo "===> [$sample_name] IDR analysis for pooled pseudoreplicates"
idr --samples \
    "$sorted_peak_00" \
    "$sorted_peak_01" \
    --input-file-type narrowPeak \
    --rank "$idr_rank" \
    --output-file "$output_dir/idr/${sample_name}_pp_idr.txt" \
    --plot \
    --log-output-file "$output_dir/idr/${sample_name}_pp_idr.log"

end_time=$(date +%s)
elapsed=$((end_time - start_time))
printf '[OK] Pooled-pseudoreplicate branch completed in %02d:%02d:%02d\n' \
    "$((elapsed / 3600))" "$(((elapsed % 3600) / 60))" "$((elapsed % 60))"
printf 'IDR output: %s\n' "$output_dir/idr/${sample_name}_pp_idr.txt"
