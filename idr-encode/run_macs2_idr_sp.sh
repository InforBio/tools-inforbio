#!/usr/bin/env bash
#
# ChIP-seq Analysis Pipeline: MACS2 Peak Calling with IDR For Self-Pseudo Replicates

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_macs2_idr_sp.sh --sample NAME --replicate-label LABEL
      --treatment BAM --control BAM --output-dir DIR
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
  --seed N                    Accepted for orchestrator compatibility [0]
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
    local collated="$tmp_dir/${label}.collated.bam"
    local header="$tmp_dir/${label}.header.sam"
    local body_prefix="$tmp_dir/${label}.collated."
    local body_0="${body_prefix}00"
    local body_1="${body_prefix}01"
    local record_count
    local records_per_half

    samtools collate -@ "$threads" -o "$collated" "$input_bam"
    samtools view -H "$collated" > "$header"

    # Preserve the original self-consistency strategy: collate the BAM and
    # divide its alignment-record stream into two equal-sized chunks.
    record_count=$(samtools view "$collated" | wc -l)
    [[ "$record_count" -ge 2 ]] || die "cannot split $input_bam: fewer than two alignment records"
    records_per_half=$(( (record_count + 1) / 2 ))
    samtools view "$collated" | split -d -l "$records_per_half" - "$body_prefix"

    [[ -s "$body_0" && -s "$body_1" ]] || die "cannot split $input_bam into two nonempty pseudoreplicates"
    { cat "$header" "$body_0"; } | samtools view -@ "$threads" -bS - > "$output_0"
    { cat "$header" "$body_1"; } | samtools view -@ "$threads" -bS - > "$output_1"
}

sample_name=""
replicate_label=""
treatment_bam=""
control_bam=""
output_dir=""
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
        --replicate-label) require_value "$@"; replicate_label="$2"; shift 2 ;;
        --treatment) require_value "$@"; treatment_bam="$2"; shift 2 ;;
        --control) require_value "$@"; control_bam="$2"; shift 2 ;;
        --output-dir) require_value "$@"; output_dir="$2"; shift 2 ;;
        --format) require_value "$@"; macs2_format="$2"; shift 2 ;;
        --genome-size) require_value "$@"; genome_size="$2"; shift 2 ;;
        --bandwidth) require_value "$@"; bandwidth="$2"; shift 2 ;;
        --mfold-low) require_value "$@"; mfold_low="$2"; shift 2 ;;
        --mfold-high) require_value "$@"; mfold_high="$2"; shift 2 ;;
        --pvalue) require_value "$@"; cutoff_type="pvalue"; cutoff_value="$2"; shift 2 ;;
        --qvalue) require_value "$@"; cutoff_type="qvalue"; cutoff_value="$2"; shift 2 ;;
        --rank) require_value "$@"; idr_rank="$2"; shift 2 ;;
        --threads) require_value "$@"; threads="$2"; shift 2 ;;
        # The original self-consistency method is collate-then-split and has no
        # random seed. Accept this option because the orchestrator supplies it.
        --seed) require_value "$@"; seed="$2"; shift 2 ;;
        --macs2-extra-arg) require_value "$@"; macs2_extra_args+=("$2"); shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done

[[ -n "$sample_name" ]] || die "--sample is required"
[[ "$sample_name" =~ ^[A-Za-z0-9._-]+$ ]] || die "--sample may contain only letters, numbers, dots, underscores, and hyphens"
[[ -n "$replicate_label" ]] || die "--replicate-label is required"
[[ "$replicate_label" =~ ^[A-Za-z0-9._-]+$ ]] || die "--replicate-label may contain only letters, numbers, dots, underscores, and hyphens"
[[ -n "$output_dir" ]] || die "--output-dir is required"
[[ -n "$treatment_bam" && -f "$treatment_bam" ]] || die "treatment BAM not found: $treatment_bam"
[[ -n "$control_bam" && -f "$control_bam" ]] || die "control BAM not found: $control_bam"
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
[[ "$seed" =~ ^-?[0-9]+$ ]] || die "--seed must be an integer"

case "$idr_rank" in
    signal.value) rank_column=7 ;;
    p.value) rank_column=8 ;;
    q.value) rank_column=9 ;;
    *) die "--rank must be p.value, q.value, or signal.value" ;;
esac

for program in samtools macs2 idr sort split wc; do
    command -v "$program" >/dev/null 2>&1 || die "$program was not found on PATH"
done

start_time=$(date +%s)
mkdir -p "$output_dir/idr" "$output_dir/tmp"
tmp_dir=$(mktemp -d "$output_dir/tmp/self-pseudoreplicates.XXXXXX")
trap 'rm -rf -- "$tmp_dir"' EXIT

comparison="${sample_name}_${replicate_label}"
echo "========================================================================"
echo "Self-consistency IDR branch"
echo "Sample: $sample_name"
echo "Replicate: $replicate_label"
echo "Start time: $(date)"
echo "========================================================================"

treatment_00="$tmp_dir/${comparison}.treatment.00.bam"
treatment_01="$tmp_dir/${comparison}.treatment.01.bam"
control_00="$tmp_dir/${comparison}.control.00.bam"
control_01="$tmp_dir/${comparison}.control.01.bam"

echo "===> [$comparison] Creating treatment and control pseudoreplicates"
split_bam "$treatment_bam" "$treatment_00" "$treatment_01" "treatment"
split_bam "$control_bam" "$control_00" "$control_01" "control"

macs2_common=(
    --format "$macs2_format"
    --gsize "$genome_size"
    --bw "$bandwidth"
    --mfold "$mfold_low" "$mfold_high"
    "--$cutoff_type" "$cutoff_value"
    --tempdir "$tmp_dir"
)
macs2_common+=("${macs2_extra_args[@]}")

echo "===> [$comparison] MACS2 peak calling on self pseudoreplicate 00"
macs2 callpeak \
    -t "$treatment_00" \
    -c "$control_00" \
    --outdir "$output_dir" \
    --name "${comparison}.selfR00" \
    "${macs2_common[@]}"

echo "===> [$comparison] MACS2 peak calling on self pseudoreplicate 01"
macs2 callpeak \
    -t "$treatment_01" \
    -c "$control_01" \
    --outdir "$output_dir" \
    --name "${comparison}.selfR01" \
    "${macs2_common[@]}"

sorted_peak_00="$tmp_dir/${comparison}.selfR00.sorted.narrowPeak"
sorted_peak_01="$tmp_dir/${comparison}.selfR01.sorted.narrowPeak"
LC_ALL=C sort -k"$rank_column","$rank_column"nr "$output_dir/${comparison}.selfR00_peaks.narrowPeak" > "$sorted_peak_00"
LC_ALL=C sort -k"$rank_column","$rank_column"nr "$output_dir/${comparison}.selfR01_peaks.narrowPeak" > "$sorted_peak_01"

echo "===> [$comparison] IDR self-consistency analysis"
idr --samples \
    "$sorted_peak_00" \
    "$sorted_peak_01" \
    --input-file-type narrowPeak \
    --rank "$idr_rank" \
    --output-file "$output_dir/idr/${comparison}_sp_idr.txt" \
    --plot \
    --log-output-file "$output_dir/idr/${comparison}_sp_idr.log"

end_time=$(date +%s)
elapsed=$((end_time - start_time))
printf '[OK] Self-consistency branch completed in %02d:%02d:%02d\n' \
    "$((elapsed / 3600))" "$(((elapsed % 3600) / 60))" "$((elapsed % 60))"
printf 'IDR output: %s\n' "$output_dir/idr/${comparison}_sp_idr.txt"
