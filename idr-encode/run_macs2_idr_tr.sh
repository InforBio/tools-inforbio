#!/usr/bin/env bash
#
# ChIP-seq Analysis Pipeline: MACS2 Peak Calling with IDR For True Replicates

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_macs2_idr_tr.sh --sample NAME --output-dir DIR
      [--treatment-1 BAM --control-1 BAM --treatment-2 BAM --control-2 BAM]
      [--peak-1 NARROWPEAK --peak-2 NARROWPEAK]
      [MACS2 and IDR options]

If both --peak-1 and --peak-2 are supplied, they are reused and MACS2 is not
run for this branch. Otherwise, all four BAM options are required.

MACS2 options (original script defaults):
  --format FORMAT             BAMPE
  --genome-size SIZE          1.87e9
  --bandwidth BP              300
  --mfold-low N               2
  --mfold-high N              50
  --pvalue P                  0.01 (mutually exclusive with --qvalue)
  --qvalue Q                  Use an FDR cutoff instead of a p-value cutoff
  --macs2-extra-arg ARG       Append one literal MACS2 argument; repeat as needed

IDR options:
  --rank METHOD               p.value, q.value, or signal.value [p.value]
  --idr-threshold P           Global IDR cutoff used to retain peaks [0.05]
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_value() {
    [[ $# -ge 2 ]] || die "option $1 requires a value"
}

sample_name=""
output_dir=""
treatment_1=""
treatment_2=""
control_1=""
control_2=""
peak_1=""
peak_2=""
macs2_format="BAMPE"
genome_size="1.87e9"
bandwidth="300"
mfold_low="2"
mfold_high="50"
cutoff_type="pvalue"
cutoff_value="0.01"
idr_rank="p.value"
idr_threshold="0.05"
macs2_extra_args=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) require_value "$@"; sample_name="$2"; shift 2 ;;
        --output-dir) require_value "$@"; output_dir="$2"; shift 2 ;;
        --treatment-1) require_value "$@"; treatment_1="$2"; shift 2 ;;
        --treatment-2) require_value "$@"; treatment_2="$2"; shift 2 ;;
        --control-1) require_value "$@"; control_1="$2"; shift 2 ;;
        --control-2) require_value "$@"; control_2="$2"; shift 2 ;;
        --peak-1) require_value "$@"; peak_1="$2"; shift 2 ;;
        --peak-2) require_value "$@"; peak_2="$2"; shift 2 ;;
        --format) require_value "$@"; macs2_format="$2"; shift 2 ;;
        --genome-size) require_value "$@"; genome_size="$2"; shift 2 ;;
        --bandwidth) require_value "$@"; bandwidth="$2"; shift 2 ;;
        --mfold-low) require_value "$@"; mfold_low="$2"; shift 2 ;;
        --mfold-high) require_value "$@"; mfold_high="$2"; shift 2 ;;
        --pvalue) require_value "$@"; cutoff_type="pvalue"; cutoff_value="$2"; shift 2 ;;
        --qvalue) require_value "$@"; cutoff_type="qvalue"; cutoff_value="$2"; shift 2 ;;
        --rank) require_value "$@"; idr_rank="$2"; shift 2 ;;
        --idr-threshold) require_value "$@"; idr_threshold="$2"; shift 2 ;;
        --macs2-extra-arg) require_value "$@"; macs2_extra_args+=("$2"); shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done

[[ -n "$sample_name" ]] || die "--sample is required"
[[ "$sample_name" =~ ^[A-Za-z0-9._-]+$ ]] || die "--sample may contain only letters, numbers, dots, underscores, and hyphens"
[[ -n "$output_dir" ]] || die "--output-dir is required"

case "$idr_rank" in
    signal.value) rank_column=7 ;;
    p.value) rank_column=8 ;;
    q.value) rank_column=9 ;;
    *) die "--rank must be p.value, q.value, or signal.value" ;;
esac
[[ "$idr_threshold" =~ ^(0[.][0-9]*[1-9][0-9]*|1([.]0+)?)$ ]] || die "--idr-threshold must be greater than 0 and no greater than 1"

for value in "$genome_size" "$bandwidth" "$mfold_low" "$mfold_high" "$cutoff_value"; do
    [[ -n "$value" ]] || die "MACS2 parameter values cannot be empty"
done

if [[ -n "$peak_1" || -n "$peak_2" ]]; then
    [[ -n "$peak_1" && -n "$peak_2" ]] || die "--peak-1 and --peak-2 must be supplied together"
    [[ -f "$peak_1" ]] || die "narrowPeak file not found: $peak_1"
    [[ -f "$peak_2" ]] || die "narrowPeak file not found: $peak_2"
    input_peak_1="$peak_1"
    input_peak_2="$peak_2"
    input_mode="BAM + narrowPeak"
else
    for value in "$treatment_1" "$control_1" "$treatment_2" "$control_2"; do
        [[ -n "$value" ]] || die "all four BAM options are required when narrowPeak files are not supplied"
        [[ -f "$value" ]] || die "BAM file not found: $value"
    done
    input_mode="BAM only"
fi

command -v idr >/dev/null 2>&1 || die "idr was not found on PATH"
if [[ "$input_mode" == "BAM only" ]]; then
    command -v macs2 >/dev/null 2>&1 || die "macs2 was not found on PATH"
fi

start_time=$(date +%s)
mkdir -p "$output_dir/idr" "$output_dir/tmp"
tmp_dir=$(mktemp -d "$output_dir/tmp/true-replicates.XXXXXX")
trap 'rm -rf -- "$tmp_dir"' EXIT

echo "========================================================================"
echo "True-replicate IDR branch"
echo "Sample: $sample_name"
echo "Input mode: $input_mode"
echo "Start time: $(date)"
echo "========================================================================"

if [[ "$input_mode" == "BAM only" ]]; then
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

    echo "===> [$sample_name] MACS2 peak calling on biological replicate 1"
    macs2 callpeak \
        -t "$treatment_1" \
        -c "$control_1" \
        --outdir "$output_dir" \
        --name "${sample_name}_rep1" \
        "${macs2_common[@]}"

    echo "===> [$sample_name] MACS2 peak calling on biological replicate 2"
    macs2 callpeak \
        -t "$treatment_2" \
        -c "$control_2" \
        --outdir "$output_dir" \
        --name "${sample_name}_rep2" \
        "${macs2_common[@]}"

    input_peak_1="$output_dir/${sample_name}_rep1_peaks.narrowPeak"
    input_peak_2="$output_dir/${sample_name}_rep2_peaks.narrowPeak"
fi

sorted_peak_1="$tmp_dir/${sample_name}_rep1.sorted.narrowPeak"
sorted_peak_2="$tmp_dir/${sample_name}_rep2.sorted.narrowPeak"
LC_ALL=C sort -k"$rank_column","$rank_column"nr "$input_peak_1" > "$sorted_peak_1"
LC_ALL=C sort -k"$rank_column","$rank_column"nr "$input_peak_2" > "$sorted_peak_2"

echo "===> [$sample_name] IDR analysis for true replicates"
idr --samples \
    "$sorted_peak_1" \
    "$sorted_peak_2" \
    --input-file-type narrowPeak \
    --rank "$idr_rank" \
    --idr-threshold "$idr_threshold" \
    --output-file "$output_dir/idr/${sample_name}_tr_idr.txt" \
    --plot \
    --log-output-file "$output_dir/idr/${sample_name}_tr_idr.log"

end_time=$(date +%s)
elapsed=$((end_time - start_time))
printf '[OK] True-replicate branch completed in %02d:%02d:%02d\n' \
    "$((elapsed / 3600))" "$(((elapsed % 3600) / 60))" "$((elapsed % 60))"
printf 'IDR output: %s\n' "$output_dir/idr/${sample_name}_tr_idr.txt"
