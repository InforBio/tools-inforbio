#!/usr/bin/env bash
# Run true-replicate, pooled-pseudoreplicate, and self-consistency IDR branches.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_idr_pipeline.sh --sample NAME --output-dir DIR
      --treatment-bam BAM --treatment-bam BAM
      --control-bam BAM --control-bam BAM
      [--narrowpeak PEAK --narrowpeak PEAK]
      [--replicate-label LABEL --replicate-label LABEL]
      [MACS2 and IDR options]

Input modes:
  BAM only
      Supply exactly two treatment BAMs and two matched control BAMs. MACS2 is
      run on the original BAMs for the true-replicate comparison.

  BAM + narrowPeak
      Additionally supply exactly two narrowPeak files in replicate order.
      They are reused for the true-replicate comparison. The BAMs are still
      required to construct pooled and self pseudoreplicates. Supply the same
      MACS2 parameters that were used to create the narrowPeak files.

MACS2 options (defaults retained from the original scripts):
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
  --threads N                 samtools threads [GALAXY_SLOTS or 1]
  --seed N                    pseudo-replicate assignment seed [0]
EOF
}

# Error reporting function. Prints all arguments passed to die as an error message.
# $* = all arguments passed to die, separated by spaces
die() {
    echo "ERROR: $*" >&2
    exit 1
}

# function to check that an option has a value. Exits with an error if not.
# $1 = option name, $2 = option value, $# = number of arguments recieved
require_value() {
    [[ $# -ge 2 ]] || die "option $1 requires a value"
}

sample_name=""
output_dir=""
treatment_bams=()
control_bams=()
narrowpeaks=()
replicate_labels=()
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

###### command line parsing ######

while [[ $# -gt 0 ]]; do # while there are still arguments to process
    case "$1" in
        --sample) require_value "$@"; sample_name="$2"; shift 2 ;; # $@ means pass all remaining command-line arguments individually.
        --output-dir) require_value "$@"; output_dir="$2"; shift 2 ;;
        --treatment-bam) require_value "$@"; treatment_bams+=("$2"); shift 2 ;;
        --control-bam) require_value "$@"; control_bams+=("$2"); shift 2 ;;
        --narrowpeak) require_value "$@"; narrowpeaks+=("$2"); shift 2 ;;
        --replicate-label) require_value "$@"; replicate_labels+=("$2"); shift 2 ;;
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

#####################################


###### parameter validation ######


[[ -n "$sample_name" ]] || die "--sample is required"
[[ "$sample_name" =~ ^[A-Za-z0-9._-]+$ ]] || die "--sample may contain only letters, numbers, dots, underscores, and hyphens"
[[ -n "$output_dir" ]] || die "--output-dir is required"
[[ ${#treatment_bams[@]} -eq 2 ]] || die "supply exactly two --treatment-bam values"
[[ ${#control_bams[@]} -eq 2 ]] || die "supply exactly two --control-bam values"
[[ ${#narrowpeaks[@]} -eq 0 || ${#narrowpeaks[@]} -eq 2 ]] || die "supply either zero or exactly two --narrowpeak values"
[[ ${#replicate_labels[@]} -eq 0 || ${#replicate_labels[@]} -eq 2 ]] || die "supply either zero or exactly two --replicate-label values"
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
[[ "$seed" =~ ^-?[0-9]+$ ]] || die "--seed must be an integer"


#####################################


###### IDR ranking validation ######

case "$idr_rank" in
    p.value|q.value|signal.value) ;;
    *) die "--rank must be p.value, q.value, or signal.value" ;;
esac

#####################################

###### MACS2 parameter validation ######

# if labels of replicates are not provided, they are labelled rep1 and rep2
if [[ ${#replicate_labels[@]} -eq 0 ]]; then
    replicate_labels=("rep1" "rep2")
fi
[[ "${replicate_labels[0]}" != "${replicate_labels[1]}" ]] || die "replicate labels must be unique"
for label in "${replicate_labels[@]}"; do
    [[ "$label" =~ ^[A-Za-z0-9._-]+$ ]] || die "replicate labels may contain only letters, numbers, dots, underscores, and hyphens"
done


###### Input file validation ######

for path in "${treatment_bams[@]}" "${control_bams[@]}" "${narrowpeaks[@]}"; do
    [[ -f "$path" ]] || die "input file not found: $path"
done


###### Locating the three branch scripts ######

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
true_script="$script_dir/run_macs2_idr_tr.sh"
pooled_script="$script_dir/run_macs2_idr_pp.sh"
self_script="$script_dir/run_macs2_idr_sp.sh"
for script in "$true_script" "$pooled_script" "$self_script"; do
    [[ -f "$script" ]] || die "branch script not found: $script"
done

mkdir -p "$output_dir"


###### Construct common argument list for all three branches ######

common_args=(
    --sample "$sample_name"
    --format "$macs2_format"
    --genome-size "$genome_size"
    --bandwidth "$bandwidth"
    --mfold-low "$mfold_low"
    --mfold-high "$mfold_high"
    "--$cutoff_type" "$cutoff_value"  ############################################
    --rank "$idr_rank"
)
for extra_arg in "${macs2_extra_args[@]}"; do
    common_args+=(--macs2-extra-arg "$extra_arg")
done


###### Run the three branches ######

# By default, the pipeline assumes that only BAM files are supplied
input_mode="bam_only"
if [[ ${#narrowpeaks[@]} -eq 2 ]]; then
    input_mode="bam_and_narrowpeak"
fi

# Record the complete, normalized run configuration for reproducibility and for
# a future Galaxy wrapper to expose as a regular tabular output.
parameter_file="$output_dir/run_parameters.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'sample\t%s\n' "$sample_name"
    printf 'input_mode\t%s\n' "$input_mode"
    printf 'treatment_bam_1\t%s\n' "${treatment_bams[0]}"
    printf 'treatment_bam_2\t%s\n' "${treatment_bams[1]}"
    printf 'control_bam_1\t%s\n' "${control_bams[0]}"
    printf 'control_bam_2\t%s\n' "${control_bams[1]}"
    if [[ "$input_mode" == "bam_and_narrowpeak" ]]; then
        printf 'narrowpeak_1\t%s\n' "${narrowpeaks[0]}"
        printf 'narrowpeak_2\t%s\n' "${narrowpeaks[1]}"
    fi
    printf 'replicate_label_1\t%s\n' "${replicate_labels[0]}"
    printf 'replicate_label_2\t%s\n' "${replicate_labels[1]}"
    printf 'macs2_format\t%s\n' "$macs2_format"
    printf 'genome_size\t%s\n' "$genome_size"
    printf 'bandwidth\t%s\n' "$bandwidth"
    printf 'mfold_low\t%s\n' "$mfold_low"
    printf 'mfold_high\t%s\n' "$mfold_high"
    printf '%s\t%s\n' "$cutoff_type" "$cutoff_value"
    printf 'idr_rank\t%s\n' "$idr_rank"
    printf 'threads\t%s\n' "$threads"
    printf 'seed\t%s\n' "$seed"
    for extra_arg in "${macs2_extra_args[@]}"; do
        printf 'macs2_extra_arg\t%s\n' "$extra_arg"
    done
} > "$parameter_file"

echo "========================================================================"
echo "Complete IDR pipeline"
echo "Sample: $sample_name"
echo "Input mode: $input_mode"
echo "Output directory: $output_dir"
echo "========================================================================"

true_args=(
    "${common_args[@]}"
    --output-dir "$output_dir/true_replicates"
    --treatment-1 "${treatment_bams[0]}"
    --control-1 "${control_bams[0]}"
    --treatment-2 "${treatment_bams[1]}"
    --control-2 "${control_bams[1]}"
)
if [[ "$input_mode" == "bam_and_narrowpeak" ]]; then
    true_args+=(--peak-1 "${narrowpeaks[0]}" --peak-2 "${narrowpeaks[1]}")
fi
bash "$true_script" "${true_args[@]}"

bash "$pooled_script" \
    "${common_args[@]}" \
    --output-dir "$output_dir/pooled_pseudoreplicates" \
    --treatment-1 "${treatment_bams[0]}" \
    --control-1 "${control_bams[0]}" \
    --treatment-2 "${treatment_bams[1]}" \
    --control-2 "${control_bams[1]}" \
    --threads "$threads" \
    --seed "$seed"

for index in 0 1; do
    label="${replicate_labels[$index]}"
    bash "$self_script" \
        "${common_args[@]}" \
        --replicate-label "$label" \
        --output-dir "$output_dir/self_consistency/$label" \
        --treatment "${treatment_bams[$index]}" \
        --control "${control_bams[$index]}" \
        --threads "$threads" \
        --seed "$((seed + index))"
done

results_file="$output_dir/idr_results.tsv"
{
    printf 'branch\tcomparison\tidr_result\n'
    printf 'true_replicates\t%s_vs_%s\t%s\n' \
        "${replicate_labels[0]}" "${replicate_labels[1]}" \
        "$output_dir/true_replicates/idr/${sample_name}_tr_idr.txt"
    printf 'pooled_pseudoreplicates\tpooled\t%s\n' \
        "$output_dir/pooled_pseudoreplicates/idr/${sample_name}_pp_idr.txt"
    for label in "${replicate_labels[@]}"; do
        printf 'self_consistency\t%s\t%s\n' "$label" \
            "$output_dir/self_consistency/$label/idr/${sample_name}_${label}_sp_idr.txt"
    done
} > "$results_file"

echo "========================================================================"
echo "[OK] All IDR branches completed"
echo "Result manifest: $results_file"
echo "Run parameters: $parameter_file"
echo "========================================================================"
