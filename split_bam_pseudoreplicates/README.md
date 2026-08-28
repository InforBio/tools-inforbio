# Create IDR pseudoreplicates — Galaxy tool

This Galaxy tool splits a BAM file into two complementary pseudoreplicates for
a ChIP-seq IDR analysis.

## Guarantees

- Every BAM alignment record is written to exactly one output.
- All records sharing the same query name remain together.
- The two mates of paired-end fragments are never separated.
- The same BAM file and seed always produce the same partitions.
- Both output BAM files are coordinate-sorted, validated, and indexed.
- A TSV report describes the resulting split.

The partition is pseudorandom and approximately 50/50. The tool does not force
the two outputs to contain exactly the same number of alignment records. For
normal ChIP-seq library sizes, the relative difference should be small.

## Requirements

- Python 3.12
- samtools 1.22.1

The versions are declared in the Galaxy tool wrapper and can be installed by a
configured Conda dependency resolver.

## Installing the tool in Galaxy

Copy the `split_bam_pseudoreplicates` directory into the local Galaxy tools
directory. Add the following entry to `tool_conf.xml`:

```xml
<section id="idr_local" name="IDR local tools">
    <tool file="split_bam_pseudoreplicates/split_bam_pseudoreplicates.xml" />
</section>
```

Restart Galaxy. The tool will appear as **Create IDR pseudoreplicates**.

## Inputs and outputs

Input:

- one individual filtered BAM file for a self-consistency analysis; or
- one BAM file produced by merging biological replicates for a pooled
  pseudoreplicate analysis;
- a non-negative seed used to make the split reproducible.

Outputs:

- pseudoreplicate 1 BAM;
- pseudoreplicate 2 BAM;
- the index for each output BAM;
- a TSV split report.

## Pooled pseudoreplicate workflow

1. Coordinate-sort the ChIP and control BAM files with **Samtools sort**.
2. Merge all ChIP BAM files with **Samtools merge**.
3. Separately merge all matched control BAM files with **Samtools merge**.
4. Run this tool once on the pooled ChIP BAM.
5. Run this tool once on the pooled control BAM.
6. Run MACS2 on ChIP pseudoreplicate 1 versus control pseudoreplicate 1.
7. Run MACS2 on ChIP pseudoreplicate 2 versus control pseudoreplicate 2.
8. Sort the two narrowPeak files and compare them with IDR.

Do not merge ChIP and control BAM files together.

## Self-consistency workflow

Run the tool independently on each biological replicate and its matched
control. Call peaks on the two corresponding ChIP/control pseudoreplicate
pairs, then compare their narrowPeak files with IDR.

## Assignment method

The tool applies a stable keyed BLAKE2b hash to the complete BAM query name.
The selected hash bit assigns every query name to partition 1 or partition 2.
Consequently, every record sharing that query name—including mates, secondary
alignments, and supplementary alignments—is assigned to the same partition.

The partitions are complementary, but their sizes are only approximately
equal. This approach must not be replaced by two independent 50% downsampling
jobs: independent samples can overlap and do not form complementary
partitions.

## Running the program outside Galaxy

The underlying program requires `samtools` in `PATH`:

```bash
python split_bam_pseudoreplicates.py \
  --input input.bam \
  --output-a pseudo1.bam --output-b pseudo2.bam \
  --index-a pseudo1.bam.bai --index-b pseudo2.bam.bai \
  --report report.tsv --seed 12345 --threads 2
```

## Testing

The repository includes a small paired-end BAM fixture and a Galaxy wrapper
test. The test checks that all 16 alignment records are retained, split 8/8 for
the test seed, and reported as complementary, sorted, and indexed outputs.

## Scientific context

This tool only creates BAM pseudoreplicates. It does not run MACS2 or IDR.
Its outputs are intended for peak calling before the resulting ranked peak
lists are submitted to IDR.

IDR reference: Li Q, Brown JB, Huang H, Bickel PJ. *Measuring reproducibility
of high-throughput experiments*. Annals of Applied Statistics. 2011.
DOI: 10.1214/11-AOAS466.
