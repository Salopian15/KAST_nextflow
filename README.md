# KAST Pipeline (Nextflow)

This repository contains a Nextflow DSL2 workflow for window-based sequence scoring, filtering, genomic feature intersection, and GO enrichment analysis.

Pipeline file:
- main.nf

## What The Pipeline Does

The workflow runs five stages:

1. Stage 1: Windowing + entropy/D2 scoring
- Splits input FASTA into windows.
- Runs KAST and entropy calculations.
- Produces score tables and a contamination PDF.

2. Stage 2: Quantile filtering
- Builds quantile-based FASTA subsets from D2 score distributions.
- Produces cumulative and non-cumulative filtered FASTA files.

3. Stage 3: Genomic coordinate intersection:
  - coords mode: skips alignment, derives BED coordinates from segment headers and intersects with GFF.

4. Stage 4: Gene/protein extraction
- Extracts CDS/protein identifiers from intersection output.

5. Stage 5: GO enrichment + plotting
- Runs GO enrichment with and without propagation.
- Produces summary tables and plots.

## Requirements

The workflow expects these tools to be available in your runtime environment:
- Nextflow (DSL2)
- Python 3
- Rscript
- bedtools
- KAST executable (kast)

Python packages used in embedded scripts include:
- pandas
- biopython
- numpy
- matplotlib
- seaborn
- goatools

## Input Files

Default filenames are configured in main.nf, but all are overridable by command-line params:
- Ref_genome.fasta
- genomic.gff
- all.txt
- go.txt
- cds.list

## Parameters

Core parameters:
- fasta: input FASTA
- window_size: window length (default 5000)
- kmer_length: KAST k-mer length (default 25)
- increment: window step size (default window_size / 2)
- ref_genome: reference FASTA
- gff: genome annotation GFF
- all_txt: GO population/background list
- go_txt: gene-to-GO association file
- cds_list: CDS to protein mapping table

Stage 3 configuration:
- stage3_mode: bwa or coords (default coords)
- reference_contig: contig name to use in coords mode when reference has multiple contigs (default null)

## Stage 3 Modes

### Mode 1: coords (default)

Use this when your stage2 FASTA headers are segment IDs like >123 and you already know window_size/increment.

Behavior:
- Reconstructs BED intervals as:
  - start = (segment_id - 1) * increment
  - end = start + window_size
- Runs bedtools intersect with your GFF.

Important:
- If your reference FASTA contains multiple contigs, set reference_contig explicitly.
- If your segment headers are not integer IDs, coords mode will fail by design.

### Mode 2: bwa

Use this for alignment-driven coordinates.

Behavior:
- Runs bwa index, bwa mem, samtools conversion/sort, bedtools bamtobed, and bedtools intersect.

## Example Runs

From the repository root:

Minimal run using defaults:
nextflow run main.nf

Explicit coords mode:
nextflow run main.nf \
  --stage3_mode coords \
  --fasta Ref_genome.fasta \
  --ref_genome Ref_genome.fasta \
  --gff genomic.gff

Coords mode with multi-contig reference:
nextflow run main.nf \
  --stage3_mode coords \
  --reference_contig chr1

Legacy alignment mode:
nextflow run main.nf \
  --stage3_mode bwa

## Outputs

Primary outputs include:
- Stage 1:
  - input.dat
  - win.input
  - input.contam.pdf
- Stage 2:
  - Filtered_sequences/*.fasta
- Stage 3:
  - *.intersect
- Stage 4:
  - *.intersect.gene
  - *.intersect.protein
- Stage 5:
  - *.go.noprop
  - *.go.prop
  - GO_results/
  - CumulativeGO_results/
  - published results under results/dist

## Notes And Assumptions

- coords mode assumes segment IDs map directly to your configured increment and window_size.
- Coordinate system follows BED conventions used by bedtools.
- Stage4 parsing expects standard GFF-like fields in bedtools intersect output.

## Troubleshooting

1. Error: Invalid --stage3_mode
- Use only bwa or coords.

2. Error in coords mode about multiple contigs
- Provide --reference_contig with a contig name present in your reference FASTA.

3. Error in coords mode about non-integer segment header
- Ensure filtered FASTA headers are integer segment IDs (for example: >123).

4. No GO output
- Verify all_txt and go_txt are valid and compatible with your protein IDs.

## Quick Command Template

Adjust paths as needed:
nextflow run main.nf \
  --fasta Ref_genome.fasta \
  --ref_genome Ref_genome.fasta \
  --gff genomic.gff \
  --all_txt all.txt \
  --go_txt go.txt \
  --cds_list cds.list \
  --window_size 5000 \
  --increment 2500 \
  --kmer_length 25 \
  --stage3_mode coords
