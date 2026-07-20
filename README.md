# KAST Pipeline (Nextflow)

Nextflow DSL2 workflow for window-based sequence scoring, quantile filtering,
genomic feature intersection, GO enrichment, and reporting. Used for
alignment-free detection of horizontally transferred regions.

Pipeline file: `main.nf` (self-contained; all Python and R is embedded)

## Getting help

```
nextflow run main.nf --help
```

Prints every parameter with its current default. Every run also logs the
resolved parameters at startup, including values that are derived rather than
set directly (`increment`, `entropy_incr`), so the log records what actually
ran.

## What the pipeline does

Seven stages:

1. **Stage 1: Windowing, KAST scoring, entropy**
   Splits the FASTA into windows, runs KAST (Manhattan distance) against the
   reference, and scans sub-windows for minimum k-mer entropy. Produces
   `<fasta>.dat` and a diagnostic PDF.

2. **Stage 2: Quantile filtering**
   Computes quantiles of the D2 score distribution and writes FASTA subsets:
   per-quantile bins, plus cumulative bins from both ends. 58 files in total.

3. **Stage 3: Coordinate intersection**
   Two modes, see below. Produces `.intersect` per bin.

4. **Stage 4: Gene and protein extraction**
   Pulls CDS `Name=` values out of the intersect output and maps them to
   protein IDs via `cds.list`.

5. **Stage 5: GO enrichment**
   Runs goatools per bin, with and without count propagation.

6. **Stage 6: Per-bin plots**
   Enrichment/purification counts and namespace abundances across bins.

7. **Stage 7: Report**
   Genome-wide score track, score distribution, GO enrichment figure,
   candidate region table, and a self-contained `report.html`.

## Requirements

Tools on `PATH`:

- Nextflow (tested on 24.10)
- Python 3
- Rscript (base R only, no extra packages)
- bedtools
- `kast`
- `bwa` and `samtools`, only if using `--stage3_mode bwa`

Python packages:

```
pip install pandas biopython matplotlib seaborn goatools
```

A virtualenv is the easiest route on macOS, since recent system and Homebrew
Pythons reject bare `pip install`:

```
python3 -m venv .venv
source .venv/bin/activate
pip install pandas biopython matplotlib seaborn goatools
```

The venv must be active in the shell that launches Nextflow; tasks inherit its
`PATH`.

### KAST on macOS

The KAST release binary is a Linux x86-64 static ELF. It will not run natively
on macOS, on Intel or Apple Silicon (Rosetta does not help; it translates
Mach-O, not ELF). Options:

- Run the pipeline on a Linux host.
- Wrap KAST in Docker. `setup_kast_and_run.sh` builds a `linux/amd64` image
  from the 1.0.1 release and installs a shim at `~/scripts/kast`, so the
  pipeline itself needs no changes. On Apple Silicon this runs under emulation,
  which is acceptable because KAST is called once per run.

### Filesystem

Nextflow's cache requires POSIX file locks. Launching from an exFAT/FAT32
volume (most USB sticks) fails with `Can't open cache DB`. Launch from a local
disk. `-w` alone does not fix this, because the cache lives in the launch
directory.

## Input files

Defaults assume all files are in the launch directory.

| File | Param | Format |
|---|---|---|
| `Ref_genome.fasta` | `--fasta`, `--ref_genome` | Genome FASTA |
| `genomic.gff` | `--gff` | GFF3, seqids matching the FASTA headers |
| `cds.list` | `--cds_list` | TSV, no header: `GENE_NAME<TAB>PROTEIN_ID` |
| `all.txt` | `--all_txt` | Population protein IDs, one per line |
| `go.txt` | `--go_txt` | TSV: `PROTEIN_ID<TAB>GO:0006810;GO:0016020` |
| `go-basic.obo` | `--obo` | GO ontology |

`go-basic.obo` is required. goatools does not download it and errors with
`COULD NOT READ(go-basic.obo)` if absent:

```
wget http://geneontology.org/ontology/go-basic.obo
```

Stage 4 walks the chain GFF `Name=` -> `cds.list` column 1 -> column 2 ->
`all.txt` -> `go.txt`. A break anywhere in that chain produces empty results
rather than an error, so it is worth checking the IDs line up before a long
run.

## Parameters

Input files:

| Param | Default |
|---|---|
| `--fasta` | `Ref_genome.fasta` |
| `--ref_genome` | `Ref_genome.fasta` |
| `--gff` | `genomic.gff` |
| `--cds_list` | `cds.list` |
| `--all_txt` | `all.txt` |
| `--go_txt` | `go.txt` |
| `--obo` | `go-basic.obo` |

Windowing and scoring:

| Param | Default | Notes |
|---|---|---|
| `--window_size` | `100000` | Window length in bp |
| `--increment` | `window_size / 2` | Step between windows |
| `--kmer_length` | `4` | See memory note below |
| `--entropy_window` | `10000` | Must be smaller than `window_size` |
| `--entropy_incr` | `entropy_window / 2` | |

Coordinate mapping and reporting:

| Param | Default | Notes |
|---|---|---|
| `--stage3_mode` | `coords` | `coords` or `bwa` |
| `--reference_contig` | `null` | Required in `coords` mode for multi-contig references |
| `--top_quantile` | `0.95` | Bin treated as the candidate set in the report |

### k-mer length and memory

KAST allocates a dense 4^k array, so memory grows fourfold per unit of k:

| k | RAM |
|---|---|
| <= 12 | < 64 MB |
| 14 | ~1 GB |
| 15 | ~4 GB |
| 16 | ~16 GB |
| >= 20 | KAST refuses to start |

The workflow rejects `--kmer_length` outside 1-16 up front. KAST's own default
is 3.

### Entropy scan

The entropy scan looks for the lowest-entropy sub-window inside each window, so
`entropy_window` must be strictly smaller than `window_size`. If they are equal
there are zero sub-windows to scan and stage 1 fails.

The entropy value is used only for colouring the stage 1 diagnostic PDF and the
report's second panel. It does not affect which regions are called as
candidates: all filtering uses the D2 score alone. Changing `--entropy_window`
produces byte-identical filtered FASTA output.

## Stage 3 modes

**coords** (default) derives intervals arithmetically from the segment headers
written by stage 2:

```
start = (segment_id - 1) * increment
end   = start + window_size
```

then runs `bedtools intersect` against the GFF. No alignment, so it is fast,
but it assumes segment IDs map directly onto the configured window geometry.
With a multi-contig reference it aborts unless `--reference_contig` is given,
since segment indices are otherwise ambiguous.

**bwa** aligns the filtered windows back to the reference with `bwa mem`, sorts,
converts to BED, and intersects. Slower, but does not assume the window
geometry.

Both modes emit a 6-column BED before intersection, so stage 4 sees the same
column offsets either way.

## Example runs

Defaults:

```
nextflow run main.nf
```

Explicit k-mer length:

```
nextflow run main.nf --kmer_length 4
```

Multi-contig reference:

```
nextflow run main.nf --reference_contig NC_008095.1
```

Alignment-based coordinates:

```
nextflow run main.nf --stage3_mode bwa -resume
```

Full template:

```
nextflow run main.nf \
  --fasta Ref_genome.fasta \
  --ref_genome Ref_genome.fasta \
  --gff genomic.gff \
  --all_txt all.txt \
  --go_txt go.txt \
  --cds_list cds.list \
  --obo go-basic.obo \
  --window_size 100000 \
  --increment 50000 \
  --kmer_length 4 \
  --stage3_mode coords \
  -with-report report.html -with-trace
```

## Outputs

Published under `results/`:

```
results/dist/                 per-bin GO enrichment tables (*.go.prop, *.go.noprop)
results/dist/GO_results/      per-bin plots
results/dist/CumulativeGO_results/
results/report/report.html    self-contained summary, open in a browser
results/report/01_genome_track.png
results/report/02_score_distribution.png
results/report/03_go_top_terms.png
results/report/candidate_regions.tsv
results/report/go_summary.tsv
results/report/summary.txt
```

Intermediate files stay in the work directory: `<fasta>.dat` (D2 score, window
number, entropy), `win.<fasta>`, `<fasta>.contam.pdf`, `Filtered_sequences/`,
`*.intersect`, `*.gene`, `*.protein`.

`candidate_regions.tsv` lists windows above the `--top_quantile` threshold with
coordinates, scores, gene counts and gene names. `go_summary.tsv` flattens
significant terms from all bins into one table.

`publishDir` uses copy mode, so re-running overwrites `results/`. Move or rename
it first if comparing parameter settings.

## Notes and assumptions

- `coords` mode assumes segment IDs map directly onto the configured
  `window_size` and `increment`.
- Stage 2 assumes a single-contig reference: the window index is per-record, so
  a multi-contig reference would need a global counter.
- Coordinates follow BED conventions.
- Stage 4 expects standard GFF fields in the intersect output, with feature type
  in column 9 and attributes in column 15.
- Stage 5 runs goatools twice per bin (116 invocations for 58 bins), each
  re-parsing the 31 MB ontology. This dominates runtime.
- Killing a run mid-stage can leave partially staged symlinks, after which
  `-resume` reports `ln: failed to create symbolic link ... File exists`. Delete
  the task directory named in the error and resume.
