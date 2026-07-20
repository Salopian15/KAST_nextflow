nextflow.enable.dsl=2

// ---- Parameters -----------------------------------------------------------
params.fasta            = "Ref_genome.fasta"
params.ref_genome       = "Ref_genome.fasta"
params.window_size      = 100000    // paper default
params.kmer_length      = 4         // KAST allocates 4^k; see --help for RAM scaling
params.increment        = null      // derived: window_size / 2 (= 50000) if unset
params.entropy_window   = 10000     // sub-window for entropy scan; MUST be < window_size
params.entropy_incr     = null      // derived: entropy_window / 2 (= 5000) if unset
params.gff              = "genomic.gff"
params.all_txt          = "all.txt"
params.go_txt           = "go.txt"
params.cds_list         = "cds.list"
params.obo              = "go-basic.obo"   // GO ontology; goatools will NOT fetch this for you
params.top_quantile     = "0.95"           // quantile bin treated as the HGT candidate set in the report
params.stage3_mode      = "coords"  // "coords" or "bwa"
params.reference_contig = null      // required for multi-contig refs in coords mode

params.help             = false

// ---- Help -----------------------------------------------------------------
def helpMessage() {
    log.info """
================================================================================
 Alignment-free HGT detection pipeline
================================================================================

 USAGE
   nextflow run main.nf [options]
   nextflow run main.nf --help

 INPUT FILES (defaults assume they sit in the launch directory)
   --fasta <file>            Genome to scan                  [${params.fasta}]
   --ref_genome <file>       Reference for KAST comparison   [${params.ref_genome}]
   --gff <file>              Annotation, GFF3                [${params.gff}]
   --cds_list <file>         TSV: GENE_NAME<TAB>PROTEIN_ID   [${params.cds_list}]
   --all_txt <file>          Population protein IDs          [${params.all_txt}]
   --go_txt <file>           TSV: PROTEIN_ID<TAB>GO:x;GO:y   [${params.go_txt}]
   --obo <file>              GO ontology (not auto-fetched)  [${params.obo}]

 WINDOWING
   --window_size <int>       Window length in bp             [${params.window_size}]
   --increment <int>         Step between windows            [window_size / 2]

 SCORING
   --kmer_length <int>       KAST k-mer length               [${params.kmer_length}]
                             KAST allocates a dense 4^k array:
                               k<=12  ~64 MB      k=14  ~1 GB
                               k=15   ~4 GB       k=16  ~16 GB
                             k>16 is impractical; k>=20 will refuse to run.
                             KAST's own default is 3.

 ENTROPY SCAN (sub-window scan inside each window; diagnostic plot only,
 it does NOT affect which regions are called as HGT candidates)
   --entropy_window <int>    Sub-window, MUST be < window_size [${params.entropy_window}]
   --entropy_incr <int>      Sub-window step                   [entropy_window / 2]

 COORDINATE MAPPING
   --stage3_mode <str>       'coords' (arithmetic, fast) or
                             'bwa' (align windows back)      [${params.stage3_mode}]
   --reference_contig <str>  Contig name; required in coords
                             mode when the reference has >1
                             contig                          [${params.reference_contig ?: 'AUTO'}]

 REPORTING
   --top_quantile <str>      Quantile bin treated as the
                             candidate set in the report     [${params.top_quantile}]

 OUTPUT
   results/dist/             Per-bin GO enrichment tables
   results/report/           report.html, figures, candidate_regions.tsv

 EXAMPLES
   nextflow run main.nf --kmer_length 4
   nextflow run main.nf --window_size 50000 --increment 25000
   nextflow run main.nf --reference_contig NC_008095.1
   nextflow run main.nf --stage3_mode bwa -resume
================================================================================
""".stripIndent()
}

if( params.help ) {
    helpMessage()
    exit 0
}

// ---- Stage 1: window, KAST Manhattan score, entropy, plot -----------------
process stage1 {
    input:
    path fasta
    val win_size
    val klen
    val incr
    val ent_win
    val ent_incr

    output:
    path "${fasta}.dat",         emit: dat_file
    path "win.${fasta}",         emit: fwin_file
    path "${fasta}.contam.pdf",  emit: pdf_file

    shell:
    '''
    #!/bin/bash
    set -eo pipefail

    FWIN="win.!{fasta}"
    OUT="!{fasta}.kast.out"

    cat << 'END_SCRIPT' > fastaWindowed3.py
#!/usr/bin/env python3
import sys
from operator import itemgetter

def processRecord(header, param, sequence, incr):
    new = header.replace(">", "")
    window = param
    if len(sequence) >= window:
        i = 0
        p1 = 0
        p2 = p1 + window
        while p2 < len(sequence):
            frag = sequence[p1:p2]
            p1 += incr
            p2 = p1 + window
            i += 1
            newHeader = f">{i}_{window}_{new}"
            print(newHeader)
            print(frag)

def main(afile, param, incr):
    sequence = ""
    header = ""
    with open(afile) as f_in:
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header != "":
                    processRecord(header, param, sequence, incr)
                header = line
                sequence = ""
            else:
                sequence += line
        if header != "":
            processRecord(header, param, sequence, incr)

def doit(argv):
    if len(argv) < 4:
        print("Chops a Fasta into windows e.g. size=2000, incr=1000")
        print("USAGE: $0 <fasta file> <windowSize> <incr>")
        sys.exit(0)
    afile = argv[1]
    param = int(argv[2])
    incr = int(argv[3])
    main(afile, param, incr)

doit(sys.argv)
END_SCRIPT

    cat << 'END_SCRIPT' > fastaEntropy3.py
#!/usr/bin/env python3
import sys
import math
from operator import itemgetter

def defineKmers(klen, nletters):
    if nletters == -4:
        letter = ["A", "C", "G", "T"]
    elif nletters == -40:
        letter = ["A", "C", "G", "U"]
        nletters = 4
    else:
        print("ERROR: wrong number of letters!", nletters)
        sys.exit(1)
    num = 4
    kmers = list(letter)
    for _ in range(klen - 1):
        new = []
        for j in range(num):
            for kmer in kmers:
                new.append(kmer + letter[j])
        kmers = list(new)
    return kmers

def countKmers(kmers, sequence, klen):
    kmerCount = {kmer: 0 for kmer in kmers}
    sequence = sequence.upper()
    maxi = len(sequence) - klen
    i = 0
    while i <= maxi:
        word = sequence[i:i + klen]
        if word in kmerCount:
            kmerCount[word] += 1
        i += 1
    return kmerCount

def calcEntropy(kmerCount, expectn, totalk):
    shannon1 = shannon2 = 0.0
    totalk = float(totalk)
    if totalk <= 0:
        return 0.0, 0.0
    expect = float(expectn) / totalk
    for word in kmerCount:
        count = float(kmerCount[word])
        if count > 0:
            kfreq = count / totalk
            shannon1 += kfreq * math.log(kfreq, 2)
            shannon2 += expect * math.log(kfreq, 2)
    return -shannon1, -shannon2

def processRecord(header, param, sequence, incr, mers, f_o):
    new = header.replace(">", "")
    window = param
    if len(sequence) < window:
        return
    p1 = 0
    p2 = p1 + window
    E = []
    while p2 < len(sequence):
        frag = sequence[p1:p2]
        p1 += incr
        p2 = p1 + window
        totH1 = totH2 = 0.0
        for kmers in mers:
            klen = len(kmers[0])
            kmerCount = countKmers(kmers, frag, klen)
            expectn = len(frag) / klen
            totalk = len(frag) - klen
            shan1, shan2 = calcEntropy(kmerCount, expectn, totalk)
            totH1 += shan1
            totH2 += shan2
        E.append([p1, totH1, totH2])
    if not E:
        return
    p1, totH1, totH2 = sorted(E, key=itemgetter(2))[0]
    ostr = str(p1) + "\\t" + str(totH1) + "\\t" + str(totH2) + "\\t" + new
    f_o.write(ostr + "\\n")

def main(afile, param, incr, f_o):
    nletters = -4
    mers = [defineKmers(k, nletters) for k in [1, 2, 3, 4, 5]]
    f_o.write("#pos\\tH1\\tH2\\tSeqID\\n")
    sequence = ""
    header = ""
    with open(afile) as f_in:
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header != "":
                    processRecord(header, param, sequence, incr, mers, f_o)
                header = line
                sequence = ""
            else:
                sequence += line
        if header != "":
            processRecord(header, param, sequence, incr, mers, f_o)

def doit(argv):
    if len(argv) < 4:
        print("Outputs lowest-entropy sub-window per sequence.")
        print("USAGE: $0 <fasta file> <window> <incr>")
        sys.exit(0)
    afile = argv[1]
    maxWin = int(argv[2])
    incr = int(argv[3])
    with open(afile + ".entropy", "w") as f_o:
        main(afile, maxWin, incr, f_o)

doit(sys.argv)
END_SCRIPT

    cat << 'END_SCRIPT' > HGT.R
args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]
myID <- args[3]

w1 <- read.csv(file=infile, sep=",", head=TRUE)

pdf(outfile)
rbPal <- colorRampPalette(c('green','red'), alpha=TRUE)

safeCut <- function(x) {
    if (length(unique(x[!is.na(x)])) < 2) {
        return(rep('black', length(x)))
    }
    rbPal(100)[as.numeric(cut(x, breaks = 100))]
}

w1$ColD <- safeCut(w1$D2score)
w1$ColE <- safeCut(w1$Entropy)

par(mfrow=c(2,1))
title1 <- paste(myID, "   K-mer stat")
title2 <- paste(myID, "   Entropy")
plot(w1$Number, w1$D2score, main=title1, xlab="Position", ylab="Stat1", pch=20, col=w1$ColE)
plot(w1$Number, w1$Entropy, main=title2, xlab="Position", ylab="Stat2", pch=20, col=w1$ColD)
invisible(dev.off())
END_SCRIPT

    command -v kast    >/dev/null 2>&1 || { echo "ERROR: 'kast' is not in PATH";    exit 1; }
    command -v python3 >/dev/null 2>&1 || { echo "ERROR: 'python3' is not in PATH"; exit 1; }
    command -v Rscript >/dev/null 2>&1 || { echo "ERROR: 'Rscript' is not in PATH"; exit 1; }

    python3 fastaWindowed3.py "!{fasta}" "!{win_size}" "!{incr}" > "$FWIN"

    kast -q "$FWIN" -r "!{fasta}" -t "manhattan" -k "!{klen}" -o "$OUT.1" -n 0 -f blastlike -nh

    python3 fastaEntropy3.py "$FWIN" "!{ent_win}" "!{ent_incr}"
    mv "$FWIN.entropy" "$OUT.3"

    echo "D2score,Number,Entropy" > "!{fasta}.dat"
    cut -f 2 "$OUT.1" | cut -d '_' -f 1 > tmp.2
    cut -f 6 "$OUT.1"                   > tmp.1
    grep -v '#' "$OUT.3" | cut -f 3     > tmp.3

    c1=$(wc -l < tmp.1); c2=$(wc -l < tmp.2); c3=$(wc -l < tmp.3)
    if [ "$c1" != "$c3" ] || [ "$c1" != "$c2" ]; then
        echo "ERROR: KAST rows ($c1) / segment rows ($c2) / entropy rows ($c3) differ — ordering mismatch"
        exit 1
    fi

    paste -d ',' tmp.1 tmp.2 tmp.3 >> "!{fasta}.dat"

    Rscript HGT.R "!{fasta}.dat" "!{fasta}.contam.pdf" "!{fasta}"
    '''
}

// ---- Stage 2: map scores back to sequences, quantile-filter ---------------
process stage2 {
    input:
    path fasta
    path dat_file
    path fwin_file

    output:
    path "Filtered_sequences/*.fasta", emit: filtered_seqs
    path "quantile_statistics.txt",    emit: quantile_stats

    shell:
    '''
    #!/bin/bash
    set -eo pipefail

    cat << 'END_SCRIPT' > kastdatWindowedfasta.py
# NOTE: assumes a single-contig reference (Number is a per-record window index).
from Bio import SeqIO
import pandas as pd
import sys, csv, os

score_file = sys.argv[1]
dna_file = sys.argv[2]
output_file = sys.argv[3]

if not os.path.exists(score_file):
    print(f"ERROR: Score file not found: {score_file}"); sys.exit(1)
if not os.path.exists(dna_file):
    print(f"ERROR: DNA file not found: {dna_file}"); sys.exit(1)

dna_list = [str(rec.seq) for rec in SeqIO.parse(dna_file, "fasta")]

rows = {"Segment": [], "D2 Score": [], "Sequence": []}
with open(score_file, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if 'Number' in row and 'D2score' in row:
            segment = int(row['Number'])
            if segment < 1 or segment > len(dna_list):
                continue
            rows['Segment'].append(segment)
            rows['D2 Score'].append(float(row['D2score']))
            rows['Sequence'].append(dna_list[segment - 1])

df = pd.DataFrame(rows)
df.to_csv(output_file, index=False)
END_SCRIPT

    cat << 'END_SCRIPT' > D2kaststats.py
import pandas as pd
import sys

csvfile = sys.argv[1]
output_file = "quantile_statistics.txt"
df = pd.read_csv(csvfile)

if 'D2score' not in df.columns:
    print("The 'D2score' column does not exist in the CSV file."); sys.exit(1)

quantiles = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
             0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
with open(output_file, 'w') as f:
    for q in quantiles:
        f.write(f"Quantile {q}: {df['D2score'].quantile(q)}\\n")
END_SCRIPT

    cat << 'END_SCRIPT' > D2_filterCumulative.py
import sys, os
import pandas as pd

if len(sys.argv) < 3:
    print("Usage: python D2_filterCumulative.py input.csv output_folder"); sys.exit(1)

csv_input = sys.argv[1]
output_folder = sys.argv[2]
os.makedirs(output_folder, exist_ok=True)

df = pd.read_csv(csv_input)

quantiles = {}
with open("quantile_statistics.txt") as fh:
    for line in fh:
        parts = line.split(":")
        quantiles[float(parts[0].split()[1])] = float(parts[1].strip())

def write_fasta(rows, path):
    with open(path, 'w') as out:
        for _, row in rows.iterrows():
            out.write(f">{row['Segment']}\\n{row['Sequence']}\\n")

qkeys = list(quantiles)
for i, quantile in enumerate(qkeys):
    lower = 0 if i == 0 else quantiles[qkeys[i - 1]]
    upper = quantiles[quantile]
    sel = df[(df["D2 Score"] >= lower) & (df["D2 Score"] <= upper)]
    write_fasta(sel, f"{output_folder}/quantile_{lower}_{quantile}.fasta")

lower = quantiles[qkeys[-1]]
write_fasta(df[df["D2 Score"] >= lower], f"{output_folder}/quantile_0.95_inf.fasta")

for quantile in qkeys:
    upper = quantiles[quantile]
    write_fasta(df[df["D2 Score"] <= upper], f"{output_folder}/Cumulative_0_{quantile}.fasta")

for quantile in qkeys:
    lower = quantiles[quantile]
    write_fasta(df[df["D2 Score"] >= lower], f"{output_folder}/Cumulative_{quantile}_inf.fasta")
END_SCRIPT

    command -v python3 >/dev/null 2>&1 || { echo "ERROR: 'python3' is not in PATH"; exit 1; }

    python3 kastdatWindowedfasta.py "!{dat_file}" "!{fwin_file}" "!{fasta}.csv"
    python3 D2kaststats.py "!{dat_file}"

    mkdir -p Filtered_sequences
    python3 D2_filterCumulative.py "!{fasta}.csv" Filtered_sequences
    '''
}

// ---- Stage 3 (bwa mode): align filtered windows, intersect with GFF -------
process stage3 {
    input:
    path fasta_file
    path ref_genome
    path gff

    output:
    path "${fasta_file.baseName}.intersect", emit: intersect_files

    shell:
    '''
    #!/bin/bash
    set -eo pipefail

    for tool in bwa samtools bedtools; do
        command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' is not in PATH"; exit 1; }
    done

    BASE="!{fasta_file.baseName}"
    bwa index "!{ref_genome}"
    bwa mem "!{ref_genome}" "!{fasta_file}" > "$BASE.sam"
    samtools view -S -b "$BASE.sam" > "$BASE.bam"
    samtools sort "$BASE.bam" -o "$BASE.sorted.bam"
    bedtools bamtobed -i "$BASE.sorted.bam" > "$BASE.sorted.bed"
    bedtools intersect -a "$BASE.sorted.bed" -b "!{gff}" -wa -wb > "$BASE.intersect"
    '''
}

// ---- Stage 3 (coords mode): map segment index -> reference coordinates -----
process stage3_coords {
    input:
    path fasta_file
    path ref_genome
    path gff
    val win_size
    val incr
    val reference_contig

    output:
    path "${fasta_file.baseName}.intersect", emit: intersect_files

    shell:
    '''
    #!/bin/bash
    set -eo pipefail

    cat << 'END_SCRIPT' > make_bed_from_segments.py
# Emits a 6-column BED (chrom,start,end,name,score,strand) so column offsets
# match the bwa/bamtobed path consumed by stage4.
import sys

if len(sys.argv) != 6:
    print("Usage: make_bed_from_segments.py <segments_fasta> <ref_fasta> <window_size> <increment> <contig|AUTO>")
    sys.exit(1)

segments_fasta = sys.argv[1]
ref_fasta = sys.argv[2]
window_size = int(sys.argv[3])
increment = int(sys.argv[4])
reference_contig = sys.argv[5]

def read_contigs(path):
    contigs = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                contigs.append(line[1:].strip().split()[0])
    return contigs

contigs = read_contigs(ref_fasta)
if not contigs:
    print("ERROR: No FASTA headers found in reference genome"); sys.exit(1)

if reference_contig == "AUTO":
    if len(contigs) != 1:
        print("ERROR: Multiple contigs in reference and --reference_contig not set")
        print("Contigs:", ",".join(contigs)); sys.exit(1)
    contig = contigs[0]
else:
    if reference_contig not in contigs:
        print(f"ERROR: reference_contig '{reference_contig}' not found in reference")
        print("Contigs:", ",".join(contigs)); sys.exit(1)
    contig = reference_contig

bed_out = f"{segments_fasta.rsplit('.', 1)[0]}.sorted.bed"
with open(segments_fasta) as in_fh, open(bed_out, "w") as out_fh:
    for line in in_fh:
        if not line.startswith(">"):
            continue
        header = line[1:].strip().split()[0]
        try:
            seg = int(header)
        except ValueError:
            print(f"ERROR: Segment header '{header}' is not an integer."); sys.exit(1)
        start = (seg - 1) * increment
        end = start + window_size
        if start < 0:
            print(f"ERROR: Negative start for segment {seg}"); sys.exit(1)
        out_fh.write(f"{contig}\\t{start}\\t{end}\\t{header}\\t.\\t+\\n")
END_SCRIPT

    command -v bedtools >/dev/null 2>&1 || { echo "ERROR: 'bedtools' is not in PATH"; exit 1; }
    command -v python3  >/dev/null 2>&1 || { echo "ERROR: 'python3' is not in PATH";  exit 1; }

    BASE="!{fasta_file.baseName}"
    python3 make_bed_from_segments.py "!{fasta_file}" "!{ref_genome}" "!{win_size}" "!{incr}" "!{reference_contig}" \
        || { echo "ERROR: make_bed_from_segments.py failed"; exit 1; }
    bedtools intersect -a "$BASE.sorted.bed" -b "!{gff}" -wa -wb > "$BASE.intersect"
    '''
}

// ---- Stage 4: extract CDS gene names, convert to protein IDs ---------------
process stage4 {
    input:
    path intersect_file
    path cds_list

    output:
    path "${intersect_file}.gene",    emit: gene_files
    path "${intersect_file}.protein", emit: protein_files

    shell:
    '''
    #!/bin/bash
    set -eo pipefail

    cat << 'END_SCRIPT' > ProteinIntersect.py
import pandas as pd
import sys, os

if len(sys.argv) != 3:
    print("Usage: ProteinIntersect.py input_file output_file"); sys.exit(1)

input_file, output_file = sys.argv[1], sys.argv[2]

# Empty intersect (window mapped nowhere) -> empty gene list, not an error.
if not os.path.exists(input_file) or os.path.getsize(input_file) == 0:
    open(output_file, 'w').close(); sys.exit(0)

try:
    df = pd.read_csv(input_file, sep='\\t', header=None)
except pd.errors.EmptyDataError:
    open(output_file, 'w').close(); sys.exit(0)

# Need at least the 6-col BED + 9-col GFF layout (type at 8, attributes at 14).
if df.shape[1] < 15:
    open(output_file, 'w').close(); sys.exit(0)

cds_df = df[df[8] == "CDS"]
unique_names = set()
for annotation in cds_df[14]:
    for field in str(annotation).split(';'):
        if field.startswith("Name="):
            unique_names.add(field.split('=', 1)[1])

with open(output_file, 'w') as outfile:
    for name in sorted(unique_names):
        outfile.write(name + "\\n")
END_SCRIPT

    cat << 'END_SCRIPT' > ProteinIDconversion.py
import pandas as pd
import sys, os

infile1, infile2, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

if not os.path.exists(infile2) or os.path.getsize(infile2) == 0:
    open(output_file, 'w').close(); sys.exit(0)

df1 = pd.read_csv(infile1, sep='\\t', header=None)
df2 = pd.read_csv(infile2, sep='\\t', header=None)
mapping = dict(zip(df1.iloc[:, 0], df1.iloc[:, 1]))

with open(output_file, 'w') as f:
    for value in df2.iloc[:, 0]:
        if value in mapping:
            f.write(str(mapping[value]) + "\\n")
END_SCRIPT

    command -v python3 >/dev/null 2>&1 || { echo "ERROR: 'python3' is not in PATH"; exit 1; }

    python3 ProteinIntersect.py "!{intersect_file}" "!{intersect_file}.gene"
    python3 ProteinIDconversion.py "!{cds_list}" "!{intersect_file}.gene" "!{intersect_file}.protein"
    '''
}

// ---- Stage 5: GO enrichment (propagated + non-propagated) ------------------
process stage5 {
    publishDir "results/dist", mode: 'copy'

    input:
    path protein_file
    path all_txt
    path go_txt
    path obo

    output:
    path "*.go.noprop", emit: noprop
    path "*.go.prop",   emit: prop

    shell:
    '''
    #!/bin/bash
    set -eo pipefail

    cat << 'END_SCRIPT' > find_enrichment1.py
#!/usr/bin/env python3
"""GO enrichment wrapper around goatools find_enrichment."""
import sys
import os.path as op
from goatools.cli.find_enrichment import GoeaCliArgs
from goatools.cli.find_enrichment import GoeaCliFnc

sys.path.insert(0, op.join(op.dirname(__file__), ".."))

def main():
    obj = GoeaCliFnc(GoeaCliArgs().args)
    results_specified = obj.get_results()
    obj.prt_results(results_specified)

if __name__ == "__main__":
    main()
END_SCRIPT

    command -v python3 >/dev/null 2>&1 || { echo "ERROR: 'python3' is not in PATH"; exit 1; }

    if [ -s "!{protein_file}" ]; then
        python3 find_enrichment1.py "!{protein_file}" "!{all_txt}" "!{go_txt}" --obo "!{obo}" --pval 0.05 --no_propagate_counts > "!{protein_file}.go.noprop"
        python3 find_enrichment1.py "!{protein_file}" "!{all_txt}" "!{go_txt}" --obo "!{obo}" --pval 0.05                     > "!{protein_file}.go.prop"
    else
        echo "# No proteins found - empty enrichment" > "!{protein_file}.go.noprop"
        echo "# No proteins found - empty enrichment" > "!{protein_file}.go.prop"
    fi
    '''
}

// ---- Stage 6: aggregate GO results into histograms -------------------------
process stage6 {
    publishDir "results/dist", mode: 'copy'

    input:
    path noprop_files
    path prop_files

    output:
    path "GO_results/*",           emit: go_plots,         optional: true
    path "CumulativeGO_results/*", emit: cumulative_plots, optional: true

    shell:
    '''
    #!/bin/bash
    set -eo pipefail

    cat << 'END_SCRIPT' > PlotGOHistograms.py
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import glob, io

def load(files):
    out = {}
    for file in files:
        with open(file) as f:
            lines = f.readlines()
        start = next((i for i, l in enumerate(lines) if l.startswith('GO\\tNS\\t')), None)
        if start is None:
            continue
        out[file] = pd.read_csv(io.StringIO(''.join(lines[start:])), sep='\\t')
    return out

nopropdfs = load([f for f in glob.glob('*.noprop') if not f.startswith('Cumulative')])
propdfs   = load([f for f in glob.glob('*.prop')   if not f.startswith('Cumulative')])

if not nopropdfs and not propdfs:
    print("No enrichment result files found - skipping histogram plots"); raise SystemExit

def counts(dfs):
    de, dp = [], []
    for name, df in dfs.items():
        try:
            q = float(name.split('_')[1])
        except (IndexError, ValueError):
            continue
        de.append({'Quantile Start': q, 'Count': df[df['enrichment'] == 'e'].shape[0]})
        dp.append({'Quantile Start': q, 'Count': df[df['enrichment'] == 'p'].shape[0]})
    return pd.DataFrame(de).sort_values('Quantile Start'), pd.DataFrame(dp).sort_values('Quantile Start')

for dfs, out, pal in [(nopropdfs, 'histogram_noprop.png', None),
                      (propdfs, 'histograms_prop.png', 'cool')]:
    if not dfs:
        continue
    de, dp = counts(dfs)
    fig, axs = plt.subplots(2)
    sns.barplot(x='Quantile Start', y='Count', data=de, ax=axs[0], color='steelblue')
    sns.barplot(x='Quantile Start', y='Count', data=dp, ax=axs[1], color='indianred')
    axs[0].set_title('e counts'); axs[1].set_title('p counts')
    fig.tight_layout(); plt.savefig(out); plt.close(fig)
END_SCRIPT

    cat << 'END_SCRIPT' > GOstatgraphs.py
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import glob, io

def load(files):
    out = {}
    for file in files:
        with open(file) as f:
            lines = f.readlines()
        start = next((i for i, l in enumerate(lines) if l.startswith('GO\\tNS\\t')), None)
        if start is None:
            continue
        out[file] = pd.read_csv(io.StringIO(''.join(lines[start:])), sep='\\t')
    return out

nopropdfs = load([f for f in glob.glob('*.noprop') if not f.startswith('Cumulative')])
propdfs   = load([f for f in glob.glob('*.prop')   if not f.startswith('Cumulative')])

if not nopropdfs and not propdfs:
    print("No enrichment result files found - skipping NS abundance plots"); raise SystemExit

data = []
def process(dfs, dtype):
    for name, df in dfs.items():
        try:
            q = float(name.split('_')[1])
        except (IndexError, ValueError):
            continue
        for enr in ['e', 'p']:
            counts = df[df['enrichment'] == enr]['NS'].value_counts(normalize=True)
            for ns in ['CC', 'MF', 'BP']:
                if ns in counts:
                    data.append({'Quantile Start': q, 'NS': ns,
                                 'Relative Abundance %': counts[ns],
                                 'Type': dtype, 'Enrichment': enr})

process(nopropdfs, 'noprop'); process(propdfs, 'prop')
df = pd.DataFrame(data)
for dtype in ['noprop', 'prop']:
    for enr in ['e', 'p']:
        sub = df[(df['Type'] == dtype) & (df['Enrichment'] == enr)]
        if sub.empty:
            continue
        plt.figure(figsize=(10, 6))
        sns.barplot(x='Quantile Start', y='Relative Abundance %', hue='NS',
                    data=sub, palette='Set2', errorbar=None)
        plt.title(f'Relative Abundances CC/MF/BP ({dtype}, {enr})')
        plt.savefig(f'RelAbundance_{dtype}_{enr}.png'); plt.close()
END_SCRIPT

    cat << 'END_SCRIPT' > PlotGOHistogramsCumulative.py
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import glob, io

def load(files):
    out = {}
    for file in files:
        with open(file) as f:
            lines = f.readlines()
        start = next((i for i, l in enumerate(lines) if l.startswith('GO\\tNS\\t')), None)
        if start is None:
            continue
        out[file] = pd.read_csv(io.StringIO(''.join(lines[start:])), sep='\\t')
    return out

def pct_from_tail(name):
    tok = name.split('_')[-1].split('.')
    return int(float(tok[0] + '.' + tok[1]) * 100)

def pct_from_field1(name):
    return int(float(name.split('_')[1]) * 100)

ctn = load(glob.glob('Cumulative_0_*.noprop'))
ctp = load(glob.glob('Cumulative_0_*.prop'))
cbn = load(glob.glob('*_inf.intersect.protein.go.noprop'))
cbp = load(glob.glob('*_inf.intersect.protein.go.prop'))

if not any([ctn, ctp, cbn, cbp]):
    print("No cumulative enrichment files - skipping cumulative plots"); raise SystemExit

def plot(dfs, keyfn, out, pal):
    if not dfs:
        return
    de, dp = [], []
    for name, df in dfs.items():
        try:
            k = keyfn(name)
        except (IndexError, ValueError):
            continue
        de.append({'Cumulative Start': k, 'Count': df[df['enrichment'] == 'e'].shape[0]})
        dp.append({'Cumulative Start': k, 'Count': df[df['enrichment'] == 'p'].shape[0]})
    de = pd.DataFrame(de).sort_values('Cumulative Start')
    dp = pd.DataFrame(dp).sort_values('Cumulative Start')
    fig, axs = plt.subplots(2, figsize=(16, 10))
    sns.barplot(x='Cumulative Start', y='Count', data=de, ax=axs[0], color='steelblue')
    sns.barplot(x='Cumulative Start', y='Count', data=dp, ax=axs[1], color='indianred')
    fig.tight_layout(); plt.savefig(out); plt.close(fig)

plot(ctn, pct_from_tail, 'histogram_cumulative_top_noprop.png', None)
plot(ctp, pct_from_tail, 'histogram_cumulative_top_prop.png', 'cool')
plot(cbn, pct_from_field1, 'histogram_cumulative_bottom_noprop.png', None)
plot(cbp, pct_from_field1, 'histogram_cumulative_bottom_prop.png', 'cool')
END_SCRIPT

    cat << 'END_SCRIPT' > GOstatgraphsCumulative.py
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import glob, io

def load(files):
    out = {}
    for file in files:
        with open(file) as f:
            lines = f.readlines()
        start = next((i for i, l in enumerate(lines) if l.startswith('GO\\tNS\\t')), None)
        if start is None:
            continue
        out[file] = pd.read_csv(io.StringIO(''.join(lines[start:])), sep='\\t')
    return out

nopropdfs = load(glob.glob('*.noprop'))
propdfs   = load(glob.glob('*.prop'))

if not nopropdfs and not propdfs:
    print("No cumulative enrichment files - skipping NS abundance plots"); raise SystemExit

data = []
def process(dfs, dtype):
    for name, df in dfs.items():
        try:
            q = float(name.split('_')[1])
        except (IndexError, ValueError):
            continue
        for enr in ['e', 'p']:
            counts = df[df['enrichment'] == enr]['NS'].value_counts(normalize=True)
            for ns in ['CC', 'MF', 'BP']:
                if ns in counts:
                    data.append({'Quantile Start': q, 'NS': ns,
                                 'Relative Abundance %': counts[ns],
                                 'Type': dtype, 'Enrichment': enr})

process(nopropdfs, 'noprop'); process(propdfs, 'prop')
df = pd.DataFrame(data)
for dtype in ['noprop', 'prop']:
    for enr in ['e', 'p']:
        sub = df[(df['Type'] == dtype) & (df['Enrichment'] == enr)]
        if sub.empty:
            continue
        plt.figure(figsize=(10, 6))
        sns.barplot(x='Quantile Start', y='Relative Abundance %', hue='NS',
                    data=sub, palette='Set2', errorbar=None)
        plt.title(f'Relative Abundances CC/MF/BP ({dtype}, {enr})')
        plt.savefig(f'CumRelAbundance_{dtype}_{enr}.png'); plt.close()
END_SCRIPT

    command -v python3 >/dev/null 2>&1 || { echo "ERROR: 'python3' is not in PATH"; exit 1; }

    mkdir -p GO_results CumulativeGO_results
    cp *.noprop *.prop GO_results/ 2>/dev/null || true
    cp Cumulative*.noprop Cumulative*.prop CumulativeGO_results/ 2>/dev/null || true

    ( cd GO_results && python3 ../PlotGOHistograms.py && python3 ../GOstatgraphs.py ) || true
    ( cd CumulativeGO_results && python3 ../PlotGOHistogramsCumulative.py && python3 ../GOstatgraphsCumulative.py ) || true
    '''
}


// ---- Stage 7: summary report (plots + tables + HTML) ----------------------
process stage7_report {
    publishDir "results", mode: 'copy'

    input:
    path dat_file
    path quantile_stats
    path intersect_files
    path noprop_files
    path prop_files
    val win_size
    val incr
    val top_q

    output:
    path "report/*", emit: report_files, optional: true

    shell:
    '''
    #!/bin/bash
    set -eo pipefail

    cat << 'END_SCRIPT' > make_report.py
#!/usr/bin/env python3
"""Summarise HGT pipeline output: genome track, score distribution,
candidate regions and GO enrichment -> report/ (PNGs, TSVs, report.html)."""
import base64
import glob
import io
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

WIN = int(sys.argv[1])
INCR = int(sys.argv[2])
TOPQ = sys.argv[3] if len(sys.argv) > 3 else "0.95"

OUT = "report"
os.makedirs(OUT, exist_ok=True)

ACCENT = "#c0392b"
BASE = "#3b6ea5"
GREY = "#8a8a8a"


def read_quantiles():
    q = {}
    if os.path.exists("quantile_statistics.txt"):
        for line in open("quantile_statistics.txt"):
            parts = line.split(":")
            if len(parts) == 2:
                try:
                    q[float(parts[0].split()[1])] = float(parts[1].strip())
                except ValueError:
                    pass
    return q


def read_dat():
    files = [f for f in glob.glob("*.dat")]
    if not files:
        return None
    df = pd.read_csv(files[0])
    if not {"D2score", "Number", "Entropy"} <= set(df.columns):
        return None
    df = df.sort_values("Number").reset_index(drop=True)
    df["start"] = (df["Number"] - 1) * INCR
    df["end"] = df["start"] + WIN
    df["mid"] = df["start"] + WIN // 2
    return df


def read_enrichment(path):
    """Return the significant-terms table from a goatools output file."""
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        lines = fh.readlines()
    start = next((i for i, l in enumerate(lines) if l.startswith("GO\\tNS\\t")), None)
    if start is None:
        return None
    try:
        df = pd.read_csv(io.StringIO("".join(lines[start:])), sep="\\t")
    except Exception:
        return None
    return df if not df.empty else None


def genes_by_segment(path):
    """segment -> sorted list of CDS gene names, from a .intersect file."""
    out = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    try:
        df = pd.read_csv(path, sep="\\t", header=None)
    except Exception:
        return out
    if df.shape[1] < 15:
        return out
    df = df[df[8] == "CDS"]
    for seg, grp in df.groupby(3):
        names = set()
        for attr in grp[14]:
            for field in str(attr).split(";"):
                if field.startswith("Name="):
                    names.add(field.split("=", 1)[1])
        try:
            out[int(seg)] = sorted(names)
        except (TypeError, ValueError):
            pass
    return out


def b64(path):
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


dat = read_dat()
quant = read_quantiles()
thresh = quant.get(float(TOPQ))
figs = []

# ---------------------------------------------------------------- genome track
if dat is not None:
    fig, axs = plt.subplots(2, 1, figsize=(13, 7), sharex=True)
    hi = dat[dat["D2score"] >= thresh] if thresh is not None else dat.iloc[0:0]
    lo = dat.drop(hi.index)

    axs[0].scatter(lo["mid"] / 1e6, lo["D2score"], s=16, color=BASE,
                   alpha=.75, label="below threshold")
    if not hi.empty:
        axs[0].scatter(hi["mid"] / 1e6, hi["D2score"], s=48, color=ACCENT,
                       edgecolor="black", linewidth=.4, zorder=3,
                       label=f"top {(1 - float(TOPQ)) * 100:.0f}% candidates")
    if thresh is not None:
        axs[0].axhline(thresh, color=ACCENT, ls="--", lw=1,
                       label=f"Q{TOPQ} = {thresh:.4g}")
    axs[0].set_ylabel("KAST distance (D2score)")
    axs[0].set_title(f"Compositional anomaly across the genome "
                     f"({WIN:,} bp windows, {INCR:,} bp step)")
    axs[0].legend(loc="upper right", fontsize=8, framealpha=.9)
    axs[0].grid(alpha=.25)

    axs[1].scatter(dat["mid"] / 1e6, dat["Entropy"], s=16, color=GREY, alpha=.8)
    if not hi.empty:
        axs[1].scatter(hi["mid"] / 1e6, hi["Entropy"], s=48, color=ACCENT,
                       edgecolor="black", linewidth=.4, zorder=3)
    axs[1].set_ylabel("min k-mer entropy")
    axs[1].set_xlabel("genome position (Mb)")
    axs[1].grid(alpha=.25)

    fig.tight_layout()
    p = f"{OUT}/01_genome_track.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    figs.append((p, "Genome-wide compositional signal. Red points are windows "
                    "above the top-quantile threshold — the HGT candidates."))

# --------------------------------------------------------- score distribution
if dat is not None:
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.hist(dat["D2score"], bins=40, color=BASE, alpha=.85, edgecolor="white")
    for q in (0.5, 0.9, float(TOPQ)):
        if q in quant:
            ax.axvline(quant[q], ls="--", lw=1.2,
                       color=ACCENT if q == float(TOPQ) else GREY)
            ax.text(quant[q], ax.get_ylim()[1] * .95, f" Q{q}", fontsize=8,
                    color=ACCENT if q == float(TOPQ) else GREY, rotation=90,
                    va="top")
    ax.set_xlabel("KAST distance (D2score)")
    ax.set_ylabel("windows")
    ax.set_title("Distribution of window scores, with quantile cut-offs")
    ax.grid(alpha=.25)
    fig.tight_layout()
    p = f"{OUT}/02_score_distribution.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    figs.append((p, "Where the quantile bins fall. The filtering thresholds are "
                    "percentiles of this distribution, not absolute cut-offs."))

# ------------------------------------------------------------ candidate table
top_intersect = f"quantile_{TOPQ}_inf.intersect"
seg_genes = genes_by_segment(top_intersect)
cand = pd.DataFrame()
if dat is not None and thresh is not None:
    cand = dat[dat["D2score"] >= thresh].copy()
    cand["n_genes"] = cand["Number"].map(lambda s: len(seg_genes.get(s, [])))
    cand["genes"] = cand["Number"].map(
        lambda s: ", ".join(seg_genes.get(s, [])[:12])
        + (" ..." if len(seg_genes.get(s, [])) > 12 else ""))
    cand = cand.sort_values("D2score", ascending=False)
    cand[["Number", "start", "end", "D2score", "Entropy", "n_genes", "genes"]].to_csv(
        f"{OUT}/candidate_regions.tsv", sep="\\t", index=False)

# --------------------------------------------------------------- GO summary
go_rows = []
for path in sorted(glob.glob("*.go.prop")) + sorted(glob.glob("*.go.noprop")):
    df = read_enrichment(path)
    if df is None or "enrichment" not in df.columns:
        continue
    label = os.path.basename(path).split(".intersect")[0]
    kind = "propagated" if path.endswith(".prop") else "non-propagated"
    for _, r in df.iterrows():
        go_rows.append({
            "bin": label, "counts": kind, "GO": r.get("GO"),
            "NS": r.get("NS"), "term": str(r.get("name", "")).strip(),
            "e/p": r.get("enrichment"),
            "study": r.get("ratio_in_study"), "pop": r.get("ratio_in_pop"),
            "fold": r.get("fold_enrichment"),
            "p_fdr_bh": r.get("p_fdr_bh"),
        })
go_df = pd.DataFrame(go_rows)
if not go_df.empty:
    go_df.to_csv(f"{OUT}/go_summary.tsv", sep="\\t", index=False)

# ------------------------------------------------- GO barplot for the top bin
top_bin = f"quantile_{TOPQ}_inf"
top_go = go_df[(go_df["bin"] == top_bin) & (go_df["counts"] == "propagated")] \\
    if not go_df.empty else pd.DataFrame()
if not top_go.empty:
    t = top_go.copy()
    t["fold"] = pd.to_numeric(t["fold"], errors="coerce")
    t = t.dropna(subset=["fold"]).sort_values("fold", ascending=True).tail(12)
    if not t.empty:
        fig, ax = plt.subplots(figsize=(9, max(3, .45 * len(t) + 1.5)))
        colors = [ACCENT if e == "e" else BASE for e in t["e/p"]]
        ax.barh([f"{g}  {n[:44]}" for g, n in zip(t["GO"], t["term"])],
                t["fold"], color=colors)
        ax.set_xlabel("fold enrichment")
        ax.set_title(f"Enriched GO terms in the top {(1 - float(TOPQ)) * 100:.0f}% "
                     f"of windows (propagated counts)")
        ax.grid(alpha=.25, axis="x")
        fig.tight_layout()
        p = f"{OUT}/03_go_top_terms.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        figs.append((p, "Functional signal in the candidate set. Red = enriched, "
                        "blue = purified."))

# ------------------------------------------------------------------- summary
n_windows = 0 if dat is None else len(dat)
n_cand = 0 if cand.empty else len(cand)
n_genes = len({g for s in (cand["Number"] if not cand.empty else [])
               for g in seg_genes.get(s, [])})
n_sig = 0 if go_df.empty else len(go_df[go_df["bin"] == top_bin])

summary = [
    ("Windows scored", f"{n_windows:,}"),
    ("Window size / step", f"{WIN:,} bp / {INCR:,} bp"),
    (f"Candidate windows (>= Q{TOPQ})", f"{n_cand:,}"),
    ("Threshold score", "n/a" if thresh is None else f"{thresh:.6g}"),
    ("Genes in candidate windows", f"{n_genes:,}"),
    ("Significant GO terms (top bin)", f"{n_sig:,}"),
    ("Quantile bins analysed",
     f"{len(set(go_df['bin'])) if not go_df.empty else 0}"),
]
with open(f"{OUT}/summary.txt", "w") as fh:
    for k, v in summary:
        fh.write(f"{k}\\t{v}\\n")

# ---------------------------------------------------------------- HTML report
def table_html(df, maxrows=40):
    if df is None or df.empty:
        return "<p class='none'>No rows.</p>"
    d = df.head(maxrows)
    th = "".join(f"<th>{c}</th>" for c in d.columns)
    tr = "".join("<tr>" + "".join(f"<td>{v}</td>" for v in row) + "</tr>"
                 for row in d.values)
    more = (f"<p class='none'>Showing {maxrows} of {len(df)} rows — "
            f"full table in the TSV.</p>" if len(df) > maxrows else "")
    return f"<table><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table>{more}"


cand_cols = ["Number", "start", "end", "D2score", "Entropy", "n_genes", "genes"]
html = [f"""<!doctype html><html><head><meta charset="utf-8">
<title>HGT pipeline report</title><style>
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
margin:0 auto;max-width:1150px;padding:32px 24px;color:#222;line-height:1.5}}
h1{{margin:0 0 4px;font-size:26px}} h2{{margin:34px 0 10px;font-size:18px;
border-bottom:2px solid #eee;padding-bottom:6px}}
.sub{{color:#777;margin-bottom:24px;font-size:14px}}
.cards{{display:flex;flex-wrap:wrap;gap:12px;margin:18px 0}}
.card{{border:1px solid #e3e3e3;border-radius:8px;padding:12px 16px;min-width:150px}}
.card .v{{font-size:22px;font-weight:600;color:{ACCENT}}}
.card .k{{font-size:12px;color:#666;text-transform:uppercase;letter-spacing:.4px}}
img{{max-width:100%;border:1px solid #eee;border-radius:6px;margin-top:8px}}
.cap{{color:#666;font-size:13px;margin:6px 0 0}}
table{{border-collapse:collapse;width:100%;font-size:12.5px;margin-top:8px}}
th,td{{border:1px solid #e6e6e6;padding:5px 8px;text-align:left}}
th{{background:#f7f7f7}} tr:nth-child(even) td{{background:#fbfbfb}}
.none{{color:#888;font-size:13px;font-style:italic}}
code{{background:#f4f4f4;padding:1px 5px;border-radius:3px}}
</style></head><body>
<h1>HGT detection report</h1>
<div class="sub">Alignment-free compositional scan &middot;
{WIN:,} bp windows, {INCR:,} bp step &middot; candidates = top
{(1 - float(TOPQ)) * 100:.0f}% of windows by KAST distance</div>
<div class="cards">"""]
for k, v in summary:
    html.append(f"<div class='card'><div class='v'>{v}</div>"
                f"<div class='k'>{k}</div></div>")
html.append("</div>")

for path, cap in figs:
    title = {"01": "Where the candidates are",
             "02": "Score distribution",
             "03": "Functional enrichment"}[os.path.basename(path)[:2]]
    html.append(f"<h2>{title}</h2><img src='data:image/png;base64,{b64(path)}'>"
                f"<p class='cap'>{cap}</p>")

html.append("<h2>Candidate regions</h2>")
html.append(table_html(cand[cand_cols] if not cand.empty else pd.DataFrame()))
html.append("<h2>Significant GO terms (top bin)</h2>")
html.append(table_html(go_df[go_df["bin"] == top_bin].drop(columns=["bin"])
                       if not go_df.empty else pd.DataFrame()))
html.append("<p class='cap'>Full tables: <code>candidate_regions.tsv</code>, "
            "<code>go_summary.tsv</code>. Per-bin enrichment output is in "
            "<code>results/dist/</code>.</p></body></html>")

with open(f"{OUT}/report.html", "w") as fh:
    fh.write("\\n".join(html))

print(f"Report written to {OUT}/report.html")
for k, v in summary:
    print(f"  {k}: {v}")

END_SCRIPT

    command -v python3 >/dev/null 2>&1 || { echo "ERROR: 'python3' is not in PATH"; exit 1; }

    python3 make_report.py "!{win_size}" "!{incr}" "!{top_q}" || {
        echo "WARNING: report generation failed; continuing"; mkdir -p report; }
    '''
}

// ---- Workflow -------------------------------------------------------------
workflow {
    def valid_stage3_modes = ["bwa", "coords"]
    if( !valid_stage3_modes.contains(params.stage3_mode) )
        error "Invalid --stage3_mode '${params.stage3_mode}'. Valid options: bwa, coords"

    int win  = params.window_size as int
    int klen = params.kmer_length as int
    int incr = (params.increment    != null) ? (params.increment    as int) : win.intdiv(2)
    int entw = params.entropy_window as int
    int enti = (params.entropy_incr != null) ? (params.entropy_incr as int) : entw.intdiv(2)

    if( entw >= win )
        error "entropy_window (${entw}) must be smaller than window_size (${win})"
    if( incr < 1 || enti < 1 )
        error "increment/entropy_incr must be >= 1"
    if( klen < 1 || klen > 16 )
        error """--kmer_length ${klen} is outside the usable range.
       KAST allocates a dense 4^k array: k=14 needs ~1 GB, k=15 ~4 GB, k=16 ~16 GB,
       and k>=20 is refused outright by KAST. Use a smaller k (KAST's default is 3)."""


    log.info """
    ------------------------------------------------------------------
     HGT pipeline - resolved parameters
    ------------------------------------------------------------------
     genome            : ${params.fasta}
     reference         : ${params.ref_genome}
     annotation        : ${params.gff}
     ontology          : ${params.obo}
     window / step     : ${win} bp / ${incr} bp
     k-mer length      : ${klen}
     entropy win/step  : ${entw} bp / ${enti} bp
     coordinate mode   : ${params.stage3_mode}
     reference contig  : ${params.reference_contig ?: 'AUTO (single-contig)'}
     candidate bin     : top ${((1 - (params.top_quantile as double)) * 100) as int}% (Q${params.top_quantile})
     launch dir        : ${workflow.launchDir}
    ------------------------------------------------------------------
    """.stripIndent()

    fasta_ch   = Channel.fromPath(params.fasta)
    ref_genome = file(params.ref_genome)
    gff_file   = file(params.gff)
    cds_list   = file(params.cds_list)
    all_txt    = file(params.all_txt)
    go_txt     = file(params.go_txt)
    obo_file   = file(params.obo)

    stage1_out = stage1(fasta_ch, win, klen, incr, entw, enti)
    stage2_out = stage2(fasta_ch, stage1_out.dat_file, stage1_out.fwin_file)

    stage3_input = stage2_out.filtered_seqs.flatten()
    def stage3_out
    if( params.stage3_mode == "bwa" )
        stage3_out = stage3(stage3_input, ref_genome, gff_file)
    else
        stage3_out = stage3_coords(stage3_input, ref_genome, gff_file, win, incr, (params.reference_contig ?: 'AUTO'))

    stage4_out = stage4(stage3_out.intersect_files, cds_list)
    stage5_out = stage5(stage4_out.protein_files, all_txt, go_txt, obo_file)
    stage6(stage5_out.noprop.collect(), stage5_out.prop.collect())

    stage7_report(
        stage1_out.dat_file,
        stage2_out.quantile_stats,
        stage3_out.intersect_files.collect(),
        stage5_out.noprop.collect(),
        stage5_out.prop.collect(),
        win, incr, params.top_quantile
    )
}

workflow.onComplete {
    log.info """
    ------------------------------------------------------------------
     ${workflow.success ? 'Completed successfully' : 'FAILED'}
     duration : ${workflow.duration}
     ${workflow.success ? "report   : ${workflow.launchDir}/results/report/report.html" : "error    : ${workflow.errorMessage ?: 'see .nextflow.log'}"}
    ------------------------------------------------------------------
    """.stripIndent()
}
