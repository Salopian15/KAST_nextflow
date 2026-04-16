nextflow.enable.dsl=2

params.fasta = "Ref_genome.fasta"
params.window_size = 5000
params.kmer_length = 25
params.increment = params.window_size / 2
params.ref_genome = "Ref_genome.fasta"
params.gff = "genomic.gff"
params.all_txt = "all.txt"
params.go_txt = "go.txt"
params.cds_list = "cds.list"
params.stage3_mode = "coords"
params.reference_contig = null

process stage1 {
    input:
    path fasta
    val win_size
    val klen
    val incr

    output:
    path "${fasta}.dat", emit: dat_file
    path "win.${fasta}", emit: fwin_file
    path "${fasta}.contam.pdf", emit: pdf_file

    script:
    '''
    #!/bin/bash
    FWIN="win.${fasta}"
    OUT="${fasta}.kast.out"

    cat << 'END_SCRIPT' > fastaWindowed3.py
#!/usr/bin/python
import os, sys, shutil, string
import subprocess
from operator import itemgetter, attrgetter

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

            if line[0] == ">":
                if header != "":
                    processRecord(header, param, sequence, incr)

                header = line
                sequence = ""
            else:
                sequence += line

        processRecord(header, param, sequence, incr)

def doit(argv):
    if len(argv) < 4:
        print("Chops a Fasta into windows e.g. size=2000, incr=1000")
        print("USAGE: \$0 <fasta file> <windowSize> <incr>")
        sys.exit(0)

    afile = argv[1]
    param = int(argv[2])
    incr = int(argv[3])

    main(afile, param, incr)

doit(sys.argv)
END_SCRIPT

    cat << 'END_SCRIPT' > fastaEntropy3.py
#!/usr/bin/python3
import os, sys, shutil, string
import subprocess
import math
from operator import itemgetter, attrgetter
import zlib
# import lzma
# lzObj = lzma.LZMACompressor()

### inner
def defineKmers(klen, nletters):

    # print(klen, nletters)
    if nletters == -4:
        letter = ["A", "C", "G", "T"]

    elif nletters == -40:
        letter = ["A", "C", "G", "U"]
        nletters = 4

    elif nletters == 20:
        letter = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
    elif nletters == 15:
        letter = ['L', 'C', 'A', 'G', 'S', 'T', 'P', 'F', 'W', 'E', 'D', 'N', 'Q', 'K', 'H']
    elif nletters == 10:
        letter = ['L', 'C', 'A', 'G', 'S', 'P', 'F', 'E', 'K', 'H']
    elif nletters == 8:
        letter = ['L', 'A', 'S', 'P', 'F', 'E', 'K', 'H']
    elif nletters == 4:
        letter = ['L', 'A', 'F', 'E']
    elif nletters == 2:
        letter = ['L', 'E']
    else:
        print("ERROR: wrong number of letters!", nletters)
        sys.exit()

    if nletters == -40 or nletters == -4:
        num = 4
    else:
        num = nletters

    ## make k-mers by starting with k-mers of size 1, and grow them to required length
    kmers = list(letter)
    kmer = ""
    current = []
    current.append(kmer)

    for i in range(klen-1):
        new = []
        for j in range(num):
            for kmer in kmers:
                kmer = kmer + letter[j]
                new.append(kmer)
        kmers = list(new)
    return kmers

def addToKmers(kmers, letter):
    i = 0
    new = []
    for kmer in kmers:
        klast = list(kmer)[-1]
        i = letter.index(klast)
        while i < len(letter):
            newk = kmer + letter[i]
            new.append(newk)
            i += 1
    return new

def countKmers(kmers, sequence, klen):
    kmerCount = {}
    for kmer in kmers:
        kmerCount[kmer] = 0
    sequence = sequence.upper()
    seqL = sequence
    max = len(seqL) - klen
    i = 0
    while i <= max:
        j = i + klen
        wordL = seqL[i:j]
        word = "".join(wordL)
        count = kmerCount.get(word, None)
        if count is not None:
            kmerCount[word] = kmerCount[word] + 1
        i = i+1
    return kmerCount



def calcEntropy( kmerCount, expectn, totalk ) :

   shannon1 = shannon2 = 0.0
   totalk = float(totalk)
   expect = float(expectn)/totalk
   for word in kmerCount.keys() :
     count = float(kmerCount[word])
     if count > 0 :
        kfreq = count/totalk
        shannon1 += kfreq * math.log(kfreq,2)
        shannon2 += expect * math.log(kfreq,2)

   shannon1 = shannon1 * -1
   shannon2 = shannon2 * -1

   return shannon1,shannon2


def processRecord(header, param, sequence, incr, mers1, mers2, mers3, mers4, mers5, f_o):

    new = header.replace(">", "")
    window = param
    if len(sequence) >= window:

        i = 0
        p1 = 0
        p2 = p1 + window
        E = []

        while p2 < len(sequence):

            frag = sequence[p1:p2]
            p1 += incr
            p2 = p1 + window
            i += 1

            totH1 = totH2 = 0.0
            kmers = []
            for k in [1, 2, 3, 4, 5]:
                if k == 1:
                    kmers = mers1
                elif k == 2:
                    kmers = mers2
                elif k == 3:
                    kmers = mers3
                elif k == 4:
                    kmers = mers4
                elif k == 5:
                    kmers = mers5
                else:
                    print("No k-mer", k)
                    continue

                klen = len(kmers[0])
                kmerCount = countKmers(kmers, frag, klen)

                expectn = len(frag) / klen
                totalk = len(frag) - klen
                shan1, shan2 = calcEntropy(kmerCount, expectn, totalk)
                totH1 += shan1
                totH2 += shan2

            E.append([p1, totH1, totH2])

        [p1, totH1, totH2] = sorted(E, key=itemgetter(2))[0]

        ostr = str(p1) + "\t" + str(totH1) + "\t" + str(totH2) + "\t" + new
        f_o.write(ostr)
        f_o.write("\n")
        print(ostr)

        return


def main(afile, param, incr, f_o):

    nletters = -4
    mers1 = defineKmers(1, nletters)
    mers2 = defineKmers(2, nletters)
    mers3 = defineKmers(3, nletters)
    mers4 = defineKmers(4, nletters)
    mers5 = defineKmers(5, nletters)

    f_in = open(afile)

    ostr = "#pos\t" + "H1\t" + "H2\t" + "SeqID"
    f_o.write(ostr)
    f_o.write("\n")

    sequence = ""
    header = ""
    for line in f_in:

        line = line.strip()

        if line[0] == ">":

            if header != "":
                processRecord(header, param, sequence, incr, mers1, mers2, mers3, mers4, mers5, f_o)

            header = line
            sequence = ""

        else:
            sequence = sequence + line

    processRecord(header, param, sequence, incr, mers1, mers2, mers3, mers4, mers5, f_o)

    f_in.close()


def doit(argv):

    if len(argv) < 2:
        print("Chops a Fasta into windows e.g. size=2000, incr=1000..")
        print("Outputs lowest entropy window only for each sequence.")
        print("USAGE: \$0 <fasta file> <max windowSize> <incr>")  # <kmer>")
        exit(0)

    afile = argv[1]
    maxWin = int(argv[2])
    incr = int(argv[3])

    nletters = -4

    ofile = afile + ".entropy"
    f_o = open(ofile, "w")
    main(afile, maxWin, incr, f_o)
    f_o.close()


doit(sys.argv)
END_SCRIPT

    cat << 'END_SCRIPT' > HGT.R
#library("scatterplot3d")

args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]
myID <- args[3]
#covLim <- args[4]


w1 <- read.csv(file=infile,sep=",",head=TRUE)
names(w1)

pdf(outfile)

#Get mean coverage of longest 10 sequences (assuming data ordered)
#covs <- w1\$Coverage[0:10]
#mcov <- mean(covs)
#covRange1 <- mcov * 2
#covRange2 <- mcov * 20


#mcov



#rbPal <- colorRampPalette(c('red','yellow','green','blue','purple'),alpha=TRUE)
rbPal <- colorRampPalette(c('green','red'),alpha=TRUE)
w1\$ColD <- rbPal(100)[as.numeric(cut(w1\$D2score,breaks = 100))]
w1\$ColE <- rbPal(100)[as.numeric(cut(w1\$Entropy,breaks = 100))]
w1\$ColN <- rbPal(100)[as.numeric(cut(w1\$Number,breaks = 100))]



#scatterplot3d(w1)
par(mfrow=c(2,1))

title1<-paste(myID,"   K-mer stat")
title2<-paste(myID,"   Entropy")

plot(w1\$Number,w1\$D2score,main=title1,xlab="Position", ylab="Stat1", pch=20, col=w1\$ColE )   #, xlim=c(0,covRange2) )
plot(w1\$Number,w1\$Entropy,main=title2,xlab="Position", ylab="Stat2", pch=20, col=w1\$ColD ) #, xlim=c(0,covRange2) )



END_SCRIPT

    command -v kast >/dev/null 2>&1 || { echo "ERROR: 'kast' is not in PATH"; exit 1; }

    python fastaWindowed3.py "${fasta}" "${win_size}" "${incr}" > "$FWIN"

    T="manhattan"
    kast -q "$FWIN" -r "${fasta}" -t "$T" -k "${klen}" -o "$OUT.1" -n 0 -f blastlike -nh

    python fastaEntropy3.py "$FWIN" "${win_size}" "${incr}" > "$OUT.3"

    echo "D2score,Number,Entropy" > "${fasta}.dat"
    cat "$OUT.1" | cut -f 2 | cut -d '_' -f 1  > tmp.2
    cat "$OUT.1" | cut -f 6 > tmp.1
    cat "$OUT.3" | grep -v '#' | cut -f 3 > tmp.3

    c1=$(wc -l < tmp.1); c3=$(wc -l < tmp.3)
    if [ "$c1" != "$c3" ]; then
        echo "ERROR: KAST output ($c1 rows) and entropy output ($c3 rows) have different row counts — possible ordering mismatch"
        exit 1
    fi

    paste -d ',' tmp.1 tmp.2 tmp.3 >> "${fasta}.dat"

    Rscript HGT.R "${fasta}.dat" "${fasta}.contam.pdf" "${fasta}"
    '''
}

process stage2 {
    input:
    path fasta
    path dat_file
    path fwin_file

    output:
    path "Filtered_sequences/*.fasta", emit: filtered_seqs

    script:
    '''
    #!/bin/bash
    cat << 'END_SCRIPT' > kastdatWindowedfasta.py
from Bio import SeqIO
import pandas as pd
import sys
import csv

score_file = sys.argv[1]  # The score file name
dna_file = sys.argv[2]  # The DNA file name
output_file = sys.argv[3]  # The output file name

import os
if not os.path.exists(score_file):
    print(f"ERROR: Score file not found: {score_file}")
    sys.exit(1)
if not os.path.exists(dna_file):
    print(f"ERROR: DNA file not found: {dna_file}")
    sys.exit(1)

sequences = {
    "Segment": [],
    "D2 Score": [],
    "Sequence": []
}

# Read the DNA sequences from the file and store them in a list
dna_list = []
with open(dna_file, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        dna_list.append(str(record.seq))

# Read the scores and corresponding segments from the score file
with open(score_file, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if 'Number' in row and 'D2score' in row:
            segment = int(row['Number'])
            d2score = float(row['D2score'])
            # Append the segment and score to the dictionary
            sequences['Segment'].append(segment)
            sequences['D2 Score'].append(d2score)
            sequences['Sequence'].append(dna_list[segment - 1])



# Create the dataframe from the dictionary
df = pd.DataFrame(sequences)
print(df.head)



# Save the dataframe to the output file
df.to_csv(output_file, index=False)
END_SCRIPT

    cat << 'END_SCRIPT' > D2kaststats.py
import pandas as pd
import numpy as np
import sys

csvfile = sys.argv[1]
output_file = "quantile_statistics.txt"

# Load the CSV file
df = pd.read_csv(csvfile)

# Ensure the 'D2score' column exists
if 'D2score' not in df.columns:
    print("The 'D2score' column does not exist in the CSV file.")
else:
    # Calculate the statistics
    quantiles = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
    quantile_stats = {q: df['D2score'].quantile(q) for q in quantiles}

    # Write the statistics to the output file
    with open(output_file, 'w') as f:
        for q, value in quantile_stats.items():
            f.write(f"Quantile {q}: {value}\n")

    print("Quantile statistics have been written to the file:", output_file)
END_SCRIPT

    cat << 'END_SCRIPT' > D2_filterCumulative.py
import sys 
import pandas as pd 
import os 
if len(sys.argv) < 2: 
    print("Usage: python D2_filter.py input.csv output_folder") 
    sys.exit() 
csv_input = sys.argv[1] 
output_folder = sys.argv[2] 
df = pd.read_csv(csv_input) 
quantiles_file = "quantile_statistics.txt" 
quantiles = {} 
if not os.path.exists(output_folder):  # Check if the output folder exists, if not, create it 
    os.makedirs(output_folder) 
with open(quantiles_file, 'r') as file: 
    for line in file: 
        parts = line.split(":") 
        quantile = float(parts[0].split()[1]) 
        value = float(parts[1].strip()) 
        quantiles[quantile] = value 
print(quantiles) 
for i, quantile in enumerate(quantiles): 
    lower = 0 if i == 0 else list(quantiles.values())[i-1] 
    upper = quantiles[quantile]
    # Modify the filtering to include the upper bound 
    score_select = df[(df["D2 Score"] >= lower) & (df["D2 Score"] <= upper)] 
    output_file = f"{output_folder}/quantile_{lower}_{quantile}.fasta"
    with open(output_file, 'w') as fasta_output:
        for index, row in score_select.iterrows():
            sequence_id = f">{row['Segment']}"
            sequence = row['Sequence']
            #print(f"Number of rows for range {lower}-{upper}: {len(score_select)}")
            fasta_output.write(f"{sequence_id}\n{sequence}\n")
    #print(f"Output file created: {output_file}")

# Handle scores that are greater than the last quantile
lower = list(quantiles.values())[-1]
upper = float('inf')
score_select = df[df["D2 Score"] >= lower]
output_file = f"{output_folder}/quantile_0.95_inf.fasta"
with open(output_file, 'w') as fasta_output:
    for index, row in score_select.iterrows():
        sequence_id = f">{row['Segment']}"
        sequence = row['Sequence']
        fasta_output.write(f"{sequence_id}\n{sequence}\n")
#print(f"Output file created: {output_file}")

# Cumulative sections

for quantile in quantiles:
    lower = 0
    upper = quantiles[quantile]
    score_select = df[(df["D2 Score"] >= lower) & (df["D2 Score"] <= upper)] 
    output_file = f"{output_folder}/Cumulative_0_{quantile}.fasta"
    with open(output_file, 'w') as fasta_output:
        for index, row in score_select.iterrows():
            sequence_id = f">{row['Segment']}"
            sequence = row['Sequence']
            fasta_output.write(f"{sequence_id}\n{sequence}\n")
    #print(f"Output file created: {output_file}")

for quantile in quantiles:
    lower = quantiles[quantile]
    upper = float('inf')
    score_select = df[(df["D2 Score"] >= lower) & (df["D2 Score"] <= upper)] 
    output_file = f"{output_folder}/Cumulative_{quantile}_inf.fasta"
    with open(output_file, 'w') as fasta_output:
        for index, row in score_select.iterrows():
            sequence_id = f">{row['Segment']}"
            sequence = row['Sequence']
            fasta_output.write(f"{sequence_id}\n{sequence}\n")
    #print(f"Output file created: {output_file}")
END_SCRIPT

    python kastdatWindowedfasta.py "${dat_file}" "${fwin_file}" "${fasta}.csv"
    python D2kaststats.py "${dat_file}"
    
    mkdir -p Filtered_sequences
    python D2_filterCumulative.py "${fasta}.csv" Filtered_sequences
    '''
}

process stage3 {
    input:
    path fasta_file
    path ref_genome
    path gff

    output:
    path "${fasta_file.baseName}.intersect", emit: intersect_files

    script:
    '''
    #!/bin/bash
    bwa index "${ref_genome}"
    bwa mem "${ref_genome}" "${fasta_file}" > "${fasta_file.baseName}.sam"
    samtools view -S -b "${fasta_file.baseName}.sam" > "${fasta_file.baseName}.bam"
    samtools sort "${fasta_file.baseName}.bam" -o "${fasta_file.baseName}.sorted.bam"
    bedtools bamtobed -i "${fasta_file.baseName}.sorted.bam" > "${fasta_file.baseName}.sorted.bed"
    bedtools intersect -a "${fasta_file.baseName}.sorted.bed" -b "${gff}" -wa -wb > "${fasta_file.baseName}.intersect"
    '''
}

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

    script:
    '''
    #!/bin/bash
    cat << 'END_SCRIPT' > make_bed_from_segments.py
import sys

if len(sys.argv) != 6:
    print("Usage: python make_bed_from_segments.py <segments_fasta> <ref_fasta> <window_size> <increment> <reference_contig_or_AUTO>")
    sys.exit(1)

segments_fasta = sys.argv[1]
ref_fasta = sys.argv[2]
window_size = int(sys.argv[3])
increment = int(sys.argv[4])
reference_contig = sys.argv[5]

def read_contigs(ref_path):
    contigs = []
    with open(ref_path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                contigs.append(line[1:].strip().split()[0])
    return contigs

contigs = read_contigs(ref_fasta)
if not contigs:
    print("ERROR: No FASTA headers found in reference genome")
    sys.exit(1)

if reference_contig == "AUTO":
    if len(contigs) != 1:
        print("ERROR: Multiple contigs found in reference and params.reference_contig is not set")
        print("Contigs:", ",".join(contigs))
        sys.exit(1)
    contig = contigs[0]
else:
    if reference_contig not in contigs:
        print(f"ERROR: reference_contig '{reference_contig}' was not found in reference genome")
        print("Contigs:", ",".join(contigs))
        sys.exit(1)
    contig = reference_contig

bed_out = f"{segments_fasta.rsplit('.', 1)[0]}.sorted.bed"
with open(segments_fasta, "r") as in_fh, open(bed_out, "w") as out_fh:
    for line in in_fh:
        if not line.startswith(">"):
            continue
        header = line[1:].strip().split()[0]
        try:
            seg = int(header)
        except ValueError:
            print(f"ERROR: Segment header '{header}' is not an integer. Expected format '>123'.")
            sys.exit(1)

        start = (seg - 1) * increment
        end = start + window_size
        if start < 0:
            print(f"ERROR: Computed negative start for segment {seg}")
            sys.exit(1)
        out_fh.write(f"{contig}\t{start}\t{end}\t{header}\n")
END_SCRIPT

    REF_CONTIG_INPUT="AUTO"
    if [ -n "${reference_contig}" ] && [ "${reference_contig}" != "null" ]; then
      REF_CONTIG_INPUT="${reference_contig}"
    fi

    python make_bed_from_segments.py "${fasta_file}" "${ref_genome}" "${win_size}" "${incr}" "$REF_CONTIG_INPUT" || { echo "ERROR: make_bed_from_segments.py failed — check segment headers and reference_contig parameter"; exit 1; }
    bedtools intersect -a "${fasta_file.baseName}.sorted.bed" -b "${gff}" -wa -wb > "${fasta_file.baseName}.intersect"
    '''
}

process stage4 {
    input:
    path intersect_file
    path cds_list

    output:
    path "${intersect_file}.gene", emit: gene_files
    path "${intersect_file}.protein", emit: protein_files

    script:
    '''
    #!/bin/bash
    cat << 'END_SCRIPT' > ProteinIntersect.py
import pandas as pd
import sys

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 3:
    print("Usage: python script.py input_file output_file")
    sys.exit(1)

# Extract command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the input file into a DataFrame
df = pd.read_csv(input_file, sep='\t', header=None)

# Filter rows where the 9th column is "CDS"
cds_df = df[df[8] == "CDS"]

# Extract unique names following "Name="
unique_names = set()
for annotation in cds_df[14]:
    annotations = annotation.split(';')
    for annotation in annotations:
        if annotation.startswith("Name="):
            name_content = annotation.split('=')[1]
            unique_names.add(name_content)

# Write unique names to the output file
with open(output_file, 'w') as outfile:
    for name_content in unique_names:
        outfile.write(name_content + '\n')

# Print message after extraction is complete
#print("Extraction complete. Unique names written to", output_file)
END_SCRIPT

    cat << 'END_SCRIPT' > ProteinIDconversion.py
import pandas as pd
import sys

# Get the input file paths from command line arguments
infile1 = sys.argv[1]
infile2 = sys.argv[2]
output_file = sys.argv[3]

# Read the first input file into a DataFrame
df1 = pd.read_csv(infile1, sep='\t', header=None)
#print(df1.head)
# Read the second input file into a DataFrame
df2 = pd.read_csv(infile2, sep='\t', header=None)
#print(df2.head)
# Create a dictionary mapping values from the first column of df1 to the second column of df1
mapping = dict(zip(df1.iloc[:, 0], df1.iloc[:, 1]))

# Open the output file to write the corresponding values
with open(output_file, 'w') as f:
    # Iterate over the values in the first column of df2
    for value in df2.iloc[:, 0]:
        # Check if the value is in the mapping
        if value in mapping:
            # Write the corresponding value from the second column of df1 to the output file
            f.write(str(mapping[value]) + '\n')
END_SCRIPT

    python ProteinIntersect.py "${intersect_file}" "${intersect_file}.gene" || true
    python ProteinIDconversion.py "${cds_list}" "${intersect_file}.gene" "${intersect_file}.protein" || true
    '''
}

process stage5 {
    publishDir "results/dist", mode: 'copy'

    input:
    path protein_file
    path all_txt
    path go_txt

    output:
    path "*.go.noprop", emit: noprop
    path "*.go.prop", emit: prop

    script:
    '''
    #!/bin/bash
    cat << 'END_SCRIPT' > find_enrichment1.py
#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
python find_enrichment.py study.file population.file gene-association.file

This program returns P-values for functional enrichment in a cluster of study
genes using Fisher's exact test, and corrected for multiple testing (including
Bonferroni, Holm, Sidak, and false discovery rate).

About significance cutoff:
--alpha: test-wise alpha; for each GO term, what significance level to apply
        (most often you don't need to change this other than 0.05 or 0.01)
--pval: experiment-wise alpha; for the entire experiment, what significance
        level to apply after Bonferroni correction
"""

__copyright__ = "Copyright (C) 2010-present, H Tang et al. All rights reserved."
__author__ = "various"

import sys
import os.path as op
from goatools.cli.find_enrichment import GoeaCliArgs
from goatools.cli.find_enrichment import GoeaCliFnc

sys.path.insert(0, op.join(op.dirname(__file__), ".."))


def main():
    """Run gene enrichment analysis."""
    # Load study, population, associations, and GoDag. Run GOEA.
    obj = GoeaCliFnc(GoeaCliArgs().args)
    # Reduce results to significant results (pval<value)
    results_specified = obj.get_results()
    # Print results in a flat list
    obj.prt_results(results_specified)
    # if obj.sections and obj.args.outfile_detail:
    #     #fout_detail = obj.args.outfile_detail if obj.args.outfile_detail else "goea_details.txt"
    #     objaart = obj.get_objaart()
    #     objaart.run("GOEA", results, sys.stdout)
    #### prt_grouped(results, objgoea, args)


if __name__ == "__main__":
    main()

# Copyright (C) 2010-present, H Tang et al. All rights reserved.

END_SCRIPT

    cat << 'END_SCRIPT' > PlotGOHistograms.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import io

# Get all .noprop and .prop files — exclude Cumulative files
noprop_files = [f for f in glob.glob('*.noprop') if not f.startswith('Cumulative')]
prop_files = [f for f in glob.glob('*.prop') if not f.startswith('Cumulative')]

propdfs = {}
nopropdfs = {}

# Read .noprop files
for file in noprop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    nopropdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

# Read .prop files
for file in prop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    propdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

# Create two subplots: one for 'e' and one for 'p'
#fig, axs = plt.subplots(2, figsize=(16,10))
fig, axs = plt.subplots(2)
# Create a DataFrame for 'e' and 'p' counts
data_e = []
data_p = []

# Iterate over all dataframes
for i, (name, df) in enumerate(nopropdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    quantile_start = float(name.split('_')[1])  # Extract the first quantile value from the filename
    data_e.append({'Quantile Start': quantile_start, 'Count': e_count, 'File': name})
    data_p.append({'Quantile Start': quantile_start, 'Count': p_count, 'File': name})

df_e = pd.DataFrame(data_e)
df_p = pd.DataFrame(data_p)

# Sort the dataframes by 'Quantile Start'
df_e = df_e.sort_values('Quantile Start')
df_p = df_p.sort_values('Quantile Start')

# 'e' histogram
sns.barplot(x='Quantile Start', y='Count', hue='File', data=df_e, ax=axs[0])

# 'p' histogram
sns.barplot(x='Quantile Start', y='Count', hue='File', data=df_p, ax=axs[1])

# Set labels and titles
axs[0].set_xlabel('Quantile Start')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts')
#axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=45, fontsize=8)  # Add rotation here
axs[0].legend().remove()

axs[1].set_xlabel('Quantile Start')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts')
#axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=45, fontsize=8)  # And here
axs[1].legend().remove()
fig.tight_layout()

# Save the figure
plt.savefig('histogram_noprop.png')

#plt.show()

# Create two subplots: one for 'e' and one for 'p'
fig, axs = plt.subplots(2)

# Create a DataFrame for 'e' and 'p' counts
data_e_prop = []
data_p_prop = []

# Iterate over all dataframes
for i, (name, df) in enumerate(propdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    quantile_start = float(name.split('_')[1])  # Extract the first quantile value from the filename
    data_e_prop.append({'Quantile Start': quantile_start, 'Count': e_count, 'File': name})
    data_p_prop.append({'Quantile Start': quantile_start, 'Count': p_count, 'File': name})

df_e_prop = pd.DataFrame(data_e_prop)
df_p_prop = pd.DataFrame(data_p_prop)

# Sort the dataframes by 'Quantile Start'
df_e_prop = df_e_prop.sort_values('Quantile Start')
df_p_prop = df_p_prop.sort_values('Quantile Start')

# 'e' histogram
sns.barplot(x='Quantile Start', y='Count', hue='File', data=df_e_prop, ax=axs[0], palette='cool')

# 'p' histogram
sns.barplot(x='Quantile Start', y='Count', hue='File', data=df_p_prop, ax=axs[1], palette='cool')

# Set labels and titles
axs[0].set_xlabel('Quantile Start')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts (prop)')
#axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=45, fontsize=8)  # Add rotation here
axs[0].legend().remove()

axs[1].set_xlabel('Quantile Start')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts (prop)')
#axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=45, fontsize=8)  # And here
axs[1].legend().remove()
fig.tight_layout()

# Save the figure
plt.savefig('histograms_prop.png')

#plt.show()

END_SCRIPT

    cat << 'END_SCRIPT' > GOstatgraphs.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import io

# Get all .noprop and .prop files — exclude Cumulative files
noprop_files = [f for f in glob.glob('*.noprop') if not f.startswith('Cumulative')]
prop_files = [f for f in glob.glob('*.prop') if not f.startswith('Cumulative')]
propdfs = {}
nopropdfs = {}

# Read .noprop files
for file in noprop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    nopropdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

# Read .prop files
for file in prop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    propdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

# Create a list for storing data
data = []

# Function to process dataframes
def process_dfs(dfs, df_type):
    for name, df in dfs.items():
        for enrichment in ['e', 'p']:
            # Filter dataframe by enrichment type
            df_filtered = df[df['enrichment'] == enrichment]
            # Count the occurrences of 'CC', 'MF', and 'BP' in the 'NS' column
            counts = df_filtered['NS'].value_counts(normalize=True)
            for ns in ['CC', 'MF', 'BP']:
                if ns in counts:
                    # Extract the first part of the filename as quantile_start
                    quantile_start = float(name.split('_')[1])
                    data.append({'Quantile Start': quantile_start, 'NS': ns, 'Relative Abundance %': counts[ns], 'Type': df_type, 'Enrichment': enrichment})

# Process nopropdfs and propdfs
process_dfs(nopropdfs, 'noprop')
process_dfs(propdfs, 'prop')

# Convert the list to a DataFrame
df = pd.DataFrame(data)

# Create separate plots for 'noprop' and 'prop', and for 'e' and 'p'
for df_type in ['noprop', 'prop']:
    for enrichment in ['e', 'p']:
        plt.figure(figsize=(10, 6))
        sns.barplot(x='Quantile Start', y='Relative Abundance %', hue='NS', data=df[(df['Type'] == df_type) & (df['Enrichment'] == enrichment)], palette='Set2', ci=None)
        plt.title(f'Relative Abundances of CC, MF, and BP ({df_type}, {enrichment})')
        plt.savefig(f'Relative Abundances of CC, MF, and BP ({df_type}, {enrichment}).png')
        #plt.show()
END_SCRIPT

    cat << 'END_SCRIPT' > PlotGOHistogramsCumulative.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import io

ctnpdfs = {}
ctpdfs = {}
cbnpdfs = {}
cbpdfs = {}


cumulative_top_noprop_files = glob.glob('Cumulative_0_*.noprop')
cumulative_top_prop_files = glob.glob('Cumulative_0_*.prop')
cumulative_bottom_noprop_files = glob.glob('*_inf.intersect.protein.go.noprop')
cumulative_bottom_prop_files = glob.glob('*_inf.intersect.protein.go.prop')



for file in cumulative_top_noprop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    ctnpdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

for file in cumulative_top_prop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    ctpdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

for file in cumulative_bottom_noprop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    cbnpdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

for file in cumulative_bottom_prop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    cbpdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')



# Create two subplots for cumulative top noprop files
fig, axs = plt.subplots(2, figsize=(16, 10))

# Create a DataFrame for 'e' and 'p' counts for cumulative top noprop files
data_e_ctn = []
data_p_ctn = []

for i, (name, df) in enumerate(ctnpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    #quantile_start = float(name.split('_')[2])
    quantile_start = float(name.split('_')[-1].split('.')[0] + '.' + name.split('_')[-1].split('.')[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_e_ctn.append({'Cumulative Start': quantile_start_percentage, 'Count': e_count, 'File': name})
    data_p_ctn.append({'Cumulative Start': quantile_start_percentage, 'Count': p_count, 'File': name})

df_e_ctn = pd.DataFrame(data_e_ctn)
df_p_ctn = pd.DataFrame(data_p_ctn)

df_e_ctn = df_e_ctn.sort_values('Cumulative Start')
df_p_ctn = df_p_ctn.sort_values('Cumulative Start')

sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_e_ctn, ax=axs[0])
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_p_ctn, ax=axs[1])

axs[0].set_xlabel('Cumulative fraction of genome assessed % -- Counting from bottom')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts (Cumulative count to Top - No Prop)')
axs[0].legend().remove()

axs[1].set_xlabel('Cumulative fraction of genome assessed % -- Counting from bottom')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts (Cumulative count to Top - No Prop)')
axs[1].legend().remove()

fig.tight_layout()
plt.savefig('histogram_cumulative_top_noprop.png')

# Create two subplots for cumulative top prop files
fig, axs = plt.subplots(2, figsize=(16, 10))

# Create a DataFrame for 'e' and 'p' counts for cumulative top prop files
data_e_ctp = []
data_p_ctp = []

for i, (name, df) in enumerate(ctpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    # quantile_start = float(name.split('_')[2]) 
    quantile_start = float(name.split('_')[-1].split('.')[0] + '.' + name.split('_')[-1].split('.')[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_e_ctp.append({'Cumulative Start': quantile_start_percentage, 'Count': e_count, 'File': name})
    data_p_ctp.append({'Cumulative Start': quantile_start_percentage, 'Count': p_count, 'File': name})


df_e_ctp = pd.DataFrame(data_e_ctp)
df_p_ctp = pd.DataFrame(data_p_ctp)

df_e_ctp = df_e_ctp.sort_values('Cumulative Start')
df_p_ctp = df_p_ctp.sort_values('Cumulative Start')

sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_e_ctp, ax=axs[0], palette='cool')
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_p_ctp, ax=axs[1], palette='cool')

axs[0].set_xlabel('Cumulative fraction of genome assessed % -- Counting from bottom')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts (Cumulative count to Top - Prop)')
axs[0].legend().remove()

axs[1].set_xlabel('Cumulative fraction of genome assessed % -- Counting from bottom')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts (Cumulative count to Top - Prop)')
axs[1].legend().remove()

fig.tight_layout()
plt.savefig('histogram_cumulative_top_prop.png')

#plt.show()
# Create two subplots for cumulative bottom noprop files
fig, axs = plt.subplots(2, figsize=(16, 10))

# Create a DataFrame for 'e' and 'p' counts for cumulative bottom noprop files
data_e_cbn = []
data_p_cbn = []

for i, (name, df) in enumerate(cbnpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    quantile_start = float(name.split('_')[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_e_cbn.append({'Cumulative Start': quantile_start_percentage, 'Count': e_count, 'File': name})
    data_p_cbn.append({'Cumulative Start': quantile_start_percentage, 'Count': p_count, 'File': name})

df_e_cbn = pd.DataFrame(data_e_cbn)
df_p_cbn = pd.DataFrame(data_p_cbn)

df_e_cbn = df_e_cbn.sort_values('Cumulative Start')
df_p_cbn = df_p_cbn.sort_values('Cumulative Start')

sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_e_cbn, ax=axs[0])
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_p_cbn, ax=axs[1])

axs[0].set_xlabel('Cumulative fraction of genome assessed % -- Counting from top')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts (Cumulative count to Bottom - No Prop)')
axs[0].legend().remove()

axs[1].set_xlabel('Cumulative fraction of genome assessed % -- Counting from top')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts (Cumulative count to Bottom - No Prop)')
axs[1].legend().remove()

fig.tight_layout()
plt.savefig('histogram_cumulative_bottom_noprop.png')

# Create two subplots for cumulative bottom prop files
fig, axs = plt.subplots(2, figsize=(16, 10))

# Create a DataFrame for 'e' and 'p' counts for cumulative bottom prop files
data_e_cbp = []
data_p_cbp = []

for i, (name, df) in enumerate(cbpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    quantile_start = float(name.split('_', 2)[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_e_cbp.append({'Cumulative Start': quantile_start_percentage, 'Count': e_count, 'File': name})
    data_p_cbp.append({'Cumulative Start': quantile_start_percentage, 'Count': p_count, 'File': name})

df_e_cbp = pd.DataFrame(data_e_cbp)
df_p_cbp = pd.DataFrame(data_p_cbp)

df_e_cbp = df_e_cbp.sort_values('Cumulative Start')
df_p_cbp = df_p_cbp.sort_values('Cumulative Start')

sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_e_cbp, ax=axs[0], palette='cool')
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_p_cbp, ax=axs[1], palette='cool')

axs[0].set_xlabel('Cumulative fraction of genome assessed % -- Counting from top')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts (Cumulative count to Bottom - Prop)')
axs[0].legend().remove()

axs[1].set_xlabel('Cumulative fraction of genome assessed % -- Counting from top')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts (Cumulative count to Bottom - Prop)')
axs[1].legend().remove()

fig.tight_layout()
plt.savefig('histogram_cumulative_bottom_prop.png')

#plt.show()

fig, axs = plt.subplots(2, figsize=(16, 10))

# Create a DataFrame for combined 'e' and 'p' counts for cumulative top noprop files
data_combined_ctn = []
for i, (name, df) in enumerate(ctnpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    combined_count = e_count + p_count
    quantile_start = float(name.split('_')[-1].split('.')[0] + '.' + name.split('_')[-1].split('.')[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_combined_ctn.append({'Cumulative Start': quantile_start_percentage, 'Count': combined_count, 'File': name})
df_combined_ctn = pd.DataFrame(data_combined_ctn)
df_combined_ctn = df_combined_ctn.sort_values('Cumulative Start')
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_combined_ctn, ax=axs[0])
axs[0].set_xlabel('Cumulative fraction of genome assessed % -- Counting from bottom')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of combined e and p counts (Cumulative count to Top - No Prop)')
axs[0].legend().remove()

# Create a DataFrame for combined 'e' and 'p' counts for cumulative top prop files
data_combined_ctp = []
for i, (name, df) in enumerate(ctpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    combined_count = e_count + p_count
    quantile_start = float(name.split('_')[-1].split('.')[0] + '.' + name.split('_')[-1].split('.')[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_combined_ctp.append({'Cumulative Start': quantile_start_percentage, 'Count': combined_count, 'File': name})
df_combined_ctp = pd.DataFrame(data_combined_ctp)
df_combined_ctp = df_combined_ctp.sort_values('Cumulative Start')
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_combined_ctp, ax=axs[1], palette='cool')
axs[1].set_xlabel('Cumulative fraction of genome assessed % -- Counting from bottom')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of combined e and p counts (Cumulative count to Top - Prop)')
axs[1].legend().remove()

fig.tight_layout()
plt.savefig('histogram_cumulative_combined_top.png')

#plt.show()

fig, axs = plt.subplots(2, figsize=(16, 10))
# Create a DataFrame for combined 'e' and 'p' counts for cumulative bottom noprop files
data_combined_cbn = []
for i, (name, df) in enumerate(cbnpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    combined_count = e_count + p_count
    quantile_start = float(name.split('_')[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_combined_cbn.append({'Cumulative Start': quantile_start_percentage, 'Count': combined_count, 'File': name})
df_combined_cbn = pd.DataFrame(data_combined_cbn)
df_combined_cbn = df_combined_cbn.sort_values('Cumulative Start')
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_combined_cbn, ax=axs[0])
axs[0].set_xlabel('Cumulative fraction of genome assessed % -- Counting from top')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of combined e and p counts (Cumulative count to Bottom - No Prop)')
axs[0].legend().remove()

# Create a DataFrame for combined 'e' and 'p' counts for cumulative bottom prop files
data_combined_cbp = []
for i, (name, df) in enumerate(cbpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    combined_count = e_count + p_count
    quantile_start = float(name.split('_', 2)[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_combined_cbp.append({'Cumulative Start': quantile_start_percentage, 'Count': combined_count, 'File': name})
df_combined_cbp = pd.DataFrame(data_combined_cbp)
df_combined_cbp = df_combined_cbp.sort_values('Cumulative Start')
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_combined_cbp, ax=axs[1], palette='cool')
axs[1].set_xlabel('Cumulative fraction of genome assessed % -- Counting from top')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of combined e and p counts (Cumulative count to Bottom - Prop)')
axs[1].legend().remove()

fig.tight_layout()
plt.savefig('histogram_cumulative_combined_bottom.png')
END_SCRIPT

    cat << 'END_SCRIPT' > GOstatgraphsCumulative.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import io

# Get all .noprop and .prop files
noprop_files = glob.glob('*.noprop')
prop_files = glob.glob('*.prop')
propdfs = {}
nopropdfs = {}

# Read .noprop files
for file in noprop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    nopropdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

# Read .prop files
for file in prop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    propdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

# Create a list for storing data
data = []

# Function to process dataframes
def process_dfs(dfs, df_type):
    for name, df in dfs.items():
        for enrichment in ['e', 'p']:
            # Filter dataframe by enrichment type
            df_filtered = df[df['enrichment'] == enrichment]
            # Count the occurrences of 'CC', 'MF', and 'BP' in the 'NS' column
            counts = df_filtered['NS'].value_counts(normalize=True)
            for ns in ['CC', 'MF', 'BP']:
                if ns in counts:
                    # Extract the first part of the filename as quantile_start
                    quantile_start = float(name.split('_')[1])
                    data.append({'Quantile Start': quantile_start, 'NS': ns, 'Relative Abundance %': counts[ns], 'Type': df_type, 'Enrichment': enrichment})

# Process nopropdfs and propdfs
process_dfs(nopropdfs, 'noprop')
process_dfs(propdfs, 'prop')

# Convert the list to a DataFrame
df = pd.DataFrame(data)

# Create separate plots for 'noprop' and 'prop', and for 'e' and 'p'
for df_type in ['noprop', 'prop']:
    for enrichment in ['e', 'p']:
        plt.figure(figsize=(10, 6))
        sns.barplot(x='Quantile Start', y='Relative Abundance %', hue='NS', data=df[(df['Type'] == df_type) & (df['Enrichment'] == enrichment)], palette='Set2', errorbar=None)
        plt.title(f'Relative Abundances of CC, MF, and BP ({df_type}, {enrichment})')
        plt.savefig(f'Relative Abundances of CC, MF, and BP ({df_type}, {enrichment}).png')
        #plt.show()
END_SCRIPT

    if [ -s "${protein_file}" ]; then
        python find_enrichment1.py "${protein_file}" "${all_txt}" "${go_txt}" --pval 0.05 --no_propagate_counts > "${protein_file}.go.noprop"
        python find_enrichment1.py "${protein_file}" "${all_txt}" "${go_txt}" --pval 0.05 > "${protein_file}.go.prop"
    else
        echo "# No proteins found — empty enrichment" > "${protein_file}.go.noprop"
        echo "# No proteins found — empty enrichment" > "${protein_file}.go.prop"
    fi
    '''
}

process stage6 {
    publishDir "results/dist", mode: 'copy'

    input:
    path noprop_files
    path prop_files

    output:
    path "GO_results/*", emit: go_plots, optional: true
    path "CumulativeGO_results/*", emit: cumulative_plots, optional: true

    script:
    '''
    #!/bin/bash
    cat << 'END_SCRIPT' > PlotGOHistograms.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import io

# Get all .noprop and .prop files — exclude Cumulative files
noprop_files = [f for f in glob.glob('*.noprop') if not f.startswith('Cumulative')]
prop_files = [f for f in glob.glob('*.prop') if not f.startswith('Cumulative')]

propdfs = {}
nopropdfs = {}

# Read .noprop files
for file in noprop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    nopropdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

# Read .prop files
for file in prop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    propdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

if not nopropdfs and not propdfs:
    print("No enrichment result files found — skipping histogram plots")
    exit(0)

fig, axs = plt.subplots(2)
data_e = []
data_p = []

for i, (name, df) in enumerate(nopropdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    quantile_start = float(name.split('_')[1])
    data_e.append({'Quantile Start': quantile_start, 'Count': e_count, 'File': name})
    data_p.append({'Quantile Start': quantile_start, 'Count': p_count, 'File': name})

df_e = pd.DataFrame(data_e)
df_p = pd.DataFrame(data_p)
df_e = df_e.sort_values('Quantile Start')
df_p = df_p.sort_values('Quantile Start')

sns.barplot(x='Quantile Start', y='Count', hue='File', data=df_e, ax=axs[0])
sns.barplot(x='Quantile Start', y='Count', hue='File', data=df_p, ax=axs[1])

axs[0].set_xlabel('Quantile Start')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts')
axs[0].legend().remove()

axs[1].set_xlabel('Quantile Start')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts')
axs[1].legend().remove()
fig.tight_layout()
plt.savefig('histogram_noprop.png')

fig, axs = plt.subplots(2)
data_e_prop = []
data_p_prop = []

for i, (name, df) in enumerate(propdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    quantile_start = float(name.split('_')[1])
    data_e_prop.append({'Quantile Start': quantile_start, 'Count': e_count, 'File': name})
    data_p_prop.append({'Quantile Start': quantile_start, 'Count': p_count, 'File': name})

df_e_prop = pd.DataFrame(data_e_prop)
df_p_prop = pd.DataFrame(data_p_prop)
df_e_prop = df_e_prop.sort_values('Quantile Start')
df_p_prop = df_p_prop.sort_values('Quantile Start')

sns.barplot(x='Quantile Start', y='Count', hue='File', data=df_e_prop, ax=axs[0], palette='cool')
sns.barplot(x='Quantile Start', y='Count', hue='File', data=df_p_prop, ax=axs[1], palette='cool')

axs[0].set_xlabel('Quantile Start')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts (prop)')
axs[0].legend().remove()

axs[1].set_xlabel('Quantile Start')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts (prop)')
axs[1].legend().remove()
fig.tight_layout()
plt.savefig('histograms_prop.png')
END_SCRIPT

    cat << 'END_SCRIPT' > GOstatgraphs.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import io

# Get all .noprop and .prop files — exclude Cumulative files
noprop_files = [f for f in glob.glob('*.noprop') if not f.startswith('Cumulative')]
prop_files = [f for f in glob.glob('*.prop') if not f.startswith('Cumulative')]
propdfs = {}
nopropdfs = {}

for file in noprop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    nopropdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

for file in prop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    propdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

if not nopropdfs and not propdfs:
    print("No enrichment result files found — skipping NS abundance plots")
    exit(0)

data = []

def process_dfs(dfs, df_type):
    for name, df in dfs.items():
        for enrichment in ['e', 'p']:
            df_filtered = df[df['enrichment'] == enrichment]
            counts = df_filtered['NS'].value_counts(normalize=True)
            for ns in ['CC', 'MF', 'BP']:
                if ns in counts:
                    quantile_start = float(name.split('_')[1])
                    data.append({'Quantile Start': quantile_start, 'NS': ns, 'Relative Abundance %': counts[ns], 'Type': df_type, 'Enrichment': enrichment})

process_dfs(nopropdfs, 'noprop')
process_dfs(propdfs, 'prop')

df = pd.DataFrame(data)

for df_type in ['noprop', 'prop']:
    for enrichment in ['e', 'p']:
        plt.figure(figsize=(10, 6))
        sns.barplot(x='Quantile Start', y='Relative Abundance %', hue='NS', data=df[(df['Type'] == df_type) & (df['Enrichment'] == enrichment)], palette='Set2', errorbar=None)
        plt.title(f'Relative Abundances of CC, MF, and BP ({df_type}, {enrichment})')
        plt.savefig(f'Relative Abundances of CC, MF, and BP ({df_type}, {enrichment}).png')
END_SCRIPT

    cat << 'END_SCRIPT' > PlotGOHistogramsCumulative.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import io

ctnpdfs = {}
ctpdfs = {}
cbnpdfs = {}
cbpdfs = {}

cumulative_top_noprop_files = glob.glob('Cumulative_0_*.noprop')
cumulative_top_prop_files = glob.glob('Cumulative_0_*.prop')
cumulative_bottom_noprop_files = glob.glob('*_inf.intersect.protein.go.noprop')
cumulative_bottom_prop_files = glob.glob('*_inf.intersect.protein.go.prop')

for file in cumulative_top_noprop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    ctnpdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

for file in cumulative_top_prop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    ctpdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

for file in cumulative_bottom_noprop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    cbnpdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

for file in cumulative_bottom_prop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    cbpdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

if not ctnpdfs and not cbnpdfs:
    print("No cumulative enrichment files found — skipping cumulative histogram plots")
    exit(0)

fig, axs = plt.subplots(2, figsize=(16, 10))
data_e_ctn = []
data_p_ctn = []

for i, (name, df) in enumerate(ctnpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    quantile_start = float(name.split('_')[-1].split('.')[0] + '.' + name.split('_')[-1].split('.')[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_e_ctn.append({'Cumulative Start': quantile_start_percentage, 'Count': e_count, 'File': name})
    data_p_ctn.append({'Cumulative Start': quantile_start_percentage, 'Count': p_count, 'File': name})

df_e_ctn = pd.DataFrame(data_e_ctn)
df_p_ctn = pd.DataFrame(data_p_ctn)
df_e_ctn = df_e_ctn.sort_values('Cumulative Start')
df_p_ctn = df_p_ctn.sort_values('Cumulative Start')

sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_e_ctn, ax=axs[0])
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_p_ctn, ax=axs[1])

axs[0].set_xlabel('Cumulative fraction of genome assessed % -- Counting from bottom')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts (Cumulative count to Top - No Prop)')
axs[0].legend().remove()

axs[1].set_xlabel('Cumulative fraction of genome assessed % -- Counting from bottom')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts (Cumulative count to Top - No Prop)')
axs[1].legend().remove()

fig.tight_layout()
plt.savefig('histogram_cumulative_top_noprop.png')

fig, axs = plt.subplots(2, figsize=(16, 10))
data_e_ctp = []
data_p_ctp = []

for i, (name, df) in enumerate(ctpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    quantile_start = float(name.split('_')[-1].split('.')[0] + '.' + name.split('_')[-1].split('.')[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_e_ctp.append({'Cumulative Start': quantile_start_percentage, 'Count': e_count, 'File': name})
    data_p_ctp.append({'Cumulative Start': quantile_start_percentage, 'Count': p_count, 'File': name})

df_e_ctp = pd.DataFrame(data_e_ctp)
df_p_ctp = pd.DataFrame(data_p_ctp)
df_e_ctp = df_e_ctp.sort_values('Cumulative Start')
df_p_ctp = df_p_ctp.sort_values('Cumulative Start')

sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_e_ctp, ax=axs[0], palette='cool')
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_p_ctp, ax=axs[1], palette='cool')

axs[0].set_xlabel('Cumulative fraction of genome assessed % -- Counting from bottom')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts (Cumulative count to Top - Prop)')
axs[0].legend().remove()

axs[1].set_xlabel('Cumulative fraction of genome assessed % -- Counting from bottom')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts (Cumulative count to Top - Prop)')
axs[1].legend().remove()

fig.tight_layout()
plt.savefig('histogram_cumulative_top_prop.png')

fig, axs = plt.subplots(2, figsize=(16, 10))
data_e_cbn = []
data_p_cbn = []

for i, (name, df) in enumerate(cbnpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    quantile_start = float(name.split('_')[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_e_cbn.append({'Cumulative Start': quantile_start_percentage, 'Count': e_count, 'File': name})
    data_p_cbn.append({'Cumulative Start': quantile_start_percentage, 'Count': p_count, 'File': name})

df_e_cbn = pd.DataFrame(data_e_cbn)
df_p_cbn = pd.DataFrame(data_p_cbn)
df_e_cbn = df_e_cbn.sort_values('Cumulative Start')
df_p_cbn = df_p_cbn.sort_values('Cumulative Start')

sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_e_cbn, ax=axs[0])
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_p_cbn, ax=axs[1])

axs[0].set_xlabel('Cumulative fraction of genome assessed % -- Counting from top')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts (Cumulative count to Bottom - No Prop)')
axs[0].legend().remove()

axs[1].set_xlabel('Cumulative fraction of genome assessed % -- Counting from top')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts (Cumulative count to Bottom - No Prop)')
axs[1].legend().remove()

fig.tight_layout()
plt.savefig('histogram_cumulative_bottom_noprop.png')

fig, axs = plt.subplots(2, figsize=(16, 10))
data_e_cbp = []
data_p_cbp = []

for i, (name, df) in enumerate(cbpdfs.items()):
    e_count = df[df['enrichment'] == 'e'].shape[0]
    p_count = df[df['enrichment'] == 'p'].shape[0]
    quantile_start = float(name.split('_', 2)[1])
    quantile_start_percentage = int(quantile_start * 100)
    data_e_cbp.append({'Cumulative Start': quantile_start_percentage, 'Count': e_count, 'File': name})
    data_p_cbp.append({'Cumulative Start': quantile_start_percentage, 'Count': p_count, 'File': name})

df_e_cbp = pd.DataFrame(data_e_cbp)
df_p_cbp = pd.DataFrame(data_p_cbp)
df_e_cbp = df_e_cbp.sort_values('Cumulative Start')
df_p_cbp = df_p_cbp.sort_values('Cumulative Start')

sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_e_cbp, ax=axs[0], palette='cool')
sns.barplot(x='Cumulative Start', y='Count', hue='File', data=df_p_cbp, ax=axs[1], palette='cool')

axs[0].set_xlabel('Cumulative fraction of genome assessed % -- Counting from top')
axs[0].set_ylabel('Counts')
axs[0].set_title('Histogram of e counts (Cumulative count to Bottom - Prop)')
axs[0].legend().remove()

axs[1].set_xlabel('Cumulative fraction of genome assessed % -- Counting from top')
axs[1].set_ylabel('Counts')
axs[1].set_title('Histogram of p counts (Cumulative count to Bottom - Prop)')
axs[1].legend().remove()

fig.tight_layout()
plt.savefig('histogram_cumulative_bottom_prop.png')
END_SCRIPT

    cat << 'END_SCRIPT' > GOstatgraphsCumulative.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import io

# Get all .noprop and .prop files in the cumulative results directory
noprop_files = glob.glob('*.noprop')
prop_files = glob.glob('*.prop')
propdfs = {}
nopropdfs = {}

for file in noprop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    nopropdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

for file in prop_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('GO\tNS\t'))
    table = ''.join(lines[start:])
    propdfs[file] = pd.read_csv(io.StringIO(table), sep='\t')

if not nopropdfs and not propdfs:
    print("No cumulative enrichment files found — skipping NS abundance plots")
    exit(0)

data = []

def process_dfs(dfs, df_type):
    for name, df in dfs.items():
        for enrichment in ['e', 'p']:
            df_filtered = df[df['enrichment'] == enrichment]
            counts = df_filtered['NS'].value_counts(normalize=True)
            for ns in ['CC', 'MF', 'BP']:
                if ns in counts:
                    quantile_start = float(name.split('_')[1])
                    data.append({'Quantile Start': quantile_start, 'NS': ns, 'Relative Abundance %': counts[ns], 'Type': df_type, 'Enrichment': enrichment})

process_dfs(nopropdfs, 'noprop')
process_dfs(propdfs, 'prop')

df = pd.DataFrame(data)

for df_type in ['noprop', 'prop']:
    for enrichment in ['e', 'p']:
        plt.figure(figsize=(10, 6))
        sns.barplot(x='Quantile Start', y='Relative Abundance %', hue='NS', data=df[(df['Type'] == df_type) & (df['Enrichment'] == enrichment)], palette='Set2', errorbar=None)
        plt.title(f'Relative Abundances of CC, MF, and BP ({df_type}, {enrichment})')
        plt.savefig(f'Relative Abundances of CC, MF, and BP ({df_type}, {enrichment}).png')
END_SCRIPT

    mkdir -p GO_results
    cp *.noprop *.prop GO_results/ 2>/dev/null || true

    mkdir -p CumulativeGO_results
    cp GO_results/Cumulative*.noprop GO_results/Cumulative*.prop CumulativeGO_results/ 2>/dev/null || true

    cd GO_results
    python ../PlotGOHistograms.py || true
    python ../GOstatgraphs.py || true
    cd ..

    cd CumulativeGO_results
    python ../PlotGOHistogramsCumulative.py || true
    python ../GOstatgraphsCumulative.py || true
    cd ..
    '''
}

workflow {
    def valid_stage3_modes = ["bwa", "coords"]
    if( !valid_stage3_modes.contains(params.stage3_mode) ) {
        error "Invalid --stage3_mode '${params.stage3_mode}'. Valid options: bwa, coords"
    }

    fasta_ch = Channel.fromPath(params.fasta)
    ref_genome = file(params.ref_genome)
    gff_file = file(params.gff)
    cds_list = file(params.cds_list)
    all_txt = file(params.all_txt)
    go_txt = file(params.go_txt)

    stage1_out = stage1(fasta_ch, params.window_size, params.kmer_length, params.increment)
    stage2_out = stage2(fasta_ch, stage1_out.dat_file, stage1_out.fwin_file)
    def stage3_input = stage2_out.filtered_seqs.flatten()
    def stage3_out
    if( params.stage3_mode == "bwa" ) {
        stage3_out = stage3(stage3_input, ref_genome, gff_file)
    }
    else {
        stage3_out = stage3_coords(stage3_input, ref_genome, gff_file, params.window_size as Integer, params.increment as Integer, params.reference_contig)
    }
    stage4_out = stage4(stage3_out.intersect_files, cds_list)
    stage5_out = stage5(stage4_out.protein_files, all_txt, go_txt)
    stage6(stage5_out.noprop.collect(), stage5_out.prop.collect())
}
