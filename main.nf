#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/aligngenomes
========================================================================================
 nf-core/aligngenomes Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/aligngenomes
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    Workflow to align a set of reads in FASTQ format against a set of reference genome(s) in FASTA format.

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/aligngenomes --reads 'reads/*.fastq.gz' --genomes 'genomes/ref1.fasta.gz,genomes.ref2.fasta.gz' -profile docker

    Arguments:
      --reads                       Path to reads in FASTQ format (must be surrounded with quotes)
      --genomes                     Path to genomes in FASTA format (must be surrounded with quotes)
      --outdir                      The output directory where the results will be saved
      --output_csv                  Name used to create the final summary CSV

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

/*
 * Create a channel for input read files and genome files
 */
reads_ch = Channel
    .from(params.reads.split(","))
    .map { it -> file(it) }
    .ifEmpty { exit 1, "params.reads was empty - no input files supplied" }

genomes_ch = Channel
    .from(params.genomes.split(","))
    .map { it -> file(it) }
    .ifEmpty { exit 1, "params.genomes was empty - no input files supplied" }


// Make sure that every read has a unique name
process correctHeaders {
  container "ubuntu:16.04"
  errorStrategy "retry"
  
  input:
  file fastq from reads_ch
  
  output:
  file "${fastq.name.replace(".gz", "").replace(".fastq", "")}.unique.headers.fastq.gz" into reads_for_alignment, reads_for_counting

  afterScript "rm *"

  """
set -e

echo Checking for input file existance
[ -s "${fastq}" ]

echo "Correcting headers"
gunzip -c "${fastq}" | \
awk '{if(NR % 4 == 1){print("@" 1 + ((NR - 1) / 4))}else{print}}' | \
gzip -c > \
"${fastq.name.replace(".gz", "").replace(".fastq", "")}.unique.headers.fastq.gz"

echo Done
  """

}

process alignGenome {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  errorStrategy 'retry'
  publishDir "${params.outdir}/bam/"

  input:
  file genome_fasta from genomes_ch
  each file(input_fastq) from reads_for_alignment
  val threads from 8
  
  output:
  set file("*.bam"), val("${genome_fasta.name}"), val("${input_fastq.name}") into bam_ch

  afterScript "rm -r *"

  """
#!/bin/bash

set -e 

echo "Processing ${genome_fasta}"

if [[ \$(gzip -t ${genome_fasta}) ]]; then

    echo "Decompressing ${genome_fasta}"
    gunzip ${genome_fasta}
    genome_fasta="${genome_fasta.name.replace(".gz\$", "")}"

else

    echo "${genome_fasta} is not compressed"
    genome_fasta="${genome_fasta}"

fi

echo "Indexing \$genome_fasta"
bwa index \$genome_fasta

echo "Aligning ${input_fastq} against \$genome_fasta"

bwa mem -t ${threads} \$genome_fasta ${input_fastq} | \
samtools view -b -F 4 -o ${input_fastq}.\$genome_fasta.bam
echo "Done aligning to \$genome_fasta"

    """

}

process alignmentStats {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  errorStrategy 'retry'
  publishDir "${params.outdir}/stats/"

  input:
  set file(bam), val(genome_fasta_name), val(input_fastq_name) from bam_ch
  
  output:
  set file("${bam}.idxstats"), file("${bam}.stats"), file("${bam}.pileup"), file("${bam}.positions"), val("${genome_fasta_name}"), val("${input_fastq_name}") into stats_ch

  afterScript "rm -r *"

  """
#!/bin/bash

set -e

echo "Processing ${bam}"

echo "Sorting alignments"
samtools sort ${bam} > ${bam}.sorted

echo "Calculating stats"
samtools stats ${bam}.sorted > ${bam}.stats

echo "Indexing alignments"
samtools index ${bam}.sorted

echo "Calculating idxstats"
samtools idxstats ${bam}.sorted > ${bam}.idxstats

echo "Formatting pileup"
samtools mpileup ${bam}.sorted > ${bam}.pileup

# Make a file with three columns, the bitwise flag, and the leftmost position, and the length of the mapped segment
echo "Formatting positions TSV"
samtools view ${bam}.sorted | awk '{print \$2 "\\t" \$4 "\\t" length(\$10)}' > ${bam}.positions

echo "Done"
  """

}

process summarizeEach {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  errorStrategy 'retry'

  input:
  set file(idxstats), file(stats), file(pileup), file(positions), val(genome_fasta_name), val(input_fastq_name) from stats_ch
  
  output:
  file "*json" into all_stats_ch

  afterScript "rm -r *"

  """
#!/usr/bin/env python3
import os
import json
import pandas as pd
from math import log as ln

def read_line(fp, prefix):
    with open(fp, 'rt') as f:
        for line in f:
            if line.startswith(prefix):
                return line.replace(prefix, '').strip(" ").strip("\\t")

def parse_flags(int_flags):
    output = dict()
    for n, flag in [
        (2048, "supplementary"),
        (1024, "duplicate"),
        (512, "fail_filter"),
        (256, "secondary"),
        (128, "last"),
        (64, "first"),
        (32, "next_rc"),
        (16, "rc"),
        (8, "next_unmapped"),
        (4, "unmapped"),
        (2, "aligned_properly"),
        (1, "multiple_segments"),
    ]:
        if int_flags >= n:
            output[flag] = True
            int_flags = int_flags - n
        else:
            output[flag] = False
    assert int_flags == 0, int_flags
    return output


def shannon_divesity(counts):
    # See https://gist.github.com/audy/783125
    
    def p(n, N):
        # Relative abundance
        if n is 0:
            return 0
        else:
            return (float(n)/N) * ln(float(n)/N)
            
    N = sum(counts)
    
    return -sum(p(n, N) for n in counts if n is not 0)


# Read in the summary of alignment positions
positions = pd.read_csv("${positions}", sep="\\t", header=None)
positions = pd.concat([positions, pd.DataFrame(map(parse_flags, positions[0]))], axis=1)
# If the read is aligned in the forward direction, use the leftmost position, otherwise use the rightmost
position_list = positions.apply(
    lambda r: r[1] + r[2] if r["rc"] else r[1],
    axis=1
).tolist()

# Calculate Shannon diversity
sdi = shannon_divesity(pd.Series(position_list).value_counts().values)

pileup = pd.read_csv("${pileup}", sep="\\t", header=None)
n_reads = int(read_line("${stats}", "SN\\treads mapped:"))
reflen = int(pd.read_csv("${idxstats}", sep="\\t", header=None).iloc[0, 1])
covlen = pileup[3].shape[0]
nerror = float(read_line("${stats}", "SN\\terror rate:").split("\\t")[0])
nbases = pileup[3].sum()

output = dict()
output["alignment_file"] = "${pileup}".replace(".pileup", "")
output["depth"] = nbases / reflen
output["coverage"] = covlen / reflen
output["error"] = nerror
output["genome_length"] = reflen
output["n_reads"] = len(position_list)
output["entropy"] = sdi
output["genome"] = "${genome_fasta_name}"
output["input_fastq"] = "${input_fastq_name}"

json_fp = "${pileup}".replace(".pileup", ".json")
assert json_fp.endswith(".json")
with open(json_fp, "wt") as fo:
    fo.write(json.dumps(output, indent=4))

assert os.path.exists(json_fp)

  """

}

// Count the number of reads
process countReads {
  container "ubuntu:16.04"
  errorStrategy 'retry'
  
  input:
  file fastq from reads_for_counting
  
  output:
  file "${fastq}.counts.csv" into counts_ch

  afterScript "rm *"

  """
#!/bin/bash

set -e

gzip -t "${fastq}"

n=\$(gunzip -c "${fastq}" | awk 'NR % 4 == 1' | wc -l)

(( \$n > 0 ))

echo "${fastq},\$n" > "${fastq}.counts.csv"

  """

}

process collectCounts {
  container "ubuntu:16.04"
  errorStrategy 'retry'

  input:
  file readcounts from counts_ch.collect()
  
  output:
  file "readcounts.csv" into readcounts_csv

  afterScript "rm -r *"

  """
#!/bin/bash

set -e

echo file,n_reads > TEMP
cat *csv >> TEMP && rm *csv && mv TEMP readcounts.csv
  """

}

process collectAll {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  publishDir "${params.outdir}/"
  errorStrategy 'retry'

  input:
  file all_jsons from all_stats_ch.collect()
  file readcounts_csv
  
  output:
  file "${params.output_csv}"

  afterScript "rm -r *"

  """
#!/usr/bin/env python3
import os
import json
import pandas as pd

readcounts = pd.read_csv("${readcounts_csv}").set_index(
    "file"
)["n_reads"].to_dict()

df = pd.DataFrame([
    json.load(open(fp, "rt"))
    for fp in os.listdir(".")
    if fp.endswith(".json")
])

df["total_reads"] = df["input_fastq"].apply(readcounts.get)
assert df["total_reads"].isnull().sum() == 0
assert (df["total_reads"] > 0).all()
df["prop_reads"] = df["n_reads"] / df["total_reads"]

df.to_csv("${params.output_csv}", index=None)

"""

}
