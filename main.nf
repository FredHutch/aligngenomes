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

    Other options:
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

/*
 */

