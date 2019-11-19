![](images/ikmb_bfx_logo.png)
# Outbreak monitoring read preprocessing

## UNDER CONSTRUCTION

## Overview

This pipeline performs preprocessing of short reads to be used in Outbreak Monitoring. 
* Read trimming and cleaning (Trimmomatic)
* Taxonomic assignment using bloom filters
* Fall-back taxonomic assignment using Pathoscope 2.0.6 (in cases where bloom filters fail)

## Executing the pipeline 

The pipeline starts from a folder (should contain PE FastQ files), generates a Fastqc report, performs adapter trimming (Nextera) and finally runs a subset of 1Mio reads against a set of bloom filters for species identification.

`nextflow -c nextflow.config run main.nf --folder /path/to/reads`

Note: This pipeline cannot be executed directly from gitlab as it requires certain files from the project to be local (the bloom filters). 

## Reading the output

The output will be gathered in the folder "output" (can be modified by using the --outdir flag), divided by library ID. 
In the topfolder you will find one subfolder per input PE library with
* Trimmed reads (using Fastp with automatic adapter recognition)
* Trimming statistics
* Bloomfilter raw outut
* Taxonomic assignment using Pathoscope (if Bloomfilters fail)

In addition, the topfolder will contain one folder per supported species with the trimmed read files. Currently, the following species are supported:
* E. coli
* S. aureus
* S. pneumoniae
* K. pneumoniae
* E. aerogenes
* E. faecalis
* E. faecium
* A. baumannii
* P. aeruginosa
* K. oxytoca
* S. marcescens
* P. mirabilis

Should your data contain any other species, the respective reads will be put in the folder "noMatch" fur manual follow up. Should the reads belong to a species not yet supported, consider adding a new bloom filter. 

Reference genome sequences can be downloaded from:
ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/

To produce a bloom filter, load BioBloom, place the genome sequence into the subfolder "filter, rename it to GENUS_SPECIES.fa and index it with Samtools using `samtools faidx GENUS_SPECIES.fa`.

Then run the Bloomfilter: biobloommaker -p GENUS_SPECIES GENUS_SPECIES.fa



