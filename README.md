# Outbreak monitoring read preprocessing

## UNDER CONSTRUCTION

## Overview

This pipeline performs preprocessing of short reads to be used in Outbreak Monitoring. 
* Read trimming and cleaning (Trimmomatic)
* Taxonomic assignment using bloom filters

## Executing the pipeline 

The pipeline starts from a folder (should contain PE FastQ files), generates a Fastqc report, performs adapter trimming (Nextera) and finally runs a subset of 1Mio reads against a set of bloom filters for species identification.

`nextflow -c nextflow.config run main.nf --folder /path/to/reads`

## Reading the output

The output will be gathered in the folder "output" (can be modified by using the --outdir flag), divided by library ID. 
In the topfolder you will find one subfolder per input PE library with
* Trimmed reads (using Trimmomatic and Nextera adapters)
* FastQC statistics
* Bloomfilter raw outut

In addition, the topfolder will contain one file per input library with a species assignment based on bloom filter mappings. Currently, the following species are supported:
* E. coli
* S. aureus
* S. pneumoniae
* K. pneumoniae
* E. faecalis
* A. baumannii

Should your data contain any other species, it will NOT be detected in this manner. Consider expanding the bloom filter set should additional species be of interest. 



