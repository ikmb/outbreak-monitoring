# Outbreak monitoring read preprocessing

## UNDER CONSTRUCTION

## Overview

This pipeline performs preprocessing of short reads to be used in Outbreak Monitoring. 
* Read trimming and cleaning (Trimmomatic)
* Taxonomic assignment using Metaphlan

## Executing the pipeline 

The pipeline starts from a folder (should contain PE FastQ files), generates a Fastqc report, performs adapter trimming (Nextera) and finally runs a subset of 1Mio reads against a set of bloom filters for species identification.

`nextflow -c nextflow.config run main.nf --folder /path/to/reads`

The output will be gathered in the folder "output" (can be modified by using the --outdir flag), divided by library ID. 


