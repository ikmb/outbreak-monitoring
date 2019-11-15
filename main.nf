/* Preprocessing pipeline for short reads to be used in Outbreak monitoring */

REPORT_SCRIPT = "$baseDir/scripts/report.rb"

BLOOMFILTER = params.bloomfilter

PATHOSCOPE_INDEX_DIR=file(params.pathoscope_index_dir)

OUTDIR = params.outdir

FOLDER=file(params.folder)

log.info "## Outbreak Monitoring Pre-processing"
log.info "Container engine:	${workflow.containerEngine}"
log.info "#####################################"
Channel
        .fromFilePairs(FOLDER + "/*_R{1,2}_001.fastq.gz", flat: true)
        .set { reads }

inputMerge = reads.groupTuple(by: 0)

process Merge {

	tag "${id}"
        publishDir("${OUTDIR}/Data/${id}")

        input:
        set id,file(forward_reads),file(reverse_reads) from inputMerge

        output:
        set val(id),file(left_merged),file(right_merged) into inputBioBloom

        script:
        left_merged = id + "_R1.fastq.gz"
        right_merged = id + "_R2.fastq.gz"

	if (forward_reads.size() > 1 && forward_reads.size() < 100) {
	        """
        	        zcat ${forward_reads.join(" ")} | gzip > $left_merged
			zcat ${reverse_reads.join(" ")} | gzip > $right_merged
	        """
	} else {
		"""	
			cp $forward_reads $left_merged
			cp $reverse_reads $right_merged
		"""
	}
}

process Bloomfilter {

  tag "${id}"
  publishDir "${OUTDIR}/Data/${id}", mode: 'copy'

  input:
  set id,file(left_reads),file(right_reads) from inputBioBloom

  output:
  set id,file(bloom),file(left_reads),file(right_reads) into outputBloomfilter

  script:

  bloom = id + "_summary.tsv"

  """
	biobloomcategorizer -p $id -e -s 0.01 -t 4 -f "$BLOOMFILTER" $left_reads $right_reads

  """

}

process resultBiobloom {

  tag "${id}"
  publishDir "${OUTDIR}/Data/${id}", mode: 'copy'

  input: 
  set id,file(bloom),file(left_reads),file(right_reads) from outputBloomfilter

  output:
  set id,file(bloomresult),file(left_reads),file(right_reads) into inputFastp

  script:

  bloomresult = id + ".species"

  """
	ruby $REPORT_SCRIPT $bloom > $bloomresult
	sleep 30
  """
}

process runFastp {

	publishDir "${OUTDIR}/${organism}", mode: 'copy'

	input:
	set val(id),val(organism_file),fastqR1,fastqR2 from inputFastp

	output:
	set val(id),val(organism),file("*val_1.fq.gz"),file("*val_2.fq.gz") into inputPathoscopeMap
   	set file(json),file(html) into fastp_logs

	script:
	organism = file(organism_file).getText().trim()

	left = file(fastqR1).getBaseName() + "_trimmed.fastq.gz"
	right = file(fastqR2).getBaseName() + "_trimmed.fastq.gz"
	json = file(fastqR1).getBaseName() + ".fastp.json"
	html = file(fastqR1).getBaseName() + ".fastp.html"

	"""
		fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
	"""

}

process runMultiQCFastq {

    tag "Generating fastq level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Fastqc", mode: 'copy'

    input:
    file('*') from fastp_logs.flatten().toList()

    output:
    file("fastq_multiqc*") into runMultiQCFastqOutput

    script:

    """
    multiqc -n fastq_multiqc *.zip *.html
    """
}

process runPathoscopeMap {

   tag "${id}"
   //publishDir "${OUTDIR}/Data/${id}/Pathoscope", mode: 'copy'

   input:
   set id,val(organism),file(left_reads),file(right_reads) from inputPathoscopeMap

   output:
   set id,file(pathoscope_sam) into inputPathoscopeId

   when:
   organism == "noMatch"

   script:
   pathoscope_sam = id + ".sam"

   """
	pathoscope2.py MAP -1 $left_reads -2 $right_reads -indexDir $PATHOSCOPE_INDEX_DIR -filterIndexPrefixes hg19_rRNA \
	-targetIndexPrefix A-Lbacteria.fa,M-Zbacteria.fa,virus.fa -outAlign $pathoscope_sam -expTag $id -numThreads 8
   """

}

process runPathoscopeId {

   tag "${id}"
   publishDir "${OUTDIR}/Data/${id}/Pathoscope", mode: 'copy'

   input:
   set id,file(samfile) from inputPathoscopeId

   output:
   set id,file(pathoscope_tsv) into outputPathoscopeId

   script:

   //pathoscope_sam = "updated_" + samfile
   pathoscope_tsv = id + "-sam-report.tsv"

   """
	pathoscope2.py ID -alignFile $samfile -fileType sam -expTag $id
   """

}
