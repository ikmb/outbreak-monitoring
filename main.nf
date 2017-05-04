/* Preprocessing pipeline for short reads to be used in Outbreak monitoring */

params.bloomfilter = "$workflow.projectDir/filter/Acinetobacter_baumannii.bf $workflow.projectDir/filter/Enterococcus_faecalis_V583.bf $workflow.projectDir/filter/Staphylococcus_aureus_NCTC8325.bf $workflow.projectDir/filter/Streptococcus_pneumoniae_R6.bf $workflow.projectDir/filter/Escherichia_coli_K12.bf $workflow.projectDir/filter/Klebsiella_pneumoniae.bf $workflow.projectDir/filter/Pseudomonas_aeruginosa_PAO1.bf $workflow.projectDir/filter/Enterobacter_cloacae.bf $workflow.projectDir/filter/Enterobacter_aerogenes.bf $workflow.projectDir/filter/Enterococcus_faecium.bf"

REPORT_SCRIPT = workflow.projectDir + "/scripts/report.rb"

TRIMMOMATIC = file(params.trimmomatic)
BLOOMFILTER = params.bloomfilter
OUTDIR=params.outdir

PATHOSCOPE_INDEX_DIR=file(params.pathoscope_index_dir)
PATHOSCOPE=file(params.pathoscope)

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen
adapters = params.adapters

FOLDER=file(params.folder)

Channel
  .fromFilePairs(FOLDER + "/*_R{1,2}_001.fastq.gz", flat: true)
  .into { inputBioBloom }

process Bloomfilter {

  tag "${id}"
  publishDir "${OUTDIR}/Data/${id}"

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
  publishDir "${OUTDIR}/Data/${id}"

  input: 
  set id,file(bloom),file(left_reads),file(right_reads) from outputBloomfilter

  output:
  set id,file(bloomresult),file(left_reads),file(right_reads) into inputTrimmomatic

  script:

  bloomresult = id + ".species"

  """
	ruby $REPORT_SCRIPT $bloom > $bloomresult
  """
}


process Trimmomatic {

   tag "${id}"
   publishDir "${OUTDIR}/${organism}", mode: 'copy'

   input:
   set id,file(organism_file),file(left_reads),file(right_reads) from inputTrimmomatic

   output:
   set id,organism,file("${id}_R1_paired.fastq.gz"), file("${id}_R2_paired.fastq.gz") into inputFastqc,inputPathoscopeMap

   script:

   organism = file("${OUTDIR}/Data/${id}/${id}.species").text.trim()   

    """
        java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar PE -threads 8 $left_reads $right_reads \
        ${id}_R1_paired.fastq.gz ${id}_1U.fastq.gz ${id}_R2_paired.fastq.gz ${id}_2U.fastq.gz \
        ILLUMINACLIP:${TRIMMOMATIC}/adapters/${adapters}:2:30:10:3:TRUE\
        LEADING:${leading} TRAILING:${trailing} SLIDINGWINDOW:${slidingwindow} MINLEN:${minlen} && sleep 5
   """

}

process Fastqc {

   tag "${id}"
   publishDir "${OUTDIR}/Data/${id}/FastQC/", mode: 'copy'

    input:
    set id, file(left_reads), file(right_reads) from inputFastqc

    output:
    set file("*.zip"), file("*.html") into outputFastqc

    script:
    """
    fastqc -t 1 -o . ${left_reads} ${right_reads}
    """

}

process runMultiQCFastq {

    tag "Generating fastq level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Fastqc"

    input:
    file('*') from outputFastqc.flatten().toList()

    output:
    file("fastq_multiqc*") into runMultiQCFastqOutput

    script:

    """
    multiqc -n fastq_multiqc *.zip *.html
    """
}

process runPathoscopeMap {

   tag "${id}"
   publishDir "${OUTDIR}/Data/${id}/Pathoscope"

   input:
   set id,organism,file(left_reads),file(right_reads) from inputPathoscopeMap

   output:
   set id,file(pathoscope_sam) into inputPathoscopeId

   when:
   organism == "noMatch"

   script:
   pathoscope_sam = id + ".sam"

   """
	$PATHOSCOPE MAP -1 $left_reads -2 $right_reads -indexDir $PATHOSCOPE_INDEX_DIR -filterIndexPrefixes hg19_rRNA \
	-targetIndexPrefix A-Lbacteria.fa,M-Zbacteria.fa,virus.fa -outAlign $pathoscope_sam -expTag $id -numThreads 8
   """

}

process runPathoscopeId {

   tag "${id}"
   publishDir "${OUTDIR}/Data/${id}/Pathoscope"

   input:
   set id,file(samfile) from inputPathoscopeId

   output:
   set id,file(pathoscope_tsv),file(pathoscope_sam) into outputPathoscopeId

   script:

   pathoscope_sam = "updated_" + samfile
   pathoscope_tsv = id + "-sam-report.tsv"

   """
	$PATHOSCOPE ID -alignFile $samfile -fileType sam -expTag $id
   """

}
