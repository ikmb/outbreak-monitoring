/* Preprocessing pipeline for short reads to be used in Outbreak monitoring */

params.bloomfilter = "$workflow.projectDir/filter/Acinetobacter_baumannii.bf $workflow.projectDir/filter/Enterococcus_faecalis_V583.bf $workflow.projectDir/filter/Staphylococcus_aureus_NCTC8325.bf $workflow.projectDir/filter/Streptococcus_pneumoniae_R6.bf $workflow.projectDir/filter/Escherichia_coli_K12.bf $workflow.projectDir/filter/Klebsiella_pneumoniae.bf $workflow.projectDir/filter/Pseudomonas_aeruginosa_PAO1.bf"

REPORT_SCRIPT = workflow.projectDir + "/scripts/report.rb"

TRIMMOMATIC = file(params.trimmomatic)
BLOOMFILTER = params.bloomfilter
OUTDIR=params.outdir

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen
adapters = params.adapters

FOLDER=file(params.folder)

Channel
  .fromFilePairs(FOLDER + "/*_R{1,2}_001.fastq.gz", flat: true)
  .into { inputBloomfilter; inputFastqc ; inputReads }

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

process Bloomfilter {

  tag "${id}"
  publishDir "${OUTDIR}/Data/${id}"

  input:
  set id,file(left_reads),file(right_reads) from inputBloomfilter

  output:
  set id,file(bloom),file(left_reads),file(right_reads) into outputBiobloom

  script:

  bloom = id + "_summary.tsv"

  """
	zcat $left_reads | head -n 4000000 | biobloomcategorizer -p $id -f "$BLOOMFILTER" -

  """

}

process resultBiobloom {

  tag "${id}"
  publishDir "${OUTDIR}/Data/${id}"

  input: 
  set id,file(bloom),file(left_reads),file(right_reads) from outputBiobloom

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
   set id,file(bloomresult),file(left_reads),file(right_reads) from inputTrimmomatic
   
   output:
   set id, file("${id}_R1_paired.fastq.gz"), file("${id}_R2_paired.fastq.gz") into outputTrimmomatic

   script:
   
   organism = file("${OUTDIR}/${id}/${id}.species").text.trim()
 
   """
	java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar PE -threads 8 $left_reads $right_reads ${id}_R1_paired.fastq.gz ${id}_1U.fastq.gz ${id}_R2_paired.fastq.gz ${id}_2U.fastq.gz ILLUMINACLIP:${TRIMMOMATIC}/adapters/${adapters}:2:30:10:3:TRUE LEADING:${leading} TRAILING:${trailing} SLIDINGWINDOW:${slidingwindow} MINLEN:${minlen}
   """

}
