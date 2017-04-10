/* Preprocessing pipeline for short reads to be used in Outbreak monitoring */

params.bloomfilter = "$workflow.projectDir/filter/Acinetobacter_baumannii.bf $workflow.projectDir/filter/Enterococcus_faecalis_V583.bf $workflow.projectDir/filter/Staphylococcus_aureus_NCTC8325.bf $workflow.projectDir/filter/Streptococcus_pneumoniae_R6.bf $workflow.projectDir/filter/Escherichia_coli_K12.bf $workflow.projectDir/filter/Klebsiella_pneumoniae.bf"

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
  .into { inputTrimmomatic; inputFastqc }

process Fastqc {

   tag "${id}"
   publishDir "${OUTDIR}/${id}/FastQC/", mode: 'copy'

    input:
    set id, file(left_reads), file(right_reads) from inputFastqc

    output:
    set file("*.zip"), file("*.html") into outputFastqc

    script:
    """
    fastqc -t 1 -o . ${left_reads} ${right_reads}
    """

}

process Trimmomatic {

   tag "${id}"
   publishDir "${OUTDIR}/${id}/trimmomatic"

   input:
   set id, file(left_reads), file(right_reads) from inputTrimmomatic

   output:
   
   set id, file("${id}_1P.fastq.gz"), file("${id}_2P.fastq.gz") into (inputBiobloom)

   script:

    """
        java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar PE -threads 8 $left_reads $right_reads ${id}_1P.fastq.gz ${id}_1U.fastq.gz ${id}_2P.fastq.gz ${id}_2U.fastq.gz ILLUMINACLIP:${TRIMMOMATIC}/adapters/${adapters}:2:30:10:3:TRUE LEADING:${leading} TRAILING:${trailing} SLIDINGWINDOW:${slidingwindow} MINLEN:${minlen}
    """

}

process Bloomfilter {

  tag "${id}"
  publishDir "${OUTDIR}/${id}"

  input:
  set id,file(left_reads),file(right_reads) from inputBiobloom

  output:
  set id,file(bloom) into outputBiobloom

  script:

  bloom = id + "_summary.tsv"

  """
	zcat $left_reads | head -n 4000000 | biobloomcategorizer -p $id -f "$BLOOMFILTER" -

  """

}

process resultBiobloom {

  tag "${id}"
  publishDir "${OUTDIR}"

  input: 
  set id,file(bloom) from outputBiobloom

  output:
  set id,file(bloomresult) into outputBloomresult

  script:

  bloomresult = id + ".species"

  """
	#!/bin/env ruby

	lines = IO.readlines('$bloom')

	header = lines.shift.split(' ')

	species = "unknown_species"
	lines.each do |line|
		elements = line.strip.split(' ')
		hit_rate = elements[4].to_f
		hits = elements[1].to_i
		organism = elements[0]
		species = organism if hits > 500000	
	end		

	f = File.new('$bloomresult','w+')
	f.puts species
	f.close
  """
			
}
