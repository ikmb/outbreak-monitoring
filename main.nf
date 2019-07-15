/* Preprocessing pipeline for short reads to be used in Outbreak monitoring */

REPORT_SCRIPT = workflow.projectDir + "/scripts/report.rb"

BLOOMFILTER = params.bloomfilter

PATHOSCOPE_INDEX_DIR=file(params.pathoscope_index_dir)

params.saveTrimmed = true

if (params.ariba.containsKey(params.ariba_db) == false) {
   exit 1, "Specified unknown ariba database, please consult the documentation for valid databases."
}

ARIBA_DB=params.ariba[params.ariba_db].database

// Trimming parameters
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

FOLDER=file(params.folder)

Channel
        .fromFilePairs(FOLDER + "/*_R{1,2}_001.fastq.gz", flat: true)
        .map { prefix, file1, file2 ->  tuple(prefix.split("_")[0..1].join("_"), file1, file2) }
        .groupTuple()
        .into { inputMerge }

process Merge {

	tag "${id}"
        publishDir("${OUTDIR}/Data/${id}")

        input:
        set id,forward_reads,reverse_reads from inputMerge

        output:
        set id,file(left_merged),file(right_merged) into inputBioBloom

        script:
        left_merged = id + "_R1.fastq.gz"
        right_merged = id + "_R2.fastq.gz"

        """
                zcat ${forward_reads.join(" ")} | gzip > $left_merged
		zcat ${reverse_reads.join(" ")} | gzip > $right_merged
        """
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
  set id,file(bloomresult),file(left_reads),file(right_reads) into inputTrimgalore

  script:

  bloomresult = id + ".species"

  """
	ruby $REPORT_SCRIPT $bloom > $bloomresult
	sleep 30
  """
}


process runTrimgalore {

   tag "${id}"
   publishDir "${OUTDIR}/${organism}", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else params.saveTrimmed ? filename : null
        }

   input:
   set val(id),val(organism_file),left,right from inputTrimgalore

   output:
   set val(id),val(organism),file("*val_1.fq.gz"),file("*val_2.fq.gz") into (inputPathoscopeMap, inputAriba)
   file "*trimming_report.txt" into trimgalore_results, trimgalore_logs 
   file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
   
   script:
   organism = file(organism_file).getText().trim()
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''

    """
    trim_galore --paired --fastqc --length 35 --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $left $right
    """

}

process runAriba {

   tag "${id}"
   publishDir "${OUTDIR}/${organism}/Ariba/${params.ariba_db}", mode: 'copy'

   input:
   set id,val(organism),file(left), file(right) from inputAriba

   output:
   set val(organism),file(report) into AribaReport

   when:
   params.antibiotics == true

   script:

   report = "report.${id}.tsv"

   """
	ariba run $ARIBA_DB $left $right out.${id}.run && cp out.${id}.run/report.tsv report.${id}.tsv
   """	

}

AribaReport
	.groupTuple(by: 0 )
	.set { AribaReportByOrganism }

process runAribaSummary {

   tag "SummarizeAriba (ALL)"
   publishDir "${OUTDIR}/${organism}/Ariba", mode: 'copy'

   input:
   set organism,reports from AribaReportByOrganism

   output:
   set file(summary),file(summary_phandango),file(summary_phandango_tre) into outputAribaSummary

   when:
   reports.size() > 1

   script:

   base = "ariba.summary"
   summary = "ariba.summary.csv"
   summary_phandango = "ariba.summary.phandango.csv"
   summary_phandango_tre = "ariba.summary.phandango.tre"

   """
        ariba summary $base ${reports.join(" ")}
   """

}

process runMultiQCFastq {

    tag "Generating fastq level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Fastqc", mode: 'copy'

    input:
    file('*') from trimgalore_fastqc_reports.flatten().toList()

    output:
    file("fastq_multiqc*") into runMultiQCFastqOutput

    script:

    """
    /ifs/data/nfs_share/ikmb_repository/software/multiqc_local/1.2/bin/multiqc -n fastq_multiqc *.zip *.html
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
	$PATHOSCOPE MAP -1 $left_reads -2 $right_reads -indexDir $PATHOSCOPE_INDEX_DIR -filterIndexPrefixes hg19_rRNA \
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
	$PATHOSCOPE ID -alignFile $samfile -fileType sam -expTag $id
   """

}
