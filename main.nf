 #! /usr/bin/env nextflow
nextflow.enable.dsl=2


process mk_directories {
  //since this is the first process we need to load singularity in order to 
  //use the singularity containers
  //beforeScript "module load $params.singularity_version"

  tag "mkDir"
  
  output:
  val 'First Process'

  script:
  """
  mkdir -p $params.outdir_fastq &&
  mkdir -p $params.outdir_trim &&
  mkdir -p $params.outdir_align &&
  mkdir -p $params.outdir_peak &&
  mkdir -p $params.singularity_folder 
  """

}

process load_modules {

  tag "load_mod"

  output:
  val 'Second Process'

  input:
  val _

  script:
  """
  module load singularity/3.7.0
  """

}

process trim {
  //errorStrategy 'ignore'
  queue = "$params.queue"
  container = "https://depot.galaxyproject.org/singularity/trim-galore:0.6.7--hdfd78af_0"

  tag "trim_galore 7"
  publishDir params.outdir_trim, mode: 'copy'

  input:
  path file1
  path file2
  val _

  output:
  //path "$params.outdir_trim*"
  //path "saida_trim"

  //trim_galore --paired -o $params.outdir_trim $params.fastq/* --fastqc
  script:
  """
  trim_galore --paired $file1 $file2 -o $params.outdir_trim
  """
}


process fastqc {
  //errorStrategy 'ignore'
  queue = "$params.queue"
  container ='https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0'

  tag "fastqc t" 
  publishDir params.outdir_fastq, mode : 'copy'
  
  input:
  path query_file
  val _

  
  script:
  """
  fastqc -o $params.outdir_fastq $query_file
  """
}

process align {
  //errorStrategy 'ignore'
  queue = "$params.queue"
  container ='https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

  tag "bwa mem" 
  publishDir params.outdir_align, mode : 'copy'
  
  input:
  path file1
  path file2
  path fileAlign
  val _

  script:
  """
  bwa mem $params.align_ref/hg19.fa $file1 $file2 | samtools view -Sb - > $params.outdir_align/test7.bam
  """
}

process sort {
  //errorStrategy 'ignore'
  queue = "$params.queue"
  container ='https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0'

  tag "sort samtools" 
  publishDir params.outdir_align, mode : 'copy'
  
  input:
  path fileAlign
  val _

  script:
  """
  samtools sort $params.outdir_align/test7.bam -o $params.outdir_align/test7.sorted.bam
  """
}


process dedup {
  queue = "$params.queue"
  container = 'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0'

  tag "dedup metric 9test" 
  publishDir params.outdir_align, mode : 'copy'

  input:
  path fileAlign
  val _
  // ASSUME_SORTED=true \\
  script:
  """
  picard MarkDuplicates \\
  I=$params.outdir_align/test7.sorted.bam \\
  O=$params.outdir_align/test7.dedup.bam \\
  REMOVE_DUPLICATES=true \\
  ASSUME_SORT_ORDER=coordinate \\
  VALIDATION_STRINGENCY=LENIENT \\
  METRICS_FILE=$params.outdir_align/MarkDuplicates.metrics.txt
	"""
}

process peak_calling{

  queue = "$params.queue"
  container = 'https://depot.galaxyproject.org/singularity/macs2:2.2.7.1--py38h4a8c8d9_3'

  tag "peak calling" 
  publishDir params.outdir_align, mode : 'copy'


  input:
  path fileAlign
  val _
  
  script:
  """
  macs2 \\
  callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -f BAMPE --extsize 146 --nomodel \\
  -t $params.outdir_align/test7.bam \\
  --outdir $params.outdir_peak \\
  -n test7
  """
}

workflow {

    def rawFiles = Channel.fromPath( "$params.fastq/*" )
    def file1 = rawFiles.first()
    def file2 = rawFiles.last()

    // this is apperantly not been used by the process but if is not passed the process doesn't read the files pointed
    def alignPath = Channel.fromPath("/data/baca/users/bjf35/cfchip_pipeline/ref_files/hg19/bwa_indices/hg19/hg19.fa.sa")

    
    //this channels is only created in order to make sure that the directories will be
    //created before the first process starts otherwise it will break the workflow
    ch1 = mk_directories()
    ch2 = load_modules(ch1)
    fastqc(rawFiles,ch2)
    trim(file1, file2,ch2)
    align(file1,file2,alignPath,ch2)
    sort(alignPath,ch2)
    dedup(alignPath,ch2)
    peak_calling(alignPath,ch2)

}

