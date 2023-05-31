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
  mkdir -p $params.outdir_align
  
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

  tag "fastqc" 
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
  //container ='https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'
  container = 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'
  tag "bwa mem 76665" 
  //publishDir params.outdir_align, mode : 'copy'
  
  input:
  path file1
  path file2
  path fileAlign
  val _

  output:
  path "test.bam"

  //bwa mem $fileAlign $file1 $file2 | samtools view -Sb - > test.bam
  script:
  """
  bwa mem $fileAlign $file1 $file2
  """
}

workflow {

    def rawFiles = Channel.fromPath( "$params.fastq/*" )
    def file1 = rawFiles.first()
    def file2 = rawFiles.last()
    def alignPath = Channel.fromPath("$params.align_ref/*")
    def fileAlign = alignPath.first()
    
    //this channels is only created in order to make sure that the directories will be
    //created before the first process starts otherwise it will break the workflow
    ch1 = mk_directories()
    ch2 = load_modules(ch1)
    fastqc(rawFiles,ch2)
    trim(file1, file2,ch2)
    align(file1,file2,fileAlign,ch2)



}

