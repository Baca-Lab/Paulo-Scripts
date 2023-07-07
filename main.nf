 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

String path_fastq,path_trim,path_align,path_lib_complex,path_peak,path_bigwig,path_bamtobed,path_avalysis,path_snp

process mk_dir {
  tag "creating dir for $sampleId"

  input:
  tuple val(sampleId), val(path),val(read1), val(read2)

  output:
  val 'First Process'

  exec:
  path_avalysis = path + "/" +  params.oudir_analysis + "/"
  
  path_fastq = path_avalysis + params.outdir_fastq
  path_trim = path_avalysis + params.outdir_trim
  path_align = path_avalysis + params.outdir_align
  path_lib_complex = path_avalysis + params.outdir_lib_complex
  path_peak = path_avalysis + params.outdir_peak
  path_bigwig = path_avalysis + params.outdir_bigwig
  path_bamtobed = path_avalysis + params.outdir_bamtobed
  path_snp = path_avalysis + params.outdir_snp

  script:
  """
  mkdir -p $path_avalysis &&
  mkdir -p $path_fastq &&
  mkdir -p $path_trim &&
  mkdir -p $path_align &&
  mkdir -p $path_lib_complex &&
  mkdir -p $path_peak &&
  mkdir -p $path_bigwig &&
  mkdir -p $path_bamtobed &&
  mkdir -p $path_snp
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



process lib_complex {
  queue = "$params.queue"
  
  //Singularity Image
  //container = "https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0"

  //Docker Image
  container = "quay.io/biocontainers/picard:2.27.4--hdfd78af_0"

  tag "4.1 - lib complexity #$sampleId"
  publishDir "$path/analysis/$params.outdir_lib_complex/", mode : 'copy'

  input:
  path sampleBam
  tuple val(sampleId), val(path),path(read1), path(read2)

  output:
  path("*.csv")

  exec:
  String strBam = sampleId + '_lib_complex.csv'

  script:
  """
  picard EstimateLibraryComplexity I=$sampleBam  O=$strBam
  """
}


process align {
  queue = "$params.queue"

  //Singularity Image
  //container ='https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

  //Docker Image
  container = 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

  tag "3 - align4 #$sampleId" 
  publishDir "$path/analysis/$params.outdir_align/", mode : 'copy'
  
  input:
  tuple path(file1),path(file2)
  tuple val(sampleId), val(path),path(read1), path(read2)
  
  output:
  path("*.bam")

  exec:
  String strBam = sampleId + '.bam'

  
  script:
  """
  bwa mem -v 0 $params.align_ref/hg19.fa $file1 $file2 | samtools view -Sb -u > $strBam
  """
}

process sort {
  //errorStrategy 'ignore'
  queue = "$params.queue"

  //Singularity Image
  //container ='https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0'

  //Docker Image
  container ='quay.io/biocontainers/samtools:1.15.1--h1170115_0'

  tag "4 - sort #$sampleId" 
  publishDir "$path/analysis/$params.outdir_align/", mode : 'copy'
  
  input:
  path sampleBam
  tuple val(sampleId), val(path),path(read1), path(read2)

  output:
  path("*.bam")

  exec:
  String strBam = sampleId + '.sorted.bam'

  script:
  """
  samtools sort $sampleBam -o $strBam
  """
}

process unique_sam {
  queue = "$params.queue"

  //Singularity Image
  //container ='https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0'

  //Docker Image
  container ='quay.io/biocontainers/samtools:1.15.1--h1170115_0'

  tag "4.5 - unique #$sampleId" 
  publishDir "$path/analysis/$params.outdir_align/", mode : 'copy'
  
  input:
  path sampleBam
  tuple val(sampleId), val(path),path(read1), path(read2)

  exec:
  String strBam = sampleId + '.unique.sorted.bam'

  output:
  path("*.bam")

  script:
  """
  samtools view -b -q 1 $sampleBam > $strBam
  """
}

process enrichment {
  queue = "$params.queue"

  //Singularity Image
  //container ='https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0'

  //Docker Image
  container ='quay.io/biocontainers/samtools:1.15.1--h1170115_0'

  tag "9.2 On/Off Enrichment #$sampleId" 
  publishDir "$path/analysis/$params.outdir_bamtobed/", mode : 'copy'

  input:
  tuple path (sampleBam), val(_)
  tuple val(sampleId), val(path), val(_), val(_)
  path (scriptFile)

  exec:
  strCSV = sampleId + '_total_enrichment.csv'

  output:
  path("*.csv")

  script:
  """
  sh $launchDir/enrichment.sh $sampleBam $params.states_ref $sampleId >> $strCSV
  """
}

process index_sam {
  queue = "$params.queue"

  //Singularity Image
  //container ='https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0'

  //Docker Image
  container ='quay.io/biocontainers/samtools:1.15.1--h1170115_0'

  tag "6 - index before peak calling #$sampleId" 
  publishDir "$path/analysis/$params.outdir_align/", mode : 'copy'
  
  input:
  tuple path (sampleBam), path(txtFile)
  tuple val(sampleId), val(path),path(read1), path(read2)

  output:
  path ('*.bai')

  script:
  """
  samtools index $sampleBam
  """
}


process dedup {
  queue = "$params.queue"

  //Singularity Image
  //container = 'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0'

  //Docker Image
  container = 'quay.io/biocontainers/picard:2.27.4--hdfd78af_0'

  tag "5 - dedup #$sampleId" 
  publishDir "$path/analysis/$params.outdir_align/", mode : 'copy'

  input:
  path sampleBam
  tuple val(sampleId), val(path),path(read1), path(read2)

  exec:
  String strBam = sampleId + '.dedup.unique.sorted.bam'
  String strTxt = sampleId + '-MarkDuplicates.metrics.txt'

  output:
  tuple path("*.bam"),path("*.txt")

  script:
  """
  picard MarkDuplicates \\
  I=$sampleBam \\
  O=$strBam \\
  REMOVE_DUPLICATES=true \\
  ASSUME_SORT_ORDER=coordinate \\
  VALIDATION_STRINGENCY=LENIENT \\
  METRICS_FILE=$strTxt
	"""
}

process fetch_chrom_sizes{

  queue = "$params.queue"
  
  //Docker Image
  container = 'quay.io/biocontainers/ucsc-fetchchromsizes:377--ha8a8165_3'

  tag "7.1 - fetch chromsizes #$sampleId" 
  publishDir "$path/analysis/$params.outdir_bigwig/", mode : 'copy'


  input:
  tuple val(sampleId), val(path),path(read1), path(read2)

  output:
  path ('*.sizes')
  
  script:
  """
  fetchChromSizes hg19 > hg19.chrom.sizes
  """
}

process peak_bed_graph{
  queue = "$params.queue"

  //Singularity Image
  //container = 'https://depot.galaxyproject.org/singularity/macs2:2.2.7.1--py38h4a8c8d9_3'

  //Docker Image
  container = 'quay.io/biocontainers/macs2:2.2.7.1--py38h4a8c8d9_3'

  tag "6 - peak calling e bed graph #$sampleId" 
  publishDir "$path/analysis/$params.outdir_peak/", mode : 'copy'

  input:
  tuple path (sampleBam), path(txtFile)
  val _
  tuple val(sampleId), val(path),path(read1), path(read2)

  output:
  tuple path ('*.*'),path ('*.bdg')
  
  script:
  """
  macs2 \\
  callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -f BAMPE --extsize 146 --nomodel \\
  -t $sampleBam \\
  -n $sampleId --bdg
  """
}


process bam_to_bed {
  queue = "$params.queue"

  //Singularity Image
  //container ='https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0'

  //Docker Image
  container ='quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0'

  tag "9.0 - Bam to Bed #$sampleId" 
  publishDir "$path/analysis/$params.outdir_bamtobed/", mode : 'copy'
  
  input:
  tuple path (sampleBam), path(txtFile)
  tuple val(sampleId), val(path),path(read1), path(read2)

  exec:
  String strBed = sampleId + '.bed'

  output:
  path ('*.bed')

  script:
  """
  bedtools bamtobed -i \\
  $sampleBam -bedpe 2> /dev/null | \\
  awk 'BEGIN{{OFS="\t";FS="\t"}}(\$1==\$4){{print \$1, \$2, \$6}}' > $strBed
  """
}

process unique_frags {
  queue = "$params.queue"

  tag "9.1 - Unique Fragments #$sampleId" 
  publishDir "$path/analysis/$params.outdir_bamtobed/", mode : 'copy'

  input:
  path (sampleBed)
  tuple val(sampleId), val(path),path(read1), path(read2)

  output:
  path ('*.csv')

  exec:
  String strCSV = sampleId + '.csv'
  
  script:
  """
  echo $sampleId && wc -l $sampleBed | cut -f1 -d' '  >> $strCSV
  """

}


process fastqc {
  //errorStrategy 'ignore'
  //debug true
  queue = "$params.queue"
  
  //Singularity Image
  //container ='https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0'
  
  //Docker Image
  container = 'quay.io/biocontainers/fastqc:0.11.9--0'

  tag "1 - fastqc #$sampleId" 
  publishDir params.outdir_fastq, mode : 'copy'
  
  input:
  tuple val(sampleId), val(path),path(read1), path(read2)
  val _

  //exec:
  //path_avalysis = path + "/" +  params.oudir_analysis + "/"
  //path_fastq = path_avalysis + params.outdir_fastq + "/"
  
  script:
  """
  fastqc -o $path_fastq $read1 $read2
  """
}

process trim {
  //errorStrategy 'ignore'
  queue = "$params.queue"
  //debug true

  //Singularity Image
  //container = "https://depot.galaxyproject.org/singularity/trim-galore:0.6.7--hdfd78af_0"

  //Docker Image
  container = "quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0"

  tag "2 - trim_galore #$sampleId"
  publishDir "$path_trim", mode: 'copy'

  input:
  tuple val(sampleId), val(path),path(read1), path(read2)
  val _

  output:
  tuple path("*1.fq.gz"),path("*2.fq.gz") 
  //path "", emit:read2
  //path "*report.txt" , emit: log

  
  //exec:
  //path_avalysis = path + "/" +  params.oudir_analysis + "/"
  //path_trim = path_avalysis + params.outdir_trim + "/"

  script:
  """
  trim_galore --paired $read1 $read2 --gzip
  """
}

process snp_fingerprint {
  queue = "$params.queue"
  
  //Docker Image
  container = 'quay.io/mcrotti1/bcftools'

  tag "8 - SNP Finger Print #$sampleId" 
  publishDir "$path/analysis/$params.outdir_snp/", mode : 'copy'

  input:
  tuple path (sampleBam), val (_)
  tuple val(sampleId), val(path),val(_), val(_)

  exec:
  String strVCFgz = sampleBam + '.vcf.gz'
  
  output:
  path ("*.vcf.gz")

  script:
  """
  samtools index $sampleBam &&
  bcftools mpileup -Ou -R $params.snps_ref -f $params.fast_align $sampleBam | bcftools call -c | bgzip > $strVCFgz
  """

}

process bedGraphToBigWig {
  queue = "$params.queue"

  //Docker Image
  container = "quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1"

  tag "7.2 - bedGraphToBigWig #$sampleId" 
  publishDir "$path/analysis/$params.outdir_bigwig/", mode : 'copy'

  input:
  path (h19Sizes)
  tuple val (otherFiles), path (bdgFiles)
  tuple val(sampleId), val(path),path(read1), path(read2)

  output:
  path ("*.bw")

  //for x in $path/analysis/$params.outdir_peak/*.bdg; do
  //     bedGraphToBigWig \$x $h19Sizes \$x.bw

  exec:
  bdgFile1 = bdgFiles.first()
  bdgFile2 = bdgFiles.last()
  bdgFile1_out = bdgFile1 + ".bw"
  bdgFile2_out = bdgFile2 + ".bw"

  script:
  """
  bedGraphToBigWig $bdgFile1 $h19Sizes $bdgFile1_out &&
  bedGraphToBigWig $bdgFile2 $h19Sizes $bdgFile2_out
  """

}

process lenght_fragment_dist_step1{
  queue = "$params.queue"

  //Singularity Image
  //container ='https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0'

  //Docker Image
  container ='quay.io/biocontainers/samtools:1.15.1--h1170115_0'

  tag "9.2 - Lenght Fragment Dist step 1 #$sampleId" 
  publishDir "$path/analysis/$params.outdir_bamtobed/", mode : 'copy'

  input:
  tuple path (sampleBam), path(_)
  tuple val(sampleId), val(path),path(_), path(_)

  output:
  path ('*.txt')

  exec:
  String strtxt = sampleId + '_fragment_lengths.txt'

  script:
  """
  samtools view $sampleBam | cut -f 9 | awk ' \$1 <= 1000 && \$1 > 0 ' > $strtxt
  """

}

process lenght_fragment_dist_step2{
  queue = "$params.queue"

  //Docker Image
  container ='pegi3s/r_data-analysis'

  tag "9.3 - Lenght Fragment Dist step 2 #$sampleId" 
  publishDir "$path/analysis/$params.outdir_bamtobed/", mode : 'copy'

  output:
  path ('*.png')

  input:
  path(fragLeng)
  tuple val(sampleId), val(path),path(_), path(_)

  exec:
  String strPNG = sampleId + '_fragDist.png' 

  //R -e 'install.packages(c("ggplot2", "reshape2"))'

  script:
  """
  Rscript $launchDir/frag_plotFragDist.R $fragLeng $strPNG $sampleId
  """

}

workflow {

    ch0 = Channel.fromPath(params.samples) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId,row.path, row.read1, row.read2) }
    chFile = Channel.fromPath("$launchDir/enrichment.sh")

    ch1 = mk_dir(ch0)
    fastqc(ch0,ch1)
    ch2 = trim(ch0,ch1)
    ch3 = align(ch2,ch0)
    ch4 = sort(ch3,ch0)
    lib_complex(ch4,ch0)
    ch5 = unique_sam(ch4,ch0)
    ch6 = dedup(ch5,ch0)
    ch7 = index_sam(ch6,ch0)
    ch8 = peak_bed_graph(ch6,ch7,ch0)
    ch9 = bam_to_bed(ch6,ch0)
    unique_frags(ch9,ch0)
    ch10 = fetch_chrom_sizes(ch0)
    snp_fingerprint(ch6,ch0)
    bedGraphToBigWig(ch10,ch8,ch0)
    enrichment(ch6,ch0,chFile)
    ch11=lenght_fragment_dist_step1(ch6,ch0)
    lenght_fragment_dist_step2(ch11,ch0)
}

