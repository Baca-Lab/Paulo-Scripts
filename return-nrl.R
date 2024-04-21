## original---Human SMARCA5 is continuously required to maintain nucleosome spacing DOI: 10.1016/j.molcel.2022.12.018
## modified code from the original--- Author: Paulo Cordeiro (paulocordeiro@gmail.com), April 8, 2024

# Suppress warning messages
options(warn = -1)

calcNuclFragLen <- function(nNucl, nuclLen = 147, linkLenMin = 10, linkLenMax = 80)
  nNucl*nuclLen + nNucl*c(linkLenMin, linkLenMax)

nrl_bed <- function(bed_file_path) {
  library(readr)
  library(GenomicRanges)
  library(AnnotationDbi)
  library(dplyr)
  library(rtracklayer)
  
  # Load the BED file
  bed_data <- import(bed_file_path, format = "bed")
  # Calculate the size of each fragment
  isize <- width(bed_data)
  
  ## ... loess fit smoothing parameter 1
  ##     a first loess fit with this parameter will be fit to the fragment length distribution
  ##     in log space, in order to capture the (low-frequency) decrease of fragment frequencies
  ##     with increasing length
  lspan <- 0.35
  
  ## ... loess fit smoothing parameter 2
  ##     a second loess fit with this parameter will be fit to the residuals of the first fit
  ##     in order to obtain a smooth curve from which the nucleosome fragment lengths will be obtained
  rspan <- 0.1
  
  ##     Remark: if the fragment length distribution differs from the one in the published
  ##     ATAC data, the lspan and rspan parameters may have to be adjusted
  
  ## ... minimal fragment length included in loess fit
  ##     (shorter fragments are mostly originating from transcription factors and will not be
  ##      informative for the nucleosomal fragment length estimation)
  minxFit <- 51
  
  ###EXAMPLE OF CALCULATED NRL using previously published loess model 
  ## ... define parameters to read fragment size (isize)
  ##     only from first reads (R1, to avoid redundancy)
  ##     of uniquely mapped proper pairs (mapq >= 11L)
  ##     to autosomes (chrs)
  
  ## ... tabulate absolute insert sizes
  ##     mx is the maximum observed fragment length
  ##     cnt[i] is the number of fragments of length i
  ##     freq[i] is the fraction of fragments of length i
  mx <- max(abs(isize))
  cnt <- tabulate(abs(isize))
  freq <- cnt / sum(cnt)
  
  nuclBinsL <- list(mono = calcNuclFragLen(1),
                    di   = calcNuclFragLen(2),
                    tri  = calcNuclFragLen(3))
  
  xs <- 1:mx
  keep <- freq > 0 & xs >= minxFit
  
  #### Memory Sensitive
  fit1 <- loess(log(freq) ~ log(xs), subset = keep, span = lspan, surface="direct")
  pfreq <- exp(predict(fit1, log(xs)))
  rfreq <- freq - pfreq
  
  fit2 <- loess(rfreq ~ xs, subset = keep, span = rspan)
  
  avgFragLength <- sapply(nuclBinsL, function(xr) {
    xrs <- seq.int(xr[1], xr[2])
    xrs[which.max(predict(fit2, xrs))]
  })
  
  estNRL <- mean(avgFragLength / seq_along(avgFragLength))
  return(estNRL)
}


# Check if a command-line argument is provided
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  # Get the string from the command-line argument
  input_file <- commandArgs(trailingOnly = TRUE)
  # Call the function to count characters
  nrl_bed(input_file)
} else {
  cat("Usage: Rscript return-nrl [string]\n")
}


