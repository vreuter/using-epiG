############################################################
###
### LIBS ect
###

library(epiG)

############################################################
###
### SETUP
###

### bam file
bam.file <- '/home/martin/Documents/epiG-scripts/GNAS.bam'
name <- "test-sample"
use_paired_reads <- TRUE

### site
chr <- "chr20"
start <- 57380000
end <- 57458000
site.name <- "GNAS"

############################################################
###
### RUN epiG
###

### epiG configuration
config <- auto_config(
  bam_file = bam.file,
  ref_file = '/home/martin/Documents/epiG-scripts/hg19_rCRSchrm_chr20.fa',
  alt_file = '/home/martin/Documents/epiG-scripts/dbsnp_135.hg19_chr20.fa',
  chr = chr,
  start = start,
  end = end,
  use_paired_reads = use_paired_reads,
  chunk_size = 150000,
  delta = 0.9)

fit <- epiG(config)

fit <- fetch.reads(fit)
fit <- fetch_ref(fit)
fit <- fetch_alt(fit)

### SAVE

save(fit, file = paste('home/martin/Documents/epiG-scripts/fit_', name,'_', site.name,'.RData', sep=''))
