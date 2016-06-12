library(epiG)

### site
chr <- "chr20"
start <- 57380000
end <- 57458000

# Create epiG configuration
config <- auto_config(
  bam_file = "raw-data-example/GNAS.bam",
  ref_file = "raw-data-example/hg19_rCRSchrm_chr20.fa",
  alt_file = "raw-data-example/dbsnp_135.hg19_chr20.fa",
  seq_type = "BSeq", # Std Bisulphite sequencing
  chr = chr,
  start = start,
  end = end)

# Run epiG
# This may take some time
g <- epiG(max_threads = 3, config = config)

# Load additional data into model
g <- fetch_reads(g) # Load reads
g <- fetch_ref(g) # Load reference genome
g <- fetch_alt(g) # Load alternative nucleotides

# Print a summary of the inferred epigenomic haplotypes
g

save(g, file = 'fitted-data-example/GNAS.RData')
