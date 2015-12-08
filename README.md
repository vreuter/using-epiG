# Using epiG

## Fitting and plotting an epi-allelic model for a WGBS sample

See script: TODO

```R
library(epiG)
```

```R
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
```

```R
fit <- fetch.reads(fit)
fit <- fetch_ref(fit)
fit <- fetch_alt(fit)
```

### Making a noise reduction plot

### Plotting the methylation pattern (Full plot)

### Plotting the methylation pattern (Tetris plot)
