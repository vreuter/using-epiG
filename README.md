# Using epiG

## Fitting and plotting an epi-allelic model for a WGBS sample

See script: TODO

```R
library(epiG)
```

```R
config <- auto_config(
  bam_file = bam.file,
  ref_file = 'hg19_rCRSchrm_chr20.fa',
  alt_file = 'dbsnp_135.hg19_chr20.fa',
  chr = chr,
  start = start,
  end = end,
  use_paired_reads = use_paired_reads,
  chunk_size = 150000,
  delta = 0.9)
```

```R
fit <- epiG(config)

fit <- fetch.reads(fit)
fit <- fetch_ref(fit)
fit <- fetch_alt(fit)
```

### Making a noise reduction plot
![alt tag](https://raw.github.com/vincent-dk/using-epiG/master/fig1.png)

```R
load('/home/martin/Documents/epiG-scripts/fit_test-sample_GNAS.RData')
nreads.threshold <- 150
```

```R
cinfo <- chain.info(fit)
```

```R
theme_set(theme_bw())

p <- ggplot(cinfo, aes(y=log(nreads), x=chain.id, fill = nreads < nreads.threshold)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(type = "qua", palette = 2) +
  geom_hline(yintercept = log(nreads.threshold), color = "black") +
  ylab("log(reads in chain)") +
  xlab("Haplotype chain") +
  theme(legend.position="none")
```

### Plotting the methylation pattern (Full plot)

### Plotting the methylation pattern (Tetris plot)
