library(epiG)
library(ggplot2)

### Load fit
load('/home/martin/Documents/epiG-scripts/fit_test-sample_GNAS.RData')
nreads.threshold <- 150

cinfo <- chain.info(fit)

theme_set(theme_bw())

p <- ggplot(cinfo, aes(y=log(nreads), x=chain.id, fill = nreads < nreads.threshold)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(type = "qua", palette = 2) +
  geom_hline(yintercept = log(nreads.threshold), color = "black") +
  ylab("log(reads in chain)") +
  xlab("Haplotype chain") +
  theme(legend.position="none")

png('/home/martin/Documents/epiG-scripts/noise_reduction_plot.png', width = 700, height = 350)
print(p)
dev.off()
