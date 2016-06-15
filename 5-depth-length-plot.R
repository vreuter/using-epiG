library(epiG)
library(ggplot2)

### Noise reduction
min.df <- 0.25
min.length <- 10000
intercept <- min.length
slope <- -min.length/min.df

### Load sample
load('fitted-data-example/GNAS.RData')

### Exstract chain info
cinfo <- chain_info(g)

### Make plot
theme_set(theme_bw())

p1 <- ggplot(cinfo) +
  geom_point(aes(x = depth.fraction, y = length, size = nreads, color = intercept + slope*depth.fraction < length)) +
  geom_abline(aes(slope = slope, intercept = intercept), alpha = .5, linetype = 2) +
  ylab("Chain length (bp)") +
  xlab("Depth fraction") +
  scale_color_brewer(type = "qua", palette = 2) +
  theme(legend.position="none")

### Save plot as png file
ggsave("depth-length-plot.png", units = "cm", width = 7, height = 3, dpi = 600, scale = 2)
