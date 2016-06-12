### Libs
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
  geom_rect(aes(ymin = depth.fraction-0.005, ymax = depth.fraction+0.005, xmin = start, xmax = end, fill = intercept + slope*depth.fraction < length), color = "black", size = 0.01) +
  ylab("Depth fraction") +
  xlab("Position") +
  scale_fill_brewer(type = "qua", palette = 2) +

  geom_vline(xintercept = c(57414161,57418015, 57424981, 57431470), color = "red", linetype = 3) +

  theme(legend.position="none",
  strip.background = element_blank(),
  panel.border = element_blank(),
  axis.title.x = element_blank(),
#  axis.text.x = element_blank(),
  axis.ticks.x= element_blank()
  )

### Save plot as png file
png('haplotype-chains-plot.png', units = "cm", width = 15, height = 6, res = 1000)
print(p1)
dev.off()
