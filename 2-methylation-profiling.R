###
### LIBS
###

library(epiG)
library(ggplot2)
library(plyr)

### Load sample
load('fitted-data-example/GNAS.RData')

###
### Meth profile plot
###

# Get position info
pos.cg <- locate_CG(g) # This only uses the ref

info <- rbind(position_info(g, pos.cg), position_info(g, pos.cg+1))

# Noise reduction
min.df <- 0.25
min.length <- 10000
intercept <- min.length
slope <- -min.length/min.df

cinfo <- chain_info(g)
cinfo <- subset(cinfo, intercept + slope*depth.fraction < length)
info <- subset(info, chain.id %in% cinfo$chain.id)

# Cut tails and head of chains

cinfo$start <- cinfo$start + 50
cinfo$end <- cinfo$end - 50

info <- ddply(info, .variables = c("chain.id"), function(x) {
  start <- cinfo$start[cinfo$chain.id==x$chain.id[1]]
  end <- cinfo$end[cinfo$chain.id==x$chain.id[1]]
  subset(x, position >= start & position <= end)
  })

info.meth <- subset(info, !is.na(methylated))

stat.df <- ddply(info.meth, .variables = c("chain.id"), function(x)
    data.frame(
      start = cinfo$start[cinfo$chain.id==x$chain.id[1]],
      end = cinfo$end[cinfo$chain.id==x$chain.id[1]],
      nreads = cinfo$nreads[cinfo$chain.id==x$chain.id[1]],
      position = median(x$position),
      meth = round(mean(x$methylated),1),
      nmeh = sum(x$methylated),
      ncpg = nrow(x),
      size = cinfo$depth.fraction[cinfo$chain.id==x$chain.id[1]]
      )
    )

rownames(stat.df) <- stat.df$chain.id

# Methylation
info.meth <- info.meth[order(info.meth$position, decreasing=TRUE),]
meth.df <- ddply(info.meth, .variables = "chain.id", function(x)
    if(nrow(x) < 3) {
      data.frame(
        start = x$position[1],
        end = x$position[nrow(x)],
        meth = mean(x$methylated),
        chain.meth = stat.df[as.character(x$chain.id[1]),"meth"],
        size = stat.df[as.character(x$chain.id[1]),"size"]
      )
    } else {
      data.frame(
        start = x$position[1:(nrow(x)-2)],
        end = x$position[3:nrow(x)],
        meth = sapply(1:(nrow(x)-2), function(i) mean(x$methylated[i:(i+2)])),
        chain.meth = stat.df[as.character(x$chain.id[1]),"meth"],
        size = stat.df[as.character(x$chain.id[1]),"size"]
      )
    }
  )

theme_set(theme_bw())

p1 <- ggplot(stat.df) +

geom_vline(xintercept = c(57414161, 57418015, 57424981, 57431470), color = "red", linetype = 3) +

geom_rect(data = meth.df, aes(ymax = chain.meth+0.04*size, ymin = chain.meth-0.04*size, xmin = start, xmax = end, fill = factor(meth)), color = NA) +
geom_rect(aes(ymin = meth-0.04*size, ymax = meth+0.04*size, xmin = start, xmax = end), size = 0.1, color = "black", fill = NA) +
scale_fill_brewer(type = "div", palette = 2) +

ylab("Methylation level") +

theme(legend.position="none",
strip.background = element_blank(),
panel.border = element_blank(),
axis.title.x = element_blank(),
axis.text.x = element_blank(),
axis.ticks.x= element_blank(),
plot.margin = unit(c(-.1,0,-.1,0), "cm")
)


### Save plot as png file
png('methylation-profile.png', units = "cm", width = 15, height = 6, res = 1000)
print(p1)
dev.off()
