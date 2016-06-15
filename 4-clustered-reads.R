###
### LIBS
###

library(epiG)
library(ggplot2)
library(grid)
library(plyr)

### Load data
start <- 57414000
end <- 57414000 + 800

load('fitted-data-example/GNAS.RData')
g <- subregion(g, start = start, end = end)

####
#### Get info
####

pos <-start(g):end(g)

### Information about reads
rinfo <- read_info(g, inc.symbols = TRUE)
rinfo <- subset(rinfo, position %in% pos)

### Estimated genotype, strand and methylation status
info <- position_info(g, pos)

### Remove small chains
cinfo <- chain_info(g)
cinfo <- subset(cinfo, nreads > 10)

info <- subset(info, chain.id %in% cinfo$chain.id)
rinfo <- subset(rinfo, chain.id %in% cinfo$chain.id)


###
### Plot data and estimated model
###

base_consitence <- function(base, strand, ref)
	(ref == base) |
	(ref == "C" & base == "T" & strand == "fwd") |
	(ref == "G" & base == "A" & strand == "rev")

is_converted <- function(base, strand, ref)  {

	tmp <- rep(NA, length(ref))
	idx <- (strand == "fwd" & ref == "C") | (strand == "rev" & ref == "G")
	tmp[idx] <-	(strand[idx] == "fwd" & base[idx] == "T") | (strand[idx] == "rev" & base[idx] == "A")

	return(tmp)
}

rinfo$name <- factor(rinfo$name)
s <- ddply(rinfo, .variables = c("name"), function(x) data.frame(start = min(x$position)))
rownames(s) <- s$name

# reads info
reads.pd <- data.frame(
	position = rinfo$position,
	read = rinfo$read.id,
	symbol = rinfo$symbol,
	y = s[rinfo$name,]$start,
	chain.id = rinfo$chain.id,
	strand = rinfo$strand,
	meth = is_converted(rinfo$symbol, rinfo$strand, rinfo$ref),
	base_error = ! base_consitence(rinfo$symbol, rinfo$strand, rinfo$ref)
	)

# order y
y.order <- sort(daply(reads.pd, .variables=c("y"), function(x) min(x$position)))
reads.pd$y <- sapply(reads.pd$y, function(i) which(y.order==i)[1])

# Fix y
reads.pd <- ddply(reads.pd, .variables=c("chain.id", "strand"), function(x) {
		x$y <- as.integer(factor(x$y));
		return(x)
		})

reads.pd$y <- (reads.pd$strand == "rev")*(-reads.pd$y) +
                  (reads.pd$strand == "fwd")*reads.pd$y

# genotype info
geno.pd <- data.frame(
	position = info$position,
	symbol = info$genotype,
	ref = info$ref,
	y = 0,
	chain.id = info$chain.id,
	meth = info$methylated,
	base_error =  info$genotype != info$ref)

### Plot

theme_set(theme_bw())

#Screen readable
pos <- locate_CG(g)


p <- ggplot(reads.pd, aes(x=position, y=y+2*(y > 0)-2*(y< 0))) +
	geom_line(aes(group = read), color = "yellowgreen", size = 0.5, alpha = 0.5) +
	geom_line(aes(group = chain.id, y = 0), color = "black", size = 0.5, alpha = 0.5) +

	geom_point(data = subset(reads.pd, base_error), color = "red", size = 0.5) +
	geom_point(data = subset(reads.pd, meth & position %in% c(pos, pos+1)), color = "black", size = 0.5) +

	geom_point(data = subset(geno.pd, meth & position %in% pos), aes(y = 0, x = position +.5), color = "black", size = 2) +
	geom_point(data = subset(geno.pd, ! meth & position %in% pos), aes(y = 0, x = position +.5), color = "black",  fill = "white", shape=21, size = 2) +

	facet_grid(chain.id ~ ., scales = "free_y", space = "free_y") +

	theme(legend.position="none", strip.background = element_blank(),
	panel.grid= element_blank(), panel.border = element_blank(), strip.text = element_blank(),
	axis.text.y = element_blank(), axis.ticks.y= element_blank(), axis.title.y = element_blank(),
	#axis.text.x = element_blank(), axis.ticks.x= element_blank(), axis.title.x = element_blank(),
	strip.text.y = element_text(colour = "black", angle = 0, size = 5, hjust = 0),
	panel.margin = unit(0.8, "cm"),
	plot.margin = unit(c(.5,-.75,0,-.75), "cm"))

  ### Save plot as png file
  ggsave("clustered-reads.png", units = "cm", width = 15, height = 6, dpi = 600, scale = 2)
