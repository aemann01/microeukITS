######################
#taxonomy barchart
######################
library(reshape2)
library(ggplot2)
dat <- read.table("microeuk_simplified.txt", sep="\t", header=T)
datmelt <- melt(dat)
pdf("taxonomy_barchart.pdf")
ggplot(datmelt, aes(fill=datmelt$taxonomy, x=datmelt$variable, y=datmelt$value)) + geom_bar(stat="identity", position="fill") + theme_classic()
dev.off()