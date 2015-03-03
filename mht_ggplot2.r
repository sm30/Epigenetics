# colored mhtplot
pval.annot <- hm450[, c('CHR','MAPINFO')]
d <- cbind(pval.annot, regr.site.pv[[2]])
# d <- as.matrix(d, length(pval.annot), 3)
names(d) <- c('chr', 'pos', 'p')

#   ord <- with(d,order(chr,pos))
#   d <- d[ord,]
d$chr[d$chr=="X"] <- "23"
d$chr[d$chr=="Y"] <- "24"
d$chr <- as.numeric(d$chr)
d <- arrange(d, chr, pos)

oldpar <- par()
par(cex=0.4)
colors <- c(rep(c("blue","red", "green"),15),"red")
mhtplot(d,control=mht.control(colors=colors,gap=1000),pch=19,srt=0)

axis(2,cex.axis=1.5)
suggestiveline <- -log10(0.05/nrow(pval.annot))
#   genomewideline <- -log10(1.8E-06)
abline(h=suggestiveline, col="blue")
#   abline(h=genomewideline, col="green")
abline(h=0)

pngname <- paste0('splines_site_mhtplot_', bhv.vars[i], '.png')
#   dev.copy(png,file=pngname,width = 960, height = 640); dev.off()

par(cex=1)


library(gap)
load("4w.rda")
ord <- with(d,order(chr,pos))
d <- d[ord,]
pdf("figures/4w.pdf",height=9,width=10)
oldpar <- par()
par(cex=0.6)
colors <- c(rep(c("blue","red"),15),"red")
mhtplot(d,control=mht.control(colors=colors,gap=1000),pch=19,srt=0)


library(gap)
# png("mhtplot.png",height=10,width=16,units="cm",res=300)
data <- with(mhtdata,cbind(chr,pos,p))
glist <- c("IRS1","SPRY2","FTO","GRIK3","SNED1","HTR1A","MARCH3","WISP3","PPP1R3B",
           "RP1L1","FDFT1","SLC39A14","GFRA1","MC4R")
hdata <- subset(mhtdata,gene%in%glist)[c("chr","pos","p","gene")]
color <- rep(c("lightgray","gray"),11)
glen <- length(glist)
hcolor <- rep("red",glen)
par(las=2, xpd=TRUE, cex.axis=1.8, cex=0.4)
ops <- mht.control(colors=color,yline=1.5,xline=3)
hops <- hmht.control(data=hdata,colors=hcolor)
mhtplot(data,ops,hops,pch=19)
# dev.off()



library(ggplot2)
# the sample data here: http://de.iplantcollaborative.org/dl/d/DDEE604B-D690-41A9-86E1-25138CDC1D9E/gwas.zip
# 1. Set Working Directory and load ggplot2

setwd('~/Projects/R_lib/epigenetics/gwas_examples')
# 2. Import p-values from previous analysis

pVals <- read.csv('ricePvals.csv', row.names=1)
# 3. Calculate -log10 of pvals, which is what will be plotted

pValsLog10 <- -log10(pVals)
# 4. Import rice map

riceMap <- read.csv('sativas413map.csv')
# 5. Bind the two together. This will making plotting using ggplot2 functions easier

rcMpPvals <- cbind(riceMap, pValsLog10)

colnames(rcMpPvals) <- c('mrk', 'chr', 'pos', 'logp')
# 6. Position is defined as position on chromosome. Need to add largest position from last marker on previous chr so that markers are ordered correctly along x-axis

chrNum=12

for (i in 1:chrNum){ ndx <- which(rcMpPvals[, 2]==i)
  lstMrk <- max(rcMpPvals[ndx, 3])
  if (i < chrNum) ndx2 <- which(rcMpPvals[, 2]==i+1)
  if (i < chrNum) rcMpPvals[ndx2, 3] <- rcMpPvals[ndx2, 3] + lstMrk
}
# 7. Use this little loop to find the midposition of every chromosome so nice chr labels can be added to the plot

bpMidVec <- vector(length=chrNum)

for (i in 1:chrNum){ndx <- which(rcMpPvals[, 2]==i)
posSub <- rcMpPvals[ndx, 3]
bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub)
}
# 8. Use qplot function in ggplot2 to create Manhattan plot. You don't need to continually create new plot objects and add components to them. This is just done here for tutorial purposes. All of these components can be added in one line.

p <- ggplot(rcMpPvals)
(p2 <- p + geom_point(aes(x=pos, y=logp, size=3.5, colour=as.factor(chr)), alpha=1/3))
(p3 <- p2 + scale_color_manual(values=rep(c('black', 'dark green'), 6)))
(p4 <- p3 + theme_bw(base_size=15))
(p5 <- p4 + theme(legend.position='none'))
(p6 <- p5 + scale_x_continuous(labels=as.character(1:chrNum), breaks=bpMidVec))
(p7 <- p6 + geom_hline(y=4.08, linetype=1, col='red', lwd=1.5))
(p8 <- p7 + ggtitle('Rice Flowering Time') + xlab('') + ylab('-log10(P)'))