# options(width = 60)
paks <- c('splines', 'IMA', 'dplyr', 'xtable', 'doParallel', 'foreach', 'knitr', 'qvalue', 'gap',
          'car', 'robust', 'rmarkdown')
lapply(paks, library, character.only=T)
nofcl <- 3

#load the mvalues, pvalue, and hm450 matrix
#hm450 is the complete annotation matrix
setwd('~/Dropbox/Projects/epigenetics/data')

# subset mvalues and hm450 to size 5k
# load("MvalueAnalysisStructures.RData")
#mvalues<-mvalues[1:5000,]
#hm450<-hm450[1:5000,]
#save.image("MvalueAnalysisStructures_5000.RData")

# # load("~/epi450k/MvalueAnalysisStructures_5000.RData")
# load("~/epi450k/MvalueAnalysisStructures.RData")
# # load("/proj/design/il450k/nest13/R/MvalueAnalysisStructures.RData")
#
# source('../data_prep.R')
# source('../model_funs.R')
#
# save.image('~/epi450k/ready_to_analyze.RData')

load('~/epi450k/ready_to_analyze.RData')


# source('../model_spec_eating.R')
# source('../model_spec_adhd_pc.R')
source('../model_spec_adhd.R')


### lm with splines, site
#initialize the column we will use to store the methylation value one site at a time
# pData$mval <- NA

cl <- makeCluster(nofcl)
registerDoParallel(cl)
regr.site.pv <- foreach (y = bhv.vars, .verbose = TRUE, .packages = c('splines','robust')) %dopar% {lm.spline(mvalues, reg.vars, y)}
stopCluster(cl)

# regr.site.pv.old <- regr.site.pv

# par(resetPar())
opar <- par(no.readonly=TRUE)
bkppar <- opar
# opar <- bkppar
# bhv.vars.short <- bhv.vars
# bhv.vars.short <- c("BASC_EXT1", "BASC_INT", "BRF_GEC1")
bhv.vars.short <- c("BASC_APHY")
# rmarkdown::render('../splines_table_site.Rmd', "html_document")
knit2html('../splines_table_site.Rmd')
# knit2pdf('../splines_table_site.Rnw')

knit2html('../splines_table_site_check.Rmd')
#

source('../gene_lists.R')
# candidate.genes <- list(food, satiety)
# names(candidate.genes) <- c('food','satiety')
candidate.genes <- list(slotkin_genes)
names(candidate.genes) <- c('slotkin_genes')
candidate.genes <- list(adhd)
names(candidate.genes) <- c('adhd')
candidate.genes <- list(brain)
names(candidate.genes) <- c('brain')


knitr::knit2html('../gene_specific.Rmd')
# knitr::knit2pdf('../pvplot_gene_specific.Rnw')


### lm with splines, region
#parallelize over regions
cl <- makeCluster(nofcl)
registerDoParallel(cl)
regr.region.pv <-
  foreach (r = 1:11, .verbose = TRUE) %:%
  foreach (y = bhv.vars, .verbose = TRUE, .packages = 'splines') %dopar% {lm.spline(mval.region[[r]], reg.vars, y)}
stopCluster(cl)

length(regr.region.pv)
names(regr.region.pv) <- names(mval.region)

knitr::knit2pdf('../splines_table_region.Rnw')

knitr::knit2pdf('../pvplot_splines.Rnw')



### mlm site
#initialize the column we will use to store the methylation value one site at a time
reg.vars$mval <- NA

#idx <- seq(1, 331525, floor(331525/12))
#create chunks for each processor
m.len <- dim(mvalues)[[1]]
chunks <- nofcl
#number of elements in each chunk
div <- ceiling(m.len/chunks)
#chunks of the vector
x <- split(seq(m.len), ceiling(seq(m.len)/div))

cl <- makeCluster(nofcl)
registerDoParallel(cl)
mlm.site.pv <- foreach (n = 1:chunks, .verbose = TRUE, .packages = 'splines') %dopar% {mlm.cpgsite(mvalues, reg.vars, x[[n]])}
stopCluster(cl)

site.pv <- mlm.site.pv
knitr::knit2pdf('../mlm_table_site.Rnw')



### mlm region
#parallelize over regions
mlm.region.pv <- foreach (r = 1:11, .verbose = TRUE) %do% {mlm.regions(mval.region[[r]], reg.vars)}

# length(regr.site.pv)
# names(regr.site.pv) <- bhv.vars
# attributes(regr.site.pv)
#
# names(regr.region.pv) <- bhv.vars
# length(regr.region.pv)

region.pv <- mlm.region.pv
names(region.pv) <- names(mval.region)

knitr::knit2pdf('../mlm_table_region.Rnw')
