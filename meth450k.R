# options(width = 60)
paks <- c('splines', 'IMA', 'dplyr', 'xtable', 'doParallel', 'foreach', 'knitr')
lapply(paks, library, character.only=T)

nofcl <- 4
cl <- makeCluster(nofcl)
registerDoParallel(cl)

#load the mvalues, pvalue, and hm450 matrix
#hm450 is the complete annotation matrix
setwd('~/epi450k')

# subset mvalues and hm450 to size 5k
#load("MvalueAnalysisStructures.RData")
#mvalues<-mvalues[1:5000,]
#hm450<-hm450[1:5000,]
#save.image("MvalueAnalysisStructures_5000.RData")

# load("MvalueAnalysisStructures_5000.RData")
load("MvalueAnalysisStructures.RData")
# load("/proj/design/il450k/nest13/R/MvalueAnalysisStructures.RData")

source('~/Projects/R_lib/epigenetics/data_prep.R')
source('~/Projects/R_lib/epigenetics/model_funs.R')


library(dplyr)
source('~/Projects/R_lib/epigenetics/model_spec_eating.R')



### lm with splines, site
#initialize the column we will use to store the methylation value one site at a time
pData$mval <- NA

regr.site.pv <- foreach (y = bhv.vars, .verbose = TRUE, .packages = 'splines') %dopar% {lm.cpgsite(mvalues, reg.vars, y)}

knitr::knit2pdf('~/Projects/R_lib/epigenetics/splines_table_site.Rnw')



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

regr.site.pv <- foreach (n = 1:chunks, .verbose = TRUE, .packages = 'splines') %dopar% {mlm.cpgsite(mvalues, reg.vars, x[[n]])}
stopCluster(cl)

knitr::knit2pdf('~/Projects/R_lib/epigenetics/mlm_table_site.Rnw')



### mlm region
#parallelize over regions
regr.region.pv <- foreach (r = 1:11, .verbose = TRUE) %do% {mlm.regions(mval.region[[r]], reg.vars)}

length(regr.site.pv)
names(regr.site.pv) <- y_vars
attributes(regr.site.pv)

names(regr.region.pv) <- y_vars
length(regr.region.pv)
names(regr.region.pv) <- names(mval.region)

knitr::knit2pdf('~/Projects/R_lib/epigenetics/mlm_table_region.Rnw')