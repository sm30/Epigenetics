# options(width = 60)
paks <- c('splines', 'IMA', 'dplyr', 'xtable', 'doParallel', 'foreach', 'knitr', 'qvalue', 'gap', 'data', 'car')
lapply(paks, library, character.only=T)
nofcl <- 4

#load the mvalues, pvalue, and hm450 matrix
#hm450 is the complete annotation matrix
setwd('~/Projects/R_lib/epigenetics/data')

# subset mvalues and hm450 to size 5k
#load("MvalueAnalysisStructures.RData")
#mvalues<-mvalues[1:5000,]
#hm450<-hm450[1:5000,]
#save.image("MvalueAnalysisStructures_5000.RData")

# load("~/epi450k/MvalueAnalysisStructures_5000.RData")
load("~/epi450k/MvalueAnalysisStructures.RData")
# load("/proj/design/il450k/nest13/R/MvalueAnalysisStructures.RData")

source('../data_prep.R')
source('../model_funs.R')


# source('../model_spec_eating.R')
# source('../model_spec_adhd_pc.R')
source('../model_spec_adhd.R')


### lm with splines, site
#initialize the column we will use to store the methylation value one site at a time
pData$mval <- NA

cl <- makeCluster(nofcl)
registerDoParallel(cl)
regr.site.pv <- foreach (y = bhv.vars, .verbose = TRUE, .packages = 'splines') %dopar% {lm.spline(mvalues, reg.vars, y)}
stopCluster(cl)

# regr.site.pv.old <- regr.site.pv

# par(resetPar())
opar <- par(no.readonly=TRUE)
bkppar <- opar
# opar <- bkppar
# bhv.vars.short <- bhv.vars
bhv.vars.short <- c("BASC_EXT1", "BASC_INT1", "BRF_GEC1")
knit2html('../splines_table_site.Rmd')
# knit2pdf('../splines_table_site.Rnw')

knit2html('../splines_table_site_check.Rmd')
#

# top candidate genes by phenotype
food <- c('AGRP','BRS3','CARTPT','CCKAR','CCKBR','FYN','GALR2','GALR3','GCG','GHRL',
          'GHSR','HCRTR1','HCRTR2','HTR2C','LEP','MC4R','MCHR1','NMUR2','NPW','NPY',
          'PMCH','PPYR1','PYY','UBR3','FTO')
brain <- c('AFF2','ALK','ALX1','BPTF','CDK5R1','CEP290','CLN5','CNTN4','CTNS','DLX2',
           'DMBX1','DSCAML1','ECE2','EGR2','FOXG1','FOXP2','GLI2','GPR56','HESX1','LHX6',
           'MAP1S','MDGA1','MYO16','NCOA6','NDUFS4','NF1','NKX2-2','NNAT','OTX2','PBX1',
           'PBX3','PBX4','PCDH18','PHGDH','PITPNM1','POU6F1','PPT1','ROBO2','SHH','SHROOM2',
           'SHROOM4','SIX3','SLIT1','SMARCA1','TBR1','UBE3A','UNC5C','UTP3','VCX3A','ZIC1','ZIC2')
satiety <- c('rs9939609','FTO','rs2867125','TMEM18','rs571312','MC4R','rs10938397','GNPDA2',
             'rs10767664','BDN','rs2815752','NEGR','rs7359397','SH2B1','rs3817334','MTCH2',
             'rs29941','KCTD15','rs543874','SEC16B','rs987237','TFAP2B','rs7138803','FAIM2',
             'rs10150332','NRXN3','rs713586','POMC','rs12444979','GPRC5B','rs2241423','MAP2K5',
             'rs1514175','TNNI3K','rs10968576','LRRN6','rs887912','FANCL','rs13078807','CADM2',
             'rs1555543','PTBP2','rs206936','NUDT3','rs9568856','OLFM4','rs9299','HOXB5',
             'rs2112347','FLJ35779','rs3797580','rs4836133','ZNF608','rs6864049','rs4929949',
             'RPL27A','rs9300093','rs3810291','TMEM160','rs7250850','rs2890652','LRP1B',
             'rs9816226','ETV5','rs13107325','SLC39A8','rs4771122','MTIF3','rs11847697',
             'PRKD1','rs2287019','QPCTL')
satiety <- unique(satiety)
rsnumber <- grep("^rs[0-9]*$", satiety)
rs <- satiety[rsnumber]
satiety2 <- satiety[-rsnumber]
candidate.genes <- list(food, satiety)
names(candidate.genes) <- c('food','satiety')


# rs <- NULL
knitr::knit2pdf('../gene_specific.Rnw')
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
