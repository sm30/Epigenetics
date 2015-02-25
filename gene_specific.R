food <- c('AGRP','BRS3','CARTPT','CCKAR','CCKBR','FYN','GALR2','GALR3','GCG','GHRL',
          'GHSR','HCRTR1','HCRTR2','HTR2C','LEP','MC4R','MCHR1','NMUR2','NPW','NPY',
          'PMCH','PPYR1','PYY','UBR3')

# names(regr.site.pv[[i]]) <- rownames(hm450)
pval.sort <- sort(regr.site.pv[[i]])
# pval.sig <- pval.sort[1:20]
pval.annot <- hm450[names(pval.sort), c(17, 18, 20, 22)]
assoc.genes <- data.frame(names(pval.sort), pval.sort, pval.annot)

food.genes <- NULL
for (gene in food) {
  pattern <- paste0(";", gene, ";|^", gene, ";|;", gene, "$")
  food.genes <- rbind(food.genes, assoc.genes[grep(pattern,assoc.genes$UCSC_RefGene_Name),])
}
food.genes <- arrange(food.genes, pval.sort)

behav.caption <- gsub("_", "-", bhv.vars[i])