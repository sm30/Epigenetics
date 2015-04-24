dim(mvalues)
dim(hm450)
dim(pData)
##annot.full contains the full annotation information
##for all the CpG probes in the 450k array
annot.full <- annot
rm(annot)
table(rownames(pData) == pData$nestid)

#load the behavioral scores
epi<-read.csv("nestsr_iversen_2.csv", as.is=TRUE)
dim(epi)
length(unique(epi$nestid)) == nrow(epi)
rownames(epi) <-epi$nestid

#delete nestid for which the methylation age doens't correspond
epi<-epi[-grep('198', rownames(epi)),]

#merge patient info (race, age, etc) with behavioral scores
table(keep<-rownames(pData) %in% rownames(epi))
#only keep pData for nestid which is in epi
pData<-pData[keep,]
mvalues<-mvalues[,keep]
#rows will have the same order as pData
#epi will only have nestid which is in pData
epi <- epi[rownames(pData),]
colnames(epi)[colnames(epi) %in% colnames(pData)]
pData<-cbind(pData,epi[,colnames(epi) != "nest_id"])

#delete rows for which we do not have mother's ADHD data
keep <- (!is.na(pData$asrs_ADHD))
pData <- pData[keep,]
mvalues <- mvalues[,keep]
dim(mvalues)
rm(keep)

#only keep those probes where SNPs are less than half a percent
table(hm450$KeepHalfPct)
table(rownames(hm450)==rownames(mvalues))
mvalues<-mvalues[hm450$KeepHalfPct,]
hm450<-hm450[hm450$KeepHalfPct,]
table(rownames(hm450)==rownames(mvalues))

#remove NESTID 198 due to possible age issues

#colnames of mvalues are the barcode
#rownames of pData are the nestid (nest id and sample id are equal)
#for duplicate samples we used the higher quality mvalues
#which leads pData$barcode not corresponding to colnames of mvalues
#set rownames of pData to colnames of mvalues which is the barcode
rownames(pData) <- colnames(mvalues)

# <<create_methy450batch_object, cache = TRUE, echo = TRUE, include = TRUE>>=

##We create a methy450batch object with 15 slots
##slots: 11 annnotated regions, mvalues, detect_p, annotation, pData

##create a pDetect matrix
##initialize detectP values to one because we don't have the detection scores for the probes
detect_p <- matrix(1, nrow = dim(mvalues)[1], ncol = dim(mvalues)[2])
colnames(detect_p) <- colnames(mvalues)

mval.matrix = as.matrix(mvalues)
detect_p = as.matrix(detect_p)
annotation = as.matrix(hm450)
groupinfo = pData

# source code from the IMA package directly adapted for the analysis
source('../ima_custom.R')

# <<print_xmethy450, cache = TRUE, echo = TRUE, include = FALSE>>=
#check the length of the regions
for (slot in slotNames(x.methy450)){
  if(grepl('Ind', slot, ignore.case = TRUE)){
    print(paste0('x.methy450@', slot, '$SID'))
    print(length(eval(parse(text = paste0('x.methy450@', slot, '$SID')))))
  }
  else{
    paste0('x.methy450@', slot)
    print(dim(eval(parse(text = paste0('x.methy450@', slot)))))
  }
}

# <<summary_statistics, cache = TRUE, echo = TRUE, include = FALSE>>=
#create a list of matrices using the 11 region level summaries from the methy450 object
#to do: remove ind and change name
mval.region = list()

#loop through slots in methy450batch object
for (s in slotNames(x.methy450)){
  if(grepl('Ind', s, ignore.case = TRUE)){
    obj <- slot(x.methy450, s)
    #delete 'Ind' at the end of each region name
    #name <- substr(s, 1, nchar(s)-3)
    mval.region[[s]] <- indexregionfunc(indexlist= obj,
                                        beta=x.methy450@bmatrix, indexmethod="median")
  }
}

length(mval.region)
names(mval.region)