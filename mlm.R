# \documentclass[12pt]{article}
#
# \usepackage{times}
# \usepackage[cm]{fullpage}
#
# \setlength{\parskip}{6pt}
# \setlength{\parindent}{0pt}
#
# \begin{document}
#
# <<get_libraries_and_set_options, cache = FALSE, echo = TRUE, include = FALSE>>=
options(width = 60)
library(splines)
library(IMA)
library(xtable)
library(doParallel)
library(foreach)
cl <- makeCluster(4)
registerDoParallel(cl)

# @
#
# <<load_data, cache = TRUE, echo = TRUE, include = TRUE>>=
#load the mvalues, pvalue, and hm450 matrix
#hm450 is the complete annotation matrix
setwd('~/epi450k')

# subset mvalues and hm450 to size 5k
#load("MvalueAnalysisStructures.RData")
#mvalues<-mvalues[1:5000,]
#hm450<-hm450[1:5000,]
#save.image("MvalueAnalysisStructures_5000.RData")

load("MvalueAnalysisStructures_5000.RData")
# load("/proj/design/il450k/nest13/R/MvalueAnalysisStructures.RData")

dim(mvalues)
dim(hm450)
dim(pData)
##annot.full contains the full annotation information
##for all the CpG probes in the 450k array
annot.full <- annot
rm(annot)
table(rownames(pData) == pData$nestid)

#load the behavioral scores
epi<-read.csv("nestsr_iversen_24OCT2014.csv", as.is=TRUE)
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
#wpi will only have nestid which is in pData
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

#create categorical variables for mother's ADHD, education, parity, and pre-pregnancy BMI
pData$asrs_ADHD_2cat <- factor(pData$asrs_ADHD, labels = 0:1, levels = 0:1)
pData$education_3cat<- factor(pData$education_3cat, labels = 1:3, levels = 1:3)
pData$parity_3cat <- factor(pData$parity_3cat, labels = 0:2, levels = 0:2)
pData$prePregBMIthreeLev <- rep(NA, nrow(pData))
pData$prePregBMIthreeLev[(!is.na(pData$BMI_LMP_kgm2)) & (pData$BMI_LMP_kgm2 < 30)] <- 0
pData$prePregBMIthreeLev[(!is.na(pData$BMI_LMP_kgm2)) & (pData$BMI_LMP_kgm2 >= 30) & (pData$BMI_LMP_kgm2 < 35)] <- 1
pData$prePregBMIthreeLev[(!is.na(pData$BMI_LMP_kgm2)) & (pData$BMI_LMP_kgm2>=35)] <- 2
pData$prePregBMIthreeLev <-factor(pData$prePregBMIthreeLev, levels = 0:2, labels = c("Lt30", "30toLt35", "Ge35"))
# @
#
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

cat(".......\nSplit the annotation file to 11
    annotated region categories\n.......\n\n", fill = TRUE)
    annot = annotation
    name = "UCSC_RefGene_Name"
    cpGsite = as.character(annot[, 1])
    genelist = strsplit(as.character(annot[, name]), ";")
    genelist[which(genelist == "character(0)")] = "NA"
    #make a list with Refgene group for each site, sometimes more than one
    #input NA for those without value
    name = "UCSC_RefGene_Group"
    refgene = strsplit(as.character(annot[, name]), ";")
    #make a list with Refgene group for each site, sometimes more than one
    refgene[which(refgene == "character(0)")] = "NA"
    listlength = lapply(refgene, length)
    listlength[listlength == 0] = 1
    col0 = rep(1:nrow(annot), listlength)
    col1 = rep(cpGsite, listlength)
    col2 = unlist(genelist)
    col3 = unlist(refgene)
    col4 = rep(as.character(annot[, "Relation_to_UCSC_CpG_Island"]),
        listlength)
    col5 = rep(as.character(annot[, "UCSC_CpG_Islands_Name"]),
        listlength)
##we rep the col0, 1, 4, 5 according to how many sites each CpG probe is linked to
##col0 are the indices, col1 are the Cpg probe names, col2 are the names of the
##associated genes, col3 are the name of the associated gene regions (exon, etc),
##col4 are CpG site region (i.e. shelf, shore, etc.),
##col5 are the chromosomal location of the CpG islands
    splitToRegionlist = function(grepname = c("TSS1500", "TSS200",
        "5'UTR", "1stExon", "Gene Body", "3'UTR")) {
        index = col3 == grepname
        col1sub = col1[index]
        col2sub = col2[index]
        temp = split(col1sub, col2sub)
        returnSID = lapply(temp, unique)
        col0sub = col0[index]
        temp = split(col0sub, col2sub)
        returnPID = lapply(temp, unique)
        return(Ind = list(SID = returnSID, PID = returnPID))
    }
##First list of summary stats indexed by associated genes
    TSS1500Ind = splitToRegionlist(grepname = "TSS1500")
    TSS200Ind = splitToRegionlist(grepname = "TSS200")
    UTR5Ind = splitToRegionlist(grepname = "5'UTR")
    EXON1Ind = splitToRegionlist(grepname = "1stExon")
    GENEBODYInd = splitToRegionlist(grepname = "Body")
    UTR3Ind = splitToRegionlist(grepname = "3'UTR")

    cat("TSS1500 region contains:", length(TSS1500Ind$SID),
        "UCSC REFGENE region \nTSS200 region contains:",
        length(TSS200Ind$SID), "UCSC REFGENE region\n5'UTR region contains:",
        length(UTR5Ind$SID), "UCSC REFGENE region\n1st Exon region contains:",
        length(EXON1Ind$SID), "UCSC REFGENE region\nGene body region contains:",
        length(GENEBODYInd$SID), "UCSC REFGENE region\n3'UTR region contains:",
        length(UTR3Ind$SID), "UCSC REFGENE region\n", fill = TRUE)

##Second list of summary stats indexed by chromosomal regions
    splitToRegionlist2 = function(grepname = c("Island", "N_Shore",
        "S_Shore", "N_Shelf", "S_Shelf")) {
        index = col4 == grepname
        col1sub = col1[index]
        col5sub = col5[index]
        temp = split(col1sub, col5sub)
        returnSID = lapply(temp, unique)
        col0sub = col0[index]
        temp = split(col0sub, col5sub)
        returnPID = lapply(temp, unique)
        return(Ind = list(SID = returnSID, PID = returnPID))
    }

    ISLANDInd = splitToRegionlist2(grepname = "Island")
    NSHOREInd = splitToRegionlist2(grepname = "N_Shore")
    SSHOREInd = splitToRegionlist2(grepname = "S_Shore")
    NSHELFInd = splitToRegionlist2(grepname = "N_Shelf")
    SSHELFInd = splitToRegionlist2(grepname = "S_Shelf")

    cat("Island region contains:", length(ISLANDInd$SID),
        "UCSC CPG ISLAND region\nN_Shore region contains",
        length(NSHOREInd$SID), "UCSC CPG ISLAND region\nS_Shore region contains",
        length(SSHOREInd$SID), "UCSC CPG ISLAND region\nN_Shelf region contains",
        length(NSHELFInd$SID), "UCSC CPG ISLAND region\nS_Shelf region contains",
        length(SSHELFInd$SID), "UCSC CPG ISLAND region\n", fill = TRUE)
    setClass("methy450batch", representation(bmatrix = "matrix",
        annot = "matrix", detectP = "matrix", groupinfo = "data.frame",
        TSS1500Ind = "list", TSS200Ind = "list", UTR5Ind = "list",
        EXON1Ind = "list", GENEBODYInd = "list", UTR3Ind = "list",
        ISLANDInd = "list", NSHOREInd = "list", SSHOREInd = "list",
        NSHELFInd = "list", SSHELFInd = "list"), where = topenv(parent.frame()))
    x.methy450 = new("methy450batch", bmatrix = as.matrix(mval.matrix),
        annot = as.matrix(annotation), detectP = as.matrix(detect_p),
        groupinfo = groupinfo, TSS1500Ind = TSS1500Ind, TSS200Ind = TSS200Ind,
        UTR5Ind = UTR5Ind, EXON1Ind = EXON1Ind, GENEBODYInd = GENEBODYInd,
        UTR3Ind = UTR3Ind, ISLANDInd = ISLANDInd, NSHOREInd = NSHOREInd,
        SSHOREInd = SSHOREInd, NSHELFInd = NSHELFInd, SSHELFInd = SSHELFInd)
    cat("\nA methy450batch class is created and the slotNames are: \n",
        slotNames(x.methy450), "\n", fill = TRUE)

# @
#
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
# @
#
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
# @
#
# <<set_vals_for_reg, cache = TRUE, echo = FALSE, include = TRUE, eval = TRUE>>=
dim(mvalues)
#create BASC_AP_HY composite score
#pData$BASC_AP_HY <- rowMeans(cbind(pData$BASC_AP, pData$BASC_HY))
#behavioral scores
#bhv.vars <- c("BASC_AP_HY" , "BASC_AP", "BASC_HY", "BRF_ISCI", "BRF_FI",
#    "BRF_EMI", "BRF_GEC", "BRF_IN", "BRF_SF", "BRF_PO", "BRF_WM", "BRF_EC")

#bhv.vars <- c("BASC_AP", "BASC_HY1", "BRF_IN1", "BRF_SF1", "BRF_PO1", "BRF_WM1", "BRF_EC1")
bhv.vars <- c("CEBQ_SRSE", "CEBQ_FR", "CEBQ_EF", "CEBQ_EOE")

#pData$BASC_HY1 <- log(10 + pData$BASC_HY)
#pData$BRF_SF1 <- log(pData$BRF_SF - 5)
#pData$BRF_WM1 <- log(pData$BRF_WM - 12)
#pData$BRF_IN1 <- log(pData$BRF_IN - 5)
#pData$BRF_EC1 <- log(pData$BRF_EC - 5)
#pData$BRF_PO1 <- log(pData$BRF_PO - 5)

reg.vars <- data.frame(pData$CEBQ_SRSE, pData$CEBQ_FR, pData$CEBQ_EF, pData$CEBQ_EOE,
    pData$age_mo_SR_Surv, pData$sex, pData$birthwt_kg,
    pData$GestAge_weeks, pData$education_3cat, pData$race_final_n, pData$parity_3cat,
    pData$mom_age_delv, pData$prePregBMIthreeLev, pData$asrs_ADHD_2cat)

colnames(reg.vars) <- c("CEBQ_SRSE", "CEBQ_FR", "CEBQ_EF", "CEBQ_EOE",
    "age_mo_SR_Surv", "sex", "birthwt_kg",
    "GestAge_weeks", "education_3cat", "race_final_n", "parity_3cat",
    "mom_age_delv", "prePregBMIthreeLev", "asrs_ADHD_2cat")
# @
#
# <<mregr_func, cache = TRUE, echo = FALSE, include = TRUE, eval = TRUE>>=
lm.cpgsite <- function(m, reg.vars, vec) {
    cntrl <- c("~reg.vars$age_mo_SR_Surv", "sex", "birthwt_kg",
    "GestAge_weeks", "education_3cat", "race_final_n", "parity_3cat",
    "mom_age_delv", "prePregBMIthreeLev", "asrs_ADHD_2cat")

    fmla <- paste('cbind(', 'reg.vars$CEBQ_SRSE, ', 'reg.vars$CEBQ_FR, ', 'reg.vars$CEBQ_EF, ',
    'reg.vars$CEBQ_EOE) ',
        paste(cntrl, collapse = "+reg.vars$"), "+bs(reg.vars$mval, degree = 3)", sep = "")

    pval <- vector(mode = "numeric", length = dim(m)[1])
    for (i in vec[1]:vec[length(vec)]) {
        reg.vars$mval <- m[i,]
        lm.out <- lm(as.formula(fmla), data = reg.vars)
        pval[i] <- anova(lm.out)["bs(reg.vars$mval, degree = 3)", "Pr(>F)"]
    }
    # lm.out <- lm(pdat$BASC_AP_HY~pdat$age_mo_SR_Surv+pdat$sex+pdat$birthwt_kg+
    #      pdat$GestAge_weeks+pdat$education_3cat+pdat$Race3+pdat$parity_3cat+
    #      pdat$mom_age_delv+pdat$prePregBMIthreeLev+pdat$asrs_ADHD_2cat+pdat$mval, data = pdat)

   return(pval)

}
# @
#
# <<site_level_mregression, cache = TRUE, echo = TRUE, eval = TRUE, include = FALSE>>=
#initialize the column we will use to store the methylation value one site at a time
reg.vars$mval <- NA

#idx <- seq(1, 331525, floor(331525/12))
#create chunks for each processor
m.len <- dim(mvalues)[[1]]
chunks <- 12
#number of elements in each chunk
div <- ceiling(m.len/chunks)
#chunks of the vector
x <- split(seq(m.len), ceiling(seq(m.len)/div))

regr.site.pv <- foreach (n = 1:chunks, .verbose = TRUE, .packages = 'splines') %dopar% {lm.cpgsite(mvalues, reg.vars, x[[n]])}
stopCluster(cl)

# @
#
# \end{document}
