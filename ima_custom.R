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