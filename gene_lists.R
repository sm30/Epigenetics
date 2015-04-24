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
# rs <- NULL

# http://adhd.psych.ac.cn/topGene.do
adhd <- ("SLC6A3
DRD4
COMT
SLC6A4
DRD5
SNAP25
DBH
MAOA
BDNF
SLC6A2
HTR1B
ADRA2A
HTR2A
TPH2
DRD2
DRD1
DRD3
CHRNA4
ADRA2C
TH
MAOB
TPH1
DDC
HTR2C")

adhd <- strsplit(adhd, "\n")
adhd <- unlist(adhd)


slotkin <- readLines('slotkin.txt')
slotkin <- gsub("  ", " ", slotkin)
slotkin_genes <- slotkin[grep("^[a-zA-Z1-9]+[ \t][0-9]", slotkin)]
slotkin_genes <- strsplit(slotkin_genes, " ")
slotkin_genes <- lapply(slotkin_genes, function(x) x[[1]])
slotkin_genes <- unlist(slotkin_genes)
slotkin_genes <- toupper(slotkin_genes)
