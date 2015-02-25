### lm with splines, site or region level
lm.spline <- function(m, reg.vars, y) {
  fmla<-paste('reg.vars$', y, '~reg.vars$',
              paste(cntrl, collapse = "+reg.vars$"),
              "+bs(reg.vars$mval, degree = 3)", sep = "")

  n.probes = dim(m)[1]
  pval <- vector(mode = "numeric", length = n.probes)
  for (i in 1:n.probes) {
    reg.vars$mval <- m[i,]
    lm.out <- lm(as.formula(fmla), data = reg.vars)
    pval[i] <- anova(lm.out)["bs(reg.vars$mval, degree = 3)", "Pr(>F)"]
  }
  return(pval)
}



### mlm with splines, site level
mlm.cpgsite <- function(m, reg.vars, vec) {
  fmla <- paste('cbind(',
                paste0('reg.vars$',bhv.vars, collapse=','), ')~reg.vars$',
                paste(cntrl, collapse = "+reg.vars$"),
                "+bs(reg.vars$mval, degree = 3)", sep = "")
  pval <- vector(mode = "numeric", length = dim(m)[1])
  for (i in vec[1]:vec[length(vec)]) {
    reg.vars$mval <- m[i,]
    lm.out <- lm(as.formula(fmla), data = reg.vars)
    pval[i] <- anova(lm.out)["bs(reg.vars$mval, degree = 3)", "Pr(>F)"]
  }
  return(pval)
}


### mlm with splines, region level

# We will get a region level summary for each gene.
# the parallelized structure loops over the # behavioral scores
# the function has 2 loops:
# over each of the 11 regions and the methylation values within each of those regions

# behav_score <- pData[y_vars]
mlm.regions <- function(m, reg.vars) {
  fmla <- paste('cbind(',
                paste0('reg.vars$',bhv.vars, collapse=','), ')~reg.vars$',
                paste(cntrl, collapse = "+reg.vars$"),
                "+bs(reg.vars$mval, degree = 3)", sep = "")
  #regress against the methylation values at each region
  #do this for each of the 11 region categories
  n.probes <- dim(m)[1]
  pval <- vector(mode= "numeric", length = n.probes)
  for (i in 1:n.probes) {
    reg.vars$mval <- m[i,]
    lm.out <- lm(as.formula(fmla), data = reg.vars)
    pval[i] <- anova(lm.out)["bs(reg.vars$mval, degree = 3)", "Pr(>F)"]
  }
  return(pval)
}



pvplot <- function(pvs,titl){
  nlogpv<-sort(-log10(pvs[!is.na(pvs)]))
  n.use<-length(nlogpv)
  expctd<-((1:n.use)-0.5)/n.use
  expctd<-sort(-log10(expctd))
  ##  lambda<-summary(lm(nlogpv[1:2000]~expctd[1:2000]))$coefficients[2,1]
  par(pty="s")
  thin<-seq(1,n.use,by=1)
  if (n.use>10000) thin<-c(seq(1,floor(0.5*n.use),by=100),seq(floor(0.5*n.use),floor(0.75*n.use),by=10),
                           seq(floor(0.75*n.use),floor(0.95*n.use),by=2),seq(floor(0.95*n.use),n.use,by=1))
  if (n.use>100000) thin<-c(seq(1,floor(0.5*n.use),by=1000),seq(floor(0.5*n.use),floor(0.75*n.use),by=100),
                            seq(floor(0.75*n.use),floor(0.95*n.use),by=10),seq(floor(0.95*n.use),n.use,by=1))
  plot(nlogpv[thin],expctd[thin],pch=16,main=titl,las=1,
       cex=1.5,xlab="observed -log10(pv)",ylab="expected -log10(pv)")
  ## text(x=0.5,y=3.5,labels=paste("lambda =",round(lambda,3),sep=" "))
  legend(x="topleft",inset=0.05,
         legend=c("95% Pointwise Interval Under Null","Expected Under Null"),
         lty=c(1,1),col=c(3,2),lwd=c(4,4),cex=0.75)
  abline(a=0,b=1,col=2,lwd=2)
  q0.025allc<-numeric(n.use)
  q0.975allc<-numeric(n.use)
  ##q0.025ballc<-numeric(n.use)
  ##q0.975ballc<-numeric(n.use)
  ##bonf.adj<-n.use
  for (i in 1:n.use){
    q0.025allc[i]<-qbeta(0.025,shape1=i,shape2=(n.use-i+1))
    q0.975allc[i]<-qbeta((1-0.025),shape1=i,shape2=(n.use-i+1))
    ##q0.025ballc[i]<-qbeta(0.025/bonf.adj,shape1=i,shape2=(n.use-i+1))
    ##q0.975ballc[i]<-qbeta((1-0.025/bonf.adj),shape1=i,shape2=(n.use-i+1))
  }
  lines(sort(-log10(q0.025allc))[thin],expctd[thin],lwd=2,col=3)
  lines(sort(-log10(q0.975allc))[thin],expctd[thin],lwd=2,col=3)
  ##lines(sort(-log10(q0.025ballc)),expctd,lwd=3,col=4)
  ##lines(sort(-log10(q0.975ballc)),expctd,lwd=3,col=4)
  return(NULL)
}

uniqueGenes <- function(vec) {
  tmp <- lapply(strsplit(vec,";"),unique)
  unlist(lapply(tmp, paste0, collapse=";")) # cat doesn't work
}
