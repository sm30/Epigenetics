### lm with splines, site level
lm.cpgsite <- function(m, reg.vars, y) {
  fmla<-paste('reg.vars$', y, '~reg.vars$',
              paste(cntrl, collapse = "+reg.vars$"), "+bs(reg.vars$mval, degree = 3)", sep = "")

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
                paste0('reg.vars$',bhv.vars, sep=','), ')~reg.vars$',
#               'reg.vars$CEBQ_SRSE, ', 'reg.vars$CEBQ_FR, ', 'reg.vars$CEBQ_EF, ', 'reg.vars$CEBQ_EOE) ',
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
  fmla <- paste('cbind(', 'reg.vars$CEBQ_SRSE, ', 'reg.vars$CEBQ_FR, ', 'reg.vars$CEBQ_EF, ',
                'reg.vars$CEBQ_EOE) ',
                paste(cntrl, collapse = "+reg.vars$"), "+bs(reg.vars$mval, degree = 3)", sep = "")
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
