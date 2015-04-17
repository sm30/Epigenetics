tf <- function(x) {
  f <- paste("Petal.Length ~ Petal.Width")
  assign('x', x, envir=.GlobalEnv) # problem with S4 with GlobalEnv
  lm.out <- lmRob(as.formula(f), data=x)
  pval <- anova(lm.out)["Petal.Width", "Pr(F)"]
  remove(x, envir=.GlobalEnv)
  return(pval)
}
tf(iris)
summary(tf(iris))

f <- paste("Petal.Length ~ Petal.Width")
lm.out <- lmRob(as.formula(f), data=iris)
pval <- anova(lm.out)["Petal.Width", "Pr(F)"]
pval

# lmRob is the issue for now

### robust lm with splines, site or region level
lmRob.spline <- function(m, reg.vars, y) {
  fmla<-paste('reg.vars$', y, "~bs(reg.vars$mval, degree = 3)", sep = "")
  print(fmla)
  n.probes = dim(m)[1]
  pval <- vector(mode = "numeric", length = n.probes)
  for (i in 1:1) {
    reg.vars$mval <- m[i,]
    lm.out <- lmRob(as.formula(fmla), data = reg.vars)
    pval[i] <- anova(lm.out)["bs(reg.vars$mval, degree = 3)", "Pr(F)"]
  }
  return(pval)
}
pData$mval <- NA

cl <- makeCluster(nofcl)
registerDoParallel(cl)
regr.site.pv <- foreach (y = bhv.vars, .verbose = TRUE, .packages = c('splines','robust', 'MASS', 'lattice', 'fit.models', 'robustbase', 'rrcov')) %dopar% {lmRob.spline(mvalues, reg.vars, y)}
stopCluster(cl)
