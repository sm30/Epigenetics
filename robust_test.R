require(foreign)
require(MASS)
library(splines)
library(sfsmisc)
library(robust)

# cdata <- read.dta("http://www.ats.ucla.edu/stat/data/crime.dta")
setwd('~/Dropbox/Projects/epigenetics/data/')
cdata <- read.dta("crime.dta")
summary(cdata)

summary(ols <- lm(crime ~ poverty + single, data = cdata))
summary(rr.huber <- rlm(crime ~ poverty + single, data = cdata))
rr.bisquare <- rlm(crime ~ poverty + single, data=cdata, psi = psi.bisquare)
summary(rr.bisquare)

summary(ols)$tTable
(ols)$pvalue
anova(ols)
anova(rr.huber)
anova(rr.bisquare)


summary(ols <- lm(crime ~ poverty + bs(single, degree=3), data = cdata))
summary(rr.huber <- rlm(crime ~ poverty + bs(single, degree=3), data = cdata))
summary(rr.bisquare <- rlm(crime ~ poverty + bs(single, degree=3), data = cdata, psi = psi.bisquare))
anova(rr.huber)
["bs(reg.vars$mval, degree = 3)", "Pr(>F)"]

f.robftest(rr.huber)

summary(rr.lmRob <- lmRob(crime ~ poverty + single, data = cdata))
summary(rr.lmrob <- lmrob(crime ~ poverty + single, data = cdata))
summary(rr.bs.lmRob <- lmRob(crime ~ poverty + bs(single, degree=3), data = cdata))
summary(rr.bs.lmrob <- lmrob(crime ~ poverty + bs(single, degree=3), data = cdata))

ols.bs <- lm(crime ~ poverty + bs(single, degree=3), data = cdata)
summary(ols)
anova(ols)
anova(ols.bs)
anova(rr.lmRob)
anova(rr.lmrob)
anova(rr.bs.lmRob)["bs(single, degree = 3)", "Pr(F)"]
