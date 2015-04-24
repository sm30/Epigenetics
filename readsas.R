# read SAS
library(sas7bdat)
library(dplyr)
library(knitr)
library(readxl)
setwd("~/Dropbox/Projects/epigenetics/data")

# alles <- read.sas7bdat('nest_merge_tb_rc.sas7bdat')
# save(alles, file='alles.RData')

load('alles.RData')
alles_names <- names(alles)
alles_names[grep('_mean$', alles_names)]


sr <- read.csv("nestsr_iversen.csv", as.is=TRUE)
# sr <- read_excel("nestsr_iversen.xlsx")
sr_names <- names(sr)
sr$age_yr_SR_Surv <- floor(sr$age_mo_SR_Surv / 12)
# srs <- sr[sr$in_450k == 1,] # includes na and mess things up. don't know why!
srs <- filter(sr, in_450k == 1)
srs <- filter(srs, !(nestid == 532 & is.na(es0_age)))
srs_bkp <- srs

# newdata <- na.omit(srs[,grep('^nestid$|^es0_pa_basc|^es0_pa_brief|es0_age', sr_names)])
# newdata <- na.omit(srs[,grep('^nestid$|^es0_pa_basc|es0_age', sr_names)])
# newdata <- na.omit(srs[,grep('^nestid$|^es0_pa_brief|es0_age', sr_names)])
# newdata <- unique(srs[,grep('^es0_pa_basc', sr_names)])
# newdata <- unique(srs[,grep('^es0_pa_brief', sr_names)])

# apply(srs[,grep('^es0_pa_basc', sr_names)], 2, function(x) sum(!is.na(x)))
# apply(srs[,grep('^es0_pa_brief', sr_names)], 2, function(x) sum(!is.na(x)))
# newdata <- mutate(newdata, update = 1)
# newdata <- merge(srs, newdata, by = 'nestid')

names(srs) <- gsub("_brief_", "_brf_", names(srs))
BRF_names <- sr_names[grep('^BRF_', sr_names)]
BASC_names <- sr_names[grep('^BASC_', sr_names)]
for (x in c(BRF_names, BASC_names)) {
  x_old <- paste0(x, '_OLD')
  y <- paste0('es0_pa_', tolower(x))
  print(y)
  srs[,x_old] <- srs[,x]
  # srs <- srs %>% mutate(x = ifelse(!is.na(es0_pa_basc_hy) & !is.na(es0_pa_brf_sf), y, x))
  f <- paste0("srs <- srs %>% mutate(", x, " = ifelse(!is.na(es0_pa_basc_hy) & !is.na(es0_pa_brf_sf), ", y, ", ", x, "))")
  print(f)
  eval(parse(text=f))

#   srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs[, y]), x] <-
#     srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs[, y]), y]
  remove(y)
}

srs <- srs %>% mutate(age_mo_SR_Surv = ifelse(!is.na(es0_pa_basc_hy) & !is.na(es0_pa_brf_sf), es0_wasi_age, age_mo_SR_Surv))

# test it worked
for (x in c(BRF_names, BASC_names)) {
  x_old <- paste0(x, '_OLD')
  y <- paste0('es0_pa_', tolower(x))
  print(cbind(srs[!is.na(srs[, y]), x],
              srs[!is.na(srs[, y]), x_old]))
#   print(cbind(srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs[, y]), x],
#     srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs[, y]), x_old]))
}

summary(srs_bkp[, grep('BASC', names(sr))])
summary(srs[, grep('BASC', names(sr))])

write.csv(srs, file = 'nestsr_iversen_2.csv')

for (x in c(BRF_names, BASC_names)) {
  x_old <- paste0(x, '_OLD')
  y <- paste0('es0_pa_', tolower(x))
  plot(srs[,x_old], srs[,y])
  plot(srs[,x_old], srs[,x])
  # plot(x, x_old, data=srs)
}

BRF_OLD_names <- paste0(BRF_names, '_OLD')
BASC_OLD_names <- paste0(BASC_names, '_OLD')
plot(srs[,BRF_OLD_names])
plot(srs[,BASC_OLD_names])
# cor(srs[,BASC_OLD_names], na.rm=TRUE)
# idx <- seq(1, length(names), floor(length(names)/4))

knit2html('~/Dropbox/Projects/R_lib/niches/splot_vars.Rmd')


srs <- srs[, -grep('es0_', names(sr_names))]

# finished









for (x in BRF_names) {
  x_old <- paste0(x, '_OLD')
  y <- paste0('es0_pa_', tolower(x))
  print(x_old)
  print(y)
  srs[,x_old] <- srs[,x]
  print("s")
  srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs[, y]), x] <-
    srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs[, y]), y]
  remove(y)
}


sr_names[grep('brf', sr_names)]
sr_names[grep('BRF', sr_names)]

update_var <- function(x) {
  y <- paste0('es0_pa_', tolower(x))
  print(y)
  srs$BASC_AG_OLD <- srs$BASC_AG
  srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs[, y]), x] <-
    srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs[, y]), y]
}
# update_var('BASC_AG')
lapply(BASC_names, function(x) update_var(x))


srs$BASC_AGA <- srs$BASC_AG
srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs$es0_pa_basc_ag), 'BASC_AGA'] <-
  srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs$es0_pa_basc_ag), 'es0_pa_basc_ag']
cbind(srs$BASC_AG, srs$BASC_AGA)

sum(!(!is.na(srs$BASC_AG)==!is.na(srs$BASC_AGA)))

sr_names <- gsub("_brief_", "_BRF_", sr_names)
sr_names <- gsub("_basc_", "_BASC_", sr_names)
sr_names[grep('age', sr_names)]

srs[!is.na(srs$es0_age) > !is.na(srs$age_yr_SR_Surv), c('es0_age', 'age_yr_SR_Surv', 'es0_pa_basc_ag', 'BASC_AG')]
srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv), c('es0_age', 'age_yr_SR_Surv', 'es0_pa_basc_ag', 'BASC_AG')]

if (!is.na(srs$es0_age > srs$age_yr_SR_Surv)) srs$BASC_AG <- srs$es0_pa_basc_ag

if (!is.na(srs$es0_age > srs$age_yr_SR_Surv))
srs$BASC_AG <- srs$es0_pa_basc_ag



srss <- mutate(srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs$es0_pa_basc_ag), ], BASC_AG = es0_pa_basc_ag)
srss[,'es0_pa_basc_ag']
srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv) & !is.na(srs$es0_pa_basc_ag),'es0_pa_basc_ag']


srs[!is.na(srs$es0_age > srs$age_yr_SR_Surv), c('es0_age', 'age_yr_SR_Surv', 'es0_pa_basc_ag', 'BASC_AG')]

srs[, c('es0_age', 'survey_age', 'age_mo_SR_Surv', 'age_yr_SR_Surv', 'es0_pa_basc_ag')]


