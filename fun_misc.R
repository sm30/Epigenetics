countby <- function(dat,byvar,freqvar) {
  aggdata <- by(dat$byvar,tail,n=1)
  aggdatadf <- do.call("rbind", as.list(aggdata))
  table(aggdatadf$freqvar)

}
countby <- function(dat,byvar,freqvar) {
  aggdata <- with(dat, by(byvar,tail,n=1))
  aggdatadf <- do.call("rbind", as.list(aggdata))
  table(aggdatadf$freqvar)
}