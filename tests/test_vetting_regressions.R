source('periodograms.R');source('periodoframe.R');source('functions.R')
# The SHO kernel must be continuous at critical damping and finite at long lags.
t <- c(0,.1,10,1e6); S0 <- 1.7; P <- 12
K <- gp_sho_cov(t,S0,log(P),log(P/(2*pi)))
stopifnot(abs(K[1,1]-S0*2*pi/P*.5)<1e-12,all(is.finite(K)))
for(q in c(.5-1e-6,.5+1e-6,.01)){
 k <- gp_sho_cov(t,S0,log(P),log(q*P/pi));stopifnot(all(is.finite(k)))
 if(abs(q-.5)<1e-4) stopifnot(max(abs(k-K))<1e-5)
}
# A constant proxy and missing proxy values must not be replaced with noise;
# simultaneous but different observations must be preserved.
f <- tempfile();writeLines(c('Time RV eRV S Ha','1 2 1 0.2 NA','1 3 1 0.2 0.1','2 4 1 0.2 0.2','2 4 1 0.2 0.2'),f)
d<-read.agatha.table(f);stopifnot(nrow(d)==3,all(d$S==.2),is.na(d$Ha[1]),sum(d$Time==1)==2)
cat('SHO boundary and data-preservation regressions passed\n')
