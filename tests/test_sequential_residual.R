source('periodograms.R')
source('periodoframe.R')
source('functions.R')
source('mcmc_func.R')
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

####################################################
## Two Keplerian signals in one data set: the second is searched in the data
## minus the MCMC (maximum-likelihood sample) solution of the first
####################################################
set.seed(21)
renew <- TRUE; Nsamp <- 1
kepler_rv <- function(t,P,K,e,omega,M0){
    M <- (M0+2*pi*t/P)%%(2*pi)
    E <- kep.mt2(M,e)
    nu <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
    K*(cos(omega+nu)+e*cos(omega))
}
P1 <- 15.3; P2 <- 41
t <- sort(runif(110,0,400))
dy <- rep(1,length(t))
y <- kepler_rv(t,P1,7,0.2,1.1,0.4)+kepler_rv(t,P2,3.5,0.1,2.0,1.5)+5+rnorm(length(t),0,dy)

b1 <- BFP(t=t,y=y,dy=dy,Nma=0,Nar=0,model.type='man',ofac=2,fmin=1/60,fmax=1/8,quantify=FALSE,renew=TRUE,progress=FALSE,Nh=2)
expect_true(abs(b1$Popt[1]-P1)<0.3,paste('BFP found',b1$Popt[1],'instead of',P1))
fit <- sigfit(per=b1,data=cbind(t,y,dy),SigType='kepler',mcf=TRUE,Niter=1000,Ncores=1,Pconv=FALSE)
####the residual handed to the next periodogram is the data minus the MCMC signal
expect_true(!is.null(fit$per$res.s) && isTRUE(all.equal(as.numeric(fit$per$res.s),as.numeric(y-fit$ysig0))),
            'sigfit() with MCMC does not pass on the data minus the MCMC signal as the residual')
expect_true(!isTRUE(all.equal(as.numeric(fit$ysig0),as.numeric(b1$ysig))),
            'the MCMC signal is identical to the periodogram optimum, so the MCMC solution was not used')
expect_true(abs(exp(fit$ParSig[['per1']])-P1)<0.3,'the MCMC did not stay on the first signal')
b2 <- BFP(t=t,y=as.numeric(fit$per$res.s),dy=dy,Nma=0,Nar=0,model.type='man',ofac=2,fmin=1/60,fmax=1/8,quantify=FALSE,renew=TRUE,progress=FALSE,Nh=2)
expect_true(abs(b2$Popt[1]-P2)<1,paste('the periodogram of the MCMC residual found',b2$Popt[1],'instead of',P2))

####the same through calc.1Dper with sequential signals and MCMC
pp <- list(ns=c('RV','Window Function'),ofac=2,frange=c(1/60,1/8),per.type='BFP',per.target='star',sequence=TRUE,
           Nmas=0,Nars=0,Inds=list(0),Nsig.max=2,per.type.seq='BFP',Niter=1000,SigType='kepler',Nh=2)
out <- calc.1Dper(Nmax.plots=50,vars='RV',per.par=pp,data=list(star=data.frame(Time=t,RV=y,eRV=dy)),Ncores=1)
pl <- out$par.list$RV
expect_true(is.matrix(pl) && all(c('P1','P2')%in%colnames(pl)),'the joint MCMC of both signals did not report P1 and P2')
expect_true(abs(pl['med','P1']-P1)<0.3 && abs(pl['med','P2']-P2)<1,
            paste('the sequential MCMC fit gave periods',pl['med','P1'],pl['med','P2'],'instead of',P1,P2))
mc <- out$model.comp$RV
expect_true(is.data.frame(mc) && nrow(mc)==2 && mc$source[2]=='MCMC' && all(mc$accepted),
            paste('the joint-MCMC model comparison gave',paste(mc$source,collapse=','),paste(round(mc$lnBF,1),collapse=',')))
expect_true(out$Nopt$RV==2,paste('the most plausible number of signals is',out$Nopt$RV,'instead of 2'))
per2 <- out$per.list$RV[,3]
expect_true(abs(out$per.list$RV[which.max(per2),1]-P2)<1,'the second periodogram does not peak at the second signal')

cat('sequential residual tests passed\n')
