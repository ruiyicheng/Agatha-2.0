source('periodograms.R')
source('periodoframe.R')
source('functions.R')
library(fields)
library(magicaxis)
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

####################################################
## The moving (2D) periodogram across data sets, periodogram types, proxies,
## windows with too little data, and two data sets combined. Each case must
## run through per2D.data() and plotMP() (the app path) without error.
####################################################
data <- list(HIP88962_PFS=read.agatha.table('data/HIP88962_PFS.vels'),
             HD361_PFS=read.agatha.table('data/HD361/HD361_PFS.vels'),
             HD361_HARPSpre=read.agatha.table('data/HD361/HD361_HARPSpre.dat'),
             HD189567=read.agatha.table('data/HD189567_DACE_HARPS03.dat'))

mp.case <- function(targets,per.type,Dtfrac,Nbin,Nma=0,Inds=0,noise.model='ARMA',gp.Prot=NA,gp.tau=NA){
    d <- data[targets]
    tsp <- diff(range(unlist(lapply(d,function(x) x[,1]))))
    pars <- list(ns=c('RV','Window Function'),ofac=1,frange=c(1/tsp,1/2),per.type=per.type,per.target=targets,
                 files=NULL,Niter=0,Nmas=rep(Nma,length(targets)),Nars=rep(0,length(targets)),
                 Inds=rep(list(Inds),length(targets)),Dt=signif(tsp*Dtfrac,3),Nbin=Nbin,
                 alpha=5,scale=TRUE,pmin.zoom=2,pmax.zoom=tsp,show.signal=TRUE,noise.model=noise.model,gp.Prot=gp.Prot,gp.tau=gp.tau)
    v <- per2D.data(vars='RV',per.par=pars,data=d)
    f <- tempfile(fileext='.pdf')
    pdf(f,8,8); plotMP(v,pars); dev.off()
    expect_true(file.info(f)$size>2000,paste('empty 2D figure for',per.type,paste(targets,collapse='+')))
    unlink(f)
    v
}

####all periodogram types on a dense and a sparse data set
for(pt in c('GLS','BGLS','GLST','LS','MLP')){
    v <- mp.case('HD189567',pt,0.3,5)
    expect_true(ncol(v$zz)==5 && nrow(v$zz)==length(v$yy) && all(is.finite(v$zz)),
                paste(pt,'moving periodogram has the wrong shape or NA on dense data'))
    v <- mp.case('HIP88962_PFS',pt,0.4,5)
    expect_true(ncol(v$zz)==5,paste(pt,'moving periodogram has the wrong number of windows on sparse data'))
}

####LS must give real powers (it used to read a field lsp() does not return)
v <- mp.case('HD361_PFS','LS',0.5,6)
expect_true(all(is.finite(v$zz)) && diff(range(v$zz))>0,'LS moving periodogram is empty or constant')

####proxies and a moving-average term in small windows: singular fits must not abort
v <- mp.case('HD361_PFS','MLP',0.5,6,Nma=1,Inds=1)
expect_true(ncol(v$zz)==6,'MLP with proxy and MA(1) did not produce all windows')

####windows with too few points are skipped, not fatal; the others are computed
v <- mp.case('HIP88962_PFS','GLS',0.1,10)
expect_true(ncol(v$zz)==10,'tiny-window case did not return one column per window')
expect_true(any(is.na(v$zz)) && any(is.finite(v$zz)),
            'tiny-window case should have some skipped and some computed windows')

####two data sets combined
for(pt in c('GLS','MLP','BFP')){
    v <- mp.case(c('HD361_PFS','HD361_HARPSpre'),pt,0.5,6)
    expect_true(ncol(v$zz)==6 && length(v$idata)==2,paste('two-set moving periodogram failed for',pt))
}

####windows skipped before the period grid is known (a sparse start) still get
####a full column, so the matrix has one column per window and lines up with
####the window times (the app once died with 'subscript out of bounds')
set.seed(12)
tsp <- 1000
tg <- c(0,sort(runif(80,0.45*tsp,tsp)))
dg <- list(gap=data.frame(Time=tg,RV=4*cos(2*pi*tg/17)+rnorm(length(tg),0,1),eRV=rep(1,length(tg))))
parsg <- list(ns=c('RV','Window Function'),ofac=1,frange=c(1/tsp,1/2),per.type='MLP',per.target='gap',files=NULL,Niter=0,
              Nmas=0,Nars=0,Inds=list(0),Dt=signif(tsp*0.3,3),Nbin=5,alpha=5,scale=TRUE,pmin.zoom=2,pmax.zoom=tsp,show.signal=TRUE,
              noise.model='ARMA',gp.Prot=NA,gp.tau=NA)
vg <- per2D.data(vars='RV',per.par=parsg,data=dg)
expect_true(is.matrix(vg$zz) && nrow(vg$zz)==length(vg$yy) && ncol(vg$zz)==5 && ncol(vg$zz)==length(vg$xx),
            paste('a sparse start broke the moving periodogram matrix:',paste(dim(vg$zz),collapse='x')))
expect_true(all(is.na(vg$zz[,1])) && any(is.finite(vg$zz[,5])),'the skipped first window is not blank or the last window is empty')
f <- tempfile(fileext='.pdf'); pdf(f,8,8); plotMP(vg,parsg); dev.off()
expect_true(file.info(f)$size>2000,'the moving periodogram with a sparse start does not draw'); unlink(f)
####adaptive windows on the same data: every window holds data, the widths differ
parsa <- parsg; parsa$adaptive <- TRUE
va <- per2D.data(vars='RV',per.par=parsa,data=dg)
expect_true(isTRUE(va$adaptive) && ncol(va$zz)==5 && all(apply(va$zz,2,function(z) any(is.finite(z)))),'adaptive windows still leave an empty window')
w <- va$tend-va$tstart
expect_true(length(w)==5 && all(w>0) && w[1]>w[5]*1.5 && all(va$xx>=va$tstart & va$xx<=va$tend),paste('adaptive window widths:',paste(round(w),collapse=',')))
expect_true(all(abs(va$yy[apply(va$zz,2,which.max)]-17)<1),'the signal is not the strongest in every adaptive window')
f <- tempfile(fileext='.pdf'); pdf(f,8,8); plotMP(va,parsa); dev.off()
expect_true(file.info(f)$size>2000,'the adaptive moving periodogram does not draw'); unlink(f)

####the moving periodogram of a clean signal peaks at the right period in every window
set.seed(2)
t <- sort(runif(200,0,400)); P <- 21.3
d1 <- data.frame(Time=t,RV=5*sin(2*pi*t/P)+rnorm(200,0,0.5),eRV=rep(0.5,200))
data$sim <- d1
v <- mp.case('sim','GLS',0.5,4)
pk <- apply(v$zz,2,function(z) v$yy[which.max(z)])
expect_true(all(abs(pk-P)<1),paste('window peaks',paste(round(pk,2),collapse=','),'do not track the injected period',P))

####GP red noise inside the windows: free hyperparameters, a fixed oscillation
####period, and two data sets whose GP is removed per set before combining
v <- mp.case('HD361_PFS','MLP',0.6,4,noise.model='GP')
expect_true(ncol(v$zz)==4 && any(is.finite(v$zz)) && grepl('_GP',v$fname),'MLP with a GP in the windows failed')
v <- mp.case('HD361_PFS','MLP',0.6,4,noise.model='GP',gp.Prot=30)
expect_true(ncol(v$zz)==4 && any(is.finite(v$zz)),'MLP with a GP and a fixed oscillation period failed')
v <- mp.case(c('HD361_PFS','HD361_HARPSpre'),'MLP',0.5,4,noise.model='GP')
expect_true(ncol(v$zz)==4 && length(v$idata)==2,'two-set moving periodogram with GP failed')
####proxies are ignored with a GP only if none were selected; with one selected they still enter
v <- mp.case('HD361_PFS','MLP',0.6,4,noise.model='GP',Inds=1)
expect_true(ncol(v$zz)==4,'MLP with GP and a proxy failed')

####################################################
## The moving periodogram of the signal-only series of a 1D fit: two data sets
## with different offsets are combined on a common zero point
####################################################
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)
renew <- TRUE; Nsamp <- 1
set.seed(8)
mk <- function(n,a,b,off){ tt <- sort(runif(n,a,b)); data.frame(Time=tt,RV=6*cos(2*pi*tt/23.7)+off+rnorm(n,0,1),eRV=rep(1,n)) }
d2 <- list(STAR_HARPS=mk(60,0,400,15),STAR_PFS=mk(50,300,800,-20))
pp1 <- list(ns=c('RV','Window Function'),ofac=2,frange=c(1/200,1/5),per.type='BFP',per.target=names(d2),sequence=FALSE,
            Nmas=c(0,0),Nars=c(0,0),Inds=list(0,0),Nsig.max=1,per.type.seq='BFP',Niter=0,SigType='circular',Nh=1)
fit1 <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=pp1,data=d2,Ncores=1)
tsp <- 800
pars <- list(ns=c('RV','Window Function'),ofac=1,frange=c(1/tsp,1/5),per.type='MLP',per.target=names(d2),files=NULL,Niter=0,
             Nmas=c(0,0),Nars=c(0,0),Inds=list(0,0),Dt=signif(tsp*0.4,3),Nbin=5,alpha=5,scale=TRUE,pmin.zoom=5,pmax.zoom=tsp,
             show.signal=TRUE,noise.model='ARMA',gp.Prot=NA,gp.tau=NA,use.fit=TRUE)
v <- per2D.data(vars='RV',per.par=pars,data=d2,fit1D=fit1)
expect_true(isTRUE(v$signal.only) && length(v$idata)==2 && length(v$t)==110,'the signal-only moving periodogram did not use the 1D fit')
ph <- fit1$phase.list$RV
expect_true(isTRUE(all.equal(v$y,as.numeric(ph[order(ph[,'t0']),'y_all']))),'the moving periodogram input is not the combined-fit series')
####the offsets (15 and -20) are gone: both sets scatter about the same level
m <- sapply(v$idata,function(x) mean(x[,2]))
expect_true(abs(m[1]-m[2])<2,paste('the offsets were not removed from the combined series:',paste(round(m,2),collapse=', ')))
expect_true(all(is.finite(v$zz)) && ncol(v$zz)==5,'the signal-only moving periodogram has the wrong shape')
####the signal is present in every window
expect_true(all(abs(v$yy[apply(v$zz,2,which.max)]-23.7)<1.5),'the injected signal is not the strongest one in every window of the signal-only series')
f <- tempfile(fileext='.pdf'); pdf(f,8,8); plotMP(v,pars); dev.off()
expect_true(file.info(f)$size>2000,'empty signal-only 2D figure'); unlink(f)
####BFP on the combined signal-only series
parsb <- pars; parsb$per.type <- 'BFP'
vb <- per2D.data(vars='RV',per.par=parsb,data=d2,fit1D=fit1)
expect_true(all(is.finite(vb$zz)) && all(abs(vb$yy[apply(vb$zz,2,which.max)]-23.7)<1.5),
            'the BFP moving periodogram of the signal-only series does not find the signal in every window')
####the 1D fit must be for the same data sets
pars1 <- pars; pars1$per.target <- 'STAR_HARPS'; pars1$Nmas <- 0; pars1$Nars <- 0; pars1$Inds <- list(0)
err <- tryCatch({per2D.data(vars='RV',per.par=pars1,data=d2,fit1D=fit1); NULL},error=function(e) conditionMessage(e))
expect_true(!is.null(err) && grepl('same data sets',err),'a 1D fit for other data sets was not rejected')

cat('moving periodogram tests passed\n')
