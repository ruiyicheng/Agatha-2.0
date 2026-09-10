source('periodograms.R')
source('periodoframe.R')
source('functions.R')
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)
Nmax.plots <- 50

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

####################################################
## A two-observable result (RV plus one activity proxy) so that the periodogram
## panels of the second observable exercise the level bookkeeping
####################################################
set.seed(9)
P <- 15.3
t <- sort(runif(80,0,300))
rv <- 4*cos(2*pi*t/P)+rnorm(80,0,0.8)
proxy <- 0.3*sin(2*pi*t/40)+rnorm(80,0,0.2)
d <- data.frame(Time=t,RV=rv,eRV=rep(0.8,80),Sindex=proxy)
renew <- TRUE; Nsamp <- 1
pp <- list(ns=c('RV','Sindex','Window Function'),ofac=1,frange=c(1/100,1/5),per.type='BFP',
           per.target='star',sequence=FALSE,Nmas=0,Nars=0,Inds=list(0),Nsig.max=1,
           per.type.seq='BFP',Niter=0,SigType='kepler',Nh=2)
out <- calc.1Dper(Nmax.plots=50,vars=c('RV','Sindex'),per.par=pp,data=list(star=d),Ncores=1)

sp <- list.single.plots(out)
expect_true(is.data.frame(sp) && nrow(sp)>=6,'list.single.plots did not enumerate the panels')
expect_true(all(c('periodogram','phase','fit','residual')%in%sp$kind),'list.single.plots is missing a panel kind')
expect_true(any(sp$ypar=='Sindex' & sp$kind=='periodogram'),'the proxy periodogram is not offered for download')

####every panel in every format must produce a non-empty file
for(k in 1:nrow(sp)) for(fmt in c('pdf','png','jpg')){
    f <- tempfile(fileext=paste0('.',fmt))
    save.single.plot(f,format=fmt,width=6,height=4.5,dpi=120,
                     plot1D.single(out,kind=sp$kind[k],ypar=sp$ypar[k],index=sp$index[k],SigType='kepler'))
    expect_true(file.exists(f) && file.info(f)$size>1000,
                paste('empty or missing',fmt,'for',sp$label[k]))
    unlink(f)
}

####the bundled multi-panel figures still draw with the new panel functions
f <- tempfile(fileext='.pdf')
pdf(f,8,8)
per1D.plot(out$per.list,out$tits,out$pers,out$levels,ylabs=out$ylabs,download=TRUE,SigType='kepler')
phase1D.plot(out$phase.list,out$sim.list,out$tits,download=TRUE,repar=FALSE)
dev.off()
expect_true(file.info(f)$size>5000,'the bundled PDF is empty')

####the periodogram panel reports the peak it annotates
####(with two harmonics the raw maximum may be the 2P alias; the annotation
####must follow the reported period)
pdf(NULL)
pk <- plot1D.single(out,'periodogram','RV',1,SigType='kepler')
dev.off()
expect_true(abs(pk$Popt-P)<0.5,paste('the annotated peak',pk$Popt,'is not the reported period',P))

####titles are human-readable, without the internal ";" coding
expect_true(!grepl(';',pretty.title(out$tits[1])) && grepl('BFP',pretty.title(out$tits[1])),
            paste('pretty.title produced',pretty.title(out$tits[1])))

####the parameter table is offered as a figure and lists the fitted parameters
expect_true(any(sp$kind=='partable' & sp$ypar=='RV'),'the parameter table is not offered for download')
expect_true(!any(sp$kind=='partable' & sp$ypar=='Window Function'),'a parameter table is offered for the window function')
tb <- par.table(out$par.list$RV)
expect_true(is.data.frame(tb) && 'P1'%in%tb$Parameter && ncol(tb)==2,'par.table of a maximum-likelihood fit is wrong')
stat <- rbind(xopt=1:3,xminus.1sig=0:2,xplus.1sig=2:4,mode=1:3,med=1:3,mean=1:3)
colnames(stat) <- c('P1','K1','e1')
tb <- par.table(stat)
expect_true(all(c('Parameter','MAP','Mean','Median','q16','q84')==colnames(tb)) && tb$q84[3]==4,
            'par.table of an MCMC summary is wrong')
pdf(NULL)
tb2 <- panel.partable(list(RV=stat),'RV',out$tits[1])
dev.off()
expect_true(identical(tb,tb2),'panel.partable does not draw the par.table')

####################################################
## Two data sets: the points are coloured by set and a legend names the sets
####################################################
set.seed(11)
t1 <- sort(runif(40,0,300)); t2 <- sort(runif(30,50,350))
mk <- function(tt,off) data.frame(Time=tt,RV=4*cos(2*pi*tt/P)+off+rnorm(length(tt),0,0.8),eRV=rep(0.8,length(tt)))
d2 <- list(harps=mk(t1,0),pfs=mk(t2,12))
pp2 <- list(ns=c('RV','Window Function'),ofac=1,frange=c(1/100,1/5),per.type='BFP',
            per.target=c('harps','pfs'),sequence=FALSE,Nmas=c(0,0),Nars=c(0,0),Inds=list(0,0),Nsig.max=1,
            per.type.seq='BFP',Niter=0,SigType='circular',Nh=1)
out2 <- calc.1Dper(Nmax.plots=50,vars='RV',per.par=pp2,data=d2,Ncores=1)
ph2 <- out2$phase.list$RV
expect_true('set0'%in%colnames(ph2),'the phase data of a multi-set fit have no set index')
expect_true(identical(attr(ph2,'sets'),c('harps','pfs')),'the set names are not attached to the phase data')
####the set index must follow the (time-sorted) points: every time of set 2 is a time of the pfs data
expect_true(all(ph2[ph2[,'set0']==2,'t0']%in%d2$pfs$Time) && sum(ph2[,'set0']==1)==40,
            'the set index does not match the data sets')
si <- set.info(ph2)
expect_true(identical(si$names,c('harps','pfs')) && length(si$col)==2 && si$col[1]!=si$col[2],
            'set.info does not describe the two sets')
expect_true(is.null(set.info(out$phase.list$RV)),'a single data set must not get a set legend')
sp2 <- list.single.plots(out2)
for(k in which(sp2$kind%in%c('phase','fit','residual','partable'))){
    f <- tempfile(fileext='.png')
    save.single.plot(f,format='png',width=6,height=4.5,dpi=120,
                     plot1D.single(out2,kind=sp2$kind[k],ypar=sp2$ypar[k],index=sp2$index[k],SigType='circular'))
    expect_true(file.exists(f) && file.info(f)$size>1000,paste('empty multi-set figure',sp2$label[k]))
    unlink(f)
}
f <- tempfile(fileext='.pdf')
pdf(f,8,8)
phase1D.plot(out2$phase.list,out2$sim.list,out2$tits,download=TRUE,repar=FALSE,par.list=out2$par.list)
dev.off()
expect_true(file.info(f)$size>5000,'the multi-set bundled PDF is empty')

####data-set labels keep the instrument only
expect_true(identical(set.label(c('HD1_PFS','HD1_HARPS')),c('PFS','HARPS')),'the star name is not dropped from the set labels')
expect_true(identical(set.label(c('HD1_PFS','HD2_PFS')),c('HD1_PFS','HD2_PFS')),'colliding instrument names must keep the full set names')
expect_true(identical(set.label(c('harps','pfs')),c('harps','pfs')),'names without a star part must be kept')
sets2 <- c('HD1_PFS','HD1_HARPS')
stat2 <- rbind(xopt=1:4,med=1:4,mean=1:4,xminus.1sig=0:3,xplus.1sig=2:5); colnames(stat2) <- c('P1','K1','gamma_HD1_PFS','sj_HD1_HARPS')
attr(stat2,'sets') <- sets2
expect_true(identical(par.table(stat2)$Parameter,c('P1','K1','gamma_PFS','sj_HARPS')),'the parameter table does not shorten the set names')

####################################################
## Sequential signals without MCMC: two injected signals, three searched; the
## model comparison accepts two and the combined model keeps only those
####################################################
set.seed(5)
t3 <- sort(runif(120,0,500))
y3 <- 5*cos(2*pi*t3/15.3)+3*sin(2*pi*t3/41)+rnorm(120,0,0.8)
d3 <- list(star=data.frame(Time=t3,RV=y3,eRV=rep(0.8,120)))
pp3 <- list(ns=c('RV','Window Function'),ofac=2,frange=c(1/100,1/5),per.type='BFP',per.target='star',sequence=TRUE,
            Nmas=0,Nars=0,Inds=list(0),Nsig.max=3,per.type.seq='BFP',Niter=0,SigType='circular',Nh=1,lnBF.min=5)
out3 <- calc.1Dper(Nmax.plots=50,vars='RV',per.par=pp3,data=d3,Ncores=1)
mc3 <- out3$model.comp$RV
expect_true(is.data.frame(mc3) && nrow(mc3)==3,'the model comparison does not list the three searched signals')
expect_true(all(mc3$accepted[1:2]) && !mc3$accepted[3],paste('the model comparison accepted',paste(mc3$accepted,collapse=','),'lnBF',paste(round(mc3$lnBF,1),collapse=',')))
expect_true(out3$Nopt$RV==2,paste('the most plausible number of signals is',out3$Nopt$RV,'instead of 2'))
expect_true(length(grep('^ysig_sig',colnames(out3$phase.list$RV)))==2,'the rejected signal is still part of the combined model')
expect_true(ncol(out3$per.list$RV)==4,'the periodogram of the rejected signal is not reported')
expect_true(any(list.single.plots(out3)$kind=='modelcomp'),'the model comparison is not offered as a figure')
f <- tempfile(fileext='.png')
save.single.plot(f,format='png',width=6,height=4.5,dpi=120,plot1D.single(out3,kind='modelcomp',ypar='RV'))
expect_true(file.info(f)$size>1000,'the model comparison figure is empty'); unlink(f)
####the bundled figures with tables spanning both columns
f <- tempfile(fileext='.pdf')
pdf(f,8,8)
combined.plot(out3$per.list,out3$phase.list,out3$sim.list,out3$tits,out3$pers,out3$levels,out3$ylabs,SigType='circular',download=TRUE,par.list=out3$par.list,model.comp=out3$model.comp)
dev.off()
expect_true(file.info(f)$size>5000,'the bundled PDF with the model comparison is empty')

####singular normal equations fall back to the nearest positive-definite matrix
####and still return a plain, nameable vector (the app once died on an S4 matrix)
lm <- matrix(c(1,1,1,1),2); v <- solve.try(lm,c(1,1))
expect_true(is.numeric(v) && !isS4(v) && length(v)==2 && all(is.finite(v)),'solve.try did not return a plain finite vector for a singular system')
names(v) <- c('a','b')
lm2 <- matrix(c(2,0.5,0.5,1),2); v2 <- solve.try(lm2,c(1,2))
expect_true(isTRUE(all.equal(as.numeric(v2),as.numeric(solve(lm2,c(1,2))))),'solve.try changed the regular solution')

####merging the parameters of sequential signals never runs off the vector,
####whatever the number of noise terms of the residual fit (the app once died
####with 'only 0's may be mixed with negative subscripts')
m <- addpar(c(P1=10,A1=1,B1=2),c(P=20,A=1,B=1,gamma=0,beta=0,sj=1,m1=0.1,w1=0.2),2)
expect_true(identical(names(m),c('P1','A1','B1','P2','A2','B2','gamma','beta','sj','m1','w1')),paste('addpar gave',paste(names(m),collapse=',')))
m2 <- addpar(c(P1=10,K1=1,e1=0.1,omega1=1,Mo1=1,gamma_a=1,gamma_b=2,beta=0,sj_a=1,sj_b=1),c(P1=20,K1=1,e1=0.2,omega1=1,Mo1=2,gamma=0,sj=1),2)
expect_true(identical(names(m2),c('P1','K1','e1','omega1','Mo1','P2','K2','e2','omega2','Mo2','gamma','sj')),paste('addpar gave',paste(names(m2),collapse=',')))
m3 <- addpar(m2,c(P1=30,K1=1,e1=0,omega1=0,Mo1=0,gamma=0,sj=1),3)
expect_true(all(c('P1','P2','P3','K3')%in%names(m3)) && sum(names(m3)=='gamma')==1,'addpar failed on the third signal')

cat('plot tests passed\n')
