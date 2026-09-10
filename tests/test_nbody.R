source('functions.R')
source('nbody.R')

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

####orbits from the RV parameters: Jupiter and Saturn from their RV signals
par <- c(P1=4332.6,K1=12.5,e1=0.049,omega1=0.26,Mo1=1,P2=10759,K2=2.76,e2=0.056,omega2=1.6,Mo2=2,gamma=0,sj=1)
o <- rv2orbits(par,Mstar=1,inc=90)
expect_true(nrow(o)==2 && abs(o$msini[1]-1)<0.02 && abs(o$msini[2]-0.3)<0.01,paste('m sin i:',paste(signif(o$msini,3),collapse=',')))
expect_true(abs(o$a[1]-5.2)<0.01 && abs(o$a[2]-9.54)<0.02,paste('a:',paste(signif(o$a,3),collapse=',')))
expect_true(abs(rv2orbits(par,1,30)$m[1]/o$m[1]-2)<1e-6,'the inclination does not scale the mass')
####circular signals are accepted too
oc <- rv2orbits(c(P1=10,A1=3,B1=4),1,90)
expect_true(!is.null(oc) && abs(oc$K-5)<1e-9 && oc$e==0,'a circular signal is not converted')

####analytical criteria
amd <- amd.stability(o,1)
expect_true(amd$AMD.stable && amd$Hill.stable && amd$Hill.separation>7,'Jupiter and Saturn are not AMD/Hill stable')
close <- rv2orbits(c(P1=365.25,K1=100,e1=0.1,omega1=0,Mo1=0,P2=420,K2=100,e2=0.1,omega2=1,Mo2=1),1,90)
amd2 <- amd.stability(close,1)
expect_true(!amd2$AMD.stable && !amd2$Hill.stable,'a close massive pair is not flagged by the analytical criteria')

####Monte Carlo samples from a 1D-like result
pl <- rbind(xopt=par,mode=par,xminus.1sig=par*0.99,xplus.1sig=par*1.01)
mc <- do.call(rbind,lapply(1:50,function(i) par*(1+rnorm(length(par),0,0.01)))); colnames(mc) <- names(par)
res <- list(par.list=list(RV=pl),mc.list=list(RV=mc))
s <- nbody.samples(res,'RV',N=5,seed=2)
expect_true(length(s$samples)==5 && grepl('MCMC',s$method) && all(sapply(s$samples,function(v) abs(v[['P1']]/par[['P1']]-1)<0.1)),'sampling from the MCMC failed')
s2 <- nbody.samples(list(par.list=list(RV=pl),mc.list=list()),'RV',N=3)
expect_true(length(s2$samples)==3 && grepl('Gaussian',s2$method),'sampling from the quantiles failed')
s3 <- nbody.samples(list(par.list=list(RV=par),mc.list=list()),'RV',N=3)
expect_true(length(s3$samples)==0,'a plain parameter vector must give the nominal solution only')

####the R integrator: a stable pair keeps its orbits, a close pair does not
r <- nbody.run(list(o),Mstar=1,tmax=1000,steps.per.orbit=100,Nout=30,engine='R')
cls <- nbody.classify(r,list(o))
expect_true(cls$stable && cls$max.da<0.02 && cls$dE<1e-3,paste('the R integrator did not keep Jupiter and Saturn stable:',cls$status,cls$max.da,cls$dE))
r2 <- nbody.run(list(close),Mstar=1,tmax=1000,steps.per.orbit=100,Nout=30,engine='R')
expect_true(!nbody.classify(r2,list(close))$stable,'the R integrator did not find the close pair unstable')
f <- tempfile(fileext='.pdf'); pdf(f,7,7); nbody.plot(r,list(o)); nbody.mc.plot(rbind(cls,cls),1000); dev.off()
expect_true(file.info(f)$size>2000,'the N-body figure is empty'); unlink(f)

####REBOUND, when available: the same two cases over a longer time
ver <- nbody.python('python3')
if(!is.null(ver)){
    rb <- nbody.run(list(o,close),Mstar=1,tmax=5e4,steps.per.orbit=20,Nout=50,engine='rebound')
    cb <- nbody.classify(rb,list(o,close))
    expect_true(cb$stable[1] && cb$max.da[1]<0.02 && cb$dE[1]<1e-4,paste('REBOUND did not keep Jupiter and Saturn stable:',cb$status[1],cb$dE[1]))
    expect_true(!cb$stable[2] && cb$t.end[2]<5e4,'REBOUND did not find the close pair unstable')
    expect_true(identical(rb$engine,'REBOUND WHFast'),'the engine label is wrong')
    cat('REBOUND',ver,'checked\n')
}else{
    cat('REBOUND not available; the Python engine was not tested\n')
}

####background jobs: launch, follow the progress, collect the result; stop a long one
st <- list(Mstar=1,tmax=2000,steps.per.orbit=100,Nout=20,escape.factor=10,encounter.hill=1,engine='R',python='python3',a.tol=0.2,method='nominal')
job <- nbody.launch(list(o),st)
expect_true(is.finite(job$pid),'the background job did not start')
res <- NULL; seen <- c()
for(k in 1:120){ s <- nbody.status(job); seen <- c(seen,s$frac); if(s$done){ res <- s$result; break }; Sys.sleep(0.5) }
expect_true(!is.null(res) && is.null(s$error) && res$class$stable[1],paste('the background job did not return a stable result:',s$error))
expect_true(any(seen>0 & seen<1) || length(seen)<=2,'the progress of the background job never advanced within the run')
job2 <- nbody.launch(list(o),modifyList(st,list(tmax=1e7)))
Sys.sleep(3)
expect_true(!nbody.status(job2)$done,'a long job finished too early')
nbody.stop(job2); Sys.sleep(1)
expect_true(!nbody.alive(job2$pid) && length(nbody.children(job2$pid))==0 && nbody.status(job2)$stopped,'the job was not stopped')

cat('N-body tests passed\n')
