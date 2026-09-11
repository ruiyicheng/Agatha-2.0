args<-commandArgs(trailingOnly=TRUE);source('periodograms.R');source('periodoframe.R');source('functions.R')
set.seed(91);t<-sort(runif(100,0,300));y<-5*cos(2*pi*t/17)+rnorm(100);dy<-runif(100,.3,.9);t<-t-min(t)
write.csv(data.frame(t=t,y=y,error=dy),file.path(args[1],'data.csv'),row.names=FALSE)
rows<-list()
for(model in c('WN','AR1','MA1','MA2')){
 p<-as.integer(model=='AR1');q<-if(model=='MA2')2 else as.integer(model=='MA1')
 v<-global.notation(t,y,dy,Nma=q,Nar=p,Indices=NULL,GP=FALSE,gp.par=rep(NA,3))
 s<-sopt(omega=2*pi/17,phi=0,Nma=q,Nar=p,Indices=NULL,data=cbind(t,y,dy),type='period',par.low=v$par.low,par.up=v$par.up,start=v$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=2)
 get<-function(n)if(n%in%names(s$par))as.numeric(s$par[n]) else NA_real_
 rows[[model]]<-data.frame(model=model,set_id='A',jitter=get('sj'),Nma=q,Nar=p,m1=get('m1'),m2=get('m2'),l1=get('l1'),tau=exp(get('logtau')),tauAR=exp(get('logtauAR')),logL=s$logL)
 write.csv(data.frame(residual=s$res),file.path(args[1],paste0(model,'.csv')),row.names=FALSE)
}
write.csv(do.call(rbind,rows),file.path(args[1],'noise.csv'),row.names=FALSE)
