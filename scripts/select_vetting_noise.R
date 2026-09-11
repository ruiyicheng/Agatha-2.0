# Fit the repository's single-instrument likelihood, without Shiny or MCMC.
args <- commandArgs(trailingOnly=TRUE)
source('periodograms.R'); source('periodoframe.R'); source('functions.R')
set.seed(7231)
input <- args[1]; output <- args[2]
instruments <- read.csv(file.path(input,'instruments.csv'))
rows <- list(); z <- 0
for(i in seq_len(nrow(instruments))){
 d <- read.csv(instruments$file[i]); t <- d$time-min(d$time); y <- d$rv-median(d$rv); dy <- d$error
 if(length(args)>2 && file.exists(args[3])){
  sub <- read.csv(args[3]); y <- y-sub$signal[sub$set_id==instruments$set_id[i]]
 }
 for(model in c('WN','AR1','MA1','MA2')){
  p <- as.integer(model=='AR1'); q <- if(model=='MA2') 2 else as.integer(model=='MA1')
  start.time <- proc.time()[3]
  msg <- ''; result <- NULL
  if(length(y)<=10 && model!='WN') msg <- 'Skipped: <=10 measurements'
  else if(length(y)<=3){
   # An instrument with 1--3 RVs cannot identify red noise or a free jitter.
   # Retain its measurements with their supplied uncertainties and an offset.
   mu <- weighted.mean(y,1/dy^2)
   result <- list(par=c(sj=0,offset=mu),logL=sum(dnorm(y,mu,dy,log=TRUE)))
   msg <- 'Sparse instrument: jitter fixed to zero; offset only'
  }
  else result <- tryCatch({
   v <- global.notation(t,y,dy,Nma=q,Nar=p,Indices=NULL,GP=FALSE,gp.par=rep(NA,3))
   sopt(omega=NA,phi=NA,Nma=q,Nar=p,Indices=NULL,data=cbind(t,y,dy),type='noise',par.low=v$par.low,par.up=v$par.up,start=v$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=3)
  },error=function(e){msg <<- conditionMessage(e);NULL})
  z<-z+1
  getpar <- function(name) if(!is.null(result) && name%in%names(result$par)) as.numeric(result$par[name]) else NA_real_
  ll <- if(is.null(result)) NA_real_ else result$logL
  k <- if(is.null(result)) NA_integer_ else length(result$par)
  if(length(y)<=3 && model=='WN' && !is.null(result)) k <- 1L
  rows[[z]] <- data.frame(set_id=instruments$set_id[i],instrument=instruments$instrument[i],n=length(y),model=model,Nar=p,Nma=q,logL=ll,k=k,BIC=-2*ll+k*log(length(y)),jitter=getpar('sj'),m1=getpar('m1'),m2=getpar('m2'),l1=getpar('l1'),tau=exp(getpar('logtau')),tauAR=exp(getpar('logtauAR')),seconds=proc.time()[3]-start.time,message=msg)
  cat(instruments$set_id[i],model,'BIC',rows[[z]]$BIC,'seconds',rows[[z]]$seconds,'\n')
 }
}
res <- do.call(rbind,rows); res$selected <- FALSE
for(id in unique(res$set_id)){
 inds <- which(res$set_id==id & is.finite(res$BIC))
 if(!length(inds)) stop(paste('No valid noise fit for',id))
 res$selected[inds[which.min(res$BIC[inds])]] <- TRUE
}
write.csv(res,output,row.names=FALSE)
