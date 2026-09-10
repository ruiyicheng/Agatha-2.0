source('periodograms.R')
source('periodoframe.R')
source('functions.R')
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)
renew <- TRUE; Nsamp <- 1

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

####################################################
## Two RV data sets with two signals; the 71 d one is mirrored by an activity
## index (H-alpha) and must be flagged, the 23.7 d one must pass every test
####################################################
set.seed(1)
mk <- function(n,a,b,off){
    tt <- sort(runif(n,a,b))
    data.frame(Time=tt,RV=6*cos(2*pi*tt/23.7)+2.5*sin(2*pi*tt/71)+off+rnorm(n,0,1),eRV=rep(1,n),
               Halpha=0.5*sin(2*pi*tt/71)+rnorm(n,0,0.3),BIS=rnorm(n,0,1))
}
d <- list(STAR_HARPS=mk(70,0,500,15),STAR_PFS=mk(60,300,900,-20))
dg <- diagnose.signals(d,names(d),noise.models=c('W','MA'),Nsig.max=3,ofac=1,frange=c(1/500,1/2),lnBF.min=5,Nwin=4,Ncores=1)
expect_true(identical(names(dg$combined),c('W','MA')),'the combined periodograms are missing a noise model')
expect_true(length(dg$individual)==2 && length(dg$proxies)==2,'the individual or proxy periodograms are missing')
expect_true(length(dg$periods)>=2 && any(abs(dg$periods-23.7)<1) && any(abs(dg$periods-71)<7),
            paste('the accepted periods are',paste(signif(dg$periods,4),collapse=', ')))
expect_true(length(dg$moving)==length(dg$periods) && isTRUE(dg$moving[[1]]$signal.only),'the moving periodograms of the signal-only series are missing')

tab <- diagnosis.table(dg)
print(tab)
expect_true(is.data.frame(tab) && nrow(tab)==length(dg$periods),'the diagnosis table has the wrong number of rows')
k1 <- which(abs(tab[['Period [d]']]-23.7)<1)[1]; k2 <- which(abs(tab[['Period [d]']]-71)<7)[1]
expect_true(tab$`Activity overlap`[k1]=='none','the 23.7 d signal is wrongly matched with an activity index')
expect_true(grepl('Halpha',tab$`Activity overlap`[k2]),'the 71 d signal is not matched with H-alpha')
expect_true(grepl('activity',tab$Verdict[k2]),paste('the verdict of the 71 d signal is',tab$Verdict[k2]))
expect_true(tab$`Robust to noise model`[k1]=='yes' && tab$`ln(BF) W`[k1]>5 && tab$`ln(BF) MA`[k1]>5,'the 23.7 d signal is not robust to the noise model')
expect_true(grepl('^candidate',tab$Verdict[k1]),paste('the verdict of the 23.7 d signal is',tab$Verdict[k1]))
tcf <- function(k){ c(as.integer(sub('/.*','',tab$`Consistent over time`[k])),as.integer(sub('.*/(\\d+).*','\\1',tab$`Consistent over time`[k]))) }
expect_true(all(is.finite(tcf(k1))) && tcf(k1)[2]==4 && tcf(k1)[1]>=3,paste('the 23.7 d signal is consistent in only',tab$`Consistent over time`[k1]))
####with the 23.7 d signal subtracted the 71 d one is seen in the windows too
expect_true(all(is.finite(tcf(k2))) && tcf(k2)[1]>=3,paste('the 71 d signal is consistent in only',tab$`Consistent over time`[k2]))
####a third accepted period within 10 per cent of the first is flagged as an alias or residual
k3 <- which(seq_len(nrow(tab))>k1 & abs(tab[['Period [d]']]-tab[['Period [d]']][k1])<0.1*tab[['Period [d]']][k1])
expect_true(length(k3)==0 || all(grepl('alias',tab$Verdict[k3]) & tab[[paste0('ln(BF) W')]][k3]<tab[[paste0('ln(BF) W')]][k1]),
            'a period close to an earlier signal is not flagged')

####the figure
pn <- diagnosis.panels(dg)
expect_true(length(pn)>=8 && all(sapply(pn,function(p) p$kind)%in%c('combined','window','individual','proxy')),'the diagnosis panels are incomplete')
f <- tempfile(fileext='.pdf')
pdf(f,10,9); n <- diagnosis.plot(dg,ncol=3,Prot=30); plotMP(dg$moving[[1]],dg$moving.par); dev.off()
expect_true(n==length(pn) && file.info(f)$size>5000,'the diagnosis figure is empty'); unlink(f)

####the GP noise model as one of the compared periodograms (one set, no extras, to keep it quick)
dgp <- diagnose.signals(d,'STAR_HARPS',noise.models='GP',Nsig.max=1,frange=c(1/500,1/2),individual=FALSE,proxies=FALSE,moving=FALSE,gp.Prot=30,Ncores=1)
expect_true(identical(names(dgp$combined),'GP') && length(dgp$periods)==1 && abs(dgp$periods-23.7)<1,
            paste('the GP diagnosis found',paste(signif(dgp$periods,4),collapse=', ')))
expect_true(grepl('GP',diagnosis.panels(dgp)[[1]]$title),'the GP panel is not labelled')
expect_true(!is.null(diagnosis.table(dgp)) && 'ln(BF) GP'%in%colnames(diagnosis.table(dgp)),'the GP column is missing from the table')

cat('signal diagnosis tests passed\n')
