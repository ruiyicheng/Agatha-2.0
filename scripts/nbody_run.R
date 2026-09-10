#!/usr/bin/env Rscript
###Runs an N-body job prepared by nbody.launch() (nbody.R) in its own process:
###  Rscript scripts/nbody_run.R config.rds jobdir
###writes progress.txt while integrating, then result.rds and 'done'
args <- commandArgs(trailingOnly=TRUE)
cfg <- readRDS(args[1]); dir <- args[2]
setwd(cfg$wd)
suppressMessages({ source('functions.R'); source('nbody.R') })
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)
s <- cfg$settings
pf <- file.path(dir,'progress.txt')
res <- tryCatch({
    run <- nbody.run(cfg$orbit.list,Mstar=s$Mstar,tmax=s$tmax,steps.per.orbit=s$steps.per.orbit,Nout=s$Nout,
                     escape.factor=s$escape.factor,encounter.hill=s$encounter.hill,engine=s$engine,python=s$python,
                     progress=function(i,n) writeLines(sprintf('%d %d 0',i,n),pf),progress.file=pf)
    cls <- nbody.classify(run,cfg$orbit.list,a.tol=s$a.tol)
    list(run=run,orbit.list=cfg$orbit.list,class=cls,method=s$method,tmax=s$tmax)
},error=function(e){ writeLines(conditionMessage(e),file.path(dir,'error.txt')); NULL })
if(!is.null(res)) saveRDS(res,file.path(dir,'result.rds'))
writeLines('done',file.path(dir,'done'))
