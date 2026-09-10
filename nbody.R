###########################################################################
####N-body stability of the orbital solutions.
####The planets come from the RV parameters of a 1D result (P, K, e, omega,
####Mo per signal, or P, A, B for circular signals) with a stellar mass and an
####inclination. Three levels of test: the analytical AMD (Laskar & Petit 2017;
####Petit, Laskar & Boue 2018) and Hill criteria; a numerical integration with
####REBOUND's WHFast (scripts/nbody_rebound.py) over e.g. 10 Myr; a short-term
####leapfrog in R when Python/REBOUND are not available. Either the nominal
####(MAP) solution or a Monte Carlo sample of the posterior is integrated.

nbody.G <- 4*pi^2                     #AU^3 Msun^-1 yr^-2
Mjup.Msun <- 9.5458e-4
nbody.script <- 'scripts/nbody_rebound.py'

nbody.python <- function(python='python3'){
###the REBOUND version seen by the interpreter, or NULL when it cannot be used
    if(is.null(python) || !nzchar(python) || is.na(Sys.which(python)) || !nzchar(Sys.which(python))) return(NULL)
    out <- tryCatch(suppressWarnings(system2(python,c('-c','"import rebound; print(rebound.__version__)"'),stdout=TRUE,stderr=TRUE)),
                    error=function(e) NULL)
    if(is.null(out) || length(out)==0) return(NULL)
    st <- attr(out,'status')
    if(!is.null(st) && st!=0) return(NULL)
    ver <- trimws(out[length(out)])
    if(!grepl('^[0-9]',ver)) return(NULL)
    ver
}

rv.signals <- function(par){
###the signals of an RV parameter vector: for every P_k the K, e, omega and Mo
###(Keplerian) or A and B (circular: K = sqrt(A^2+B^2), e = 0, Mo = -atan2(B,A))
    par <- unlist(par)
    nm <- names(par)
    ks <- sort(as.integer(sub('^P','',grep('^P[0-9]+$',nm,value=TRUE))))
    if(length(ks)==0) return(NULL)
    rows <- list()
    for(k in ks){
        g <- function(x) if(paste0(x,k)%in%nm) as.numeric(par[[paste0(x,k)]]) else NA
        P <- g('P'); K <- g('K'); e <- g('e'); om <- g('omega'); Mo <- g('Mo')
        if(is.na(K)){
            A <- g('A'); B <- g('B')
            if(is.na(A) || is.na(B)) next
            K <- sqrt(A^2+B^2); e <- 0; om <- 0; Mo <- -atan2(B,A)
        }
        if(is.na(e)) e <- 0
        if(is.na(om)) om <- 0
        if(is.na(Mo)) Mo <- 0
        if(!is.finite(P) || P<=0 || !is.finite(K)) next
        rows[[length(rows)+1]] <- data.frame(signal=k,P=P,K=abs(K),e=min(max(e,0),0.99),omega=om,Mo=Mo)
    }
    if(length(rows)==0) return(NULL)
    do.call(rbind,rows)
}

rv2orbits <- function(par,Mstar=1,inc=90){
###orbital elements of the planets: a [AU], m sin i [Mjup], m [Msun] for the
###inclination inc [deg] (coplanar system), the planet's argument of
###periastron (the RV omega is the star's, omega_p = omega + pi) and the mean
###anomaly Mo at the reference epoch (the first observation)
    sig <- rv.signals(par)
    if(is.null(sig)) return(NULL)
    sini <- max(sin(inc*pi/180),1e-3)
    msini <- 4.919e-3*sig$K*sqrt(1-sig$e^2)*sig$P^(1/3)*Mstar^(2/3)
    m <- msini/sini*Mjup.Msun
    a <- ((Mstar+m)*(sig$P/365.25)^2)^(1/3)
    data.frame(signal=sig$signal,P=sig$P,K=sig$K,e=sig$e,omega=sig$omega,Mo=sig$Mo%%(2*pi),a=a,msini=msini,m=m,
               omega.p=(sig$omega+pi)%%(2*pi),inc=inc)
}

amd.stability <- function(orb,Mstar=1){
###AMD-stability of every adjacent pair (Laskar & Petit 2017) with the Hill
###limit of Petit, Laskar & Boue (2018), for coplanar orbits: the relative AMD
###of the whole system C (in units of the outer planet's Lambda) against the
###critical values for collision and for Hill stability, and the separation in
###mutual Hill radii (stable for two planets when above 2 sqrt(3))
    if(is.null(orb) || nrow(orb)<2) return(NULL)
    o <- orb[order(orb$a),]
    Lam <- o$m*sqrt(nbody.G*Mstar*o$a)
    AMD <- sum(Lam*(1-sqrt(1-o$e^2)))
    rows <- list()
    for(k in 1:(nrow(o)-1)){
        m1 <- o$m[k]; m2 <- o$m[k+1]; a1 <- o$a[k]; a2 <- o$a[k+1]
        alpha <- a1/a2; gam <- m1/m2; eps <- (m1+m2)/Mstar
        C <- AMD/Lam[k+1]
###collision: the least relative AMD on the collision curve alpha e1 + e2 = 1
        f <- function(e1){ e2 <- 1-alpha*e1; gam*sqrt(alpha)*(1-sqrt(1-e1^2))+(1-sqrt(1-e2^2)) }
        Cc.coll <- optimize(f,c(0,min(1,1/alpha)-1e-9))$objective
        Cc.hill <- gam*sqrt(alpha)+1-(1+gam)^1.5*sqrt(alpha/(gam+alpha)*(1+3^(4/3)*eps^(2/3)*gam/(1+gam)^2))
        Cc <- min(Cc.coll,Cc.hill,na.rm=TRUE)
        rh <- ((m1+m2)/(3*Mstar))^(1/3)*(a1+a2)/2
        rows[[k]] <- data.frame(pair=paste0(o$signal[k],'-',o$signal[k+1]),a.inner=a1,a.outer=a2,relative.AMD=C,
                                critical.collision=Cc.coll,critical.Hill=Cc.hill,AMD.ratio=C/Cc,
                                Hill.separation=(a2-a1)/rh,AMD.stable=C<Cc,Hill.stable=(a2-a1)/rh>2*sqrt(3),
                                stringsAsFactors=FALSE)
    }
    do.call(rbind,rows)
}

nbody.samples <- function(res1D,ypar='RV',N=0,seed=NULL){
###the nominal RV solution (MAP, or the mode/maximum likelihood) and N Monte
###Carlo draws of the signal parameters: rows of the MCMC sample when it holds
###every signal, else Gaussian draws from the 1-sigma quantiles of the
###parameter table, else none
    pl <- res1D$par.list[[ypar]]
    if(is.null(pl)) stop('no fitted parameters for ',ypar)
    if(is.matrix(pl)){
        r <- if('xopt'%in%rownames(pl) && all(is.finite(pl['xopt',]))) 'xopt' else if('mode'%in%rownames(pl)) 'mode' else rownames(pl)[1]
        nominal <- setNames(as.numeric(pl[r,]),colnames(pl))
    }else{
        nominal <- unlist(pl)
    }
    sig <- rv.signals(nominal)
    if(is.null(sig)) stop('no Keplerian or circular signal among the fitted parameters')
    need <- unlist(lapply(sig$signal,function(k) paste0(c('P','K','e','omega','Mo','A','B'),k)))
    need <- need[need%in%names(nominal)]
    out <- list(nominal=nominal,samples=list(),method='nominal solution only')
    if(N>0){
        if(!is.null(seed)) set.seed(seed)
        mc <- res1D$mc.list[[ypar]]
        if(is.matrix(mc) && all(need%in%colnames(mc)) && nrow(mc)>=2){
            rows <- sample(nrow(mc),N,replace=nrow(mc)<N)
            out$samples <- lapply(rows,function(i){ v <- nominal; v[need] <- as.numeric(mc[i,need]); v })
            out$method <- paste0(N,' draws from the MCMC posterior (',nrow(mc),' samples)')
        }else if(is.matrix(pl) && all(c('xminus.1sig','xplus.1sig')%in%rownames(pl))){
            sd <- setNames(as.numeric(pl['xplus.1sig',need]-pl['xminus.1sig',need])/2,need)
            sd[!is.finite(sd)] <- 0
            ie <- grep('^e[0-9]+$',need); iP <- grep('^P[0-9]+$',need)
            out$samples <- lapply(1:N,function(i){
                v <- nominal
                v[need] <- nominal[need]+rnorm(length(need),0,sd)
                if(length(ie)>0) v[need[ie]] <- pmin(pmax(v[need[ie]],0),0.95)
                if(length(iP)>0) v[need[iP]] <- pmax(v[need[iP]],1e-3)
                v
            })
            out$method <- paste0(N,' Gaussian draws from the 1-sigma quantiles')
        }else{
            out$method <- 'nominal solution only (no MCMC sample or uncertainties available)'
        }
    }
    out
}

####numerical integration ######################################################
nbody.config <- function(orbit.list,Mstar,tmax,steps.per.orbit,Nout,escape.factor,encounter.hill,file,progress.file=NULL,offset=0,total=length(orbit.list)){
    con <- file(file,'w'); on.exit(close(con))
    writeLines(c(paste('Mstar',Mstar),paste('tmax',tmax),paste('steps_per_orbit',steps.per.orbit),paste('Nout',Nout),
                 paste('escape_factor',escape.factor),paste('encounter_hill',encounter.hill),
                 paste('system_offset',offset),paste('system_total',total),
                 if(!is.null(progress.file)) paste('progress_file',progress.file)),con)
    for(i in seq_along(orbit.list)){
        o <- orbit.list[[i]]
        writeLines(paste('system',i),con)
        writeLines(sprintf('planet %.10g %.10g %.10g %.10g %.10g %.10g',o$m,o$a,o$e,o$omega.p,o$Mo,0),con)
    }
}

nbody.rebound <- function(orbit.list,Mstar=1,tmax=1e7,steps.per.orbit=20,Nout=200,escape.factor=10,encounter.hill=1,
                          python='python3',progress=NULL,progress.file=NULL){
###WHFast integration of every system through scripts/nbody_rebound.py, one
###call per system so that progress can be reported; the driver writes its
###progress within a system to progress.file when given
    script <- nbody.script
    if(!file.exists(script)) stop('the REBOUND driver ',script,' is missing')
    summ <- list(); tracks <- list()
    for(i in seq_along(orbit.list)){
        if(is.function(progress)) progress(i,length(orbit.list))
        cfg <- tempfile(fileext='.txt'); orb <- tempfile(fileext='.txt'); sm <- tempfile(fileext='.txt')
        nbody.config(orbit.list[i],Mstar,tmax,steps.per.orbit,Nout,escape.factor,encounter.hill,cfg,progress.file=progress.file,offset=i-1,total=length(orbit.list))
        log <- suppressWarnings(system2(python,c(script,cfg,orb,sm),stdout=TRUE,stderr=TRUE))
        if(!file.exists(sm) || file.info(sm)$size==0){
            stop('the REBOUND integration failed: ',paste(tail(log,3),collapse=' '))
        }
        s <- read.table(sm,header=TRUE,stringsAsFactors=FALSE); s$system <- i
        tr <- read.table(orb,header=TRUE); tr$system <- i
        summ[[i]] <- s; tracks[[i]] <- tr
        unlink(c(cfg,orb,sm))
    }
    list(summary=do.call(rbind,summ),tracks=do.call(rbind,tracks),engine='REBOUND WHFast')
}

####a short-term leapfrog in R (no Python needed) ###############################
kepler.E <- function(M,e){
###eccentric anomaly by Newton's method
    E <- if(e<0.8) M else pi
    for(k in 1:50){
        d <- (E-e*sin(E)-M)/(1-e*cos(E))
        E <- E-d
        if(all(abs(d)<1e-12)) break
    }
    E
}

kep2cart <- function(a,e,inc,Omega,omega,M,mu){
###heliocentric position and velocity from the elements (angles in radians)
    E <- kepler.E(M%%(2*pi),e)
    nu <- 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2))
    r <- a*(1-e*cos(E))
    h <- sqrt(mu*a*(1-e^2))
    xp <- r*cos(nu); yp <- r*sin(nu)
    vxp <- -mu/h*sin(nu); vyp <- mu/h*(e+cos(nu))
    cO <- cos(Omega); sO <- sin(Omega); ci <- cos(inc); si <- sin(inc); cw <- cos(omega); sw <- sin(omega)
    R <- rbind(c(cO*cw-sO*sw*ci,-cO*sw-sO*cw*ci,sO*si),
               c(sO*cw+cO*sw*ci,-sO*sw+cO*cw*ci,-cO*si),
               c(sw*si,cw*si,ci))
    list(r=as.numeric(R%*%c(xp,yp,0)),v=as.numeric(R%*%c(vxp,vyp,0)))
}

cart2kep <- function(r,v,mu){
###semi-major axis and eccentricity from a heliocentric state vector
    rr <- sqrt(sum(r^2)); vv2 <- sum(v^2)
    a <- 1/(2/rr-vv2/mu)
    h <- c(r[2]*v[3]-r[3]*v[2],r[3]*v[1]-r[1]*v[3],r[1]*v[2]-r[2]*v[1])
    ev <- (c(v[2]*h[3]-v[3]*h[2],v[3]*h[1]-v[1]*h[3],v[1]*h[2]-v[2]*h[1]))/mu-r/rr
    c(a=a,e=sqrt(sum(ev^2)))
}

nbody.leapfrog <- function(orb,Mstar=1,tmax=1e3,steps.per.orbit=100,Nout=200,escape.factor=10,encounter.hill=1,progress.file=NULL,label=c(1,1)){
###kick-drift-kick leapfrog of the star and its planets in barycentric
###coordinates, with the planets' a and e (heliocentric) at log-spaced times;
###stops at an escape or a close encounter, like the REBOUND driver
    n <- nrow(orb)
    m <- c(Mstar,orb$m)
    x <- matrix(0,n+1,3); v <- matrix(0,n+1,3)
    for(k in 1:n){
        s <- kep2cart(orb$a[k],orb$e[k],0,0,orb$omega.p[k],orb$Mo[k],nbody.G*(Mstar+orb$m[k]))
        x[k+1,] <- s$r; v[k+1,] <- s$v
    }
    x <- sweep(x,2,colSums(x*m)/sum(m)); v <- sweep(v,2,colSums(v*m)/sum(m))
    Pmin <- min(orb$P)/365.25
    dt <- Pmin/steps.per.orbit
    amax0 <- max(orb$a)
    o <- orb[order(orb$a),]
    rh <- if(n>1) min(((o$m[-n]+o$m[-1])/(3*Mstar))^(1/3)*(o$a[-n]+o$a[-1])/2) else Inf
    accel <- function(x){
        acc <- matrix(0,n+1,3)
        for(i in 1:(n+1)){
            d <- sweep(x,2,x[i,])
            d2 <- rowSums(d^2); d2[i] <- Inf
            acc[i,] <- colSums(d*(nbody.G*m/d2^1.5))
        }
        acc
    }
    t1 <- max(10*Pmin,tmax/1e4)
    times <- if(tmax>t1) c(0,t1*(tmax/t1)^((0:(Nout-1))/(Nout-1))) else c(0,tmax)
    E0 <- 0.5*sum(m*rowSums(v^2))-sum(sapply(1:n,function(i) sum(nbody.G*m[i]*m[(i+1):(n+1)]/sqrt(rowSums(sweep(x[(i+1):(n+1),,drop=FALSE],2,x[i,])^2)))))
    energy <- function(){ 0.5*sum(m*rowSums(v^2))-sum(sapply(1:n,function(i) sum(nbody.G*m[i]*m[(i+1):(n+1)]/sqrt(rowSums(sweep(x[(i+1):(n+1),,drop=FALSE],2,x[i,])^2))))) }
    elements <- function(){
        t(sapply(1:n,function(k) cart2kep(x[k+1,]-x[1,],v[k+1,]-v[1,],nbody.G*(Mstar+m[k+1]))))
    }
    tr <- list(); status <- 'stable'; t.end <- tmax; t <- 0
    a <- accel(x)
    for(j in seq_along(times)){
        while(t<times[j]-1e-12){
            h <- min(dt,times[j]-t)
            v <- v+0.5*h*a; x <- x+h*v; a <- accel(x); v <- v+0.5*h*a
            t <- t+h
        }
        el <- elements()
        if(!is.null(progress.file)) try(writeLines(sprintf('%d %d %.4f',label[1],label[2],t/tmax),progress.file),silent=TRUE)
        tr[[j]] <- data.frame(system=1,t=t,planet=1:n,a=el[,'a'],e=el[,'e'])
        rr <- sqrt(rowSums(sweep(x[-1,,drop=FALSE],2,x[1,])^2))
        if(any(!is.finite(el)) || any(el[,'a']<=0) || any(el[,'e']>=1) || any(rr>escape.factor*amax0)){ status <- 'escape'; t.end <- t; break }
        if(n>1){
            dmin <- min(dist(x[-1,,drop=FALSE]))
            if(dmin<encounter.hill*rh){ status <- 'encounter'; t.end <- t; break }
        }
    }
    dE <- abs((energy()-E0)/E0)
    list(summary=data.frame(system=1,status=status,t_end=t.end,dE=dE,stringsAsFactors=FALSE),tracks=do.call(rbind,tr),engine='R leapfrog')
}

nbody.run <- function(orbit.list,Mstar=1,tmax=1e7,steps.per.orbit=20,Nout=200,escape.factor=10,encounter.hill=1,
                      engine=c('rebound','R'),python='python3',progress=NULL,progress.file=NULL){
###integrate every system of orbit.list (data frames from rv2orbits) and
###return the summary (status, t_end, energy error) and the a(t), e(t) tracks
    engine <- match.arg(engine)
    if(engine=='rebound'){
        return(nbody.rebound(orbit.list,Mstar,tmax,steps.per.orbit,Nout,escape.factor,encounter.hill,python=python,progress=progress,progress.file=progress.file))
    }
    summ <- list(); tracks <- list()
    for(i in seq_along(orbit.list)){
        if(is.function(progress)) progress(i,length(orbit.list))
        r <- nbody.leapfrog(orbit.list[[i]],Mstar,tmax,steps.per.orbit,Nout,escape.factor,encounter.hill,progress.file=progress.file,label=c(i,length(orbit.list)))
        r$summary$system <- i; r$tracks$system <- i
        summ[[i]] <- r$summary; tracks[[i]] <- r$tracks
    }
    list(summary=do.call(rbind,summ),tracks=do.call(rbind,tracks),engine='R leapfrog')
}

nbody.classify <- function(run,orbit.list,a.tol=0.2,e.max=0.9){
###one row per integrated system: the integrator's verdict, refined by the
###drift of the semi-major axes (|a/a0 - 1| > a.tol) and the eccentricities
###(e > e.max) along the tracks; unstable systems report the time of the event
    s <- run$summary
    out <- list()
    for(i in seq_len(nrow(s))){
        tr <- run$tracks[run$tracks$system==s$system[i],]
        a0 <- orbit.list[[s$system[i]]]$a
        da <- abs(tr$a/a0[tr$planet]-1)
        bad <- which(da>a.tol | tr$e>e.max)
        status <- s$status[i]; t.end <- s$t_end[i]
        if(status=='stable' && length(bad)>0){ status <- 'unstable (drift)'; t.end <- min(tr$t[bad]) }
        out[[i]] <- data.frame(system=s$system[i],status=status,t.end=t.end,max.da=max(da,na.rm=TRUE),max.e=max(tr$e,na.rm=TRUE),
                               dE=s$dE[i],stable=status=='stable',stringsAsFactors=FALSE)
    }
    do.call(rbind,out)
}

nbody.plot <- function(run,orbit.list,system=1,main='Nominal solution'){
###a(t) and e(t) of one integrated system, one colour per planet
    tr <- run$tracks[run$tracks$system==system,]
    orb <- orbit.list[[system]]
    n <- nrow(orb)
    cols <- set.palette[(1:n-1)%%length(set.palette)+1]
    op <- par(mfrow=c(2,1),mar=c(4.2,5,2.5,1.2),mgp=c(2.8,0.7,0),las=1); on.exit(par(op))
    tt <- pmax(tr$t,min(tr$t[tr$t>0]))
    plot(tt,tr$a,type='n',log='xy',xlab='Time [yr]',ylab='a [AU]',main=main)
    for(k in 1:n){ s <- tr$planet==k; lines(tt[s],tr$a[s],col=cols[k],lwd=1.5) }
    legend('topleft',bg='white',box.col='white',legend=paste0('signal ',orb$signal,' (P = ',signif(orb$P,4),' d)'),col=cols,lwd=1.5,cex=0.85)
    plot(tt,tr$e,type='n',log='x',xlab='Time [yr]',ylab='e',ylim=c(0,max(0.1,max(tr$e,na.rm=TRUE)*1.1)))
    for(k in 1:n){ s <- tr$planet==k; lines(tt[s],tr$e[s],col=cols[k],lwd=1.5) }
    invisible()
}

nbody.mc.plot <- function(cls,tmax,main='Monte Carlo sample'){
###the fraction of the sampled systems still stable as a function of time
    if(is.null(cls) || nrow(cls)==0) return(invisible())
    te <- ifelse(cls$stable,Inf,cls$t.end)
    tt <- sort(unique(c(te[is.finite(te)],tmax)))
    frac <- sapply(tt,function(x) mean(te>x))
    t0 <- min(c(tt[tt>0],tmax))/10
    op <- par(mar=c(4.2,5,2.5,1.2),mgp=c(2.8,0.7,0),las=1); on.exit(par(op))
    plot(c(t0,tt),c(1,frac),type='s',log='x',xlab='Time [yr]',ylab='Fraction of stable systems',ylim=c(0,1.05),main=main,lwd=2,col=pub.col$peak)
    text(t0,0.05,paste0(sum(cls$stable),' of ',nrow(cls),' systems stable to ',format(tmax,big.mark=','),' yr'),adj=c(0,0),col=pub.col$peak)
    invisible()
}

####background jobs, so that the app stays responsive and a run can be stopped
nbody.launch <- function(orbit.list,settings,wd=getwd()){
###start scripts/nbody_run.R in a separate R process on the systems of
###orbit.list with the integration settings (a list: Mstar, tmax,
###steps.per.orbit, Nout, escape.factor, encounter.hill, engine, python,
###a.tol); returns the job (its directory and process id)
    dir <- tempfile('nbody_'); dir.create(dir)
    saveRDS(list(orbit.list=orbit.list,settings=settings,wd=normalizePath(wd)),file.path(dir,'config.rds'))
###the whole group is backgrounded with its own redirections, so that the
###process id can be read without waiting for the job
    cmd <- sprintf('{ cd %s && %s %s %s %s ; } > %s 2>&1 < /dev/null & echo $!',shQuote(normalizePath(wd)),shQuote(file.path(R.home('bin'),'Rscript')),
                   shQuote('scripts/nbody_run.R'),shQuote(file.path(dir,'config.rds')),shQuote(dir),shQuote(file.path(dir,'log.txt')))
    pid <- suppressWarnings(as.integer(system(cmd,intern=TRUE)[1]))
    list(dir=dir,pid=pid,started=Sys.time(),n=length(orbit.list))
}

nbody.status <- function(job){
###progress of a job: system i of n and the fraction of that system done;
###done=TRUE with the result once scripts/nbody_run.R has finished; an error
###message when it failed
    out <- list(done=FALSE,i=0,n=job$n,frac=0,result=NULL,error=NULL,stopped=file.exists(file.path(job$dir,'stopped')))
    pf <- file.path(job$dir,'progress.txt')
    if(file.exists(pf)){
        p <- tryCatch(scan(pf,quiet=TRUE),error=function(e) NULL)
        if(length(p)>=3){ out$i <- p[1]; out$n <- p[2]; out$frac <- p[3] }
    }
    rf <- file.path(job$dir,'result.rds')
    if(file.exists(file.path(job$dir,'done'))){
        out$done <- TRUE
        if(file.exists(rf)) out$result <- tryCatch(readRDS(rf),error=function(e) NULL)
        ef <- file.path(job$dir,'error.txt')
        if(file.exists(ef)) out$error <- paste(readLines(ef),collapse=' ')
    }else if(!out$stopped && !is.na(job$pid) && !nbody.alive(job$pid)){
###the process is gone without reporting: treat as failed and show the log
        out$done <- TRUE
        lf <- file.path(job$dir,'log.txt')
        out$error <- paste('the integration process ended unexpectedly:',if(file.exists(lf)) paste(tail(readLines(lf),3),collapse=' ') else '')
    }
    out
}

nbody.alive <- function(pid){
    if(is.na(pid)) return(FALSE)
    r <- suppressWarnings(system(sprintf('kill -0 %d 2>/dev/null',pid),intern=FALSE))
    r==0
}

nbody.children <- function(pid){
###all descendants of a process (the runner's Rscript and its Python child)
    kids <- suppressWarnings(as.integer(system(sprintf('pgrep -P %d 2>/dev/null',pid),intern=TRUE)))
    kids <- kids[is.finite(kids)]
    c(kids,unlist(lapply(kids,nbody.children)))
}

nbody.stop <- function(job){
###kill the runner, its Python child and the launching shell, and mark the job as stopped
    if(is.null(job) || is.na(job$pid)) return(invisible(FALSE))
    pids <- unique(c(nbody.children(job$pid),job$pid))
    for(p in rev(pids)) suppressWarnings(system(sprintf('kill %d 2>/dev/null',p),intern=FALSE))
    writeLines(format(Sys.time()),file.path(job$dir,'stopped'))
    invisible(TRUE)
}
