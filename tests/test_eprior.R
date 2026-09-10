source('periodograms.R')
source('periodoframe.R')
source('functions.R')
source('mcmc_func.R')

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

####################################################
## The eccentricity prior descriptions
####################################################
d <- eprior.check(NULL)
expect_true(d$type=='beta' && abs(d$a-0.867)<1e-9 && abs(d$b-3.03)<1e-9 && d$sigma==0.1,'the default prior is not the Kipping beta')
expect_true(eprior.check(list(type='Uniform'))$type=='uniform','the uniform prior is not recognised')
expect_true(eprior.check(list(type='halfgauss',sigma=-1))$sigma==0.1,'an invalid sigma is not replaced by the default')
expect_true(eprior.check(list(type='nonsense'))$type=='beta','an unknown prior type is not replaced by the default')

####each density integrates to one over [0,1)
for(ep in list(list(type='beta'),list(type='beta',a=2,b=5),list(type='halfgauss',sigma=0.1),list(type='halfgauss',sigma=0.3),list(type='uniform'))){
    f <- function(x) sapply(x,function(e) exp(eprior.log(e,ep)))
    tot <- integrate(f,0,1)$value
    tol <- if(ep$type=='halfgauss') 2*pnorm(1,0,ep$sigma,lower.tail=FALSE)+1e-3 else 1e-2
    expect_true(abs(tot-1)<tol,paste('the',eprior.describe(ep),'prior integrates to',tot))
}
expect_true(eprior.log(0.5,list(type='uniform'))==0,'the uniform prior is not flat')
expect_true(abs(eprior.log(0.5,list(type='halfgauss',sigma=0.1))-log(2*dnorm(0.5,0,0.1)))<1e-12,'the half-Gaussian prior is wrong')
expect_true(abs(eprior.log(0.3,list(type='beta'))-dbeta(0.3,0.867,3.03,log=TRUE))<1e-12,'the beta prior is wrong')
expect_true(is.finite(eprior.log(0,list(type='beta'))),'the beta prior is not finite at e=0')
expect_true(eprior.log(1.2)==-Inf && eprior.log(-0.1)==-Inf,'eccentricities outside [0,1) are not excluded')
expect_true(eprior.log(c(0.1,0.2),list(type='halfgauss'))==eprior.log(0.1,list(type='halfgauss'))+eprior.log(0.2,list(type='halfgauss')),
            'the prior of several eccentricities is not the sum')

####################################################
## prior.func() of the single-set MCMC follows the chosen prior
####################################################
par.min <- c(per1=0,K1=0,e1=0,omega1=0,Mo1=0)
par.max <- c(per1=1,K1=10,e1=1,omega1=2*pi,Mo1=2*pi)
Esd <- 0.1
ins <- 'none'; out <- list(none=list(noise=list(nqp=c(0,0,0)))); phi.min <- -1; phi.max <- 1
pars <- c(per1=0.5,K1=5,e1=0.4,omega1=1,Mo1=1)
e.prior <- list(type='uniform')
p.uni <- prior.func(pars)
e.prior <- list(type='halfgauss',sigma=0.1)
p.hg <- prior.func(pars)
e.prior <- list(type='beta')
p.beta <- prior.func(pars)
expect_true(abs((p.hg-p.uni)-log(2*dnorm(0.4,0,0.1)))<1e-9,'prior.func does not apply the half-Gaussian eccentricity prior')
expect_true(abs((p.beta-p.uni)-dbeta(0.4,0.867,3.03,log=TRUE))<1e-9,'prior.func does not apply the beta eccentricity prior')
rm(e.prior)
expect_true(abs(prior.func(pars)-p.hg)<1e-9,'without e.prior the historical half-Gaussian of scale Esd must apply')

cat('eccentricity prior tests passed\n')
