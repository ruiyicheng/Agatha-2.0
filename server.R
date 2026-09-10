library(shiny)
library(magicaxis)
library(ramify)
library(utils)
options(warn=-2)
#library(ggplot2)
# use the below options code if you wish to increase the file input limit, in this example file input limit is increased from 5MB to 9MB
# options(shiny.maxRequestSize = 9*1024^2)
options(shiny.maxRequestSize=30*1024^2)
Nmax.plots <- 50
count0 <- 0
instruments <- c('HARPS','SOHPIE','HARPN','AAT','KECK','APF','PFS')
tol <- 1e-16
progress <- TRUE
#basis <- 'linear'
#renew <- FALSE
renew <- TRUE
Nsamp <- 2
#Nsamp <- 1
basis <- 'natural'
#trend <- FALSE
##load functions
source('periodoframe.R')
source("periodograms.R")
source('functions.R',local=TRUE)
source('mcmc_func.R')
###like functions.R, sourced into the server's environment so that its plot
###functions see the shared palette and style helpers
source('nbody.R',local=TRUE)

data.files <- list.files(path='data', pattern='\\.(dat|vels|rv)$', full.name=FALSE)

shinyServer(function(input, output, session){
####select from list
    output$about <- renderUI({HTML(paste("
<html>
<head>
<style>
p {
    text-indent: 25px;
}
</style>
</head>
<body>
<br />
<p>Agatha is the name of my wife's most favorite crime novelist, Agatha Christie. Similar to the investigations of various crimes in the detective novels, the Agatha algorithm is to find the weak signals embedded in correlated noise.

This web app is based on the code in GitHub: <a href='https://github.com/phillippro/agatha'>https://github.com/phillippro/agatha</a>. If you use this web app in your work, please cite 'Feng F., Tuomi M., Jones H. R. A., 2017, Agatha: disentangling periodic signals from correlated noise in a periodogram framework, MNRAS in press'. This paper is put on <a href='https://arxiv.org/abs/1705.03089'>arxiv</a>.</p>

<p>Agatha has the following features:</p>
<ul>
  <li>Fit the time-correlated noise using the moving average model</li>
  <li>Compare noise models to select the Goldilocks noise model (Feng et al. 2016, MNRAS, 461, 2440; available <a href='https://arxiv.org/abs/1606.05196'>here</a>)</li>
  <li>Optimize the frequency-dependent linear trend simultaneously with sinusoids and noise components</li>
  <li>Account for wavelength-dependent noise by fitting a set of linear functions of the difference between radial velocities (or other wavelength-dependent proxies for non-RV data sets) measured at different wavelengths</li>
  <li>Assess the significance of signals using the BIC-estimated Bayes factor</li>
  <li>Produce the so-called \"moving/2D periodogram\" to visualize the change of signals time thus visually testing the consistency of signals.</li>
</ul>

<p>
Agatha is based on the Bayes factor periodogram (BFP) and the marginalized likelihood periodogram (MLP). The BFP is calculated by maximizing the likelihood of a combination of sinusoids, linear functions of time and noisy proxies, and the moving average model. The Bayes factor for a given frequency is derived from the maximum likelihood by approximating the Bayes factor using the Bayes Information Criterion (BIC). The MLP is calculated by marginalizing the likelihood over the amplitudes of sinusoids and the parameters in a linear function of time. Before calculating MLP, the best-fitted noise model is subtracted from the data.
</p>

<p>
The BFP and MLP can be compared with the Lomb-Scargle periodogram (LS), the generalized LS (GLS), the GLS with floating trend (GLST) and the Bayesian GLS (BGLS). All periodograms can be computed for the sub-dataset within a moving time window to form 2D periodograms, which are also called 'moving periodograms'. Moving periodograms are used to check the consistency of signals in time. The user should adjust the 'visualization parameters' to optimize the visualization of signals in moving periodograms.
</p>
</body>
</html>
"))})

  output$files <- renderUI({
    if(is.null(input$uptype)) return()
    if(input$uptype=='list'){
        selectizeInput('target','Select data files from the list',
                  choices=gsub('\\..+','',gsub('_TERRA.+','',data.files)),multiple=TRUE)
    }else if(input$uptype=='upload'){
        fileInput('files', 'Choose files', multiple=TRUE)
#selectizeInput('Nf','Number of files to upload',choices=1:10,selected=1,multiple=FALSE)
    }
  })

    output$uptext <- renderUI({
        if(is.null(input$uptype)) return()
        if(input$uptype=='upload'){
            helpText("You can upload one or more plain-text time series files. Filenames should be 'star_instrument.fmt' where 'fmt' can be any plain-text extension. The first columns should be observation time, observable/RV, and optional uncertainty; two-column files get unit uncertainties. Additional columns are treated as noise proxies. After upload, select multiple data sets in the 1D Periodogram tab to fit them simultaneously.")
        }
    })

 output$nI.max <- renderUI({
      if(is.null(input$proxy.type)) return()
      if(input$proxy.type=='cum'){
          selectInput("ni.max",'Maximum number of noise proxies',choices = 0:NI.max()[input$comp.target],selected = min(3,NI.max()[input$comp.target]))
      }
  })

    NI.max <- reactive({
        if(is.null(data())) return()
        val <- c()
        for(i in 1:length(data())){
            val <- c(val,ncol(data()[[i]])-3)
        }
        names(val) <- names(data())
        return(val)
    })

    output$per.type.seq <- renderUI({
        if(is.null(input$sequence)) return()
        if(input$sequence) selectInput("per.type.seq",'Periodogram used to find additional signals',choices=input$per.type,selected=NULL,multiple=FALSE)
    })

    output$sequential <- renderUI({
        if(is.null(Ntarget()) | is.null(input$signal.type)) return()
        if(input$signal.type!='stochastic'){
            checkboxInput('sequence','Find additional signals sequentially',value=FALSE)
        }
    })

    output$Nsig.max <- renderUI({
        if(is.null(input$sequence)) return()
        if(input$sequence){
            tagList(
                sliderInput("Nsig.max", "Maximum number of signals", min = 2, max = 10,value=2,step=1),
                numericInput('lnBF.min','ln(BF) a further signal must exceed to be accepted',value=5,min=0,max=100,step=0.5),
                helpText('Signals are added one at a time and compared with the model without them (BIC-estimated ln(BF); from the joint MCMC of all signals for a single data set when MCMC is on, otherwise from the periodogram). The search stops at the first signal below the threshold; the most plausible number of signals is the number accepted, at most the maximum.')
            )
        }
    })

    output$noise.model <- renderUI({
        if(is.null(input$per.type)) return()
        if(any(input$per.type=='MLP'|input$per.type=='BFP')){
            radioButtons('noise.model','Red noise model',
                         c('ARMA (per data set)'='ARMA','Gaussian process (quasi-periodic SHO kernel, shared by all data sets)'='GP'),selected='ARMA')
        }
    })

    output$gp.par <- renderUI({
        if(is.null(input$noise.model) || input$noise.model!='GP') return()
        tagList(
            numericInput('gp.Prot','GP oscillation period [time unit of the data] (quasi-periodicity of the correlated noise, e.g. a stellar rotation period from photometry; leave empty to fit it)',value=NA,min=0),
            numericInput('gp.tau','GP damping time scale [time unit of the data] (how long correlations stay coherent, e.g. an active-region lifetime; leave empty to fit it)',value=NA,min=0),
            helpText('The GP amplitude is always fitted. If the damping time scale is much shorter than the oscillation period the kernel is overdamped and the period is not constrained. With the Stochastic signal type the oscillation period is scanned and a fixed value is ignored.')
        )
    })

    output$nar <- renderUI({
    if(is.null(input$per.target) | is.null(input$per.type)) return()
    if(!is.null(input$noise.model) && input$noise.model=='GP') return()
    if(any(input$per.type=='MLP'|input$per.type=='BFP')){
        lapply(1:Ntarget(), function(i){
            selectizeInput(paste0("Nar",i),paste('Number of AR components for',input$per.target[i]),choices = 0:10,selected = 0,multiple=FALSE)}
	)
        }
    })

  output$nma <- renderUI({
    if(is.null(input$per.target) | is.null(input$per.type)) return()
    if(!is.null(input$noise.model) && input$noise.model=='GP') return()
    if(any(input$per.type=='MLP'|input$per.type=='BFP')){
        lapply(1:Ntarget(), function(i){
            selectizeInput(paste0("Nma",i),paste('Number of MA components for',input$per.target[i]),choices = 0:10,selected = 0,multiple=FALSE)}
	)
        }
  })

    output$noise.model2 <- renderUI({
        if(is.null(input$per.type2)) return()
        if(any(input$per.type2=='MLP'|input$per.type2=='BFP')){
            radioButtons('noise.model2','Red noise model',
                         c('ARMA (per data set)'='ARMA','Gaussian process (quasi-periodic SHO kernel, shared by all data sets)'='GP'),selected='ARMA')
        }
    })

    output$gp.par2 <- renderUI({
        if(is.null(input$noise.model2) || input$noise.model2!='GP') return()
        tagList(
            numericInput('gp.Prot2','GP oscillation period [time unit of the data] (leave empty to fit it in every window)',value=NA,min=0),
            numericInput('gp.tau2','GP damping time scale [time unit of the data] (leave empty to fit it)',value=NA,min=0),
            helpText('The GP hyperparameters are fitted inside each moving window; fixing the oscillation period (e.g. from photometry) makes the windows comparable and much faster. BFP with a GP is slow: each window is a full scan.')
        )
    })

    output$nar2 <- renderUI({
    if(is.null(input$per.target2) | is.null(input$per.type2)) return()
    if(!is.null(input$noise.model2) && input$noise.model2=='GP') return()
    if(any(input$per.type2=='MLP'|input$per.type2=='BFP')){
        lapply(1:Ntarget2(), function(i){
            selectizeInput(paste0("Nar2.",i),paste('Number of AR components for',input$per.target2[i]),choices = 0:10,selected = 0,multiple=FALSE)}
	)
    }
    })

    output$nma2 <- renderUI({
    if(is.null(input$per.target2) | is.null(input$per.type2)) return()
    if(!is.null(input$noise.model2) && input$noise.model2=='GP') return()
    if(any(input$per.type2=='MLP'|input$per.type2=='BFP')){
        lapply(1:Ntarget2(), function(i){
            selectizeInput(paste0("Nma2.",i),paste('Number of MA components for',input$per.target2[i]),choices = 0:10,selected = 0,multiple=FALSE)}
	)
        }
    })

    output$per.target <- renderUI({
        if(is.null(data())) return()
        selected <- if(length(names(data()))>1) names(data()) else names(data())[1]
        selectizeInput("per.target",'Data sets',choices=names(data()),selected=selected,multiple=TRUE)
    })

    output$signal <- renderUI({
        if(is.null(input$per.type)) return()
        if(!is.null(Ntarget()) && Ntarget()>1){
            if(any(input$per.type=='BFP')){
                radioButtons("signal.type",'Signal type',c("Circular shared signal"="circular","Keplerian shared signal"='kepler','Stochastic'='stochastic'))
            }else{
                radioButtons("signal.type",'Signal type',c("Circular shared signal"="circular","Keplerian shared signal"='kepler'))
            }
        }else if(any(input$per.type=='BFP')){
            radioButtons("signal.type",'Signal type',c("Circular"="circular","Keplerian"='kepler','Stochastic'='stochastic'))
        }else{
            radioButtons("signal.type",'Signal type',c("Circular"="circular","Keplerian"='kepler'))
        }
    })

    output$nh <- renderUI({
        if(is.null(input$per.type)) return()
        if(any(input$per.type=='MLP'|input$per.type=='BFP')){
            selectizeInput('Nh','Number of harmonics of the signal (2 or more fits eccentric orbits)',choices=1:4,selected=2,multiple=FALSE)
        }
    })

    output$mcf <- renderUI({
        if(is.null(input$per.type)) return()
#        if(!is.null(data())){
        if(any(input$per.type=='BFP')){
            choices <- c(0,100,1000,10000,100000,1000000)
            tagList(
                selectizeInput('Niter','MCMC sample size', choices=choices,selected=0,multiple=FALSE),
###eccentricity prior of the MCMC; only meaningful for Keplerian signals
                conditionalPanel("input.Niter > 0 && input['signal.type'] == 'kepler'",
                    radioButtons('e.prior','Eccentricity prior of the MCMC',
                                 c('Beta distribution (Kipping 2013)'='beta','Uniform'='uniform','Half-Gaussian'='halfgauss'),selected='beta'),
                    conditionalPanel("input['e.prior'] == 'beta'",
                        fluidRow(column(6,numericInput('e.beta.a','Beta shape a',value=0.867,min=0.01,max=20,step=0.01)),
                                 column(6,numericInput('e.beta.b','Beta shape b',value=3.03,min=0.01,max=20,step=0.01))),
                        helpText('Beta(a, b) with a=0.867 and b=3.03 is the eccentricity distribution of RV planets found by Kipping (2013, MNRAS 434, L51).')),
                    conditionalPanel("input['e.prior'] == 'halfgauss'",
                        numericInput('e.sigma','Sigma of the half-Gaussian',value=0.1,min=0.001,max=1,step=0.01),
                        helpText('Half-Gaussian on e >= 0: the density of a Gaussian of zero mean and the given sigma, folded onto positive eccentricities.'))
                )
            )
        }
    })

    Ntarget <- reactive({
        if(is.null(input$per.target)) return()
        length(input$per.target)
    })

    output$per.target2 <- renderUI({
        if(is.null(data())) return()
        selectizeInput("per.target2",'Data sets',choices=names(data()),selected=names(data())[1],multiple=TRUE)
    })

    Ntarget2 <- reactive({
        length(input$per.target2)
    })

    output$per.type <- renderUI({
        if(is.null(Ntarget())) return()
        if(Ntarget()>1){
            selectInput("per.type",'Periodogram type',
                        choices=c('BFP','MLP'),selected="BFP",multiple=FALSE)
        }else{
            selectInput("per.type",'Periodogram type',
                        choices=c('BFP','MLP','GLST','BGLS','GLS','LS'),selected="BFP",multiple=TRUE)
        }
    })

    output$per.type2 <- renderUI({
        if(is.null(Ntarget2())) return()
        if(Ntarget2()>1){
###one type at a time for combined sets (the combined series is analysed white)
            selectInput("per.type2",'Periodogram type',choices=c('BFP','MLP','GLST','BGLS','GLS','LS'),selected="MLP",multiple=FALSE)
        }else{
            selectInput("per.type2",'Periodogram type',
                        choices=c('BFP','MLP','GLST','BGLS','GLS','LS'),selected="MLP",multiple=TRUE)
        }
    })

    output$text2D <- renderText({
        "<font color=\"DarkSlateGray\"><b>To make 2D periodograms, the time series should be properly sampled and each time window should contain at least a few data points, e.g. 100 data points over a time span beyond 200 time units. </b></font>"
    })

    output$Inds <- renderUI({
        if(is.null(input$per.type) | is.null(data()) | is.null(Ntarget()) | is.null(NI.max())) return()
#        if(all(NI.max()==0)) return()
        if(!any(input$per.type=='BFP' | input$per.type=='MLP')) return()
        lapply(1:Ntarget(),function(i){
            selectInput(paste0('Inds',i),paste('Noise proxies for',input$per.target[i]),choices = 0:NI.max()[input$per.target[i]],selected = 0,multiple=TRUE)
        })
    })

    output$Inds2 <- renderUI({
        if(is.null(data()) | is.null(input$per.target2)) return()
        lapply(1:Ntarget2(),function(i){
            selectInput(paste0("Inds2.",i),paste('Noise proxies for',input$per.target2[i]),choices = 0:NI.max()[input$per.target2[i]],selected = 0,multiple=TRUE)
        })
    })

    output$prange2 <- renderUI({
        if(is.null(input$Dt) | is.null(data()) | is.null(input$per.target2)) return()
        Dt <- signif(2*tspan()*as.numeric(input$Dt),3)
        fmin <- 1/Dt
        logpmin <- -2
        logpmax <- signif(log10(Dt),2)
        sliderInput("prange2","Period range in base-10 log scale",min = logpmin,max = logpmax,value = c(0.1,logpmax),step=0.1)
    })

    output$proxy.text <- renderUI({
      helpText("If 'cumulative' is selected, the noise proxies would be arranged in decreasing order of the Pearson correlation coefficients between proxies and RVs. Then proxies would be compared cumulatively from the basic number up to the maximum number of proxies. If 'group' is selected, the proxies would not be rearranged, and would be compared in groups, which are determined by the basic number of proxies and group division numbers. For example, if the basic number is 4 and division numbers are 6,11,19, models with proxies of {1-4}, {1-6}, {1-4, 7-11}, {1-4, 12-19} would be compared. If 'manual' is selected, the user should manually input the groups of proxies for comparison. Note that '0' means no proxy. ")
  })

    output$comp.target <- renderUI({
        if(is.null(data())) return()
        selectizeInput("comp.target",'Data sets',choices=names(data()),selected=names(data())[1],multiple=FALSE)
    })

  output$proxy.type <- renderUI({
      if(is.null(NI.max())) return()
      if(any(NI.max()[input$comp.target]>0)){
          radioButtons("proxy.type",'Type of proxy comparison',c("Cumulative"="cum","Group"='group','Manual'='man'))
      }else{
          output$warn <- renderText({'No indices available for comparison!'})
          verbatimTextOutput("warn")
      }
  })

  output$nI.basic <- renderUI({
      if(is.null(input$proxy.type)) return()
      if(NI.max()[input$comp.target]>0 & input$proxy.type!='man'){
#input$proxy.type=='cum'
          selectInput("NI0",'Basic number of noise proxies',choices = 0:NI.max()[input$comp.target],selected = 0)#NI.max())
      }
  })

    output$Nman <- renderUI({
        if(is.null(input$proxy.type)) return()
        if(input$proxy.type=='man' & NI.max()[input$comp.target]>0){
#            sliderInput("Nman",'Number of proxy groups',min=1,max=NI.max(),value=1,step=1)
            selectizeInput("Nman",'Number of proxy groups',choices=1:6,selected=1,multiple=FALSE)
        }
    })

    output$nI.man <- renderUI({
        if(is.null(input$proxy.type) | is.null(input$Nman)) return()
        if(input$proxy.type=='man' & NI.max()[input$comp.target]>0){
            lapply(1:as.integer(input$Nman), function(i) {
                selectInput(paste0("NI.man",i),paste0('proxy group ',i),
                            choices=0:NI.max()[input$comp.target],multiple = TRUE)
            })
        }
    })

  output$nI.comp <- renderUI({
    if(is.null(input$proxy.type) | is.null(input$NI0)) return()
    if(input$proxy.type=='group' & NI.max()[input$comp.target]>as.integer(input$NI0)){
      selectInput("NI.group",'Group division numbers',
                  choices=(as.integer(input$NI0)+1):NI.max()[input$comp.target],multiple = TRUE )
    }
  })

  target.list <- reactive({
    f1 <- gsub('_TERRA.+','',list.files('data',full.names=FALSE))
    gsub('.dat','',f1)
  })

  target <- reactive({
    if(is.null(input$target) & is.null(input$files)) return()
    if(input$uptype=='list'){
      return(input$target)
    }else if(input$uptype=='upload'){
        if(is.null(input$files) || nrow(input$files)==0) return()
        return(make.unique(vapply(input$files$name,agatha.dataset.name,character(1)),sep='_'))
    }
  })

  instr <- reactive({
      gsub('.+_','',target())
  })

  # added "session" because updateSelectInput requires it
  data <- eventReactive(input$show,{
    ins <- instr()
    #input$show
    if(input$uptype=='upload'){
        if(!is.null(input$files) && nrow(input$files)>0){
            tmp <- list(NA)
            df <- rep(tmp,nrow(input$files))
            for(i in 1:nrow(input$files)){
                data.path <- input$files[[i,'datapath']]
                df[[i]] <- read.agatha.table(data.path,center.rv=FALSE)
            }
            names(df) <- target()
        }
    }else if(input$uptype=='list'){
        if(is.null(input$target)) return()
        tmp <- list(NA)
        df <- rep(tmp,length(input$target))
        names(df) <- input$target
        for(i in 1:length(input$target)){
            target <- input$target[i]
            star <- gsub('_.+','',target)
            dir  <- 'data/'
            ind <- grep(target,data.files)
            file <- data.files[ind[1]]
            f0 <- paste0(dir,file)
            if(!file.exists(f0)){
                f0 <- paste0(dir,target,'.dat')
            }
            df[[i]] <- read.agatha.table(f0,center.rv=FALSE)
        }
    }
    return(df)
})

    observeEvent(input$show,{
        lapply(1:length(data()),function(j)
            output[[paste0('data.out',j)]] <- downloadHandler(
                filename = function() {
                    f1 <- gsub(" ",'_',Sys.time())
                    f2 <- gsub(":",'-',f1)
                    paste('data_',names(data())[j], f2, '.txt', sep='')
                },
                content = function(file) {
                    write.table(data()[[j]], file,quote=FALSE,row.names=FALSE)#FALSE,col.names=FALSE
                }
            )
            )
    })


    observeEvent(input$show,{
        output$download.data <- renderUI({
                lapply(1:length(data()),function(j){
                    downloadButton(paste0('data.out',j), paste('Download',names(data())[j]))
                })
            })
    })

    observeEvent(input$show,{
        output$tab <- renderUI({
            if(is.null(data())) return()
                isolate({
                    tabs <- lapply(1:length(target()),function(i){
                        output[[paste0('f',target()[i])]] <- renderDataTable(data()[[i]])
                        tabPanel(target()[i],dataTableOutput(paste0('f',target()[i])))
                    })
                    do.call(tabsetPanel, tabs)
                })
        })
    })

###variable names
  ns <- reactive({
    nam <- c()
    for(i in 1:length(instr())){
      names <- colnames(data()[[i]])
      names <- names[-c(1,3)]#e.g. 'Time' and 'eRV' for RV data
      names <- c(names,'Window Function')
#      nam <- c(nam,paste(names(data())[i],names,sep=':'))
      nam <- c(nam,names)
    }
    return(unique(nam))
  })

###variable names dependent on per.target
    ns.1D <- reactive({
        if(is.null(input$per.target)) return()
        tar <- input$per.target
        names.by.target <- lapply(tar,function(target){
            names <- colnames(data()[[target]])
            names[-c(1,3)]#e.g. 'Time' and 'eRV' for RV data
        })
        if(length(names.by.target)>1){
            nam <- Reduce(intersect,names.by.target)
        }else{
            nam <- names.by.target[[1]]
        }
        nam <- c(nam,'Window Function')
        return(unique(nam))
    })

    ns.2D <- reactive({
        if(is.null(input$per.target2)) return()
        nam <- c()
        if(length(input$per.target2)>0){
            tar <- input$per.target2
        }
        for(i in 1:length(input$per.target2)){
            names <- colnames(data()[[tar[i]]])
            names <- names[-c(1,3)]#e.g. 'Time' and 'eRV' for RV data
            names <- c(names,'Window Function')
            nam <- c(nam,names)
        }
        return(unique(nam))
    })

  ns.wt <- reactive({
      nam <- c()
      lab <- c()
      for(i in 1:length(data())){
        labs <- names <- colnames(data()[[i]])
        labs[grep('RV',names)] <- 'RV [m/s]'
        labs[grep('Time',names)] <- 'Time [JD-2400000]'
        labs[!grepl(paste(names[1:3],collapse='|'),names)] <- paste('Normalized',labs[!grepl(paste(names[1:3],collapse='|'),names)])
        nam <- c(nam,paste(names(data())[i],names,sep=':'))
        lab <- c(lab,paste(names(data())[i],labs,sep=':'))
      }
    return(list(name=nam,label=lab))
  })

    output$scatter.target <- renderUI({
        if(is.null(data())) return()
        selectizeInput("scatter.target",'Data sets',choices=names(data()),selected=names(data())[1],multiple=FALSE)
    })

    output$xs <- renderUI({
        if(is.null(data()) | is.null(input$scatter.target)) return()
        names <- ns.wt()$name
        nam <- names[grep(input$scatter.target,names)]
        if(length(nam)==0){
            nam <- names[grep(gsub('\\+','.+',input$scatter.target),names)]
        }
        selectizeInput("xs", "Choose x axis",
                       choices  = nam,
                       selected = nam[1],multiple=FALSE)
    })


    output$ys <- renderUI({
        if(is.null(data()) | is.null(input$scatter.target)) return()
        names <- ns.wt()$name
        nam <- names[grepl(input$scatter.target,names)]
        if(length(nam)==0){
            nam <- names[grep(gsub('\\+','.+',input$scatter.target),names)]
        }
        selectizeInput("ys", "Choose y axis",
                       choices  = nam,
                       selected = nam[2],multiple=FALSE)
    })

    scatterInput <- function(){
        i <- 1
        tar <- gsub(':.+','',input$xs[i])
        vars <- colnames(data()[[tar]])
        instrument <- gsub(':.+','',input$xs[i])
        indx <- which(input$xs==ns.wt()$name)
        x <- gsub('.+:','',ns.wt()$label[indx])
        indy <- which(input$ys==ns.wt()$name)
        y <- gsub('.+:','',ns.wt()$label[indy])
        varx <- data()[[instrument]][,gsub('.+:','',input$xs[i])]
        vary <- data()[[instrument]][,gsub('.+:','',input$ys[i])]
        names <- vars[1:3]
        if(!grepl(paste0(names,collapse='|'),input$xs[i])) varx <- scale(varx)
        if(!grepl(paste0(names,collapse='|'),input$ys[i])) vary <- scale(vary)
        plot(varx,vary,xlab=x,ylab=y,pch=20,cex=0.5)
        ey <- data()[[tar]][,3]
        xname <- gsub('.+:','',input$xs[i])
        yname <- gsub('.+:','',input$ys[i])
        if(xname==vars[1] & yname==vars[2] & mean(ey)>0.01*mean(vary)){
            arrows(varx,vary-ey,varx,vary+ey,length=0.03,angle=90,code=3)
        }
    }

    observeEvent(input$scatter,{
        output$sca <- renderPlot({
            isolate({
                par(mfrow=c(length(input$xs),1),cex=1,cex.axis=1.5,cex.lab=1.5,mar=c(5,5,1,1))
                scatterInput()
            })
        })
    })

    observeEvent(input$scatter,{
        output$scatter <- renderUI({
            height <- 400*ceiling(length(input$xs)/2)
            plotOutput("sca", width = "400px", height = height)
        })
    })

    output$download.scatter <- downloadHandler(
        filename = function() {
            f1 <- gsub(" ",'_',Sys.time())
            f2 <- gsub(":",'-',f1)
            paste0("scatter_",input$scatter.target,'_',f2,".pdf")
        },
        content = function(file) {
            pdf(file,4,4)
            par(mar=c(5,5,1,1))
            scatterInput()
            dev.off()
        })

    observeEvent(input$scatter,{
        output$download.scatter.button <- renderUI({
            downloadButton('download.scatter', 'Download scatter plot')
        })
    })

    output$var <- renderUI({
        if(is.null(input$per.type)) return()
        if(!any(grepl('MLP',input$per.type)) & !any(grepl('BFP',input$per.type))){
            selectInput("yvar", "Choose observables", choices  = ns.1D(),selected = ns.1D()[1],multiple=TRUE)
        }else{
            selectInput("yvar", "Choose observables",
                        choices  = ns.1D()[ns.1D()!='Window Function'],
                        selected = ns.1D()[1],multiple=TRUE)
        }
    })

  output$var2 <- renderUI({
      if(is.null(input$per.type2)) return()
      selectInput("yvar2", "Choose observables",
                        choices  = ns()[ns.2D()!='Window Function'],
                        selected = ns()[1],multiple=FALSE)
  })

  output$helpvar <- renderUI({
    if(is.null(input$yvar)) return()
    helpText("If the BFP is selected, only 'RV' is available for selection. The meaning of variables are as follows: 'all'--the periodograms of all variables,
            'RVs'--periodograms of RVs, 'Indices'-- periodograms of Indices,
             'Instrument:Variable'--individual variables")
  })

  periodogram.var <- reactive({
      if(is.null(input$yvar)) return()
    vars <- input$yvar[input$yvar!='all' & input$yvar!='RVs' & input$yvar!='Indices']
    return(unique(vars))
})

  periodogram.var2 <-  reactive({
    if(is.null(input$yvar2)) return()
    vars <- c()
    if(any(input$yvar2!='Indices' & input$yvar2!='all' & input$yvar2!='RVs')){
        vars <- c(vars,input$yvar2[input$yvar2!='all' & input$yvar2!='RVs' & input$yvar2!='Indices'])
    }
    return(unique(vars))
  })

    prange <- reactive({
        if(is.null(input$prange)) return()
        as.numeric(10^input$prange)
    })

    prange2 <- reactive({
        if(is.null(input$prange2)) return()
        as.numeric(10^input$prange2)
    })

  per.par <- reactive({
###inputs that a hidden or not-yet-rendered control leaves NULL fall back to a
###default rather than silently breaking the calculation (e.g. the AR/MA
###selectors are hidden when the GP noise model is chosen)
      iget <- function(name,default){
          v <- input[[name]]
          if(is.null(v) || length(v)==0) default else v
      }
      vals <- list(ns=ns(),ofac=input$ofac,frange=1/prange()[2:1],per.type=input$per.type,per.target=input$per.target,SigType=input$signal.type,files=input$files,SigType=input$signal.type,sequence=isTRUE(input$sequence))
      if(any(input$per.type=='MLP' | input$per.type=='BFP')){
          Nars <- Nmas <- c()#consider multiple orders
          Inds <- list()
          for(i in 1:Ntarget()){
              inds <- as.integer(iget(paste0('Inds',i),0))
              if(all(inds==0)){
                  Inds <- c(Inds,list(inds))
              }else{
                  Inds <- c(Inds,list(inds[inds!=0]))
              }
              Nmas <- c(Nmas,as.integer(iget(paste0('Nma',i),0)))
              Nars <- c(Nars,as.integer(iget(paste0('Nar',i),0)))
          }
#          vals <- c(vals,Nmas=list(Nmas),Inds=list(Inds))
          vals <- c(vals,Nmas=list(Nmas),Nars=list(Nars),Inds=list(Inds))
      }else{
          vals <- c(vals,Nmas=0,Nars=0,Inds=0)
      }
      if(vals$SigType=='stochastic'){
          vals$Nsig.max <- 1
	  vals$per.type.seq <- vals$per.type
      }else if(!vals$sequence){
          vals$Nsig.max <- 1
	  vals$per.type.seq <- vals$per.type
      }else{
	  vals$Nsig.max <- as.integer(iget('Nsig.max',1))
	  vals$per.type.seq <- iget('per.type.seq',vals$per.type[1])
      }
      if(any(input$per.type=='BFP')){
          vals$Niter <- as.numeric(iget('Niter',0))
      }else{
          vals$Niter <- 0
      }
      vals$lnBF.min <- suppressWarnings(as.numeric(iget('lnBF.min',5)))
      vals$e.prior <- list(type=iget('e.prior','beta'),
                           a=suppressWarnings(as.numeric(iget('e.beta.a',0.867))),
                           b=suppressWarnings(as.numeric(iget('e.beta.b',3.03))),
                           sigma=suppressWarnings(as.numeric(iget('e.sigma',0.1))))
      vals$Nh <- as.integer(iget('Nh',1))
      vals$noise.model <- if(is.null(input$noise.model)) 'ARMA' else input$noise.model
      vals$gp.Prot <- if(is.null(input$gp.Prot)) NA else suppressWarnings(as.numeric(input$gp.Prot))
      vals$gp.tau <- if(is.null(input$gp.tau)) NA else suppressWarnings(as.numeric(input$gp.tau))
      return(vals)
  })

  per.par2 <- reactive({
      vals <- list(ns=ns(),ofac=input$ofac2,frange=1/prange2()[2:1],per.type=input$per.type2,per.target=input$per.target2,files=input$files,Niter=as.numeric(input$Niter))#SigType=input$signal.type
      if(any(input$per.type2=='MLP'|input$per.type2=='BFP')){
          Nmas <- c()
          Nars <- c()
          Inds <- list()
          iget2 <- function(name,default){
              v <- input[[name]]
              if(is.null(v) || length(v)==0) default else v
          }
          for(i in 1:Ntarget2()){
              inds <- as.integer(iget2(paste0('Inds2.',i),0))
              if(all(inds==0)){
                  Inds <- c(Inds,list(inds))
              }else{
                  Inds <- c(Inds,list(inds[inds!=0]))
              }
              Nmas <- c(Nmas,as.integer(iget2(paste0('Nma2.',i),0)))
              Nars <- c(Nars,as.integer(iget2(paste0('Nar2.',i),0)))
          }
          vals <- c(vals,Nmas=list(Nmas),Nars=list(Nars),Inds=list(Inds))
      }else{
          vals <- c(vals,Nmas=0,Nars=0,Inds=0)
      }
      vals <- c(vals,Dt=signif(tspan()*as.numeric(input$Dt),3),Nbin=as.integer(input$Nbin),alpha=as.integer(input$alpha),scale=input$scale,pmin.zoom=input$range.zoom[1],pmax.zoom=input$range.zoom[2],show.signal=input$show.signal)
      vals$use.fit <- isTRUE(input$use.fit)
      vals$adaptive <- isTRUE(input$adaptive2)
      vals$noise.model <- if(is.null(input$noise.model2)) 'ARMA' else input$noise.model2
      vals$gp.Prot <- if(is.null(input$gp.Prot2)) NA else suppressWarnings(as.numeric(input$gp.Prot2))
      vals$gp.tau <- if(is.null(input$gp.tau2)) NA else suppressWarnings(as.numeric(input$gp.tau2))
      return(vals)
  })

    model.selection <- eventReactive(input$compare,{
                                        #      instrument <- instr()[input]#gsub(':.+','',ns()[1])
        tab <- data()[[input$comp.target]]
        if(!is.null(input$NI0)){
            Nbasic <- as.integer(input$NI0)
        }else{
            Nbasic <- 0
        }
        if(!is.null(input$proxy.type)){
            if(input$proxy.type=='group'){
                groups <- input$NI.group
                proxy.type <- 'group'
                ni <- NI.max()[input$comp.target]
            }else if(input$proxy.type=='man'){
                groups <- list()
                cat('names(input)=',names(input),'\n')
                for(i in 1:as.integer(input$Nman)){
                    inds <- as.integer(input[[paste0('NI.man',i)]])
                    if(!all(inds==0)){
                        inds <- inds[inds>0]
                    }
                    groups[[i]] <- inds
                }
                cat('names(groups)=',names(groups),'\n')
                cat('length(groups)=',length(groups),'\n')
                proxy.type <- 'man'
                ni <- NI.max()[input$comp.target]
            }else{
                groups <- NULL
                proxy.type <- 'cum'
                ni <- as.integer(input$ni.max)
            }
        }else{
            groups <- NULL
            proxy.type <- 'cum'
            ni <- 0
        }
        out <- calcBF(data=tab,Nbasic=Nbasic,
                      proxy.type=proxy.type,
                      Nma.max=as.integer(input$Nma.max),
		      Nar.max=as.integer(input$Nar.max),
                      groups=groups,Nproxy=ni,progress=TRUE)

        logBF <- out$lnBF
        row.names <- c()
        for(j in 1:length(out$Inds)){
            if(all(out$Inds[[j]]==0)){
                row.names <- c(row.names,'no proxy')
            }else{
                row.names <- c(row.names,paste0('proxies: ',paste(out$Inds[[j]],collapse=',')))
            }
        }
        rownames(logBF) <- row.names

        logBF.download <- logBF
#        colnames(logBF.download) <- paste0('MA',out$Nmas)

        rnames <- c()
        for(j in 1:nrow(logBF)){
            rnames <- c(rnames,paste0('proxy',paste(out$Inds[[j]],collapse='-')))
        }
#        rownames(logBF.download) <- NULL#rnames
        return(list(logBF=logBF,out=out,logBF.download=logBF.download))
    })

  observeEvent(input$compare,{
      output$BFtab <- renderUI({
          output$table <- renderTable({model.selection()$logBF},digits=1,caption = "Logarithmic BIC-estimated Bayes factor",rownames=TRUE,colnames=TRUE,
                                      caption.placement = getOption("xtable.caption.placement", "top"),
                                      caption.width = getOption("xtable.caption.width", NULL))
          tableOutput('table')
      })
  })

  output$download.logBF <- downloadHandler(
      filename = function(){
          f1 <- gsub(" ",'_',Sys.time())
          f2 <- gsub(":",'-',f1)
#          paste('logBF_',input$comp.target,'_', f2, '.txt', sep='')
          paste('logBF_',input$comp.target,'_', f2, '.csv', sep='')
      },
      content = function(file) {
          write.csv2(round(model.selection()$logBF.download,digit=1), file,quote=FALSE)
      }
  )

  output$download.logBF.table <- renderUI({
      if(is.null(model.selection())) return()
#      downloadLink('download.logBF', 'Download the Bayes Factor table')
      downloadButton('download.logBF', 'Download the Bayes Factor table')
  })

  output$optNoise <- renderUI({
      if(is.null(input$compare)) return()
      if(input$compare>0){
          output$noise.opt <- renderText({
              Nma.opt <- model.selection()$out$Nma.opt
              Nar.opt <- model.selection()$out$Nar.opt
              Inds.opt <- model.selection()$out$Inds.opt
#              Inds.opt <- model.selection()$out$NI.opt
              if(Nma.opt==0 & Nar.opt==0){
                  t1 <- 'white noise'
              }else{
#                  t1 <- paste0('MA(',Nma.opt,')')
                  t1 <- paste0('ARMA(',Nar.opt,',',Nma.opt,')')
              }
              if(all(Inds.opt==0)){
                  t2 <- 'Optimal proxies: no proxy'
              }else{
                  t2 <- paste0('Optimal proxies: ',paste(Inds.opt,collapse=','))
              }
              text1 <- paste0('Optimal noise model: ',t1)
              text2 <- t2
              HTML(paste(text1, text2, sep = '<br/>'))
          })
          htmlOutput('noise.opt')
      }
  })

output$color <- renderUI({
    if(is.null(MP.data()) | is.null(Ntarget2())) return()
    if(Ntarget2()>1){
            ts <- c()
            for(j in 1:Ntarget2()){
                cols <- c('black','red','blue','green','orange','brown','cyan','pink')
                ts <- c(ts,paste0(cols[j],': Noise-subtracted ',input$per.target2[j]))
            }
            out <- paste(ts,collapse='<br/><br/>')
            h5(HTML(paste0('<br/>',out)))
#        htmlOutput('encode')
    }
})

  Nper <- eventReactive(input$plot1D,{
    if(is.null(input$per.type)) return()
    Nvar <- length(periodogram.var())
    Nplots <- 0
    pars <- per.par()
    Nplots <- Nplots+length(input$per.type)
    Nplots <- max(1,Nplots)*Nvar
    return(Nplots)
  })

  Nper2 <- eventReactive(input$plot2D,{
    if(is.null(input$per.type2)) return()
    Nvar <- length(periodogram.var2())
    Nplots <- 0
    pars <- per.par2()
    Nplots <- Nplots+length(input$per.type2)
    Nplots <- max(1,Nplots)*Nvar
    return(Nplots)
  })

  tvper <- reactive({
    if(is.null(data())) return()
    logic <- c()
    for(j in 1:length(data())){
      trv <- data()[[j]][,1]
      if(length(trv)>100 & (max(trv)-min(trv))>1000){
        logic <- c(logic,TRUE)
      }else{
        logic <- c(logic,FALSE)
      }
    }
    return(logic)
  })

    tspan <- reactive({
        if(is.null(data()) | is.null(input$per.target2)) return()
        ts.min <- ts.max <- c()
        for(i in 1:length(input$per.target2)){
            tmp <- data()[[input$per.target2[i]]][,1]
            ts.min <- c(ts.min,min(tmp))
            ts.max <- c(ts.max,max(tmp))
        }
        tmin <- min(ts.min)
        tmax <- max(ts.max)
        dt <- tmax-tmin
        return(dt)
    })

  output$Dt <- renderUI({
#      cat('input$plot2D=',input$plot2D,'\n')
      if(!is.null(data())){
          sliderInput('Dt','Moving time window [in unit of the whole time span]',
                  min=0.01,max=0.99,value=0.5,step=0.01)
#          sliderInput("Dt", "Moving time window", min = 100, max = ,value=min(1000,round(tmax-tmin)),step=100)
      }
  })

   output$textDt <- renderUI({
        if(is.null(input$Dt)) return()
        helpText(paste0("The time window is ",signif(tspan()*as.numeric(input$Dt),3)," time unit. The user should adjust it to guarantee the existence of a few data points in the time window for each moving step."))
    })

  output$Nbin <- renderUI({
      if(!is.null(data())){
          selectizeInput('Nbin','Number of moving steps',
                  choices=c(2,5,10,20,50,100,200,500),selected=10,multiple=FALSE)
#          sliderInput("Nbin", "Number of moving steps", min = 5, max = 500,value=10)
      }
  })

  output$alpha <- renderUI({
      if(!is.null(data())){
          sliderInput('alpha','Truncate the color bar to optimize visualization', min = 0, max = 10,value=5,step=0.1)
      }
  })

  output$zoom <- renderUI({
      if(is.null(input$prange2)) return()
      pr <- 10^as.numeric(input$prange2)
      pmin <- pr[1]-pr[1]%%0.1
      pmax <- pr[2]-pr[2]%%0.1
      plow <- pmin
      pup <- pmin+0.1*(pmax-pmin)
      pup <- pup-pup%%0.1
      sliderInput('range.zoom','Zoom-in period range', min = pmin, max = pmax,value=c(plow,pup),step=0.1)
  })

###get BFP power spectrum
    data1D <- eventReactive(input$plot1D,{
        stack <- NULL
        tryCatch(withCallingHandlers({
###use calc.1Dper() from functions.R to calculate periodogram
            calc.1Dper(Nmax.plots, periodogram.var(),per.par(),data())
        },error=function(e){
            stack <<- sys.calls()
        }),error=function(e){
###a silent failure looks like a dead button: tell the user what went wrong and
###write a full diagnostic dump so the failure can be reproduced exactly
            dump.file <- file.path(getwd(),'agatha_last_error.txt')
            try({
                con <- file(dump.file,'w')
                writeLines(c(paste('time:',format(Sys.time())),
                             paste('git:',tryCatch(system('git rev-parse --short HEAD',intern=TRUE),error=function(x) 'unknown')),
                             paste('error:',conditionMessage(e)),
                             '--- call stack ---',
                             if(is.null(stack)) '(unavailable)' else unlist(lapply(stack,function(cl) paste(deparse(cl)[1]))),
                             '--- observables ---',
                             paste(deparse(periodogram.var()),collapse='\n'),
                             '--- settings (per.par) ---',
                             paste(deparse(per.par()),collapse='\n')),con)
                close(con)
            },silent=TRUE)
            showNotification(paste('Periodogram calculation failed:',conditionMessage(e),
                                   '- diagnostics written to',dump.file),type='error',duration=NULL)
            NULL
        })
    })

    save.data <- function(li,f,format='normal'){
        fs <- c()
        if(length(li)==1){
            if(format=='normal'){
                write.table(li[[1]], file=f,quote=FALSE,row.names=FALSE)
            }else{
                tab <- li[[1]]
                if(is.null(dim(tab))){
                    write.table(t(tab), f,quote=FALSE,row.names=FALSE)
                }else{
                    write.csv(tab, f,quote=FALSE,row.names=TRUE)
                }
            }
            fs <- c(fs,f)
        }else{
            for(j in 1:length(li)){
                ff <- paste0(names(li)[j],'_',f)
                if(format=='normal'){
                    write.table(li[[j]], file=ff,quote=FALSE,row.names=FALSE)
                }else{
                    tab <- li[[j]]
                    if(is.null(dim(tab))){
                        write.table(t(tab), file=ff,quote=FALSE,row.names=FALSE)
                    }else{
                        write.csv(tab, file=ff,quote=FALSE,row.names=TRUE)
                    }
                }
                fs <- c(fs,ff)
            }
        }
        fs
    }

    output$all.data <- downloadHandler(
        filename = function() {
            paste0(data1D()$fname,'.zip')
        },
        content = function(file) {
        tmpdir <- tempdir()
        setwd(tempdir())
        fs <- c()
        fs <- c(fs,save.data(data1D()$per.list,paste0(data1D()$fname,'_Periodogram.txt')))
        fs <- c(fs,save.data(data1D()$phase.list,paste0(data1D()$fname,'_PhaseData.txt')))
        fs <- c(fs,save.data(data1D()$mc.list,paste0(data1D()$fname,'_MCposterior.txt')))
        fs <- c(fs,save.data(data1D()$sim.list,paste0(data1D()$fname,'_SimFit.txt')))
        fs <- c(fs,save.data(data1D()$par.list,paste0(data1D()$fname,'_OptPar.csv'),format='csv'))
#        save(list=ls(all=TRUE),file='test0.Robj')
        zip(zipfile=file, files=fs)
      },
      contentType = "application/zip"
    )

    output$phase1D.data <- downloadHandler(
        filename = function() {
            if(length(data1D()$phase.list)>1){
                paste0(data1D()$fname,'_PhaseData.tar')
            }else{
                paste0(data1D()$fname,'_PhaseData.txt')
            }
        },
        content = function(file) {
            save.data(data1D()$phase.list,file)
        }
    )

    output$sim1D.data <- downloadHandler(
        filename = function() {
            if(length(data1D()$sim.list)>1){
                paste0(data1D()$fname,'_SimFit.tar')
            }else{
                paste0(data1D()$fname,'_SimFit.txt')
            }
        },
        content = function(file) {
            save.data(data1D()$sim.list,file)
        }
    )

    output$par1D.data <- downloadHandler(
        filename = function() {
            if(length(data1D()$par.list)>1){
                paste0(data1D()$fname,'_OptPar.tar')
            }else{
                paste0(data1D()$fname,'_OptPar.txt')
            }
        },
        content = function(file) {
            save.data(data1D()$par.list,file,format='csv')
        }
    )

    output$download.all.data <- renderUI({
        if(is.null(data1D())) return()
        downloadButton('all.data', 'Download data behind the figures')
    })

    output$download.per1D.data <- renderUI({
        if(is.null(data1D())) return()
        downloadButton('per1D.data', 'Download data for periodograms')
    })

    output$download.phase1D.data <- renderUI({
        if(is.null(data1D())) return()
        downloadButton('phase1D.data', 'Download data for phase plots')
    })

    output$download.sim1D.data <- renderUI({
        if(is.null(data1D())) return()
        downloadButton('sim1D.data', 'Download data for model prediction')
    })

    output$download.par1D.data <- renderUI({
        if(is.null(data1D())) return()
        downloadButton('par1D.data', 'Download optimal parameter values')
    })

    output$per1D.figure <- downloadHandler(
        filename = function() {
            paste0(data1D()$fname,'_Periodogram.pdf')
        },
      content = function(file) {
        pdf(file,8,8)
        per1D.plot(data1D()$per.list,data1D()$tits,data1D()$pers,data1D()$levels,ylabs=data1D()$ylabs,download=TRUE,SigType=input$signal.type,par.list=data1D()$par.list)
        phase1D.plot(data1D()$phase.list,data1D()$sim.list,data1D()$tits,download=TRUE,repar=FALSE,par.list=data1D()$par.list,model.comp=data1D()$model.comp)
        dev.off()
      })

    output$phase1D.figure <- downloadHandler(
        filename = function() {
            paste0(data1D()$fname,'_Phase.pdf')
        },
      content = function(file){
        pdf(file,8,8)
        phase1D.plot(data1D()$phase.list,data1D()$sim.list,data1D()$tits,download=TRUE,par.list=data1D()$par.list,model.comp=data1D()$model.comp)
        dev.off()
      })


###individual, publication-quality figures
    single.plots <- reactive({
        if(is.null(data1D())) return(NULL)
        list.single.plots(data1D())
    })

    output$plot.single <- renderUI({
        if(is.null(single.plots())) return()
        tagList(
            h5('Download an individual figure'),
            selectizeInput('single.plot','Figure',choices=single.plots()$label,multiple=FALSE),
            radioButtons('single.format','Format',choices=c('PDF (vector)'='pdf','PNG'='png','JPG'='jpg'),selected='pdf',inline=TRUE),
            fluidRow(column(4,numericInput('single.width','Width [in]',value=6,min=2,max=20,step=0.5)),
                     column(4,numericInput('single.height','Height [in]',value=4.5,min=2,max=20,step=0.5)),
                     column(4,numericInput('single.dpi','DPI',value=300,min=72,max=1200,step=50))),
            downloadButton('single.figure','Download this figure')
        )
    })

    output$single.figure <- downloadHandler(
        filename = function(){
            sp <- single.plots(); k <- which(sp$label==input$single.plot)[1]
            tag <- gsub('[^A-Za-z0-9]+','_',sp$label[k])
            paste0(data1D()$fname,'_',tag,'.',input$single.format)
        },
        content = function(file){
            sp <- single.plots(); k <- which(sp$label==input$single.plot)[1]
            save.single.plot(file,format=input$single.format,width=input$single.width,height=input$single.height,dpi=input$single.dpi,
                             plot1D.single(data1D(),kind=sp$kind[k],ypar=sp$ypar[k],index=sp$index[k],SigType=input$signal.type))
        })

    output$download.per1D.plot <- renderUI({
        if(is.null(data1D())) return()
        downloadButton('per1D.figure', 'Download all figures (PDF)')
#        if(input$down.type=='all'){
#            downloadButton('per1D.figure', 'Download periodograms')
#        }else{
#            downloadButton('per1D.single', 'Download periodograms')
#        }
    })


    output$download.phase1D.plot <- renderUI({
        if(is.null(data1D())) return()
        downloadButton('phase1D.figure', 'Download model fit and residual')
    })

    output$help.per1D <- renderUI({
        if(is.null(data1D())) return()
        helpText("The column names are 'P' and 'type:Observable:power', where 'name' is the periodogram type, 'P' is period, and 'power' is the periodogram power which could be logarithmic marginalized likelihood (logML; for MLP and BGLS) or Bayes factor (logBF; for BFP) or power (for other periodograms). ")
    })

    output$combined <- renderPlot({
        if(is.null(data1D())) return()
        combined.plot(data1D()$per.list,data1D()$phase.list,data1D()$sim.list,data1D()$tits,data1D()$pers,data1D()$levels,data1D()$ylabs,SigType=input$signal.type,par.list=data1D()$par.list,model.comp=data1D()$model.comp)
    })

    output$plot.1Dcombined <- renderUI({
        plotOutput("combined", width = "750px", height = 400*ceiling(Nmax.plots/2))
    })

###the signal-only series of the 1D fit (offsets, trend and red noise removed)
###as the input of the moving periodogram, once a 1D fit exists for the same sets
    fit1D.sets <- function(){
        d <- tryCatch(data1D(),error=function(e) NULL)
        if(is.null(d) || length(d$phase.list)==0) return(NULL)
        sets <- attr(d$phase.list[[1]],'sets')
        if(is.null(sets)) return(NULL)
        sets
    }
    output$use.fit2 <- renderUI({
        if(is.null(input$per.target2)) return()
        sets <- fit1D.sets()
        if(is.null(sets)){
            return(helpText('Compute a 1D periodogram (with or without MCMC) for these data sets to enable the moving periodogram of the signal-only data: the sets combined after removing the offsets, the trend and the red noise of that fit, which extends the baseline for testing the time consistency of long-period signals.'))
        }
        ok <- setequal(sets,input$per.target2)
        tagList(
            checkboxInput('use.fit','Use the signal-only data of the 1D fit (offsets, trend and red noise removed)',value=FALSE),
            helpText(if(ok){
                'The data sets are combined after subtracting the per-set offsets, the trend and the red-noise model of the 1D fit (the MCMC solution when MCMC was run), so only the potential signals and white noise remain. The noise settings above are then ignored.'
            }else{
                paste0('The 1D fit was made for ',paste(sets,collapse=', '),'; select the same data sets here to use it.')
            })
        )
    })

  MP.data <- eventReactive(input$data.update,{
      tryCatch(per2D.data(periodogram.var2(),per.par2(),data(),fit1D=if(isTRUE(input$use.fit)) data1D() else NULL),
               error=function(e){
                   showNotification(paste('2D periodogram calculation failed:',conditionMessage(e)),type='error',duration=NULL)
                   NULL
               })
  })

    output$MP.data <- downloadHandler(
        filename = function() {
#            f1 <- gsub(" ",'_',Sys.time())
#            f2 <- gsub(":",'-',f1)
#            paste('periodogram2D_', f2, '.txt', sep='')
            paste0(MP.data()$fname,'.txt')
        },
        content = function(file) {
            tmp <- MP.data()
            if(nrow(tmp$zz)==length(tmp$xx)){
                tab <- cbind(tmp$xx,tmp$zz)
                tab <- t(rbind(c(NA,tmp$yy),tab))
            }else{
                tab <- rbind(tmp$xx,tmp$zz)
                tab <- cbind(c(NA,tmp$yy),tab)
            }
            write.table(tab, file,quote=FALSE,row.names=FALSE,col.names=FALSE)#FALSE,col.names=FALSE
        }
    )

    output$download.MP.data <- renderUI({
        if(is.null(MP.data())) return()
        downloadButton('MP.data', 'Download data of 2D periodogram')
    })

###################################################################
####N-body stability tab
    nbody.avail <- reactive({
        input$nbody.python
        nbody.python(if(is.null(input$nbody.python) || !nzchar(input$nbody.python)) 'python3' else input$nbody.python)
    })

    output$nbody.engine <- renderUI({
        ver <- nbody.avail()
        tagList(
            radioButtons('nbody.engine','Integrator',
                         choices=if(!is.null(ver)) c('REBOUND WHFast (Python, long-term)'='rebound','Leapfrog in R (short-term)'='R') else c('Leapfrog in R (short-term)'='R'),
                         selected=if(!is.null(ver)) 'rebound' else 'R'),
            helpText(if(!is.null(ver)) paste0('REBOUND ',ver,' found. Integrations of ~10 Myr are only practical with it.') else
                     'REBOUND was not found for this Python interpreter (install with: python3 -m pip install --user "rebound<5"); only the short-term R integrator is available.')
        )
    })

    nbody.orbits <- reactive({
###the planets of the current 1D solution (RV) for the chosen stellar mass and inclination
        d <- tryCatch(data1D(),error=function(e) NULL)
        if(is.null(d) || is.null(d$par.list) || length(d$par.list)==0) return(NULL)
        ypar <- names(d$par.list)[1]
        smp <- tryCatch(nbody.samples(d,ypar,N=0),error=function(e) NULL)
        if(is.null(smp)) return(NULL)
        list(ypar=ypar,nominal=smp$nominal,orb=rv2orbits(smp$nominal,Mstar=as.numeric(input$nbody.Mstar),inc=as.numeric(input$nbody.inc)))
    })

    output$nbody.source <- renderUI({
        o <- nbody.orbits()
        if(is.null(o) || is.null(o$orb)) return(helpText('Compute a 1D periodogram with a Keplerian (or circular) signal first; its fitted parameters (the MAP of the MCMC when it ran) define the planets.'))
        n <- nrow(o$orb)
        helpText(paste0(n,' planet',if(n>1) 's' else '',' from the 1D fit of ',o$ypar,': P = ',paste(signif(o$orb$P,4),collapse=', '),' d; ',
                        'a = ',paste(signif(o$orb$a,3),collapse=', '),' AU; m sin i = ',paste(signif(o$orb$msini,3),collapse=', '),' Mjup. ',
                        'Planet masses follow from the inclination (coplanar orbits).'))
    })

    output$nbody.estimate <- renderUI({
        o <- nbody.orbits()
        if(is.null(o) || is.null(o$orb)) return()
        Pmin <- min(o$orb$P)/365.25
        steps <- as.numeric(input$nbody.tmax)/(Pmin/as.numeric(input$nbody.steps))
        nsys <- 1+(if(identical(input$nbody.mode,'mc')) as.integer(input$nbody.N) else 0)
        rate <- if(identical(input$nbody.engine,'rebound')) 5e6 else 5e3
        helpText(paste0('About ',format(signif(steps,2),big.mark=',',scientific=TRUE),' steps per system (',nsys,' system',if(nsys>1) 's' else '','); roughly ',
                        signif(steps*nsys/rate/60,2),' minutes at ',format(rate,big.mark=','),' steps per second.'))
    })

    output$nbody.table <- renderTable({
        o <- nbody.orbits(); if(is.null(o) || is.null(o$orb)) return(NULL)
        tab <- o$orb[,c('signal','P','K','e','omega','Mo','a','msini','m')]
        colnames(tab) <- c('Signal','P [d]','K [m/s]','e','omega [rad]','Mo [rad]','a [AU]','m sin i [Mjup]','m [Msun]')
        tab
    },digits=4,caption='Planets of the orbital solution',caption.placement='top')

    output$nbody.amd <- renderTable({
        o <- nbody.orbits(); if(is.null(o) || is.null(o$orb) || nrow(o$orb)<2) return(NULL)
        tab <- amd.stability(o$orb,Mstar=as.numeric(input$nbody.Mstar))
        tab$AMD.stable <- ifelse(tab$AMD.stable,'yes','no'); tab$Hill.stable <- ifelse(tab$Hill.stable,'yes','no')
        colnames(tab) <- c('Pair','a inner [AU]','a outer [AU]','Relative AMD','Critical (collision)','Critical (Hill)','AMD ratio','Separation [mutual Hill radii]','AMD-stable','Hill-stable (> 2 sqrt 3)')
        tab
    },digits=3,caption='Analytical criteria for adjacent pairs (nominal orbits)',caption.placement='top')

###the integration runs in a separate process (nbody.launch), so the app
###stays responsive, the progress is read from the job's progress file, and
###the run can be stopped
    nbody.rv <- reactiveValues(job=NULL,result=NULL,message=NULL,bar=NULL)

    observeEvent(input$nbody.go,{
        if(!is.null(nbody.rv$job)){ showNotification('An integration is already running; stop it first.',type='warning'); return() }
        o <- nbody.orbits()
        if(is.null(o) || is.null(o$orb)){ showNotification('No orbital solution: compute a 1D periodogram first.',type='error'); return() }
        d <- data1D()
        Mstar <- as.numeric(input$nbody.Mstar); inc <- as.numeric(input$nbody.inc)
        N <- if(identical(input$nbody.mode,'mc')) as.integer(input$nbody.N) else 0
        smp <- nbody.samples(d,o$ypar,N=N,seed=1)
        pars <- c(list(smp$nominal),smp$samples)
        orbit.list <- lapply(pars,function(p) rv2orbits(p,Mstar=Mstar,inc=inc))
        settings <- list(Mstar=Mstar,tmax=as.numeric(input$nbody.tmax),steps.per.orbit=as.numeric(input$nbody.steps),
                         Nout=as.integer(input$nbody.Nout),escape.factor=as.numeric(input$nbody.escape),encounter.hill=as.numeric(input$nbody.hill),
                         engine=if(is.null(input$nbody.engine)) 'R' else input$nbody.engine,
                         python=if(nzchar(input$nbody.python)) input$nbody.python else 'python3',
                         a.tol=as.numeric(input$nbody.atol),method=smp$method)
        nbody.rv$result <- NULL; nbody.rv$message <- NULL
        job <- tryCatch(nbody.launch(orbit.list,settings),error=function(e){ showNotification(paste('could not start the integration:',conditionMessage(e)),type='error'); NULL })
        if(is.null(job)) return()
        bar <- shiny::Progress$new(); bar$set(message='Integrating the orbits',value=0,detail='starting')
        nbody.rv$bar <- bar
        nbody.rv$job <- job
    })

    observeEvent(input$nbody.stop,{
        job <- nbody.rv$job
        if(is.null(job)){ showNotification('No integration is running.',type='message'); return() }
        nbody.stop(job)
        if(!is.null(nbody.rv$bar)){ nbody.rv$bar$close(); nbody.rv$bar <- NULL }
        nbody.rv$message <- paste0('Integration stopped by the user after ',format(round(as.numeric(difftime(Sys.time(),job$started,units='secs')))),' s.')
        nbody.rv$job <- NULL
    })

    observe({
        job <- nbody.rv$job
        if(is.null(job)) return()
        invalidateLater(1000)
        st <- nbody.status(job)
        if(!st$done){
            if(!is.null(nbody.rv$bar)){
                frac <- if(st$n>0) (max(st$i,1)-1+st$frac)/st$n else 0
                nbody.rv$bar$set(value=frac,detail=paste0('system ',max(st$i,1),' of ',st$n,': ',round(100*st$frac),'% of the integration time; ',
                                                          format(round(as.numeric(difftime(Sys.time(),job$started,units='secs')))),' s elapsed'))
            }
            return()
        }
        if(!is.null(nbody.rv$bar)){ nbody.rv$bar$close(); nbody.rv$bar <- NULL }
        if(!is.null(st$error) || is.null(st$result)){
            nbody.rv$message <- paste('N-body integration failed:',if(is.null(st$error)) 'no result was written' else st$error)
            showNotification(nbody.rv$message,type='error',duration=NULL)
        }else{
            nbody.rv$result <- st$result
            nbody.rv$message <- paste0('Finished in ',format(round(as.numeric(difftime(Sys.time(),job$started,units='secs')))),' s.')
        }
        nbody.rv$job <- NULL
    })

    nbody.result <- reactive({ nbody.rv$result })

    output$nbody.progress <- renderUI({
        if(!is.null(nbody.rv$job)) return(helpText(HTML('<b>Integration running</b> - see the progress bar; press Stop to abort.')))
        if(!is.null(nbody.rv$message)) helpText(nbody.rv$message)
    })

    output$nbody.summary <- renderUI({
        r <- nbody.result(); if(is.null(r)) return()
        cls <- r$class
        nom <- cls[1,]
        tagList(h4(paste0('Nominal solution: ',nom$status,if(!nom$stable) paste0(' at ',signif(nom$t.end,3),' yr') else paste0(' over ',format(r$tmax,big.mark=','),' yr'))),
                helpText(paste0(r$run$engine,'; ',r$method,'. Energy error of the nominal run: ',signif(nom$dE,2),
                                '. A system is unstable at an escape, a close encounter within the chosen number of mutual Hill radii, a drift of a semi-major axis beyond the tolerance, or e > 0.9.')),
                if(nrow(cls)>1) h4(paste0(sum(cls$stable[-1]),' of ',nrow(cls)-1,' sampled systems stable')))
    })

    output$nbody.results <- renderTable({
        r <- nbody.result(); if(is.null(r)) return(NULL)
        cls <- r$class
        cls$system <- c('nominal',if(nrow(cls)>1) paste('sample',seq_len(nrow(cls)-1)))
        cls$stable <- ifelse(cls$stable,'yes','no')
        colnames(cls) <- c('System','Status','End time [yr]','max |da/a|','max e','Energy error','Stable')
        cls
    },digits=3,caption='Integrated systems',caption.placement='top')

    output$nbody.plot <- renderPlot({
        r <- nbody.result(); if(is.null(r)) return()
        nbody.plot(r$run,r$orbit.list,system=1,main=paste0('Nominal solution (',r$run$engine,')'))
    })

    output$nbody.mcplot <- renderPlot({
        r <- nbody.result(); if(is.null(r) || nrow(r$class)<2) return()
        nbody.mc.plot(r$class[-1,],r$tmax,main=r$method)
    })

    output$plot.nbody <- renderUI({
        r <- nbody.result(); if(is.null(r)) return()
        tagList(plotOutput('nbody.plot',width='700px',height='650px'),
                if(nrow(r$class)>1) plotOutput('nbody.mcplot',width='700px',height='400px'))
    })

    output$nbody.download <- downloadHandler(
        filename = function() 'nbody_stability.csv',
        content = function(file){ r <- nbody.result(); write.csv(r$class,file,row.names=FALSE) })
    output$nbody.tracks <- downloadHandler(
        filename = function() 'nbody_tracks.csv',
        content = function(file){ r <- nbody.result(); write.csv(r$run$tracks,file,row.names=FALSE) })
    output$download.nbody <- renderUI({
        if(is.null(nbody.result())) return()
        tagList(downloadButton('nbody.download','Download the stability table (CSV)'),downloadButton('nbody.tracks','Download a(t), e(t) of every run (CSV)'))
    })

###################################################################
####Signal diagnosis tab
    output$per.target3 <- renderUI({
        if(is.null(data())) return(helpText('Upload or select the data sets first.'))
        selectizeInput('per.target3','RV data sets to diagnose',choices=names(data()),selected=names(data()),multiple=TRUE)
    })

    diag <- eventReactive(input$diagnose,{
        if(is.null(data()) || is.null(input$per.target3) || length(input$per.target3)==0){
            showNotification('Select at least one data set.',type='error'); return(NULL)
        }
        sets <- input$per.target3
        noise <- if(is.null(input$noise3)) 'W' else input$noise3
        nsteps <- length(noise)+(if(isTRUE(input$individual3) && length(sets)>1) length(sets) else 0)+
                  (if(isTRUE(input$proxies3)) length(sets) else 0)+(if(isTRUE(input$moving3)) 1 else 0)
        prange <- 10^input$prange3
        tryCatch(
            withProgress(message='Diagnosing the signals',value=0,{
                diagnose.signals(data(),sets,noise.models=noise,Nsig.max=as.integer(input$Nsig3),ofac=as.numeric(input$ofac3),
                                 frange=c(1/prange[2],1/prange[1]),lnBF.min=as.numeric(input$lnBF.min3),SigType=input$SigType3,
                                 Nh=if(input$SigType3=='kepler') 2 else 1,Ncores=1,Nwin=as.integer(input$Nwin3),
                                 individual=isTRUE(input$individual3),proxies=isTRUE(input$proxies3),moving=isTRUE(input$moving3),
                                 gp.Prot=suppressWarnings(as.numeric(input$gp.Prot3)),gp.tau=suppressWarnings(as.numeric(input$gp.tau3)),
                                 adaptive=isTRUE(input$adaptive3),
                                 progress=function(msg) incProgress(1/max(nsteps,1),detail=msg))
            }),
            error=function(e){
                showNotification(paste('Signal diagnosis failed:',conditionMessage(e)),type='error',duration=NULL)
                NULL
            })
    })

    output$diag.text <- renderUI({
        d <- diag()
        if(is.null(d)) return()
        n <- length(d$periods)
        tagList(h4(paste0('Signals in ',paste(set.label(d$sets),collapse=' + '),' (',d$noise.models[1],' noise): ',n,' accepted')),
                helpText(paste0('Accepted periods: ',if(n>0) paste(signif(d$periods,5),collapse=', ') else 'none',
                                ' d. The table applies the criteria of Feng et al. (2020): ln(BF) > ',d$lnBF.min,
                                ' under every noise model, no activity index peaking within 10 per cent of the period, and the signal above the threshold in most windows of the moving periodogram (the time-consistency test is not applied to periods longer than half the time span).')))
    })

    output$diag.summary <- renderTable({
        d <- diag(); if(is.null(d)) return(NULL)
        diagnosis.table(d)
    },digits=1,caption='Diagnosis of the accepted signals',caption.placement='top')

    output$diag.plot <- renderPlot({
        d <- diag(); if(is.null(d)) return()
        diagnosis.plot(d,ncol=as.integer(input$ncol3),Prot=suppressWarnings(as.numeric(input$Prot3)))
    })

    output$plot.diag <- renderUI({
        d <- diag(); if(is.null(d)) return()
        nc <- as.integer(input$ncol3); lay <- diagnosis.layout(d,nc)
        plotOutput('diag.plot',width='900px',height=paste0(round(240*sum(lay$heights)),'px'))
    })

    observeEvent(diag(),{
        d <- diag(); if(is.null(d) || length(d$moving)==0) return()
        for(k in seq_along(d$moving)){
            local({
                kk <- k
                output[[paste0('diag.mp',kk)]] <- renderPlot({
                    if(is.null(d$moving[[kk]])) return()
                    plotMP(d$moving[[kk]],d$moving.par)
                })
            })
        }
    })

    output$plot.diag.mp <- renderUI({
        d <- diag(); if(is.null(d) || length(d$moving)==0) return()
        items <- lapply(seq_along(d$moving),function(k){
            if(is.null(d$moving[[k]])) return(NULL)
            tagList(h5(paste0('Signal ',k,' (P = ',signif(d$periods[k],5),' d)')),plotOutput(paste0('diag.mp',k),width='600px',height='600px'))
        })
        tagList(h4('Moving periodograms of the signal-only combined series'),
                helpText('For each accepted signal: the data sets combined after removing the per-set offsets, trend and red noise of the combined fit, and the other accepted signals subtracted; a real signal stays at the same period in every window.'),
                items)
    })

    output$diag.figure <- downloadHandler(
        filename = function() paste0(paste(diag()$sets,collapse='_'),'_diagnosis.pdf'),
        content = function(file){
            d <- diag(); nc <- as.integer(input$ncol3)
            lay <- diagnosis.layout(d,nc)
            pdf(file,width=3.3*nc,height=2.9*sum(lay$heights))
            diagnosis.plot(d,ncol=nc,Prot=suppressWarnings(as.numeric(input$Prot3)))
            for(mv in d$moving) if(!is.null(mv)) plotMP(mv,d$moving.par)
            dev.off()
        })

    output$diag.table <- downloadHandler(
        filename = function() paste0(paste(diag()$sets,collapse='_'),'_diagnosis.csv'),
        content = function(file) write.csv(diagnosis.table(diag()),file,row.names=FALSE)
    )

    output$download.diag <- renderUI({
        if(is.null(diag())) return()
        tagList(downloadButton('diag.figure','Download the diagnosis figure (PDF)'),
                downloadButton('diag.table','Download the diagnosis table (CSV)'))
    })

    output$MP.series <- downloadHandler(
        filename = function() paste0(MP.data()$fname,'_series.txt'),
        content = function(file){
            v <- MP.data()
            tab <- data.frame(Time=v$t,y=v$y,ey=v$dy,set=if(is.null(v$set)) rep(paste(v$ypar),length(v$t)) else v$set)
            colnames(tab)[2:3] <- c(v$ypar,paste0('e',v$ypar))
            write.table(tab,file,quote=FALSE,row.names=FALSE)
        }
    )

    output$download.MP.series <- renderUI({
        if(is.null(MP.data()) || !isTRUE(MP.data()$signal.only)) return()
        downloadButton('MP.series', 'Download the combined signal-only series')
    })

  observeEvent(input$plot2D,{
      output$per2 <- renderPlot({
          isolate({
              plotMP(MP.data(),per.par2())
          })
      })
  })

    observeEvent(input$plot2D,{
        output$plot.2Dper <- renderUI({
            plotOutput("per2", width = "600px", height = "600px")
        })
    })

    output$per2D.figure <- downloadHandler(
        filename = function() {
#            f1 <- gsub(" ",'_',Sys.time())
#            f2 <- gsub(":",'-',f1)
#            paste('periodogram2D_', f2, '.pdf', sep='')
            paste0(MP.data()$fname,'_scale',input$scale,'_Dt',signif(tspan()*as.numeric(input$Dt),3),'d.pdf')
        },
        content = function(file) {
            pdf(file,8,8)
            plotMP(MP.data(),per.par2())
            dev.off()
        })

    output$download.per2D.plot <- renderUI({
        if(is.null(MP.data())) return()
        downloadButton('per2D.figure', 'Download 2D periodogram')
    })

})
