library(shiny)
shinyUI(fluidPage(
    includeCSS("my_style.css"),
#    titlePanel(HTML("<font color='blue'>Needle</font>")),
    titlePanel(HTML("Agatha")),
    h4(HTML("<em>Disentangling periodic signals from correlated noise in a periodogram framework</em>")),
    tabsetPanel(
        tabPanel("About Agatha",uiOutput('about')
                 ),
        tabPanel("Choose File",
                 sidebarLayout(
                     sidebarPanel(
                                        # h4('Upload Type'),
                                        #  checkboxGroupInput('uptype','Upload Type',choices=c('Select from list','Upload files'),  selected=NULL),
                         radioButtons('uptype','Upload Type',
                                      choices=c('Select from the list'='list','Upload files'='upload'),  selected='list'),
#                         uiOutput('select.type'),
                         uiOutput('files'),
                         uiOutput('uptext'),
                         actionButton('show','upload and show data'),
                         uiOutput('download.data')
                                        # added interface for uploading data from
                                        # http://shiny.rstudio.com/gallery/file-upload.html
                                        #       tags$br()
                     ),
                     mainPanel(
                         uiOutput('tab')
                     )
                 )
                 ),
        tabPanel("Scatter Plot",
                 pageWithSidebar(
                     headerPanel(''),
                     sidebarPanel(
                         helpText(HTML("<b>Stand-alone tab: it only needs the data sets loaded in 'Choose File'.</b>")),
                         uiOutput("scatter.target"),
                         uiOutput("xs"),
                         uiOutput("ys"),
                         actionButton('scatter', 'show scatter plot'),
                         uiOutput('download.scatter.button')
                     ),
                     mainPanel(
                                        #plotOutput('scatter')
                         uiOutput("scatter")
                     )
                 )
                 ),
        tabPanel("Model Comparison",
                 pageWithSidebar(
                     headerPanel(''),
                     sidebarPanel(
######noise model comparison
                         helpText(HTML("<b>Stand-alone tab: it only needs the data sets loaded in 'Choose File'. Its result is not used by the other tabs.</b>")),
                         uiOutput('comp.target'),
                         sliderInput("Nar.max",'Maximum number of AR components',min=0,max=10,value=1,step=1),
                         sliderInput("Nma.max",'Maximum number of MA components',min=0,max=10,value=1,step=1),
                         uiOutput('proxy.type'),
                         uiOutput('proxy.text'),
                         uiOutput('nI.basic'),
                         uiOutput('nI.max'),
                         helpText("The number of proxies is counted from the fourth column of the data."),
                         uiOutput('Nman'),
                         uiOutput('nI.man'),
                         uiOutput('nI.comp'),
                                        # uiOutput('nI.group'),
                                        # uiOutput('mlp'),
                         actionButton('compare', 'compare noise models')
                     ,width = 6),
                     mainPanel(
                         tags$style(type="text/css", "#file_progress { max-width: 200px; }"),
                         uiOutput('BFtab'),
                         uiOutput('optNoise'),
                         uiOutput('download.logBF.table')
                     ,width=6)
                 )
                 ),
####calculate 1D periodograms
        tabPanel("1D Periodogram",
                 pageWithSidebar(
                     headerPanel(''),
                     sidebarPanel(
                         helpText(HTML("<b>Needs only the data sets loaded in 'Choose File'. Its result (the fitted signals, the sequential model comparison and the MCMC when enabled) is the input of the 'N-body Stability' tab and of the signal-only option of the '2D Periodogram' tab, so run it with the noise model, signal type and number of signals you consider final.</b>")),
                         uiOutput('per.target'),
                         helpText("If more than one data set is selected, BFP and MLP fit a shared circular signal with one RV offset per data set."),
                         uiOutput('per.type'),
                         helpText("There could be errors in the calculation of the BFP if the data is small (e.g. less than 20 data points) or not well sampled. "),
                         uiOutput('noise.model'),
                         uiOutput('gp.par'),
                         uiOutput('nar'),
                         uiOutput('nma'),
                         uiOutput('nh'),
                         uiOutput('Inds'),
                         uiOutput('signal'),
                         uiOutput('mcf'),
                         sliderInput("prange","Period range in base-10 log scale",min = -2,max = 6,value = c(0.1,4),step=0.1),
                         sliderInput("ofac", "Oversampling factor", min = 0, max = 10, value=1,step=0.1),
#                         radioButtons("signal.type",'Signal type',c("Circular"="circular","Keplerian"='kepler','Stochastic'='stochastic')),
                                        # "Empty inputs" - they will be updated after the data is uploaded
                         helpText("If the BFP is selected, only 'RV' is available for the following observable selection."),
                         uiOutput("var"),
#                         helpText("If the BFP is selected, the periodograms are only calculated for RVs. The meaning of variables are as follows: 'all'--all observables, 'RVs'--RVs, 'Proxies'-- noise proxies,'Instrument:Variable'--individual observables"),
                                        #uiOutput('helpvar'),
                                        #uiOutput('tv'),
                         uiOutput('sequential'),
                         uiOutput('per.type.seq'),
                         uiOutput('Nsig.max'),
                         actionButton('plot1D', 'Calculate periodograms'),
#                         radioButtons('down.type','Download plots',choices=c('All plots'='all','Individual plot'='individual'),  selected='individual'),
#                         uiOutput('select.type'),
                         uiOutput('plot.single'),
                         uiOutput('download.per1D.plot'),
#                         uiOutput('download.phase1D.plot'),
                         helpText("The users are encouraged to make their own periodogram figures and phase plots by downloading and using the relevant data."),
                         uiOutput('download.all.data'),
                     ),
                     mainPanel(
                         uiOutput("plot.1Dcombined"),
                     )
                 )
                 ),
####calculate 2D periodograms
        tabPanel("2D Periodogram",
                 pageWithSidebar(
                     headerPanel(''),
                     sidebarPanel(
                         helpText(HTML("<b>Works on its own with the raw data and the noise settings below. Only the option 'Use the signal-only data of the 1D fit' requires a '1D Periodogram' result for the same data sets.</b>")),
                         uiOutput('per.target2'),
                         helpText("If more than one data set is selected, the data sets are combined after subtracting the best-fitted noise components of each set, and the combined series is analysed with the chosen periodogram type (BFP or MLP are recommended)."),
                         uiOutput('use.fit2'),
                         uiOutput('per.type2'),
                         htmlOutput('text2D'),
#                         tags$br(),
                         uiOutput('nar2'),
                         uiOutput('noise.model2'),
                         uiOutput('gp.par2'),
                         uiOutput('nma2'),
                         uiOutput('Inds2'),
                         sliderInput("ofac2", "Oversampling factor", min = 0.2, max = 20, value=1,step=0.2),
#                         selectInput("yvar", "Choose observables",
#                                     choices  = 'RVs',
#                                     selected = 'RVs',multiple=FALSE),
                         uiOutput('var2'),
#                         helpText("'all': the periodograms of all variables;
#                             'RVs': periodograms of RVs;
#                            'Indices': periodograms of Indices;
#                             'Instrument:Variable': individual variables"),
                                        #uiOutput('helpvar'),
                                       #uiOutput('tv'),
                         uiOutput("Dt"),
                         uiOutput('textDt'),
                         checkboxInput('adaptive2','Adaptive windows (a fixed number of points per window)',value=FALSE),
                         helpText("With adaptive windows each window holds the number of points a window of the chosen width contains on average, so irregularly sampled data never leave a window empty: sparse epochs get wider windows and dense epochs narrower ones. The period grid stays that of the chosen width."),
                         uiOutput('prange2'),
                         br(),
                         uiOutput("Nbin"),
                         helpText("The above parameters are called 'calculating parameters', which are used for calculate the moving periodogram.","The following parameters are called 'visualization parameters', and are set to optimize the visulization of signals."),
                         uiOutput("alpha"),
                         uiOutput('zoom'),
                         checkboxInput('scale','Normalize power',value=TRUE),
                         checkboxInput('show.signal','Show significant signals',value=TRUE),
                         helpText("If you change the calculating parameters, click both 'calculate' and 'plot' to show the 2D periodogram.", "If you only change the visualization parameters, only click 'plot' to show the periodogram."),
                         actionButton('data.update', 'calculate'),
                         actionButton('plot2D', 'plot'),
                         uiOutput('download.per2D.plot'),
                         helpText("The users are encouraged to make their own plot of moving periodogram by downloading and using the relevant data. The first row is the centers of time windows. The first column is the periods, and the rest data is the matrix of periodogram powers."),
                         uiOutput('download.MP.data'),
                         uiOutput('download.MP.series')
                     ,width=6),
                     mainPanel(
#                         plotOutput("per2", width = "750px", height = 400)
                         uiOutput("plot.2Dper"),
                         htmlOutput("color"),width=6
                     )
                 )
                 ),
####overall diagnosis of the signals in several RV data sets
        tabPanel("Signal Diagnosis",
                 pageWithSidebar(
                     headerPanel(''),
                     sidebarPanel(
                         helpText(HTML("<b>Self-contained: it needs only the loaded data sets and runs its own sequential searches, per-set and activity-index periodograms and moving periodograms. It does not use, and is not used by, the other tabs; the planets integrated in 'N-body Stability' come from the '1D Periodogram' run, so if the diagnosis rejects a signal, re-run the 1D tab with matching settings before the N-body test.</b>")),
                         uiOutput('per.target3'),
                         helpText("Bayes factor periodograms (BFP) of the combined RV data with the signals subtracted in sequence, for each chosen noise model; of every data set on its own; and of the activity indices and the window function - the diagnostic figure of Feng et al. (2020, ApJS 250, 29). A Keplerian signal should be significant, robust to the noise model, absent from the activity indices, and consistent over time (moving periodogram of the signal-only combined series)."),
                         checkboxGroupInput('noise3','Noise models',choices=c('White'='W','MA(1)'='MA','AR(1)'='AR','GP'='GP'),selected=c('W','MA','AR'),inline=TRUE),
                         helpText("The first chosen model is also used for the individual data sets and defines the reported signals. GP is the shared quasi-periodic (SHO) Gaussian process; it is much slower than the ARMA models, and with a free oscillation period it can absorb a real planetary signal, so fix the period to the stellar rotation period when it is known."),
                         conditionalPanel("input.noise3 && input.noise3.indexOf('GP') > -1",
                             fluidRow(column(6,numericInput('gp.Prot3','GP oscillation period [day] (blank: free)',value=NA,min=0)),
                                      column(6,numericInput('gp.tau3','GP damping time scale [day] (blank: free)',value=NA,min=0)))),
                         sliderInput('Nsig3','Maximum number of signals',min=1,max=6,value=3,step=1),
                         numericInput('lnBF.min3','ln(BF) a signal must exceed',value=5,min=0,max=100,step=0.5),
                         radioButtons('SigType3','Signal type',c('Circular'='circular','Keplerian'='kepler'),selected='circular',inline=TRUE),
                         sliderInput("prange3","Period range in base-10 log scale",min=-1,max=5,value=c(0.3,4),step=0.1),
                         sliderInput("ofac3","Oversampling factor",min=0.2,max=10,value=1,step=0.2),
                         numericInput('Prot3','Rotation period [day] (optional, drawn as a green dotted line)',value=NA,min=0),
                         checkboxInput('individual3','Periodograms of the individual data sets',value=TRUE),
                         checkboxInput('proxies3','Periodograms of the activity indices',value=TRUE),
                         checkboxInput('moving3','Moving periodogram of the signal-only combined series',value=TRUE),
                         sliderInput('Nwin3','Windows of the moving periodogram (width: half the time span)',min=3,max=12,value=5,step=1),
                         checkboxInput('adaptive3','Adaptive windows (a fixed number of points per window, never empty)',value=TRUE),
                         sliderInput('ncol3','Panels per row',min=2,max=4,value=3,step=1),
                         actionButton('diagnose','Diagnose signals'),
                         uiOutput('download.diag')
                     ,width=4),
                     mainPanel(
                         uiOutput('diag.text'),
                         tableOutput('diag.summary'),
                         uiOutput('plot.diag'),
                         uiOutput('plot.diag.mp')
                     ,width=8)
                 )
                 ),
####N-body stability of the orbital solution
        tabPanel("N-body Stability",
                 pageWithSidebar(
                     headerPanel(''),
                     sidebarPanel(
                         helpText(HTML("<b>Requires a '1D Periodogram' result: the planets are the fitted signals of that run (the MAP of the MCMC when MCMC was enabled, otherwise the maximum-likelihood fit), and the Monte Carlo option draws from its MCMC samples. Run the 1D tab first, with a Keplerian signal type and MCMC on for posterior sampling.</b>")),
                         uiOutput('nbody.source'),
                         fluidRow(column(6,numericInput('nbody.Mstar','Stellar mass [Msun]',value=1,min=0.01,step=0.01)),
                                  column(6,numericInput('nbody.inc','Orbital inclination [deg]',value=90,min=1,max=90,step=1))),
                         helpText("Planet masses are m sin i / sin i for coplanar orbits; the semi-major axes follow from Kepler's third law."),
                         uiOutput('nbody.engine'),
                         textInput('nbody.python','Python interpreter',value='python3'),
                         numericInput('nbody.tmax','Integration time [yr]',value=1e7,min=1),
                         fluidRow(column(6,numericInput('nbody.steps','Steps per innermost orbit',value=20,min=5,step=5)),
                                  column(6,numericInput('nbody.Nout','Outputs (log-spaced)',value=200,min=10,step=10))),
                         uiOutput('nbody.estimate'),
                         radioButtons('nbody.mode','Orbits to integrate',c('Nominal solution only'='nominal','Nominal plus a Monte Carlo sample'='mc'),selected='nominal'),
                         conditionalPanel("input['nbody.mode'] == 'mc'",
                             numericInput('nbody.N','Number of sampled systems',value=20,min=1,max=1000,step=1),
                             helpText('The samples are rows of the MCMC posterior when MCMC ran, otherwise Gaussian draws from the 1-sigma quantiles.')),
                         fluidRow(column(4,numericInput('nbody.atol','Allowed |da/a|',value=0.2,min=0.01,step=0.05)),
                                  column(4,numericInput('nbody.escape','Escape at [x max a]',value=10,min=2,step=1)),
                                  column(4,numericInput('nbody.hill','Encounter [mutual Hill radii]',value=1,min=0,step=0.5))),
                         actionButton('nbody.go','Integrate'),
                         actionButton('nbody.stop','Stop'),
                         uiOutput('nbody.progress'),
                         helpText('The integration runs in a separate process: the progress bar shows the system being integrated and the fraction of its integration time done, and Stop aborts it at any moment.'),
                         uiOutput('download.nbody')
                     ,width=4),
                     mainPanel(
                         tableOutput('nbody.table'),
                         tableOutput('nbody.amd'),
                         uiOutput('nbody.summary'),
                         tableOutput('nbody.results'),
                         uiOutput('plot.nbody')
                     ,width=8)
                 )
                 )
    )
)
        )
