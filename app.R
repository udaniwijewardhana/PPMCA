library(shiny)
library(DT)
library(plyr)
library(ggplot2)
library(dplyr)
library(lattice)
library(gridExtra)
library(leaflet)
library(INLA)
library(breakpoint)
library(cumSeg)
library(changepoint)
library(bcp)

### Shiny user interface ###

ui <- fluidPage(

titlePanel(strong("Single Species Persistence Models and Changepoint Analysis", titleWidth = 350)),
hr(),

div(style="display: inline-block;vertical-align:top; width: 300px;", fileInput("file1", "Choose data CSV File", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
div(style="display: inline-block;vertical-align:top; width: 300px;", fileInput("file2", "Choose explanatory variables CSV File", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
div(style="display: inline-block;vertical-align:top; width: 300px;", selectInput("datanorm", "Count data normalization:", choices=c("rnorm", "stand", "none"), selected = "none")),
div(style="display: inline-block;vertical-align:top; width: 300px;", selectInput("prednorm", "Predictors normalization:", choices=c("rnorm", "stand", "none"), selected = "none")),

tags$style(type="text/css",
           ".shiny-output-error { visibility: hidden; }",
           ".shiny-output-error:before { visibility: hidden; }"
),

tabsetPanel(

tabPanel("Auunual Data and Predictors",
         fluidRow(style = "margin-top: 25px;",
                  column(7, p(tags$b('Annual Count Data', style = "font-size: 150%; font-family:Helvetica; color:#4c4c4c; text-align:left;"))),
                  column(5, p(tags$b('Annual Predictor Data', style = "font-size: 150%; font-family:Helvetica; color:#4c4c4c; text-align:left;")))),
         fluidRow(
         column(7, DT::dataTableOutput("contents")),
         column(5, DT::dataTableOutput("predictors")))
),

tabPanel("Persistence Proability Model",
          fluidRow(style = "margin-top: 25px;",
          column(6, p(tags$b('Persistence Probabilities for each year', style = "font-size: 150%; font-family:Helvetica; color:#4c4c4c; text-align:left;"))),
          column(6, p(tags$b('Persistence Probabilities for the most recent year', style = "font-size: 150%; font-family:Helvetica; color:#4c4c4c; text-align:left;")))),
          fluidRow(
          column(6, DTOutput("persistence1")),
          column(6, DTOutput("persistence2")))
),

tabPanel("Abundance Plots",
         fluidRow(style = "margin-top: 25px;",
                  column(5, p(tags$b('Abundance Plot', style = "font-size: 150%; font-family:Helvetica; color:#4c4c4c; text-align:left;"))),
                  column(7, p(tags$b('Location Map', style = "font-size: 150%; font-family:Helvetica; color:#4c4c4c; text-align:left;")))),
         fluidRow(column(5, plotOutput("plot1", "Abundance Plot", width = "100%", height = "400px")),
                  column(7, leafletOutput("plot2", width = "100%", height = "400px")))
),

tabPanel("Single-species Joint Model",
         sidebarLayout(
         sidebarPanel(
             sliderInput(inputId = "offset", label = "offsets for automatic boundaries",
                         value = c(0.1, 0.3), min = 0.05, max = 2),
             sliderInput(inputId = "max.edge", label = "max.edge",
                         value = c(0.05, 0.5), min = 0.01, max = 2),
             sliderInput(inputId = "cutoff", label = "cutoff",
                         value = 0.01, min = 0, max = 0.2),
             sliderInput(inputId = "min.angle", label = "min.angle",
                         value = c(21, 30), min = 1, max = 35),
             numericInput("convex", "convex for boundary", 0.5, min = NA, max = NA, step = 0.0001),
             actionButton("mesh", "mesh"),
             helpText("Zero-inflation model distribution is Binomial.
                             Applicable count models are 'Negative Binomial',
                             'Zero Inflated Negative Binomial' and 'Negative Binomial Hurdle'."),
             selectInput("distribution", "Distribution of count model:",
                         choices=c("Negative Binomial", "Zeroinflated Negative Binomial",
                                   "Negative Binomial Hurdle"), selected = "Negative Binomial Hurdle"),
             helpText("Posterior plots are only applicable for spatial or spatio-temporal models."),
             selectInput("speffect", "spatial random effect model:",
                         choices=c("spde", "iid"), selected = "iid"),
             selectInput("tempeffect", "temporal random effect model:",
                         choices=c("ar1", "iid", "rw1", "rw2"), selected = "ar1"),
             selectInput("int", "interaction between predictor variables:", choices = c("0","1"), selected = "0"),
             selectInput("term", "variables which have interactions:", choices=c("1 and 2","1 and 3", "1 and 4", "1 and 5",
                                                                                 "2 and 3", "2 and 4", "2 and 5",
                                                                                 "3 and 4", "3 and 5", "4 and 5", "all"), selected = ""),
             hr(), hr(), hr(), hr(), hr(), hr(),
             actionButton("summary", "Summary")
             ),
          mainPanel(fluidRow(column(4, tags$h3("Mesh:"), uiOutput("tab")),
                             column(8, tags$h3("Posterior Plots:"))),
                    fluidRow(column(4, plotOutput("mesh", "mesh", width = "80%", height = "400px")),
                             column(4, plotOutput("posteriormPlot")),
                             column(4, plotOutput("posteriorsdPlot"))),
                    tags$h3("Summary results of species joint model:"),
                    fluidRow(column(12, div(style='height:700px; overflow: scroll',
                             verbatimTextOutput("summary")))))
)),

tabPanel("Changepoint Analysis",
         sidebarLayout(
         sidebarPanel(div(style='height:800px; overflow: scroll',
         selectInput("changepoint", "changepoint method", choices=c("changepoint", "breakpoint", "cumSeg", "bcp"), selected = "changepoint"),
         conditionalPanel(
           condition = "input.changepoint == 'changepoint'",
             tags$h3("Algorithms of changepoint package:"),
             tags$ol(
             tags$li("cpt.mean(data,penalty='MBIC',pen.value=0,method='AMOC',
                      Q=5,test.stat='Normal',class=TRUE,param.estimates=TRUE,
                      minseglen=1)"),
             tags$li("cpt.var(data,penalty='MBIC',pen.value=0,know.mean=FALSE,
                      mu=NA,method='AMOC',Q=5,test.stat='Normal',
                      class=TRUE,param.estimates=TRUE,minseglen=2)"),
             tags$li("cpt.meanvar(data,penalty='MBIC',pen.value=0,method='AMOC',
                      Q=5,test.stat='Normal',class=TRUE,param.estimates=TRUE,
                      shape=1,minseglen=2)")),
           selectInput("changes", "type of changes to identify", choices=c("mean", "variance", "mean and variance"), selected = "mean"),
           selectInput("penalty", "penalty", choices=c("BIC", "MBIC", "AIC"), selected = "AIC"),
           selectInput("method", "method", choices=c("AMOC", "PELT", "SegNeigh", "BinSeg"), selected = "PELT"),
           numericInput("Q", "maximum number of changepoints to search", 5, min = 0, max = NA, step = 1),
           helpText("The maximum number of changepoints to search for using the BinSeg method.
                     The maximum number of segments (number of changepoints + 1) to search for using the SegNeigh method."),
           selectInput("test.stat", "test statistic", choices=c("Normal", "Gamma", "Exponential", "Poisson"), selected = "Normal"),
           helpText("Gamma, Exponential and Poisson distribution only applicable for mean and variance function."),
           selectInput("class", "class", choices=c("TRUE", "FALSE"), selected = "TRUE"),
           selectInput("param.estimates", "param.estimates", choices=c("TRUE", "FALSE"), selected = "TRUE"),
           numericInput("shape", "shape parameter for Gamma distribution", 1, min = -100, max = NA, step = 0.001),
           numericInput("minseglen", "minseglen", 1, min = 1, max = NA, step = 1),
           helpText("Positive integer giving the minimum segment length (no. of observations between
                      changes), default is the minimum allowed by theory. For cpt.var and cpt.meanvar,
                      the minimum minseglen is 2 and for cpt.mean, the minimum minseglen is 1.")),

           conditionalPanel(
           condition = "input.changepoint == 'breakpoint'",
             tags$h3("Algorithms of breakpoint package:"),
             tags$ol(
             tags$li("CE.NB(data, Nmax = 10, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8,
                               distyp = 1, penalty = 'BIC', parallel = FALSE)"),
             tags$li("CE.NB.Init(data, init.locs, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8,
                               distyp = 1, penalty = 'BIC', var.init = 1e+05, parallel = FALSE)"),
             tags$li("CE.Normal.Init.Mean(data, init.locs, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8,
                               b = 0.8, distyp = 1, penalty = 'BIC', var.init = 1e+05, parallel = FALSE)"),
             tags$li("CE.Normal.Init.MeanVar(data, init.locs, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8,
                               b = 0.8, distyp = 1, penalty = 'BIC', var.init = 1e+05, parallel = FALSE)"),
             tags$li("CE.Normal.Mean(data, Nmax = 10, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8,
                               distyp = 1, penalty = 'BIC', parallel = FALSE)"),
             tags$li("CE.Normal.MeanVar(data, Nmax = 10, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8,
                               b = 0.8, distyp = 1, penalty = 'BIC', parallel = FALSE)"),
             tags$li("CE.ZINB(data, Nmax = 10, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8,
                               distyp = 1, penalty = 'BIC', parallel = FALSE)"),
             tags$li("CE.ZINB.Init(data, init.locs, eps = 0.01, rho = 0.05, M = 200, h = 5, a = 0.8, b = 0.8,
                               distyp = 1, penalty = 'BIC', var.init = 1e+05, parallel = FALSE)")),
           selectInput("algorithm", "breakpoint function", choices=c("CE.NB", "CE.NB.Init",
                                                                     "CE.Normal.Init.Mean", "CE.Normal.Init.MeanVar", "CE.Normal.Mean", "CE.Normal.MeanVar",
                                                                     "CE.ZINB", "CE.ZINB.Init"), selected = "CE.NB"),
           numericInput("Nmax", "maximum number of breakpoints to search", 10, min = 0, max = NA, step = 1),
           numericInput("eps", "the cut-off value for the stopping criterion in the CE method", 0.01, min = 0, max = 100, step = 0.00001),
           numericInput("rho", "the fraction which is used to obtain the best performing set of sample solutions", 0.05, min = 0, max = NA, step = 0.00001),
           numericInput("M", "sample size to be used in simulating the locations of break-points", 200, min = 0, max = NA, step = 1),
           numericInput("h", "minimum aberration width", 5, min = 0, max = NA, step = 1),
           numericInput("a", "a smoothing parameter value", 0.8, min = 0, max = NA, step = 0.00001),
           helpText("It is used in the four parameter beta distribution to smooth both shape parameters.
                      When simulating from the truncated normal distribution,
                      this value is used to smooth the estimates of the mean values."),
           numericInput("b", "a smoothing parameter value", 0.8, min = 0, max = NA, step = 0.00001),
           helpText("It is used in the truncated normal distribution to smooth the estimates of the standard deviation."),
           selectInput("distyp", "distribution to simulate break-point locations", choices=c("1", "2"), selected = "1"),
           helpText("Options: 1 = four parameter beta distribution, 2 = truncated normal distribution"),
           selectInput("penalty", "penalty", choices=c("BIC", "AIC"), selected = "BIC"),
           numericInput("var.init", "Initial variance value to facilitate the search process", 100000, min = 0, max = NA, step = 0.00001),
           selectInput("parallel", "parallel", choices=c("TRUE", "FALSE"), selected = "FALSE"),
           textInput("init.locs", "Initial break-point locations - enter a vector (comma delimited) - e.g. '0,1,2'", NULL)),

           conditionalPanel(
           condition = "input.changepoint == 'cumSeg'",
             tags$h3("Algorithms of cumSeg package:"),
             tags$li("fit.control(toll = 0.001, it.max = 5, display = FALSE, last = TRUE,
                             maxit.glm = 25, h = 1, stop.if.error = FALSE)"),
             tags$li("sel.control(display = FALSE, type = c('bic', 'mdl', 'rss'), S = 1,
                             Cn = log(log(n)), alg = c('stepwise', 'lasso'), edf.psi = TRUE)"),
             tags$li("jumpoints(y, x, k = min(30, round(length(y)/10)), output = '2',
                             psi = NULL, round = TRUE, control = fit.control(),
                             selection = sel.control())"),
           hr(),
           numericInput("k", "k = the starting number of changepoints", 2, min = 0, max = NA, step = 1),
           helpText("It should be quite larger than the supposed number of (true) changepoints.
                    This argument is ignored if starting values of the
                    changepoints are specified via psi."),
           selectInput("output", "output", choices=c(1,2,3), selected = 1),
           textInput("psi", "psi = starting values for the changepoints - enter a vector (comma delimited) - e.g. '0,1,2'", NULL),
           helpText("When psi=NULL (default), k quantiles are assumed."),
           selectInput("round", "round", choices=c("TRUE", "FALSE"), selected = "TRUE"),
           numericInput("toll", "toll = positive convergence tolerance", 0.001, min = 0, max = NA, step = 0.0001),
           numericInput("it.max", "it.max = integer giving the maximal number of iterations", 5, min = 0, max = NA, step = 1),
           selectInput("display", "display", choices=c("TRUE", "FALSE"), selected = "FALSE"),
           selectInput("last", "last", choices=c("TRUE", "FALSE"), selected = "TRUE"),
           numericInput("maxit.glm", "maxit.glm", 25, min = 0, max = NA, step = 1),
           numericInput("h", "h", 1, min = 0, max = NA, step = 1),
           selectInput("stop.if.error", "stop.if.error", choices=c("TRUE", "FALSE"), selected = "FALSE"),
           helpText("logical indicating if the algorithm should stop when one or more estimated
                     changepoints do not assume admissible values."),
           selectInput("type", "type", choices=c("bic", "mdl", "rss"), selected = "bic"),
           numericInput("S", "S", 1, min = 0, max = NA, step = 1),
           helpText("If type = 'rss' the optimal model is selected when the residual sum of squares
                     decreases by the threshold S."),
           numericInput("Cn", "Cn", 1, min = 0, max = NA, step = 0.0000001),
           helpText("If type= 'bic' a character string (as a function of 'n') to specify to generalized
                     BIC. If Cn=1 the standard BIC is used."),
             selectInput("alg", "alg", choices=c("stepwise", "lasso"), selected = "stepwise"),
             selectInput("edf.psi", "edf.psi", choices=c("TRUE", "FALSE"), selected = "TRUE")),

           conditionalPanel(
           condition = "input.changepoint == 'bcp'",
             tags$h3("Algorithm of bcp package:"),
             tags$li("bcp(y, x = NULL, id = NULL, adj = NULL, w0 = NULL, p0 = 0.2,
                       d = 10, burnin = 50, mcmc = 500, return.mcmc = FALSE,
                       boundaryType = 'node', p1 = 1, freqAPP = 20, nreg = -1)"),
           hr(),
           numericInput("adjk", "k = the number of neighbors assumed for a typical vertex in adj list", 4, min = 4, max = 8, step = 4),
           helpText("An adjacency list. Indexing the observations from 1 to n, the i-th
                    element of the list is a vector of indices (offset by 1) for nodes that share an
                    edge with node i. Generates an adjacency list for a n node by m node grid,
                    assuming a maximum of k neighbors. k must be always 4 or 8. n and m are always >= 2.
                    Adjacency list algorithm = makeAdjGrid(n,m,k)"),
           selectInput("pred", "predictors", choices=c("Yes", "No"), selected = "No"),
           helpText("w0 is applicable when predictors are available."),
           textInput("w0", "w0 - enter a vector (comma delimited) - e.g. '0.2,0.5,0.7'", NULL),
           numericInput("p0", "p0", 0.2, min = 0, max = 1, step = 0.01),
           helpText("w0 and p0 are optional values which are between 0 and 1
                    for Barry and Haritgan's hyperparameters; these
                    default to the value 0.2, which has been found to work well."),
           numericInput("p1", "p1 = The proportion of Active Pixel Passes run that are the actual
                        Active Pixel Passes", 1, min = 0, max = 1, step = 0.01),
           numericInput("d", "d", 10, min = 0, max = NA, step = 1),
           helpText("a positive number only used for linear regression change point models.
                    Lower d means higher chance of fitting the full linear model
                    (instead of the intercept-only model)."),
           numericInput("burnin", "burnin", 50, min = 0, max = NA, step = 1),
           numericInput("mcmc", "mcmc", 500, min = 0, max = NA, step = 1),
           selectInput("return.mcmc", "return.mcmc", choices=c("TRUE", "FALSE"), selected = "FALSE"),
           numericInput("freqAPP", "freqAPP", 20, min = 0, max = NA, step = 1),
           selectInput("boundaryType", "boundaryType", choices=c("node", "edge"), selected = "node"),
           numericInput("nreg", "nreg", -1, min = 2, max = NA, step = 1),
           helpText("only applicable for regression; related to parameter d describing the
                    minimum number of observations needed in a block to allow for fitting a regression
                    model. Defaults to 2*number of predictors.")
         ))),

         mainPanel(
           tags$h3("Summary results of changepoint method:"),
           fluidRow(
           column(12, div(style='height:336px; overflow: scroll',
                  verbatimTextOutput("changepoint")))),
           hr(),
           tags$h3("Changepoint Plot:"),
           fluidRow(column(8, plotOutput("changepointPlot"))))
))))

### Shiny server ###

server <- function(input, output) {

# Input data csv file

filedata1 <- reactive({
    inFile1 <- input$file1
    if (is.null(inFile1)){
    return(NULL)}

    x <- as.data.frame(read.csv(inFile1$datapath))
    y <- aggregate(Count ~ Locality + Year, x, FUN = sum)
    y <- y[order(y$Locality, y$Year),]
    y$ID <- paste(y$Locality, y$Year, sep = "-")

    nLoc = length(unique(y$Locality))
    nyear = length(unique(y$Year))
    Locality <- as.character(rep(unique(y$Locality), each = nyear))
    Year <- rep(unique(x$Year), nLoc)
    Y <- as.data.frame(cbind(Locality, Year))
    Y$ID <- paste(Locality, Year, sep = "-")

    Final <- merge(Y,y,"ID",all.x=TRUE)
    names(Final)[names(Final) == "Year.x"] <- "Year"
    names(Final)[names(Final) == "Locality.x"] <- "Locality"

    y <- data.frame(Final$ID, Final$Locality, Final$Year, Final$Count)
    names(y)[names(y) == c("Final.ID", "Final.Locality", "Final.Year", "Final.Count")] <- c("ID", "Locality", "Year", "Count")

    y$Count[is.na(y$Count)] <- -1

    y$Occ <- ifelse(y$Count>=0,1,0)
    y$T <- rep(1:nyear, nLoc)
    y$N <- y$Occ
    y$t <- NA

    a <- split(y, y$Locality)
    for (j in 1:nLoc) { a[[j]]$N[1] <- ifelse((a[[j]]$Occ[1]==1), 1,0)}
    for (j in 1:nLoc){a[[j]]$t[1] = 0}
    for (j in 1:nLoc) {{for(i in 2:nyear) { a[[j]]$N[i] <- if(a[[j]]$Occ[i]==1) { (1 + a[[j]]$N[i-1]) } else { (a[[j]]$N[i-1]) }}}}
    for (j in 1:nLoc) {for(i in 2:nyear) { a[[j]]$t[i] <- (a[[j]]$N[i-1]) }}
    y$N <- unlist(lapply(seq_along(1:nLoc), function(x) as.numeric(a[[x]]$N)))
    y$t <- unlist(lapply(seq_along(1:nLoc), function(x) as.numeric(a[[x]]$t)))
    y$P <- mapply(function(N, T, t) (1/(1+(((T/t)^(N-1)-1)/(N-1)))), y$N,y$T,y$t)
    y$Conclusion <- with(y, ifelse(y$P==1, "Persisting", ifelse(y$P < 0.5 & y$P >= 0, "More likely to extinct", "More likely to extant")))
    y$Count[y$Count == -1] <- NA
    y$Count0 <- y$Count
    y$Count0[is.na(y$Count0)] <- 0

    Final <- merge(Y,y,"ID",all.x=TRUE)
    names(Final)[names(Final) == "Year.x"] <- "Year"
    names(Final)[names(Final) == "Locality.x"] <- "Locality"
    y <- Final[,!(names(Final) %in% c("Locality.y", "Year.y"))]
    y$yearNo <- rep(1:nyear, nLoc)
    y$P <- format(round(y$P, 5), nsmall = 5)
    y
})

# Input predictors csv file

filedata2 <- reactive({
    inFile1 <- input$file1
    if (is.null(inFile1)){
    return(NULL)}

    x <- as.data.frame(read.csv(inFile1$datapath))

    dat <- aggregate(Count ~ Locality + Latitude + Longitude + Year, x, FUN = sum)
    dat <- dat[order(dat$Locality, dat$Year),]
    dat$ID <- paste(dat$Locality, dat$Year, sep = "-")

    nLoc = length(unique(dat$Locality))
    nyear = length(unique(dat$Year))

    Locality <- as.character(rep(unique(dat$Locality[order(dat$Locality)]), each = nyear))
    Year <- rep(unique(dat$Year[order(dat$Year)]), nLoc)
    Latitude <- dat$Latitude
    Longitude <- dat$Longitude
    Latitude <- rep((unique(Latitude))[order(unique(Locality))], each = nyear)
    Longitude <- rep((unique(Longitude))[order(unique(Locality))], each = nyear)

    Y <- data.frame(Locality, Latitude, Longitude, Year)
    Y$ID <- paste(Y$Locality, Y$Year, sep = "-")

    Final <- merge(Y,dat,"ID",all.x=TRUE)
    names(Final)[names(Final) == "Year.x"] <- "Year"
    names(Final)[names(Final) == "Locality.x"] <- "Locality"
    names(Final)[names(Final) == "Latitude.x"] <- "Latitude"
    names(Final)[names(Final) == "Longitude.x"] <- "Longitude"
    Final <- Final[,!(names(Final) %in% c("Locality.y", "Year.y", "Latitude.y", "Longitude.y"))]
    return(Final)
})

# Create INLA mesh

mesh <- reactive({
  m <- filedata2()
  coords <- cbind(m$Latitude, m$Longitude)
  bnd = inla.nonconvex.hull(coords, convex = input$convex)

    out <- INLA::inla.mesh.2d(
    loc = coords,
    boundary = bnd,
    max.edge = input$max.edge,
    min.angle = rev(input$min.angle),
    cutoff = input$cutoff,
    offset = input$offset)
    return(out)

})

mm <- eventReactive(input$mesh, {
      plot(mesh())
      #points(coords, col = "red", pch = 16)
})

# Output of mesh

output$mesh <- renderPlot({
  return(mm())
})

filedata3 <- reactive({
  if (is.null(input$file2)){
  return(NULL)}
  df <- filedata2()
  inFile2 <- input$file2
  pred <- as.data.frame(read.csv(inFile2$datapath))

  n = nrow(df)
  count = df$Count
  nLoc = length(unique(df$Locality))
  nyear = length(unique(df$Year))

  pp = lapply(seq_along(1:ncol(pred)), function(x) rep((pred[,x]), nLoc))

  if(input$prednorm == "rnorm"){p <- lapply(seq_along(1:ncol(pred)), function(x) rnorm(pp[[x]]))
  } else if (input$prednorm == "stand"){p <- lapply(seq_along(1:ncol(pred)), function(x) scale(pp[[x]]))
  } else {p <- lapply(seq_along(1:ncol(pred)), function(x) pp[[x]])}

  pred.z = lapply(seq_along(1:ncol(pred)), function(x) c(p[[x]], rep(NA, n)))
  pred.y = lapply(seq_along(1:ncol(pred)), function(x) c(rep(NA, n), p[[x]]))
  pred.z <- data.frame(matrix(unlist(pred.z), nrow=n*2),stringsAsFactors=FALSE)
  pred.y <- data.frame(matrix(unlist(pred.y), nrow=n*2),stringsAsFactors=FALSE)

  return(cbind(pred.z, pred.y))
})

filedata4 <- reactive({
  if (is.null(filedata2())){
  return(NULL)}
  df <- filedata2()
  count = df$Count
  n.y = n.z = n = nrow(df)
  z <- rep(0,n)
  y = count

  if(input$distribution == "Negative Binomial Hurdle"){
    y[y == 0] <- NA
  } else {y = count}

  if(input$distribution == "Negative Binomial Hurdle"){
    z[which(count > 0)]<-1 } else {z[which(count > 0 | count == 0)]<-1}

  Y = matrix(rep(NA, 4 * n), ncol = 2)
  Y[1:n, 1] = z
  Y[(n + 1):(2*n), 2] = y
  return(Y)
})

# Create joint zero inflation model

fitsummary <- reactive({
    if (is.null(filedata3())){

    df1 <- filedata2()
    df2 <- filedata4()
    nLoc = length(unique(df1$Locality))
    nyear = length(unique(df1$Year))
    Year = rep(1:nyear, nLoc)
    coords <- cbind(df1$Latitude, df1$Longitude)
    n = nrow(df1)

    mu.z = rep(1:0, c(n, n))
    mu.y = rep(0:1, c(n, n))

    coords$hhid = mesh()$idx$loc
    S = rep(coords$hhid, 2)
    year = rep(Year, 2)

    spde = INLA::inla.spde2.matern(mesh(), constr = TRUE)

    sp <- if(input$speffect == "spde"){spde}else{"iid"}
    tempeffect <- as.character(if(input$tempeffect == "ar1"){"ar1"
    }else if(input$tempeffect == "rw1"){"rw1"
    }else if(input$tempeffect == "rw1"){"rw2"
    }else{"iid"})

    Y = df2
    data = list(Y = Y, mu.z = mu.z, mu.y = mu.y)

    formula = Y ~ 0 + mu.z + mu.y + f(S, model = sp) + f(year, model = tempeffect)

    distribution <- as.character(if(input$distribution == "Negative Binomial"){"nbinomial"
    } else if(input$distribution == "Zeroinflated Negative Binomial") {"zeroinflatednbinomial1"
    } else {"zeroinflatednbinomial0"})

    model <- INLA::inla(formula, data = data, family = c("binomial", distribution),
                        control.family = list(list(link = "logit"), list(link = "log")),
                        control.compute = list(dic = TRUE,cpo = TRUE, po = TRUE))
    return(model)

    } else {
    df1 <- filedata2()
    df2 <- filedata3()
    df3 <- filedata4()
    nLoc = length(unique(df1$Locality))
    nyear = length(unique(df1$Year))
    Year = rep(1:nyear, nLoc)
    coords <- cbind(df1$Latitude, df1$Longitude)
    n = nrow(df1)

    pred.z <- df2[,c(1:(ncol(df2)/2))]
    pred.y <- df2[,c(((ncol(df2)/2)+1):ncol(df2))]

    p.z = matrix(NA, ncol = 5, nrow = n*2)
    p.y = matrix(NA, ncol = 5, nrow = n*2)
    for( i in 1:ncol(pred.z)){p.z[,i] = pred.z[,i]}
    for( i in 1:ncol(pred.y)){p.y[,i] = pred.y[,i]}

    p.z1 = p.z[,1]
    p.z2 = p.z[,2]
    p.z3 = p.z[,3]
    p.z4 = p.z[,4]
    p.z5 = p.z[,5]

    p.y1 = p.y[,1]
    p.y2 = p.y[,2]
    p.y3 = p.y[,3]
    p.y4 = p.y[,4]
    p.y5 = p.y[,5]

    npred = (ncol(df2)/2)

    mu.z = rep(1:0, c(n, n))
    mu.y = rep(0:1, c(n, n))

    coords$hhid = mesh()$idx$loc
    S = rep(coords$hhid, 2)
    year = rep(Year, 2)

    spde = INLA::inla.spde2.matern(mesh(), constr = TRUE)

    sp <- if(input$speffect == "spde"){spde}else{"iid"}
    tempeffect <- as.character(if(input$tempeffect == "ar1"){"ar1"
    }else if(input$tempeffect == "rw1"){"rw1"
    }else if(input$tempeffect == "rw1"){"rw2"
    }else{"iid"})

    int <- as.integer(if(input$int == "1"){1
    }else {0})

    if(npred == 2 & int == 1){
      if(input$term == "1 and 2"){term <- as.integer(1)
      }else{return("error")}
    } else if(npred == 3 & int == 1){
      if(input$term == "1 and 2"){term <- as.integer(1)
      }else if(input$term == "1 and 3"){term <- as.integer(2)
      }else if(input$term == "2 and 3"){term <- as.integer(5)
      }else if(input$term == "all"){term <- as.character("all")
      }else{return("error")}
    } else if(npred == 4 & int == 1){
      if(input$term == "1 and 2"){term <- as.integer(1)
      }else if(input$term == "1 and 3"){term <- as.integer(2)
      }else if(input$term == "1 and 4"){term <- as.integer(3)
      }else if(input$term == "2 and 3"){term <- as.integer(5)
      }else if(input$term == "2 and 4"){term <- as.integer(6)
      }else if(input$term == "3 and 4"){term <- as.integer(8)
      }else if(input$term == "all"){term <- as.character("all")
      }else{return("error")}
    } else if(npred == 5 & int == 1){
      if(input$term == "1 and 2"){term <- as.integer(1)
      }else if(input$term == "1 and 3"){term <- as.integer(2)
      }else if(input$term == "1 and 4"){term <- as.integer(3)
      }else if(input$term == "1 and 5"){term <- as.integer(4)
      }else if(input$term == "2 and 3"){term <- as.integer(5)
      }else if(input$term == "2 and 4"){term <- as.integer(6)
      }else if(input$term == "2 and 5"){term <- as.integer(7)
      }else if(input$term == "3 and 4"){term <- as.integer(8)
      }else if(input$term == "3 and 5"){term <- as.integer(9)
      }else if(input$term == "4 and 5"){term <- as.integer(10)
      }else if(input$term == "all"){term <- as.character("all")
      }else{return("error")}
    }

    Y = df3
    data = list(Y = Y, mu.z = mu.z, mu.y = mu.y,
                p.z1 = p.z1, p.y1 = p.y1, p.z2 = p.z2, p.y2 = p.y2,
                p.z3 = p.z3, p.y3 = p.y3, p.z4 = p.z4, p.y4 = p.y4,
                p.z5 = p.z5, p.y5 = p.y5)

    if(npred == 1){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.y1 + f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 2 & int == 0){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.y1 + p.y2 + f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 2 & int == 1){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 + p.y1 * p.y2 + f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 3 & int == 0){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 + p.y1 + p.y2 + p.y3 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 3 & int == 1 & term == 1){
      formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 + p.z3 + p.y1 * p.y2 + p.y3 + f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 3 & int == 1 & term == 2){
      formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z3 + p.z2 + p.y1 * p.y3 + p.y2 + f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 3 & int == 1 & term == 5){
      formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z3 + p.y1 + p.y2 * p.y3 + f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 3 & int == 1 & term == "all"){
      formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 * p.z3 + p.y1 * p.y2 * p.y3 + f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 4 & int == 0){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 + p.z4 + p.y1 + p.y2 + p.y3 + p.y4 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 4 & int == 1 & term == 1){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 + p.z3 + p.z4 + p.y1 * p.y2 + p.y3 + p.y4 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 4 & int == 0 & term == 2){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z3 + p.z2 + p.z4 + p.y1 * p.y3 + p.y2 + p.y4 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 4 & int == 0 & term == 3){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z4 + p.z2 + p.z3 + p.y1 * p.y4 + p.y2 + p.y3 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 4 & int == 0 & term == 5){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z3 + p.z4 + p.y1 + p.y2 * p.y3 + p.y4 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 4 & int == 0 & term == 6){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z4 + p.z3 + p.y1 + p.y2 * p.y4 + p.y3 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 4 & int == 0 & term == 8){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 * p.z4 + p.y1 + p.y2 + p.y3 * p.y4 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 4 & int == 0 & term == "all"){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 * p.z3 * p.z4 + p.y1 * p.y2 * p.y3 * p.y4 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 0){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 + p.z4 + p.z5 + p.y1 + p.y2 + p.y3 + p.y4 + p.y5 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 1 & term == 1){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 + p.z3 + p.z4 + p.z5 + p.y1 * p.y2 + p.y3 + p.y4 + p.y5 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 0 & term == 2){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z3 + p.z2 + p.z4 + p.z5 + p.y1 * p.y3 + p.y2 + p.y4 + p.y5 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 0 & term == 3){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z4 + p.z2 + p.z3 + p.z5 + p.y1 * p.y4 + p.y2 + p.y3 + p.y5 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 0 & term == 4){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z5 + p.z2 + p.z3 + p.z4 + p.y1 * p.y5 + p.y2 + p.y3 + p.y4 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 0 & term == 5){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z3 + p.z4 + p.z5 + p.y1 + p.y2 * p.y3 + p.y4 + p.y5 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 0 & term == 6){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z4 + p.z3 + p.z5 + p.y1 + p.y2 * p.y4 + p.y3 + p.y5 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 0 & term == 7){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 * p.z5 + p.z3 + p.z4 + p.y1 + p.y2 * p.y5 + p.y3 + p.y4 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 0 & term == 8){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 * p.z4 + p.z5 + p.y1 + p.y2 + p.y3 * p.y4 + p.y5 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 0 & term == 9){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 * p.z5 + p.z4 + p.y1 + p.y2 + p.y3 * p.y5 + p.y4 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 0 & term == 10){formula = Y ~ 0 + mu.z + mu.y + p.z1 + p.z2 + p.z3 + p.z4 * p.z5 + p.y1 + p.y2 + p.y3 + p.y4 * p.y5 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else if(npred == 5 & int == 1 & term == "all"){formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 * p.z3 * p.z4 * p.z5 + p.y1 * p.y2 * p.y3 * p.y4 * p.z5 +
      f(S, model = sp) + f(year, model = tempeffect)
    } else {formula = Y ~ 0 + mu.z + mu.y + p.z1 * p.z2 * p.z3 * p.z4 * p.z5 + p.y1 * p.y2 * p.y3 * p.y4 * p.z5 +
      f(S, model = sp) + f(year, model = tempeffect)
    }

    distribution <- as.character(if(input$distribution == "Negative Binomial"){"nbinomial"
    } else if(input$distribution == "Zeroinflated Negative Binomial") {"zeroinflatednbinomial1"
    } else{"zeroinflatednbinomial0"})

    model <- INLA::inla(formula, data = data, family = c("binomial", distribution),
                        control.family = list(list(link = "logit"), list(link = "log")),
                        control.compute = list(dic = TRUE,cpo = TRUE, po = TRUE))
    return(model)}
})

proj <- reactive({
  if (is.null(fitsummary())){
  return(NULL)}
  proj <- INLA::inla.mesh.projector(mesh(), projection = "longlat", dims = c(150, 100))
  return(proj)
})

# Summary output of joint zero inflation model

fitsum <- eventReactive(input$summary, {
  fitsummary()
})

filedata5 <- reactive({
  inFile1 <- input$file1
  if (is.null(inFile1)){
  return(NULL)}

  x <- as.data.frame(read.csv(inFile1$datapath))
  data <- aggregate(Count ~ Year, x, FUN = sum, stringsAsFactors = FALSE)
  data <- data[order(data$Year), ]
  data$Count0 <- data$Count
  data$Count0[is.na(data$Count0)] <- 0

  if(input$datanorm == "rnorm"){data$Count0 <- rnorm(data$Count0)
  } else if (input$datanorm == "stand"){data$Count0 <- scale(data$Count0)
  } else {data$Count0 <- data$Count0}

  return(data)
})

filedata6 <- reactive({
  inFile2 <- input$file2
  if (is.null(inFile2)){
  return(NULL)}

  y <- as.data.frame(read.csv(inFile2$datapath))

  for(i in 1:ncol(y)){
    if(input$prednorm == "rnorm"){y[,i] <- rnorm(y[,i])
    } else if (input$prednorm == "stand"){y[,i] <- scale(y[,i])
    } else {y[,i] <- y[,i]}}

  return(y)
})

fit <- reactive({

  fit <- fit.control(toll = input$toll, it.max = input$it.max, display = input$display, last = input$last,
                     maxit.glm = input$maxit.glm, h = input$h, stop.if.error = input$stop.if.error)
  return(fit)
})

sel <- reactive({

  sel <- sel.control(display = input$display, type = input$type, S = input$S,
                     Cn = input$Cn, alg = input$alg, edf.psi = input$edf.psi)
  return(sel)
})

psi <- reactive({

  data <- filedata5()
  x <- 1:length(data$Count0)

  if(input$psi == ""){
    psi = quantile(x, prob= seq(0,1,l=input$k+2)[-c(1,input$k+2)], names=FALSE)
    return(psi)
  } else {
    psi = as.numeric(unlist(strsplit(input$psi,",")))
    return(psi)
  }
})

w0 <- reactive({

  inFile2 <- input$file2
  if (is.null(inFile2)){
    return(NULL)}

  predictors <- filedata6()
  w = as.numeric(unlist(strsplit(input$w0,",")))
  if (any(w > 1) | any(w < 0)){stop("Each element in w0 must be between 0 and 1.")}

  if(input$pred == "No"){return(NULL)}
  if(input$pred == "Yes"){
    if(length(w) == 0){
      w <- NULL
    } else {
      if(length(w) == ncol(predictors)){
        w <- w
      } else if(length(w) > ncol(predictors)){
        stop("Length of w0 is greater than the number of predictors.")
      } else if(length(w) == 1){
        w <- rep(w, ncol(predictors))
        print("Assume you wanted each error-to-signal ratio to be iid from U(0, w0).")
      } else {
        stop("Length of w0 is less than the number of predictors.")
      }
    }}
  return(w)
})

adjx <- reactive({

  inFile1 <- input$file1
  inFile2 <- input$file2
  n = nrow(filedata5())

  if (is.null(inFile1)){return(NULL)}
  if (is.null(inFile2)){m <- n} else {m <- ncol(filedata6())}

  adj = makeAdjGrid(n, m, input$adjk)
  return(adj)
})

# Changepoint analysis

changepoint <- reactive({

data <- filedata5()
inFile2 <- input$file2
w0 <- w0()
adj <- adjx()

if(input$changepoint == "changepoint"){

    if(input$method == "Gamma"){

    if(input$changes == "mean"){
      cpt = cpt.mean(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                     test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                     shape=input$shape, minseglen=input$minseglen)
    } else if(input$changes == "variance"){
      cpt = cpt.var(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                    test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                    shape=input$shape, minseglen=input$minseglen)
    } else{
      cpt = cpt.meanvar(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                        test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                        shape=input$shape, minseglen=input$minseglen)
    }} else{

      if(input$changes == "mean"){
        cpt = cpt.mean(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                       test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                       minseglen=input$minseglen)
      } else if(input$changes == "variance"){
        cpt = cpt.var(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                      test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                      minseglen=input$minseglen)
      } else{
        cpt = cpt.meanvar(data$Count0, penalty=input$penalty, pen.value=input$pen.value, method=input$method, Q=input$Q,
                          test.stat=input$test.stat, class=input$class, param.estimates=input$param.estimates,
                          minseglen=input$minseglen)
    }}
return(cpt)

} else if(input$changepoint == "breakpoint"){

  if(input$algorithm == "CE.NB"){
    bp = CE.NB(as.data.frame(data$Count0), Nmax = input$Nmax, eps = input$eps, rho = input$rho, M = input$M, h = input$h,
               a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty, parallel = input$parallel)
  } else if(input$algorithm == "CE.NB.Init"){
    bp = CE.NB.Init(as.data.frame(data$Count0), as.numeric(unlist(strsplit(input$init.locs,","))),
                    eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                    a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty,
                    var.init = input$var.init, parallel = input$parallel)
  } else if(input$algorithm == "CE.NB.Init.Mean"){
    bp = CE.Normal.Init.Mean(as.data.frame(data$Count0), as.numeric(unlist(strsplit(input$init.locs,","))),
                             eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                             a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty,
                             var.init = input$var.init, parallel = input$parallel)
  } else if(input$algorithm == "CE.NB.Init.MeanVar"){
    bp = CE.Normal.Init.MeanVar(as.data.frame(data$Count0), as.numeric(unlist(strsplit(input$init.locs,","))),
                                eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                                a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty,
                                var.init = input$var.init, parallel = input$parallel)
  } else if(input$algorithm == "CE.Normal.Mean"){
    bp = CE.Normal.Mean(as.data.frame(data$Count0), Nmax = input$Nmax, eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                        a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty, parallel = input$parallel)
  } else if(input$algorithm == "CE.Normal.MeanVar"){
    bp = CE.Normal.MeanVar(as.data.frame(data$Count0), Nmax = input$Nmax, eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                           a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty, parallel = input$parallel)
  } else if(input$algorithm == "CE.ZINB"){
    bp = CE.ZINB(as.data.frame(data$Count0), Nmax = input$Nmax, eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                 a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty, parallel = input$parallel)
  } else {
    bp = CE.ZINB.Init(as.data.frame(data$Count0), as.numeric(unlist(strsplit(input$init.locs,","))),
                      eps = input$eps, rho = input$rho, M = input$M, h = input$h,
                      a = input$a, b = input$b, distyp = input$distyp, penalty = input$penalty,
                      var.init = input$var.init, parallel = input$parallel)
  }
return(bp)

} else if(input$changepoint == "cumSeg"){

  x <- 1:length(data$Count0)

  cumSeg = jumpoints(as.matrix(data$Count0), x, output = input$output,
                     k = input$k, round = input$round,
                     psi = psi(),
                     control = fit(), selection = sel())
return(cumSeg)
} else {

  if(input$pred == "Yes" && !is.null(inFile2)){predictors <- as.matrix(filedata6())
  } else if(input$pred == "No" && !is.null(inFile2)){predictors <- NULL
  } else {predictors <- NULL}

  bcp = bcp::bcp(y = as.matrix(data$Count0), x = predictors, id = NULL,
                 adj = NULL, w0 = NULL,
                 p0 = input$p0, d = input$d, burnin = input$burnin, mcmc = input$mcmc,
                 return.mcmc = input$return.mcmc, boundaryType = input$boundaryType, p1 = input$p1,
                 freqAPP = input$freqApp, nreg = input$nreg)
return(bcp)
}
})

# Output of data table

output$contents <- DT::renderDataTable({
    req(input$file1)
    df <- read.csv(input$file1$datapath)
    df1 <- DT::datatable(filedata2())
    if (is.null(df)) return(NULL)
    df1
})

# Output of predictor table

output$predictors <- DT::renderDataTable({
  req(input$file2)
  df <- read.csv(input$file2$datapath)
  df1 <- DT::datatable(filedata6())
  if (is.null(df)) return(NULL)
  df1
})

# Output of abundance plot

output$plot1 <- renderPlot({
    abundanceplot <- ggplot(data=filedata1(), aes(x=Year, y=Count0, group=Locality)) +
                     xlab("Year") + ylab("Abundance") +
                     geom_line(aes(color=Locality))+ geom_point(aes(color=Locality))
    abundanceplot
})

# Output of species distribution map

output$plot2 <- renderLeaflet({
    df <- filedata2()
    map <- leaflet(df) %>% addTiles() %>% addMarkers(~Longitude, ~Latitude, label = ~as.character(Locality), clusterOptions = markerClusterOptions())
    map
})

# Persistence table output for each year

output$persistence1 <- renderDT({
    df <- filedata1()
    df <- df[,c(2:5,9,10)]
    if (is.null(df)) return(NULL)
    return(df)
})

# Persistence table output for most recent years

output$persistence2 <- renderDT({
    df <- filedata1()
    nyear = length(unique(df$Year))
    df <- df[which(df$yearNo==nyear),]
    df <- df[,c(2:5,9,10)]
    if (is.null(df)) return(NULL)
    return(df)
})

output$summary <- renderPrint({
  return(summary(fitsum()))
})

# Create posterior mean plot

pmPlot <- reactive({
  if(input$speffect == "iid"){return(NULL)
  } else {
mp <- levelplot(row.values=proj()$x, column.values=proj()$y,
                inla.mesh.project(proj(), fitsummary()$summary.random$S$mean),
                xla='Latitude', yla='Longitude',
                main='posterior mean plot', contour=TRUE,
                xlim=range(proj()$x), ylim=range(proj()$y))
      return(mp)}
})

# Create posterior standard deviation plot

psdPlot <- reactive({
  if(input$speffect == "iid"){return(NULL)
  } else {
sdp <- levelplot(row.values=proj()$x, column.values=proj()$y,
                 inla.mesh.project(proj(), fitsummary()$summary.random$S$sd),
                 xla='Latitude', yla='Longitude',
                 main='posterior standard deviation plot', contour=TRUE,
                 xlim=range(proj()$x), ylim=range(proj()$y))
       return(sdp)}
})

# Output of posterior standard deviation plot

output$posteriormPlot <- renderPlot({
  pmPlot()
})

# Output of posterior standard deviation plot

output$posteriorsdPlot <- renderPlot({
  psdPlot()
})

# Output of changepoint method summary results

output$changepoint <- renderPrint({
  changepoint()
})

# Output of changepoint method plot

output$changepointPlot <- renderPlot({
  if(input$changepoint == "breakpoint"){
  data <- filedata1()
  profilePlot(changepoint(), as.data.frame(data$Count0))
  } else {
  plot(changepoint())
  }
})

url <- a("Definition", href="https://rdrr.io/github/andrewzm/INLA/man/inla.mesh.2d.html")
output$tab <- renderUI({
   tagList("URL link:", url)
})

}

shinyApp(ui, server)
