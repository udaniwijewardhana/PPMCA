# PPMCA

## Persistence Probability Models and Changepoint Analysis

Users can use this app to identify the significant explanatory variables and the trends of the species during a span of years. The app would use truncated count models with log linear regression and zero inflation probability models with logistic regression. This app enables users to identify the significant changes in abundance level using Bayesian changepoint methods. Species persistence models in both spatial and temporal scales are used to identify the significant predictors for species count data while considering zero inflation models and detect significant changes in species abundance level using changepoint analysis. 

Udani Wijewardhana (udaniwijewardhana@gmail.com)1, Denny Meyer (dmeyer@swin.edu.au)1, and Madawa Jayawardana (madawa.jayawardana@petermac.org)1,2,3

1. Department of Statistics, Data Science and Epidemiology, Swinburne University of Technology, Hawthorn 3122 Victoria, Australia
2. Peter MacCallum Cancer Centre, Melbourne 3000 Victoria, Australia
3. Sir Peter MacCallum Department of Oncology, The University of Melbourne, Parkville 3010 Victoria, Australia

- Use a basic Persistence Probability Model to identify the persistence probabilities for each location during the most recent years. 
- Use single species distribution models considering both abundance (count model) and occupancy (zero inflation probability model) to identify the significant factors for this abundance /occurrence using INLA. This app can only take geostatistical data into account.
- Use four major changepoint methods to identify some significant changes of the species abundance level using raw data.

This application includes five main parts namely 'Annual data and Predictors', 'Persistence Probability Model', 'Abundance plots' with a distribution map, 'Single-species Joint Model' and 'Changepoint Analysis'. Figure 1 illustrates the interface of the PPMCA application. This app is only applicable for annual data. Main dataframe should have the columns Locality, Latitude, Longitude, Year and Count. These names are case sensitive. Annual predictors should enter as a separate dataframe. Possible maximum number of predictors is 5. These predictor coefficients display as p.y1, p.y2, ... for count model and p.z1, p.z2, ... for zero inflation probability model. The order is same as in the .csv file of predictors.

[![Video Tutorial](https://img.youtube.com/vi/cX_uNPx7UOQ/0.jpg)](https://www.youtube.com/watch?v=cX_uNPx7UOQ&feature=youtu.be)

Video 1. App tour.

![github-small](https://user-images.githubusercontent.com/55830031/80335620-27102f00-8898-11ea-89f7-4995a0407db9.png)

Figure 1. Persistence Probability Models and Changepoint Analysis interface.
                                            
The user can use the "HPdata.csv" (with "Predictors.csv") file to test the application. In this main file there are 5 columns with following column headings. The "Predictors.csv" file includes annual predictors (temp, rainfall and population).    
 
1. "Locality" - Detected Location 
2. "Latitude" - Latitude coordinate of the Location 
3. "Longitude" - Longitude coordinate of the Location 
4. "Year" - Detected Year
5. "Count" - Species count

First the user needs to upload the data csv file and explanatory variables csv file (optional) into the application and then select whether normalize the data and predictors or not to analyse. The explanatory variable csv file should be a table with annual data. 
The 'Annual data and Predictors' includes two data tables which are annual count data table and predictors data. 'Persistence Probability Model' also includes two data tables which give the basic persistence probabilities for each year and for the most recent year. The third window includes plots to visualise the abundance and the location distribution map. The 'Single-species Joint Model' gives a combination of count model and zero-inflation model results with spatial and temporal scale using INLA. This window also visualizes posterior mean and standard deviation plots for spatial and spatio-temporal models. This shows INLA mesh, summary results of the model and posterior plots (when applicable). Here, users can change the mesh parameters as well as spatial and temporal effect model in order to find the most suitable model. The 'Changepoint Analysis' window provides significant changes in trend using four Bayesian changepoint methods changepoint, breakpoint, cumSeg and bcp. This shows the summary of changepoint method results and the profile plot. Users could also refer to papers published by authors in order to find out further about these methods. Figure 2 shows a set of screenshots of the panels.

![github-small](https://user-images.githubusercontent.com/55830031/81125121-a2e63780-8f7a-11ea-9ac5-13e7afde8fc2.jpg)

Figure 2. Set of screenshots of app panels

1. Persistence Probability Model
2. Single-species Joint Model
3. Changepoint Analysis

## Persistence Probability Model

To obtain persistence probabilities use could consider annual occurrence data for the time period. The Bayesian formulation of Solow’s (1993) equation was used to calculate the posterior probability (p) of the species being extant (Ree and McCarthy, 2005).

                                               p = 1/(1+(1+[(T/t)^(N-1)-1])/((N-1)))

Assumed that the prior probability of the species being extant in the last year is 0.5. This prior probability is the probability of the species being extant in the last year of recording prior to considering the count data. In this formula, N is the number of years in which the species was recorded between year 0 and year T and the t is the year when the species was last recorded. The probability of persistence (p) is a score between 0 and 1.0 which is assessed in each area. The probability of persistence (p = 1) means the species  certain to be persistent at the end of the recording period, (p ≥ 0.5) means it is more likely to be extant in that area rather than extinct, and (p < 0.5) means it is more likely to be extinct than extant.

## Single-species Joint Model

To estimate the persistence, user can consider INLA models with zero-inflation probability model which means joint models. Joint temporal, spatial and spatio-temporal models are especially developed for data which has access zeros. Access zeros are divided into two types such as structural zeros and random zeros. Structural zeros refer to zero responses by those subjects whose count response will always be zero and random zeros that occur to subjects whose count response can be greater than zero but appear to be zero due to sampling variability (He et al., 2014). The responses of species data came from two different distributions such as occurrence and abundance which have two models for each response that are affected by spatial and temporal common factors. Therefore, it is better to use joint model. Negative Binomial (NB) model is more in line for zero inflation species data which resolve the oversdispersion issue as well. Therefore in this window users can get the summary outputs for Joint NB model (Cameron and Trivedi, 2013), Joint Hurdle NB model (Cragg, 1971; Mullahy, 1986) or Joint Zero Inflated NB (ZINB) model (Cameron and Trivedi, 1998) which are the most common zero inflation joint models to identify the significance of the predictors and to identify the persistency. 

A Bayesian hierarchical modelling approach is used here to allow us to conveniently account for structures in parameter uncertainty and potential dependence such as spatial and temporal structures. With a Bayesian approach, a joint posterior distribution is obtained for the process and parameters of the given data. This is conducted using Integrated Nested Laplace Approximation (INLA). INLA is popular as an approximation tool for fitting Bayesian models. INLA is an alternative robust method for the traditional Markov Chain Monte Carlo (MCMC) based Bayesian analyses (Paul et al. 2010). The key advantages of INLA are the ease with which complex models can be created and modified, without the need to write complex code, and the speed at which inference can be done even for spatial problems with hundreds of thousands of observations (Sarul, 2015). In spatio-temporal settings, it is often assumed that the covariance is separable in space and time, and thus, the temporal structure may be modelled using an auto-regressive process (Cressie, 2011). Fitted INLA with log-linear regression for the positive counts (truncated count model) and a logistic regression for the zero counts using a Bernoulli distribution in spatial and temporal scales. When we are dealing with process defined over a continuous domain, one can express a large class of random fields as solution of stochastic partial differential equations (SPDEs). In R-INLA this solution is approximated using high dimensional basis representation with simple local basis function. These basis functions are defined over a triangulation of the domain; this triangulation is the mesh (https://haakonbakka.bitbucket.io/btopic126.html). This app can only create non-convex meshes. The general hierarchical model can be described in three phases such as data model, process model and parameter model. Adapted these INLA models using R-INLA (Martins, Simpson, Lindgren and Rue, 2013).  

## Changepoint Analysis

User can use four Bayesian changepoint methods to identify the significant changes of the abundance level using raw data. Changepoint methods can be used only for a single location annual data which estimates using raw data while replacing missing values by zero. Therefore, the input data has taken as a single location to estimate changepoints. Predictors can only use with bcp method.

### changepoint package

The changepoint package implements various mainstream and specialised changepoint methods for finding single and multiple changepoints within data which includes many popular non-parametric and frequentist methods (Killick, Haynes and Eckley, 2016). 

### breakpoint package

The breakpoint package implements variants of the Cross-Entropy (CE) method to estimate both the number and the corresponding locations
of break-points in biological sequences of continuous and discrete measurements. The proposed method primarily built to detect multiple
break-points in genomic sequences. However, it can be easily extended and applied to other problems (Priyadarshana and Sofronov, 2016). 

### cumSeg package

The cumSeg package (Muggeo, 2010) estimates the of number and location of change points in mean-shift (piecewise constant) models which is useful to model genomic sequences of continuous measurements. The algorithm first estimates the highest number of change points using the efficient ‘segmented’ algorithm of Muggeo (2003) and then select some of them using a generalized BIC criterion by applying the lar’s algorithm of Efron et al. (2004) (Muggeo, 2010). 

### bcp package

The bcp package provides an implementation of the Barry and Hartigan (1993) product partition model for the normal errors change point problem using Markov Chain Monte Carlo. It also extends the methodology to regression models on a connected graph (Wang and Emerson, 2015) and allows estimation of change point models with multivariate responses (Erdman and Emerson, 2007).

## Acknowledgements 

The authors gratefully acknowledge the Ministry of Higher Education in Australia for providing a Research Training Program Stipend (RTPS) scholarship grant as funding for this research.

## Installation Instructions

### To install and explore the application in your R, you could type the following on your r console. 
                                         
> shiny::runGitHub( "PPMCA", "uwijewardhana") 
                                   
- User can access the standard R-repository to download and install package R-INLA by http://www.r-inla.org/download.
- User should install the other attached packages listed in sessionInfo() to explore the application in R console.                          
### You can also install the "PPMCA" package and run the app on your r console following the steps below.

- Download the https://github.com/uwijewardhana/PPMCA/ repository zip folder and extract.
- Open the extracted folder and open PPMCA.Rproj.
- Then run the below code in PPMCA.Rproj:
                        
> devtools::install_github("uwijewardhana/PPMCA")
> library(PPMCA)
> shiny::runApp(appDir = "./app.R")

## Session Information

> sessionInfo()

![github-small](https://user-images.githubusercontent.com/55830031/80561889-4b4f4580-8a29-11ea-81a1-023d9334023a.png)

## Manual 

[PPMCA_manual](https://github.com/cran/SpatialEpiApp/files/4613335/PPMCA_manual.pdf) 
