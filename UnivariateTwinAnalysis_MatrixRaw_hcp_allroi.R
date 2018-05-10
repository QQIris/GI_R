
# -----------------------------------------------------------------------
# Program: UnivariateTwinAnalysis_MatrixRaw.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: ACE
# DataType: Twin
# Field: Human Behavior Genetics
#
# Purpose: 
#      Univariate Twin Analysis model to estimate causes of variation
#      Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.06	added Model, Data & Field metadata
#      Hermine Maes -- 2014.11.04 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

# Load Data
alldata <- read.csv("E:/coding/R/data/func/hp200_s2_level2_MSMALL.feat/face-avg_t.csv")
results <- matrix(NA,1,3)
for(j in 1:3){
    
    ffa_l <- alldata[c(1,j)]
    mzid = read.csv("./data/MZ_subjectID.csv")
    dzid = read.csv("./data/DZ_subjectID.csv")
    t1data <- c()
    t2data <- c()
    id <- c()
    
    get_twin_data <- function(id_data,target_data){
        for(i in 1:nrow(id_data)){ #循环遍历
            sinFamily <- id_data[i,] 
            t1 <- sinFamily[[2]]
            t2 <- sinFamily[[3]]
            id <- append(id,as.character(sinFamily[[1]]))
            if((t1 %in% as.matrix(target_data)) & (t2 %in% as.matrix(target_data))){
                t1data <- append(t1data,(target_data[which(target_data==t1),2]))
                t2data <- append(t2data,(target_data[which(target_data==t2),2]))#如果双生子的id都有面孔激活数据，则提取
            }
            else{
                t1data <- append(t1data,NA)#如果双生子id中有一个没有面孔激活数据，则两者都设为NA
                t2data <- append(t2data,NA)
            }
        }
        twin_data <- data.frame(t1=t1data,t2=t2data)
    }
    
    mzData <- get_twin_data(mzid,ffa_l)
    dzData <- get_twin_data(dzid,ffa_l)
    
    # Select Variables for Analysis
    Vars      <- 't'
    nv        <- 1       # number of variables
    ntv       <- nv*2    # number of total variables
    selVars   <- paste(Vars,c(rep(1,nv),rep(2,nv)),sep="")   #c('bmi1','bmi2')
    
    
    # Generate Descriptive Statistics
    colMeans(mzData,na.rm=TRUE)
    colMeans(dzData,na.rm=TRUE)
    cov(mzData,use="complete")
    cov(dzData,use="complete")
    # Prepare Data
    # -----------------------------------------------------------------------------
    
    require(OpenMx)
    
    # Set Starting Values
    svMe      <- 2.3      # start value for means
    svPa      <- .5      # start value for path coefficients (sqrt(variance/#ofpaths))
    
    # ACE Model
    # Matrices declared to store a, d, and e Path Coefficients
    pathA     <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                           free=TRUE, values=svPa, label="a11", name="a" ) 
    pathC     <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                           free=TRUE, values=svPa, label="c11", name="c" )
    pathE     <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                           free=TRUE, values=svPa, label="e11", name="e" )
    
    # Matrices generated to hold A, C, and E computed Variance Components
    covA      <- mxAlgebra( expression=a %*% t(a), name="A" )
    covC      <- mxAlgebra( expression=c %*% t(c), name="C" ) 
    covE      <- mxAlgebra( expression=e %*% t(e), name="E" )
    
    # Algebra to compute total variances
    covP      <- mxAlgebra( expression=A+C+E, name="V" )
    
    # Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
    meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, 
                           free=TRUE, values=svMe, label="mean", name="expMean" )
    covMZ     <- mxAlgebra( expression=rbind( cbind(V, A+C), 
                                              cbind(A+C, V)), name="expCovMZ" )
    covDZ     <- mxAlgebra( expression=rbind( cbind(V, 0.5%x%A+C), 
                                              cbind(0.5%x%A+C , V)), name="expCovDZ" )
    
    # Data objects for Multiple Groups
    dataMZ    <- mxData( observed=mzData, type="raw" )
    dataDZ    <- mxData( observed=dzData, type="raw" )
    
    # Objective objects for Multiple Groups
    expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMean", 
                                      dimnames=selVars )
    expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMean", 
                                      dimnames=selVars )
    funML     <- mxFitFunctionML()
    
    # Combine Groups
    pars      <- list( pathA, pathC, pathE, covA, covC, covE, covP )
    modelMZ   <- mxModel( pars, meanG, covMZ, dataMZ, expMZ, funML, name="MZ" )
    modelDZ   <- mxModel( pars, meanG, covDZ, dataDZ, expDZ, funML, name="DZ" )
    fitML     <- mxFitFunctionMultigroup(c("MZ.fitfunction","DZ.fitfunction") )
    twinACEModel  <- mxModel( "ACE", pars, modelMZ, modelDZ, fitML )
    
    # Run Model
    twinACEFit   <- mxRun(twinACEModel, intervals=T)
    twinACESum   <- summary(twinACEFit)
    twinACESum
    # Fit ACE Model with RawData and Matrices Input
    # -----------------------------------------------------------------------------
    
    twinACEFit <- mxRun(twinACEModel)
    # Run ACE Model
    # -----------------------------------------------------------------------------
    
    # Generate ACE Model Output
    estMean   <- mxEval(mean, twinACEFit)             # expected mean
    estCovMZ  <- mxEval(expCovMZ, twinACEFit$MZ)      # expected covariance matrix for MZ's
    estCovDZ  <- mxEval(expCovDZ, twinACEFit$DZ)      # expected covariance matrix for DZ's
    estVA     <- mxEval(a*a, twinACEFit)              # additive genetic variance, a^2
    estVC     <- mxEval(c*c, twinACEFit)              # dominance variance, d^2
    estVE     <- mxEval(e*e, twinACEFit)              # unique environmental variance, e^2
    estVP     <- (estVA+estVC+estVE)                  # total variance
    estPropVA <- estVA/estVP                          # standardized additive genetic variance
    estPropVC <- estVC/estVP                          # standardized dominance variance
    estPropVE <- estVE/estVP                          # standardized unique environmental variance
    estACE    <- rbind(cbind(estVA,estVC,estVE),      # table of estimates
                       cbind(estPropVA,estPropVC,estPropVE))
    LL_ACE    <- mxEval(objective, twinACEFit)        # likelihood of ADE model
    # Get Model Output
    # -----------------------------------------------------------------------------
    
    
    # Change Model
    twinAEModel   <- mxModel( twinACEFit, name="AE" )
    twinAEModel   <- omxSetParameters( twinAEModel, labels="c11", free=FALSE, values=0 )
    twinAEFit     <- mxRun(twinAEModel)
    # Run AE Model
    # -----------------------------------------------------------------------------
    
    # Generate AE Model Output
    estVA     <- mxEval(a*a, twinAEFit)               # additive genetic variance, a^2
    estVE     <- mxEval(e*e, twinAEFit)               # unique environmental variance, e^2
    estVP     <- (estVA+estVE)                    # total variance
    estPropVA <- estVA/estVP                      # standardized additive genetic variance
    estPropVE <- estVE/estVP                      # standardized unique environmental variance
    estAE     <- rbind(cbind(estVA,estVE),        # table of estimates
                       cbind(estPropVA,estPropVE))
    LL_AE     <- mxEval(objective, twinAEFit)         # likelihood of AE model
    LRT_ACE_AE <- LL_AE - LL_ACE
    # Get Model Output

    estACE
    estAE
    LRT_ACE_AE
    # Print relevant output
    # -----------------------------------------------------------------------------
    
}
