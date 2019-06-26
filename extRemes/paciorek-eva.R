#######################################
#### Setup Fort precipitation data ####
#######################################

data(Fort, package = 'extRemes')
years <- min(Fort$year):max(Fort$year)
nYears <- length(years)
threshold <- 1.5 # inches; might want to transform everything to cm
ord <- order(Fort$year, Fort$month, Fort$day) 
Fort <- Fort[ord, ]
ind <- Fort$Prec > threshold
index <- seq_len(nrow(Fort))
FortExc <- Fort[ind, ]
head(FortExc)

## Plot annual maxes to guide interpretation

FortAnnualMax <- aggregate(Prec ~ year, data = Fort, max)
plot(FortAnnualMax$year, FortAnnualMax$Prec, xlab = 'year', ylab = 'maximum of daily precip (inches)')

##########################################
#### Stationary POT fit to daily data ####
##########################################

library(climextRemes)

out <- fit_pot(FortExc$Prec, threshold = threshold, nBlocks = nYears, 
               returnPeriod = 20, returnValue = 3.5,
               getParams = TRUE, bootSE = FALSE)
names(out)

out$returnValue ## 3.32 for a 20-year return period
exp(out$logReturnProb)  # probability of 0.04 of exceeding 3.5 inches in a given year
exp(out$logReturnPeriod)  # 25-year return period for a return value of 3.5 inches
## Confidence intervals from basic delta-method-based standard errors
out$returnValue + c(-2, 2) * out$se_returnValue  ## return value
exp(out$logReturnPeriod + c(-2, 2) * out$se_logReturnPeriod)  ## return period

## Confidence intervals from bootstrap-based standard errors

out <- fit_pot(FortExc$Prec, threshold = threshold, nBlocks = nYears, 
               blockIndex = FortExc$year, firstBlock = min(Fort$year),
               returnPeriod = 20, returnValue = 3.5,
               getParams = TRUE, bootSE = TRUE)
## Note: 'blockIndex' allows matching of observations to blocks for bootstrap resampling

out$returnValue + c(-2, 2) * out$se_returnValue_boot ## return value
exp(out$logReturnPeriod + c(-2, 2) * out$se_logReturnPeriod_boot)  ## return period

########################################################
#### Fit with GEV location parameter linear in year ####
########################################################

out <- fit_pot(FortExc$Prec, threshold = threshold,
               x = data.frame(years = years), locationFun = ~years,
               nBlocks = nYears, blockIndex = FortExc$year, firstBlock = min(Fort$year),
               returnPeriod = 20, returnValue = 3.5,
               getParams = TRUE, bootSE = FALSE,
               xNew = data.frame(years = range(Fort$year)))
## Note: 'blockIndex' allows matching of observations to covariate ('x') values

out$returnValue
exp(out$logReturnPeriod)

################################
#### Fit after declustering ####
################################

out <- fit_pot(FortExc$Prec, threshold = threshold,
               x = data.frame(years = years), locationFun = ~years,
               nBlocks = nYears, blockIndex = FortExc$year, firstBlock = min(Fort$year),
               returnPeriod = 20, returnValue = 3.5,
               getParams = TRUE, bootSE = FALSE,
               xNew = data.frame(years = range(Fort$year)),
               index = index[ind], declustering = 'noruns')
## Note: 'index' allows determination of sequential exceedances

out$returnValue
exp(out$logReturnPeriod)


####################################
#### Account for missing values ####
####################################

## There are no missing values in the Fort dataset, but for illustration, assume that
## the number missing in each year is as follows:
nMissing <- rpois(length(unique(Fort$year)), 2)
propMissing <- nMissing / 365 ## doesn't account for leap years

out <- fit_pot(FortExc$Prec, threshold = threshold, nBlocks = nYears, 
               returnPeriod = 20, returnValue = 3.5,
               getParams = TRUE, bootSE = FALSE,
               proportionMissing = propMissing)

out$returnValue
exp(out$logReturnPeriod)

####################################################################
#### Report contrasts for return values, periods, probabilities ####
####################################################################

out <- fit_pot(FortExc$Prec, threshold = threshold,
               x = data.frame(years = years), locationFun = ~years,
               nBlocks = nYears, blockIndex = FortExc$year, firstBlock = min(Fort$year), 
               returnPeriod = 20, returnValue = 3.5,
               getParams = TRUE,
               xNew = data.frame(years = max(Fort$year)),    
               xContrast = data.frame(years = min(Fort$year)))
## Values corresponding to covariates in 'xNew' contrasted against values for covariates in 'xContrast'
## Here we compare 1999 to 1900

## Change in return value 1900 to 1999:
out$returnValueDiff
out$returnValueDiff + c(-2, 2) * out$se_returnValueDiff
## Relative risk: ratio of return probabilities 1999 to 1900
exp(out$logReturnProbDiff)
exp(out$logReturnProbDiff + c(-2, 2) * out$se_logReturnProbDiff)
    
##################################################
#### Analysis of seasonal total precipitation ####
##################################################

summerData <- Fort[Fort$month %in% 6:8, ]  # June, July, August precipitation
FortSummer <- aggregate(Prec ~ year, data = summerData, sum)
summerThreshold <- quantile(FortSummer$Prec, 0.8)  # not-so-extreme given limited obs when aggregated
FortSummerExc <- FortSummer[FortSummer$Prec > summerThreshold, ]
FortSummerExc

out <- fit_pot(FortSummerExc$Prec, threshold = summerThreshold,
               x = data.frame(years = years), locationFun = ~years, 
               nBlocks = nYears, blockIndex = FortSummerExc$year, firstBlock = min(Fort$year),
               returnPeriod = 20, returnValue = 10,
               getParams = TRUE, bootSE = FALSE,
               xNew = data.frame(years = max(Fort$year)),
               xContrast = data.frame(years = min(Fort$year)))
## Note: each year (single observation) treated as a block, so return probability
## can be interpreted as probability of exceeding a value in a single year.

exp(out$logReturnProb)  ## return probability for 10 inches summer precipitation for 1999
## Relative risk: ratio of return probabilities 1999 to 1900
exp(out$logReturnProbDiff)  
exp(out$logReturnProbDiff + c(-2, 2) * out$se_logReturnProbDiff) 

#####################################
#### Analysis of replicated data ####
#####################################

## CAM5 atmosphere-only ensemble under all forcings (factual scenario)
## for Texas-Oklahoma area, March-August precipitation
## Interest lies in changes in extreme low precipitation, prompted by 2011 drought
## Here we analyze changes from 1960 to 2013 using 50 ensemble members

library(ncdf4)
nc_all <- nc_open('pr_MarAug_LBNL_CAM5-1-1degree_All-Hist_est1_v2-0_196001-201312.nc')
all <- ncvar_get(nc_all, 'pr', start=c(1,1), count = c(-1,-1))

## Preprocessing to extract 50 ensemble members available for full 1960-2013 period
avail <- c(1:10, 36:50, 61:70, 86:100)
startYear <- 1960
endYear <- 2013
fullYears <- startYear:endYear
nYears <- length(fullYears)
baseline <- 1961:2010
baseIndices <- baseline - startYear + 1

## 54 years (rows) with 50 ensemble members (columns)
all <- all[ , avail]

## Calculate (relative) anomalies
mn <- mean(all[baseIndices, ])
allAnom <- all / mn

## Set up replicated data in format required
nReplicates <- ncol(allAnom)         # 50 replicates (ensemble members)
threshold <- quantile(allAnom, .02)   # 2th percentile of distribution
allAnomVec <- c(allAnom)  ## string out data in a column-wise vector
blockIndex <- rep(fullYears, nReplicates)
replicateIndex <- rep(1:nReplicates, each = nYears)
sub <- which(allAnomVec < threshold)

out <- fit_pot(allAnom[sub], threshold = threshold,
               x = data.frame(years = fullYears), locationFun = ~years, 
               nBlocks = nYears, blockIndex = blockIndex[sub], firstBlock = startYear,
               replicateIndex = replicateIndex[sub], nReplicates = nReplicates,
               returnPeriod = 20, returnValue = 0.6, 
               xNew = data.frame(years = endYear),
               xContrast = data.frame(years = startYear),
               getParams = TRUE, bootSE = FALSE,
               upperTail = FALSE)
## Note: linear model for location parameter
## Note: 'replicateIndex' and 'nReplicates' allow for correct handling of model ensembles
## Note: 'upperTail = FALSE' for analysis of lower tail extremes

exp(out$logReturnPeriod)   ## Return period for 2013
exp(out$logReturnPeriod + c(-2, 2) * out$se_logReturnPeriod)
exp(out$logReturnProbDiff)  ## Ratio of return probabilities for end year compared to begin year
exp(out$logReturnProbDiff + c(-2, 2) * out$se_logReturnProbDiff)  ## confidence interval


#######################################
#### Model-based event attribution ####
#######################################

## Binomial-based risk ratios using 400-member ensemble for each of
## all and natural forcings for 2011 to do event attribution for
## Texas/Oklahoma summer 2011 drought (March - August).

nc_all <- nc_open('pr_MarAug_LBNL_CAM5-1-1degree_All-Hist_est1_v2-0_196001-201312.nc')
nc_nat <- nc_open('pr_MarAug_LBNL_CAM5-1-1degree_Nat-Hist_CMIP5-est1_v2-0_196001-201312.nc')
nc_obs <- nc_open('pr_MarAug_CRU-TS-3-22_observed_196001-201312.nc')

all <- ncvar_get(nc_all, 'pr', start=c(1,1), count = c(-1,-1))
nat <- ncvar_get(nc_nat, 'pr', start=c(1,1), count = c(-1,-1))
obs <- ncvar_get(nc_obs, 'pr', start=c(1,1), count = c(-1,-1))

avail <- c(1:10, 36:50, 61:70, 86:100)

startYr = 1960
endYr = 2013
fullYrs = startYr:endYr
baseline = 1961:2010
baseIndices <- baseline - startYr + 1

## Calculate anomalies based on baseline period, 50 member full ensemble
mnModel = mean(all[baseIndices, avail])
mnObs = mean(obs[baseIndices])

allAnom = all / mnModel
natAnom = nat / mnModel  # relative to all forcings baseline
obsAnom = obs / mnObs

eventYear <-  52 # 2011
level <- 0.9  # 90% confidence intervals
thresholdPerc <- 0.1

allAnom <- allAnom[eventYear, ]
natAnom <- natAnom[eventYear, ]

## Actual event of 0.4 has no exceedances in either scenario, so use slightly less extreme definition
event <- 0.5

yA <- sum(allAnom < event)
yN <- sum(natAnom < event)
n <- length(allAnom)

result <- calc_riskRatio_binom(y = c(yA, yN), n = rep(n, 2),
                                    ciType = c('koopman', 'lrt'),
                                    ciLevel = level, lrtControl = list(bounds = c(0.01, 500)))
result$riskRatio
result$ci_riskRatio_lrt  ## likelihood-ratio based interval
result$ci_riskRatio_koopman  ## Koopman-based interval

## Consider even less extreme event of 0.6
event <- 0.6

yA <- sum(allAnom < event)
yN <- sum(natAnom < event)
n <- length(allAnom)

result <- calc_riskRatio_binom(y = c(yA, yN), n = rep(n, 2),
                                    ciType = c('koopman', 'lrt'),
                                    ciLevel = level, lrtControl = list(bounds = c(0.01, 500)))
result$riskRatio
result$ci_riskRatio_lrt  ## likelihood-ratio based interval
result$ci_riskRatio_koopman  ## Koopman-based interval

