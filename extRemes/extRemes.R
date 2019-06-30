##### TUTORIAL: A short introduction to the extRemes package #####

#### Part 1: Maxima and GEV fitting ####

#Samples of size 100 of maxima of standard normal distributed samples
Zmax <- matrix( rnorm( 100 * 1000 ), 1000, 100 )
dim( Zmax )
Zmax <- apply( Zmax, 2, max )
class(Zmax)
hist(Zmax)

#load extRemes package
library( extRemes )
#fit a GEV distribution the the maxima data in Zmax
fit <- fevd( Zmax )
fit

#output for the fitted model 
plot( fit )
ci( fit, type = "parameter")
distill( fit )
strip( fit )

#real data: Fort Collins annual maximum temperatures in inches
data("ftcanmax")
plot(ftcanmax, type="l", lwd=2, ylab = "Annual Maximum Precipitation (inches/100)\nFort Collins, Colorado", col = "darkblue")
fit <- fevd( Prec, data = ftcanmax, units = "inches/100")
fit
plot( fit )
  
#with year as covariate for location
fit2 <- fevd( Prec, data = ftcanmax, 
              location.fun = ~Year, units = "inches/100" )
lr.test( fit, fit2 )

#### Part 2: Threshold exceedances ####

#GP fitting
data("Fort")
fit <- fevd( Prec, data = Fort, threshold = 0.395, type = "GP", 	units = "inches")
plot(fit)
plot( fit, type = "trace")

#### Part 3: Frequency of extremes (counting extreme events, declustering) ####

data(FCwx)
tempGT95 <- c(aggregate(FCwx$MxT, 
                        by = list(FCwx$Year), 
                        function(x) sum(x > 95, na.rm = TRUE))$x)
yr <- unique(FCwx$Year)
plot(yr, tempGT95, 
     type = "h", 
     col = "darkblue", 
     xlab = "Year",
     ylab = 
       "Number of Days with Max. Daily Temp. > 95 deg. F")
#check if there is overdispersion in the counts of exceedances per year 
#(which would be a sign of temporal clustering of threshold exceedances)
fpois( tempGT95 )
#yes!
#check if there are temporal trends in the number of counts
fit <- glm(tempGT95~yr, family = poisson())
summary(fit)
#yes, slight but significant increase over years

#fit the univariate point process model to threshold exceedances
fit <- fevd( Prec, Fort, 
             threshold = 0.395, type = "PP", 	units = "inches")
plot( fit )

#### Part 4: Threshold selection ####

#threshold selection using GP fits:
data("Denversp")
threshrange.plot( x = Denversp$Prec,
                  r = c( 0.1, 0.95 ) )
#threshold selection using point process fits:
data("Tphap")
threshrange.plot( x = 
                    Tphap$MaxT, 
                  r = c( 106, 110 ), type = "PP",
                  nint = 30 )

#### Part 5: Declustering ####
extremalindex( Tphap$MaxT, threshold = 105 )
y <- decluster( Tphap$MaxT, 
                threshold = 105, 
                r = 2 )
y
plot( y )
extremalindex( y, threshold = 105 )

#fit PP model to declustered data
Tphap2 <- Tphap
Tphap2$MaxT.dc <- c( y )
fit <- fevd( MaxT.dc, threshold = 105, 
             data = Tphap2, type = "PP",
             time.units = "62/year", 
             units = "deg F")
plot(fit)

#### Part 6: Bivariate Tail dependence analysis ####
z <- runif( 100, -1, 0 )     
w <- -1*(1 + z)     
taildep( z, w, u = 0.8 )
taildep.test( z, w )

