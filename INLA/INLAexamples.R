library(INLA)
library(fields) #for quilt.plot

setwd("")
load("DATA_TRAINING_SUBSET.RData")

#fix the number of CPU cores that INLA is allowed to use
#(by default, all available cores)
nthreads=2

# check the data: 
dim(anom.training.sub)
n=nrow(anom.training.sub) #number of observation days
m=ncol(anom.training.sub) #number of sites
head(time.sub)
tail(time.sub)
#check if no days are missing?
table(diff(time.sub)) #time difference is 2 days (48 hours) for 3 cases 
#in our modeling, we assume regular 1-day time steps and do not take account of these very slightly irregular time steps

length(month.sub)
length(year.sub)
table(year.sub)
length(day.sub)
table((day.sub))


#part of the values is missing:
mean(is.na(anom.training.sub))
#all sites have some missing values:
min(apply(is.na(anom.training.sub),2,mean))
#some time steps have no missing values:
min(apply(is.na(anom.training.sub),1,mean))

hist(anom.training.sub)
mean(anom.training.sub,na.rm=TRUE)

#any trends over years?
year.mean=apply(anom.training.sub,1,mean,na.rm=TRUE)
boxplot(year.mean~year.sub)
tmp=year.sub-min(year.sub)
#maybe there is some time trend?
summary(lm(year.mean~tmp))

#some summary maps:
#mean:
tmp=apply(anom.training.sub,2,mean,na.rm=TRUE)
quilt.plot(loc.sub,tmp,nx=50,ny=50)
#standard deviation:
tmp=apply(anom.training.sub>1,2,sd,na.rm=TRUE)
quilt.plot(loc.sub,tmp,nx=50,ny=50)
#threshold exceedance probabilities:
thresh=quantile(anom.training.sub,0.95,na.rm=TRUE)
tmp=apply(anom.training.sub>thresh,2,mean,na.rm=TRUE)
quilt.plot(loc.sub,tmp,nx=50,ny=50)


#### Some auxiliary function(s) ####

plot_random_effect=function(effect){
  ylim=range(effect$mean,effect$`0.025quant`,effect$`0.975quant`)
  plot(effect$ID,effect$mean,type="l",lwd=2,ylim=ylim)
  lines(effect$ID,effect$`0.025quant`,lwd=2,col="blue",lty=2)
  lines(effect$ID,effect$`0.975quant`,lwd=2,col="blue",lty=2)
}


#### PART 1: Modeling of a time series of threshold exceedances ####

#choose the site that we want to model
idsite=512 
obs=anom.training.sub[,idsite]
plot(obs,type="l")
mean(is.na(obs))

#fix threshold, and generate exceedance indicators and positive exceedances
prob=0.95
thresh=quantile(obs,prob,na.rm=TRUE)
excess=pmax(obs-thresh,0)
isexceed=as.numeric(obs>thresh)

#Our estimation procedure consists of two steps here. We construct and fit two separate models, the first for the probability of exceedances over the threshold, the second for the positive exceedances over the threshold.
#Notice that we could also fit the two models simultaneously with INLA, using two likelihood types (binomial and GP) in the same INLA model, with random effects that are in common or have correlation between the two submodels. However, for simplicity, we here fit the models for exceedance probability and for the size of the excess separately. 

#we generate a data frame with the response y and covariates to be included in the model:
#to check for a linear time trend, we include a time covariate normalized to the interval [0,1]
daysFromStart=as.numeric((time.sub-time.sub[1]))/86400
xtime=daysFromStart/max(daysFromStart)

# Model 1: latent AR1 time series model ####

# Step 1: Modeling the exceedance probability ####

#generate the data frame to be used with INLA:
dataDF=data.frame(y=isexceed,month=month.sub,year=year.sub,xtime=xtime,idx=1:n)

#To learn about the formulation of the AR1 model and the RW1 model in INLA, we can consult the documentation:
#(these commands open a pdf file)
#inla.doc("ar1")
#inla.doc("rw1")

#We will use Penalized Complexity Priors for the hyperparameters:
hyper.prec=list(prior="pc.prec",param=c(10,.05),initial=1)
hyper.cor=list(prior="pccor1",param=c(.75,.5),initial=log(0.75/(1-0.75)))

#write the model formula in the usual R style, but with INLA-specific f(...)-terms for random effects:
form=y~xtime+f(idx,model="ar1",hyper=list(theta1=hyper.prec,theta2=hyper.cor))
#form=y~xtime+f(month,model="rw1",cyclic=TRUE,hyper=list(theta=hyper.prec))+f(idx,model="ar1",hyper=list(theta1=hyper.prec,theta2=hyper.cor))

#We here keep the missing response values (NA) during the estimation. inla will automatically calculate the predictions (=fitted values) for the missing responses. By default, an identity link function is used for missing values. By specifying "link=1" in the control.predictor-argument, we will use the correct link function (i.e., the link function of the first "family" - here, we have only one likelihood family "binomial").
  
#running the following estimation with inla takes around 1 minute: 
fit=inla(form,
         data=dataDF,
         family="binomial",
         Ntrials=rep(1,n),
         control.predictor=list(compute=TRUE,link=1),
         control.fixed=list(prec=1,mean.intercept=log((1-prob)/prob),prec.intercept=1),
         #control.inla=list(strategy="simplified.laplace",int.strategy="ccd"),
         control.compute=list(cpo=TRUE),
         verbose=FALSE,
         num.threads=nthreads
)
summary(fit)
fit1.p=fit
#We obtain a significant time trends since the 95% credible interval of xtime does not cover 0.

#Plot data and fitted values:
plot(dataDF$y,pch=4)
lines(1:n,fit$summary.fitted.values$mean,col="blue",lwd=2)

#Plot the AR1 random effect:
plot_random_effect(fit$summary.random$idx)

#Step 2: Modeling the size of the excess ####

# Here we remove all data without observed positive exceedance for fitting the generalized Pareto distribution.

idx2keep=which(!is.na(excess) & excess > 0)
n.exc=length(idx2keep)
n.exc

dataDF=data.frame(y=excess[idx2keep],month=month.sub[idx2keep],year=year.sub[idx2keep],xtime=xtime[idx2keep],idx=idx2keep)

form=y~xtime+f(idx,model="ar1",hyper=list(theta1=hyper.prec,theta2=hyper.cor))

# The link function is designed model a certain (fixed) quantile level of the generalized Pareto distribution.
# We here use the median, and we will specify this value in the control.family-argument of inla(...).
qprob=0.5

# You may consult inla.doc("loggamma") to see how the approximate PC prior (exponential) can be implemented.
#Since the data look rather light-tailed here, we use strong shrinkage of the model towards shape=0 by setting the penalization rate to 10.
hyper.shape=list(prior="loggamma", param=c(1,10),initial=log(.1))
#running the following estimation with inla takes around 2 seconds:
fit=inla(form,
         data=dataDF,
         family="gp",
         control.family=list(list(control.link=list(quantile=qprob),hyper=list(theta=hyper.shape))),
         control.predictor=list(compute=TRUE,link=1),
         control.fixed=list(prec=1,prec.intercept=1),
         #control.inla=list(strategy="simplified.laplace",int.strategy="ccd"),
         control.compute=list(cpo=TRUE),
         verbose=FALSE,
         num.threads=nthreads
)
summary(fit)
fit1.gp=fit

#Plot data against fitted GP median:
plot(fit$summary.fitted.values$mean,dataDF$y,pch=19,cex=.5,asp=1,xlab="fitted median",ylab="observed exceedance")
#Naturally, the fitted median values in our Bayesian model are relatively smooth, such that the points are not aligned along the diagonal. 

#Plot the AR1 random effect:
plot_random_effect(fit$summary.random$idx)
#(here with line segments over periods with NA values)

#plot of posterior means of medians without non-exceeding time steps:
plot(idx2keep,fit$summary.random$idx$mean,pch=19,cex=.5,xlab="day",ylab="estimated latent AR 1 effect")


# Model 2: use the observation in t-1 for prediction of t ####

# Step 1: Modeling the exceedance probability ####

obspreced=c(NA,obs[1:(n-1)])
#we have to remove the first observation where we do not have the covariat obspreced:
dataDF=data.frame(y=isexceed[-1],obspreced=obspreced[-1],month=month.sub[-1],year=year.sub[-1],xtime=xtime[-1],idx=2:n)

form=y~xtime+obspreced

#The following estimation runs very quickly (logistic regression without random effect).
fit=inla(form,
         data=dataDF,
         family="binomial",
         Ntrials=rep(1,n-1),
         control.predictor=list(compute=TRUE,link=1),
         control.fixed=list(prec=1,mean.intercept=log((1-prob)/prob),prec.intercept=1),
         #control.inla=list(strategy="simplified.laplace",int.strategy="ccd"),
         control.compute=list(cpo=TRUE),
         verbose=FALSE,
         num.threads=nthreads
)
summary(fit)
fit2.p=fit

#Plot data and fitted values:
plot(dataDF$y,pch=4)
lines(1:(n-1),fit$summary.fitted.values$mean,col="blue",lwd=2)

#Step 2: Model the positive exceedance value ####

obspreced=c(NA,obs[1:(n-1)])
#we have to remove the first observation where we do not have the covariat obspreced:
dataDF=data.frame(y=excess[idx2keep],obspreced=obspreced[idx2keep],month=month.sub[idx2keep],year=year.sub[idx2keep],xtime=xtime[idx2keep],idx=idx2keep)

form=y~xtime+obspreced
#running the following estimation with inla takes around 1 second:
fit=inla(form,
         data=dataDF,
         family="gp",
         control.family=list(list(control.link=list(quantile=qprob),hyper=list(theta=hyper.shape))),
         control.predictor=list(compute=TRUE,link=1),
         control.fixed=list(prec=1,prec.intercept=1),
         #control.inla=list(strategy="simplified.laplace",int.strategy="ccd"),
         control.compute=list(cpo=TRUE),
         verbose=FALSE,
         num.threads=nthreads
)
summary(fit)
fit2.gp=fit


# Compare the fitted values of the models 1 and 2 : ####
plot(fit1.p$summary.fitted.values$mean[-1],fit2.p$summary.fitted.values$mean,asp=1,pch=19,cex=.2)
abline(0,1,col="blue",lwd=2)

plot(fit1.gp$summary.fitted.values$mean,fit2.gp$summary.fitted.values$mean,asp=1)
abline(0,1,col="blue",lwd=2)


#### PART 2: Spatial modeling ####

# Model 3 : A spatial model ####

# To estimate spatial trends in the model, we here implement the estimation of a model with a spatial effect. The spatial effect is not replicated in time, which means that our focus is on trends and not on dependence.
# We here fit a Matérn correlation function with fixed range parameter. Notice that the SPDE approach implemented in INLA provides a much more flexible modeling framework around Matérn-based spatial correlation models, but implementing the model is slightly more technical. 

#Data subsetting: choose the first n.sub time steps, and m.sub locations at random. 
n.sub=3*365
m.sub=100
id.n.sub=1:n.sub
set.seed(1)
id.m.sub=sample(1:m,m.sub)
loc.sub.sub=loc.sub[id.m.sub,]
#transform to metric distances (km)
tmp=SpatialPoints(loc.sub.sub,proj4string=CRS("+init=epsg:4326"))
loc.sub.sub=spTransform(tmp,CRS("+init=epsg:32639"))@coords/10^3
colnames(loc.sub.sub)=c("x","y")
range(loc.sub.sub[,1])
range(loc.sub.sub[,2])
plot(loc.sub.sub,asp=1)


datamat=anom.training.sub[id.n.sub,id.m.sub]
obs=as.numeric(datamat)
mean(is.na(obs))
id.day=rep(id.m.sub,m.sub)
id.loc=rep(id.m.sub,each=n.sub)

#fix threshold, and generate exceedance indicators and positive exceedances
prob=0.9
thresh=quantile(obs,prob,na.rm=TRUE)
excess=pmax(obs-thresh,0)
isexceed=as.numeric(obs>thresh)

daysFromStart=as.numeric((time.sub[id.n.sub]-time.sub[1]))/86400
xtime=daysFromStart/max(daysFromStart)

#Consider a Matérn correlation with fixed correlation range here:
maternrange=150
dist=dist(loc.sub.sub)
kappa=sqrt(8)/maternrange
nu=1 
cor = as.matrix(2^(1-nu)*(kappa*dist)^nu*besselK(dist*kappa,nu)/gamma(nu))
diag(cor)=1
Q0 = solve(cor)

# Step 1: Modeling the exceedance probability ####

dataDF=data.frame(y=isexceed,month=rep(month.sub[id.n.sub],m.sub),xtime=rep(xtime,m.sub),id.loc=rep(1:m.sub,each=n.sub))

hyper.prec=list(prior="pc.prec", param=c(1, 0.5),initial=log(1))

form=y~xtime+f(id.loc,model="generic0",Cmatrix=Q0, hyper=list(theta=hyper.prec))+f(month,model="rw1",cyclic=TRUE,constr=TRUE,hyper=list(theta=hyper.prec))

#The following estimation may take a bit more time due to the inclusion of the spatial latent effect: 
fit=inla(form,
         data=dataDF,
         family="binomial",
         Ntrials=rep(1,n.sub*m.sub),
         control.predictor=list(compute=TRUE,link=1),
         control.fixed=list(prec=1,mean.intercept=log((1-prob)/prob),prec.intercept=1),
         control.inla=list(strategy="gaussian",int.strategy="eb"),
         #control.inla=list(strategy="simplified.laplace",int.strategy="ccd"),
         control.compute=list(cpo=TRUE),
         verbose=FALSE,
         num.threads=nthreads
)
summary(fit)
fit3.p=fit

plot_random_effect(fit$summary.random$month)

range(fit$summary.random$id.loc$mean)
plot(loc.sub.sub,cex=.5*exp(.5*fit$summary.random$id.loc$mean/mean(fit$summary.random$id.loc$mean)))

#Step 2: Model the positive exceedance value ####

idx2keep=which(!is.na(excess) & excess>0)
n.exc=length(idx2keep)
n.exc

dataDF=data.frame(y=excess[idx2keep],month=rep(month.sub[id.n.sub],m.sub)[idx2keep],xtime=rep(xtime,m.sub)[idx2keep],id.loc=rep(1:m.sub,each=n.sub)[idx2keep])

hyper.prec=list(prior="pc.prec", param=c(1, 0.5),initial=log(1))

form=y~xtime+f(id.loc,model="generic0",Cmatrix=Q0, hyper=list(theta=hyper.prec))+f(month,model="rw1",cyclic=TRUE,constr=TRUE,hyper=list(theta=hyper.prec))

fit=inla(form,
         data=dataDF,
         family="gp",
         control.family=list(list(control.link=list(quantile=qprob),hyper=list(theta=hyper.shape))),
         control.predictor=list(compute=TRUE,link=1),
         control.fixed=list(prec=1,prec.intercept=1),
         control.inla=list(strategy="simplified.laplace",int.strategy="ccd"),
         control.compute=list(cpo=TRUE),
         verbose=FALSE,
         num.threads=nthreads
)
summary(fit)
fit3.gp=fit

plot_random_effect(fit$summary.random$month)

range(fit$summary.random$id.loc$mean)
plot(loc.sub.sub,cex=.5*exp(.5*fit$summary.random$id.loc$mean/mean(fit$summary.random$id.loc$mean)))
