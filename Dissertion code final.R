install.packages("data.table")
library(data.table)
install.packages("tidyverse")
library(tidyverse)
install.packages("xts")
library(xts)
install.packages("dplyr")
library(dplyr)
install.packages("zoo")
library(zoo)
install.packages("moveHMM")
library(moveHMM)
install.packages("momentuHMM")
library(momentuHMM)
#install.packages("CircStats")
#library(CircStats)
library(ggplot2)
install.packages("parallel")
library(parallel)
#library(circular)
#library(stats)

ncores <- detectCores() - 1

#### Importing data

#setwd("C:/Users/Mike Goodman/Documents/Sussex/Dissertation/Marsh Harrier text files/")
setwd("C:/Users/Owner/Documents/")
rawdata <- read.csv(file='MH indiv 5325667, GPS 2017, 1 Jan - 31 Dec.csv')
head(rawdata)

df_req_fields <- subset(rawdata,select=c("individual.local.identifier","timestamp","location.long","location.lat","heading","height.above.msl"))
head(df_req_fields)

df_req_fields$timestamp<- as.POSIXct(df_req_fields$timestamp,format="%Y-%m-%d %H:%M:%OS")

#### Standardising data ####

# dataframe of year in 5-min intervals
reg_times<-seq(from=as.POSIXct("2017-01-01 00:00:00",tz="CET"),
               to=as.POSIXct("2018-1-1 00:00:00",tz="CET"),by=60*5)
reg_times<-as.data.frame(reg_times)
names(reg_times)[1]<-"timestamp"

empty_cols<-c("individual.local.identifier","location.long","location.lat","ground.speed","heading","height.above.msl")
reg_times[,empty_cols]<-NA
reg_times <- reg_times[,c("individual.local.identifier","timestamp","location.long","location.lat","heading","height.above.msl")]
head(reg_times)

comb_df<-bind_rows(df_req_fields,reg_times)
comb_df<-comb_df[order(comb_df$timestamp),]

comb_z<-na.locf(comb_df, fromLast=FALSE)
head(comb_z,50)

tail(comb_z,50)
raw_times<-list(df_req_fields$timestamp)

reg_data <- filter(comb_z, comb_z$timestamp %in% reg_times$timestamp)
head(reg_data,50)

names(reg_data)[names(reg_data) == "location.long"] <- "x" 
names(reg_data)[names(reg_data) == "location.lat"] <- "y"
names(reg_data)[names(reg_data) == "individual.local.identifier"] <- "ID"
names(reg_data)[names(reg_data) == "height.above.msl"] <- "height"

state_names<-c("rest","forage","migrate")
prepd_data<-prepData(data=reg_data, type="LL",coordNames=c("x","y")) #, covNames=c("height")

# empirical histograms

dev.new(width=5, height=4, unit="in")
# hist(prepd_data$step, freq=FALSE, col="#FFFF66",
#      main="Histogram of Step Lengths with Fitted Density Function", xlab="Step Length", ylab="Density",
#      breaks=2000)

dev.new(width=5, height=4, unit="in")
hist(prepd_data$angle, freq=FALSE, col="#FFFF66",
     main="Histogram of Turning Angles with Fitted Density Function", xlab="Step Length", ylab="Density",
     xlim=c(-pi,pi), ylim=c(0,0.3),breaks=30)
lines(density(prepd_data$angle[!is.na(prepd_data$angle)]),lwd=2,col="blue")


#### von Mises model

# defining start values
mu01 <- c(0.01,10,100)
sigma01 <- c(5,10,25)
zeromass01 <- c(0.1,0.1,0.1)
stepPar01 <- c(mu01, sigma01, zeromass01)
angleMean01 <- c(0,0,0)
conc01 <- c(0.05,2,10)
anglePar01 <-c(angleMean01, conc01) 

models01 <- fitHMM(data=prepd_data, nbStates=3, dist = list(step="weibull", angle = "vm"),
                   Par0 = list(step = stepPar01, angle = anglePar01),
                   estAngleMean=list(angle=TRUE), stateNames = state_names)

#### main von Mises plot
vM_range <- seq(-pi,pi,(pi--pi)/(length(models01$data$angle)-1))
vM_data_rest<-dvm(vM_range, mu=models01$mle$angle["mean",1],kappa=models01$mle$angle["concentration",1])
vM_data_forage<-dvm(vM_range, mu=models01$mle$angle["mean",2],kappa=models01$mle$angle["concentration",2])
vM_data_migrate<-dvm(vM_range, mu=models01$mle$angle["mean",3],kappa=models01$mle$angle["concentration",3])                   #control.circular=list(units="radians"))

dev.new(width=5, height=4, unit="in")
plot(vM_range, vM_data_rest,
     xlab="Turning Angle (\u03b8)",ylab="Density",main="Fitted Von Mises PDFs for All Three States", cex.main=1.1,
     xlim=c(-pi,pi), ylim=c(0,max(vM_data_migrate)), type="l", lty=1, lwd=0.001,col="forestgreen")
lines(vM_range, vM_data_forage, lty=1,lwd=0.5, col="navyblue")
lines(vM_range, vM_data_migrate, lty=1, lwd=0.00000001, col="red")
# legend text
vM_legend_rest <- paste("rest: \U03bc=", round(models01$mle$angle["mean",1],4), ", \u03ba=", round(models01$mle$angle["concentration",1],4),sep="")
vM_legend_forage <- paste("forage: \U03bc=", round(models01$mle$angle["mean",2],4), ", \u03ba=", round(models01$mle$angle["concentration",2],4),sep="")
vM_legend_migrate <- paste("migrate: \U03bc=", round(models01$mle$angle["mean",3],4), ", \u03ba=", round(models01$mle$angle["concentration",3],4),sep="")
vM_legend_text <- c(vM_legend_rest, vM_legend_forage, vM_legend_migrate)
op <- par(cex=0.8)
vM_legend<-legend(x="topright",legend=vM_legend_text,lwd=c(1,1,1),col=c("forestgreen","navyblue","red"),)


#### Weibull plots

Weibull_range <- seq(0,max(models01$mle$step)+0.1,(max(models01$mle$step)+0.1)/(length(models01$data$step)-1))
Wb_data_rest<-dweibull(Weibull_range, shape=models01$mle$step["shape",1],scale=models01$mle$step["scale",1])
Wb_data_forage<-dweibull(Weibull_range, shape=models01$mle$step["shape",2],scale=models01$mle$step["scale",2])
Wb_data_migrate<-dweibull(Weibull_range, shape=models01$mle$step["shape",3],scale=models01$mle$step["scale",3])                   #control.circular=list(units="radians"))

# weibull rest plot
dev.new(width=5, height=4, unit="in")
plot(Weibull_range, Wb_data_rest,
     xlab="Step Length",ylab="Density",main="Fitted Weibull PDF for Rest State", cex.main=1.1,
     xlim=c(0,0.05), #ylim=c(0,max(Wb_data3),
     type="l", lty=1, lwd=0.001,col="forestgreen")
Weibull_legend_rest <- paste("\u03ba=", round(models01$mle$step["shape",1],4), ", \u03bb=", round(models01$mle$step["scale",1],4),sep="")
op <- par(cex=1.0)
wb_rest_legend <- legend(x="topright",legend=Weibull_legend_rest, lwd=1,col="forestgreen")

# Weibull forage plot
dev.new(width=5, height=4, unit="in")
plot(Weibull_range, Wb_data_forage,
     xlab="Step Length",ylab="Density",main="Fitted Weibull PDF for Forage State", cex.main=1.1,
     xlim=c(0,1), #ylim=c(0,max(Wb_data3),
     type="l", lty=1, lwd=0.001,col="navyblue")
Weibull_legend_forage <- paste("\u03ba=", round(models01$mle$step["shape",2],4), ", \u03bb=", round(models01$mle$step["scale",2],4),sep="")
op <- par(cex=1.0)
wb_rest_legend <- legend(x="topright",legend=Weibull_legend_forage, lwd=1,col="navyblue")

# Weibull migrate plot
dev.new(width=5, height=4, unit="in")
plot(Weibull_range, Wb_data_migrate,
     xlab="Step Length",ylab="Density",main="Fitted Weibull PDF for Migrate State", cex.main=1.1,
     xlim=c(0,50), #ylim=c(0,max(Wb_data3),
     type="l", lty=1, lwd=0.001,col="red")
Weibull_legend_migrate <- paste("\u03ba=", round(models01$mle$step["shape",3],4), ", \u03bb=", round(models01$mle$step["scale",3],4),sep="")
op <- par(cex=1.0)
wb_rest_migrate <- legend(x="topright",legend=Weibull_legend_migrate, lwd=1,col="red")

## Maximum Likelihood Estimators
Wb_shape_rest  <- models01$mle$step["shape",1]
Wb_shape_forage <- models01$mle$step["shape",2]
Wb_shape_migrate <- models01$mle$step["shape",3]
Wb_scale_rest <- models01$mle$step["scale",1]
Wb_scale_forage <- models01$mle$step["scale",2]
Wb_scale_migrate <- models01$mle$step["scale",3]

# caculating means from MLEs
Weibull_mean <- function(scale_param, shape_param){
  wb_mean <- scale_param*gamma(1+1/shape_param)
  return(wb_mean)
}
StepMean_rest <- Weibull_mean(Wb_scale_rest, Wb_shape_rest)
StepMean_forage <- Weibull_mean(Wb_scale_forage, Wb_shape_forage)
StepMean_migrate <- Weibull_mean(Wb_scale_migrate, Wb_shape_migrate)

StepMean_rest
StepMean_forage
StepMean_migrate


#### Wrapped Cauchy model ##

m2_dist <- list(step="weibull", angle="wrpcauchy")
#defining start values
mu02 <- c(0.1,0.5,1)
sigma02 <- c(0.5,0.1,0.05)
zeromass02 <- c(0.1,0.1,0.1)
conc02 <- c(0.2,0.5,0.9)
stepPar02 <- c(mu02, sigma02, zeromass02)
angleMean02 <- c(0,0,0)
anglePar02 <- c(angleMean02, conc02)
Par0_m2 <-list(step=stepPar02, angle=anglePar02)

models02 <- fitHMM(data=prepd_data, nbStates=3, dist = m2_dist,
                   Par0 = Par0_m2, 
                   estAngleMean=list(angle=TRUE), stateNames = state_names)

WC_range <- seq(-pi,pi,(pi--pi)/(length(models02$data$angle)-1))
WC_data_rest<-dwrpcauchy(WC_range, mu=models02$mle$angle["mean",1],rho=models02$mle$angle["concentration",1])
WC_data_forage<-dwrpcauchy(WC_range, mu=models02$mle$angle["mean",2],rho=models02$mle$angle["concentration",2])
WC_data_migrate<-dwrpcauchy(WC_range, mu=models02$mle$angle["mean",3],rho=models02$mle$angle["concentration",3])                   #control.circular=list(units="radians"))

dev.new(width=5, height=4, unit="in")
plot(WC_range, WC_data_rest,
     xlab="Turning Angle (\u03b8)",ylab="Density",main="Fitted Wrapped Cauchy PDFs for All Three States",
     xlim=c(-pi,pi), ylim=c(0,max(WC_data_migrate)), type="l", lty=1, lwd=0.001,col="forestgreen")
lines(WC_range, WC_data_forage, lty=1,lwd=0.5, col="navyblue")
lines(WC_range, WC_data_migrate, lty=1, lwd=0.00000001, col="red")
# legend text
WC_legend_rest <- paste("rest: \u03bc=", round(models02$mle$angle["mean",1],4), ", \u03c1=", round(models02$mle$angle["concentration",1],4),sep="")
WC_legend_forage <- paste("forage: \u03bc=", round(models02$mle$angle["mean",2],4), ", \u03c1=", round(models02$mle$angle["concentration",2],4),sep="")
WC_legend_migrate <- paste("migrate: \u03bc=", round(models02$mle$angle["mean",3],4), ", \u03c1=", round(models02$mle$angle["concentration",3],4),sep="")
WC_legend_text <- c(WC_legend_rest, WC_legend_forage, WC_legend_migrate)
op <- par(cex=0.9)
WC_legend<-legend(x="top",legend=WC_legend_text,lwd=c(1,1,1),col=c("forestgreen","navyblue","red"),)


###### Confidence Intervals

CIreal(models01, alpha=0.95, covs=NULL, parms=NULL)
CIreal(models02, alpha=0.95, covs=NULL, parms=NULL)

# dev.new(width=5, height=4, unit="in")
# plot(models01$data$step, ask=TRUE, sepStates=TRUE, ylim=c(0,10), plotCI=TRUE)

# dev.new(width=5, height=4, unit="in")
# plot(models02$data$step, ask=TRUE, sepStates=TRUE, ylim=c(0,10), plotCI=TRUE)

# model selection``

AIC(models01,models02)

## pseudo-residual plots

## pseudo-residuals - model 1

# qq - rest
plot.new() #dev.new(width=50, height=30, unit="in")
qplot(qweibull(ppoints(length(models01$data$step)),
               shape = models01$mle$step["shape",1],
               scale=models01$mle$step["scale",1]),
      models01$data$step) +
  labs(title="Rest Step Length", x="Theoretical", y="Sample") +
  theme(plot.title=element_text(hjust=0.5)) +
 ylim(0,10)
qqline(models01$data$step,col="forestgreen",lwd=1)

# qq forage
plot.new()
#dev.new(width=50, height=30, unit="in")
qplot(qweibull(ppoints(length(models01$data$step)),
               shape = models01$mle$step["shape",2],
               scale=models01$mle$step["scale",2]),
      models01$data$step) +
  labs(title="Forage Step Length", x="Theoretical", y="Sample") +
  theme(plot.title=element_text(hjust=0.5)) +
  ylim(0,20) +
  geom_point(alpha=0.5)
qqline(models01$data$step,col="navyblue",lwd=1)

# qq - migrate
plot.new()
#dev.new(width=50, height=30, unit="in")
qplot(qweibull(ppoints(length(models01$data$step)),
               shape = models01$mle$step["shape",3],
               scale=models01$mle$step["scale",3]),
      models01$data$step) +
  labs(title="Migrate Step Length", x="Theoretical", y="Sample") +
  theme(plot.title=element_text(hjust=0.5)) + 
  ylim(0,50)
qqline(models01$data$step,col="red",lwd=1)


dev.new(width=50, height=30, unit="in")
plotPR(models01, lag.max=24, ncores=ncores)

dev.new(width=100, height=30, unit="in")
plotPR(models02, lag.max=NULL, ncores=ncores)
