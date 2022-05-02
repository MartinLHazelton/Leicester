
# Load packages

# devtools::install_github("MartinLHazelton/LinInvCount") # Only run when package needs updating
# devtools::install_github("MartinLHazelton/transportation") # Only run when package needs updating

library(LinInvCount)
library(transportation)
library(extraDistr)
library(lpSolve)
source("functions.r")

# Setting up the network

A <- as.matrix(read.table(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/A.txt"))
Afull <- as.matrix(read.table(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/Afull.txt"))
monitored.links <- scan(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/monitored-links.txt")
D <- scan(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/D.txt") 
O <- scan(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/O.txt") 
links <- read.table(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/links.txt",header=T) 
nodes <- read.table(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/nodes.txt",header=T)
priormean <- read.table(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/priormean.txt",header=T)
counts <- as.matrix(read.table(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/counts.txt",row.names=1,header=T))

theta <- priormean$Prior.Mean
ODpair <- factor(O+D/1000)
OD <- as.numeric(ODpair)
routes.per.OD <- table(OD)
Alpha <- links$Cost0  # Convert from 0.01 hours to minutes
Beta <- links$Capacity
r <- ncol(A)
n <- nrow(A)
d <- length(theta)

ndays <- ncol(counts)


set.seed(2022)

theta.prior <- cbind(priormean$Prior.Mean,rep(1,d)) 
theta.prior[theta.prior==0] <- 0.1
theta.ini <- theta.prior[,1]
prior.mean.logit.par <- 0.02
prior.sd.logit.par <- prior.mean.logit.par/6

lambda.additional <- list()
lambda.additional$OD <- OD
lambda.additional$Afull <- Afull
lambda.additional$Alpha <- Alpha
lambda.additional$Beta <- Beta
lambda.additional$routes.per.OD <- routes.per.OD
lambda.additional$BLOCK <- TRUE
lambda.additional$theta.prior <- theta.prior
lambda.additional$logit.par.prior <- c(prior.mean.logit.par,prior.sd.logit.par)
lambda.additional$delta.logit.par <- prior.sd.logit.par/2
other.pars <- prior.mean.logit.par
lambda.additional$FIX <- F
NDRAWS <- 100000
BURNIN <- 2000
THIN <- 20
lambda.ini <- SUE(ODdemand=theta.ini,ODpair=OD,A=Afull,Alpha=Alpha,Beta=Beta,theta=prior.mean.logit.par,tune.d=2)$x
x.ini.1 <- calculate.x.ini(round(rowMeans(counts)),A,lambda.ini)
x.ini <- matrix(x.ini.1,nrow=r,ncol=ndays,byrow=F)

lambda.additional$gamma <- 0
res.time.0 <- system.time(leic.0 <- Xlambdasampler(y=counts,A=A,lambda.updater=lu,lambda.ini=lambda.ini,Method="MH",x.ini=x.ini,lambda.additional=lambda.additional,other.pars=other.pars,ndraws=NDRAWS,burnin=BURNIN,THIN=THIN))

lambda.additional$gamma <- 4
res.time.4 <- system.time(leic.4 <- Xlambdasampler(y=counts,A=A,lambda.updater=lu,lambda.ini=lambda.ini,Method="MH",x.ini=x.ini,lambda.additional=lambda.additional,other.pars=other.pars,ndraws=NDRAWS,burnin=BURNIN,THIN=THIN))

lambda.additional$gamma <- 16
res.time.16 <- system.time(leic.16 <- Xlambdasampler(y=counts,A=A,lambda.updater=lu,lambda.ini=lambda.ini,Method="MH",x.ini=x.ini,lambda.additional=lambda.additional,other.pars=other.pars,ndraws=NDRAWS,burnin=BURNIN,THIN=THIN))

lambda.additional$gamma <- 100
res.time.100 <- system.time(leic.100 <- Xlambdasampler(y=counts,A=A,lambda.updater=lu,lambda.ini=lambda.ini,Method="MH",x.ini=x.ini,lambda.additional=lambda.additional,other.pars=other.pars,ndraws=NDRAWS,burnin=BURNIN,THIN=THIN))


theta.summary.0
theta.summary.4
theta.summary.16
theta.summary.100

p.summary.0 <- p.summary(leic.0$LAMBDA)
p.summary.4 <- p.summary(leic.4$LAMBDA)
p.summary.16 <- p.summary(leic.16$LAMBDA)
p.summary.100 <- p.summary(leic.100$LAMBDA)

mean.theta.df <- data.frame(gamma0=theta.summary.0$Mean,gamma4=theta.summary.4$Mean,gamma16=theta.summary.16$Mean,gamma100=theta.summary.100$Mean)
lo.theta.df <- data.frame(gamma0=theta.summary.0$Lower,gamma4=theta.summary.4$Lower,gamma16=theta.summary.16$Lower,gamma100=theta.summary.100$Lower)
hi.theta.df <- data.frame(gamma0=theta.summary.0$Upper,gamma4=theta.summary.4$Upper,gamma16=theta.summary.16$Upper,gamma100=theta.summary.100$Upper)

mean.p.df <- data.frame(gamma0=p.summary.0$Mean,gamma4=p.summary.4$Mean,gamma16=p.summary.16$Mean,gamma100=p.summary.100$Mean)

psi.summary.0 <- psi.summary(leic.0$LAMBDA,omega=leic.0$OTHER.PARS,OD=OD,Afull=Afull,Alpha=Alpha,Beta=Beta,routes.per.OD=routes.per.OD)
psi.summary.4 <- psi.summary(leic.4$LAMBDA,omega=leic.4$OTHER.PARS,OD,Afull,Alpha,Beta,routes.per.OD)
psi.summary.16 <- psi.summary(leic.16$LAMBDA,omega=leic.16$OTHER.PARS,OD,Afull,Alpha,Beta,routes.per.OD)
psi.summary.100 <- psi.summary(leic.100$LAMBDA,omega=leic.100$OTHER.PARS,OD,Afull,Alpha,Beta,routes.per.OD)

mean.psi.df <- data.frame(gamma0=psi.summary.0$Mean,gamma4=psi.summary.4$Mean,gamma16=psi.summary.16$Mean,gamma100=psi.summary.100$Mean)


