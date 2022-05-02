# Lambda updater

lu <- function(x,lambda,lambda.additional,other.pars=numeric(0)){
	if (is.vector(x)) x <- matrix(x,ncol=1)
	OD <- lambda.additional$OD
	Afull <- lambda.additional$Afull
	a <- lambda.additional$theta.prior[,1]
	b <- lambda.additional$theta.prior[,2]
	FIX <- lambda.additional$FIX
 	BLOCK <- lambda.additional$BLOCK
	if (is.null(FIX)) FIX <- F
	if (is.null(BLOCK)) BLOCK <- T
	prior.mean.logit.par <- lambda.additional$logit.par.prior[1]
	prior.sd.logit.par <- lambda.additional$logit.par.prior[2]
	Alpha <- lambda.additional$Alpha
	Beta <- lambda.additional$Beta
	gamma <- lambda.additional$gamma
	routes.per.OD <- lambda.additional$routes.per.OD
	if (length(other.pars)==0) { logit.par <- lambda.additional$logit.par } else { logit.par <- other.pars } 
	d <- length(a)
	T <- ncol(x)
	z <- as.matrix(aggregate(x,list(OD=OD),sum)[,-1])
	xbar <- apply(x,1,mean)
	zbar <- tapply(xbar,OD,sum)
	theta <- c(tapply(lambda,OD,sum))
	p <- lambda/rep(theta,routes.per.OD)
	psi <- lambda.additional$psi
	if (is.null(psi)) {
		psi <- TrafficAssignment(ODdemand=theta,ODpair=OD,A=Afull,Alpha=Alpha,Beta=Beta,logit.par=logit.par,x.ini=lambda,routes.per.OD=routes.per.OD,FIX=FIX)
	}
	if (BLOCK){
		theta.cand <- rgamma(d,zbar*T+a,rate=T+b)
		psi.cand <- TrafficAssignment(ODdemand=theta.cand,ODpair=OD,A=Afull,Alpha=Alpha,Beta=Beta,logit.par=logit.par,x.ini=psi*rep(theta.cand,routes.per.OD),routes.per.OD=routes.per.OD,FIX=FIX)
		if (!FIX) {AR <- exp( sum(gamma*(psi.cand-psi)*log(p)-lgamma(1+gamma*psi.cand)+lgamma(1+gamma*psi)) )} else {AR <- 1}  # Ratio for Dirichlet parameters 
		if (!is.na(AR)){
			if (runif(1) < AR){
				theta <- theta.cand
				psi <- psi.cand
			}
		}
		for (i in 1:d){
			if (sum(OD==i)>1) { 
				p[OD==i] <- rdirichlet(1,1+T*xbar[OD==i]+gamma*psi[OD==i])
			} else { p[OD==i] = 1 }
		}	
		lambda <- p*rep(theta,routes.per.OD)
	}
	if (!BLOCK){
		for (i in 1:d){
			theta.cand <- theta
			theta.cand[i] <- rgamma(1,zbar[i]*T+a,rate=T+b)
			psi.cand <- TrafficAssignment(ODdemand=theta.cand,ODpair=OD,A=Afull,Alpha=Alpha,Beta=Beta,logit.par=logit.par,x.ini=psi*rep(theta.cand,routes.per.OD),routes.per.OD=routes.per.OD,FIX=FIX)
			if (!FIX) {AR <- exp( sum(gamma*(psi.cand-psi)*log(p)-lgamma(1+gamma*psi.cand)+lgamma(1+gamma*psi)) )} else {AR <- 1}  # Ratio for Dirichlet parameters 
			if (!is.na(AR)){
				if (runif(1) < AR){
					theta <- theta.cand
					psi <- psi.cand
				}
			}
			if (sum(OD==i)>1) { 
				p[OD==i] <- rdirichlet(1,1+T*xbar[OD==i]+gamma*psi[OD==i])
			} else { p[OD==i] = 1 }
			lambda <- p*rep(theta,routes.per.OD)
		}
	}
	if (FIX) logit.par <- rnorm(1,prior.mean.logit.par,sd=prior.sd.logit.par)
	if (!FIX){
		logit.par.cand <- rnorm(1,logit.par,sd=lambda.additional$delta.logit.par)
		if (logit.par.cand > 0){
			psi.cand <- TrafficAssignment(ODdemand=theta,ODpair=OD,A=Afull,Alpha=Alpha,Beta=Beta,logit.par=logit.par.cand,x.ini=psi*rep(theta,routes.per.OD),routes.per.OD=routes.per.OD,FIX=FIX)
			AR <- exp( dnorm(logit.par.cand,prior.mean.logit.par,sd=prior.sd.logit.par,log=T) - dnorm(logit.par,prior.mean.logit.par,sd=prior.sd.logit.par,log=T) +sum(gamma*(psi.cand-psi)*log(p)-lgamma(1+gamma*psi.cand)+lgamma(1+gamma*psi)) ) 
			if (!is.na(AR)){
				if (runif(1) < AR){
					logit.par <- logit.par.cand
					psi <- psi.cand
				}
			}
		}
	}
	lambda.additional$psi <- psi
	list(lambda=lambda,lambda.additional=lambda.additional,other.pars=logit.par)
}


TrafficAssignment <- function(ODdemand,ODpair,A,Alpha,Beta,logit.par,x.ini=NULL,routes.per.OD,FIX=FALSE){
	if (FIX) { res <- sue.fixed } else
	res <- SUE(ODdemand=ODdemand,ODpair=ODpair,A=A,Alpha=Alpha,Beta=Beta,theta=logit.par,x.ini=x.ini,tune.d=2)$x/rep(ODdemand,routes.per.OD)
	res
}

simulations <- function(A,theta,gamma,logit.par,lambda.additional,NDRAWS=10000,BURNIN=2000,THIN=10,SEED=2022){
	set.seed(SEED)
	OD <- lambda.additional$OD 
	Afull <- lambda.additional$Afull 
	theta.prior <- lambda.additional$theta.prior 
	Alpha <- lambda.additional$Alpha
	Beta <- lambda.additional$Beta
	routes.per.OD <- lambda.additional$routes.per.OD
	delta.logit.par <- lambda.additional$delta.logit.par
	logit.par.prior <- lambda.additional$logit.par.prior 
	lambda.additional$FIX <- F
	lambda.additional$gamma <- gamma
	prior.mean.log.par <- logit.par.prior[1]
	other.pars <- prior.mean.logit.par
	lambda <- SUE(ODdemand=theta,ODpair=OD,A=Afull,Alpha=Alpha,Beta=Beta,theta=logit.par,tune.d=2)$x
	x <- rpois(ndays*r,lambda)
	y <- A%*%matrix(x,ncol=ndays,byrow=F)
	theta.ini <- theta.prior[,1]/theta.prior[,2]
	lambda.ini <- SUE(ODdemand=theta.ini,ODpair=OD,A=Afull,Alpha=Alpha,Beta=Beta,theta=prior.mean.logit.par,tune.d=2)$x
	x.ini.1 <- calculate.x.ini(round(rowMeans(y)),A,lambda.ini)
	x.ini <- matrix(x.ini.1,nrow=r,ncol=ndays,byrow=F)
	res.time <- system.time(res <- Xlambdasampler(y=y,A=A,lambda.updater=lu,lambda.ini=lambda.ini,Method="MH",x.ini=x.ini,lambda.additional=lambda.additional,other.pars=other.pars,ndraws=NDRAWS,burnin=BURNIN,THIN=THIN))
	lambda.additional$FIX <- T
	res2.time <- system.time(res2 <- Xlambdasampler(y=y,A=A,lambda.updater=lu,lambda.ini=lambda.ini,Method="MH",x.ini=x.ini,lambda.additional=lambda.additional,other.pars=other.pars,ndraws=NDRAWS,burnin=BURNIN,THIN=THIN))
	list(res=res,res2=res2,res.time=res.time,res2.time=res2.time)		
}

		

calculate.x.ini <- function(y,A,lambda.ini,WEIGHTS=-1){
	if (!is.integer(y)) y <- as.integer(y)
	r <- length(lambda.ini)
	n <- length(y)
	if (any(WEIGHTS<0)) WEIGHTS <- c(rep(1,r),rep(100,n))
	X <- rbind(diag(r),A)
	yy <- c(lambda.ini,y)
	lambda.hat <- coef(lm(yy ~ -1 + X, weights=WEIGHTS))
	x1 <- floor(lambda.hat)
	x1[x1<0] <- 0
	while (any(y-A%*%x1 < 0)) x1 <- floor(x1*0.99)
	y2 <- y-A%*%x1
	z <- lp("min", objective.in = rep(1,r), const.mat = A, const.dir = rep("=", nrow(A)), const.rhs = y2, all.int = T)$solution
	x1+z
}

q975 <- function(x){quantile(x,p=0.975)}
q025 <- function(x){quantile(x,p=0.025)}

plot.results <- function(theta,RES,RESfix,BURNIN=0,THIN=1){
	theta.res <- as.matrix(aggregate(RES,list(OD=OD),sum))[,-1]
	theta.fix <- as.matrix(aggregate(RESfix,list(OD=OD),sum))[,-1]
	theta.res <- theta.res[,-(1:round(BURNIN/THIN))]
	theta.fix <- theta.fix[,-(1:round(BURNIN/THIN))]
	theta.res.hi <- apply(theta.res,1,q975)
	theta.res.lo <- apply(theta.res,1,q025)
	theta.fix.hi <- apply(theta.fix,1,q975)
	theta.fix.lo <- apply(theta.fix,1,q025)
	par(mfrow=c(1,2))
	plot(theta,rowMeans(theta.res))
	abline(1,1)
	segments(theta,theta.res.lo,theta,theta.res.hi)
	plot(theta,rowMeans(theta.fix))
	abline(1,1)
	segments(theta,theta.fix.lo,theta,theta.res.hi)
	par(mfrow=c(1,1))
}

trace.plot <- function(RES,INDEX){
	theta.res <- as.matrix(aggregate(RES,list(OD=OD),sum))[INDEX,-1]
	par(mfrow=c(length(INDEX),1))
	for (i in 1:length(INDEX)){
		plot(c(theta.res[i,]),type="l")
	}
}

results.summary <- function(RES,theta){
	theta.res <- as.matrix(aggregate(RES,list(OD=OD),sum))[,-1]
	theta.res.hi <- apply(theta.res,1,q975)
	theta.res.lo <- apply(theta.res,1,q025)
	theta.res.mean <- apply(theta.res,1,mean)
	MSE <- mean((theta.res.mean-theta)^2)
	coverage <- 100*mean((theta-theta.res.lo)*(theta-theta.res.hi) < 0)
	list(MSE=MSE,coverage=coverage)
}

theta.summary <- function(RES){
	theta.res <- as.matrix(aggregate(RES,list(OD=OD),sum))[,-1]
	theta.res.hi <- apply(theta.res,1,q975)
	theta.res.lo <- apply(theta.res,1,q025)
	theta.res.mean <- apply(theta.res,1,mean)
	list(Mean=theta.res.mean,Lower=theta.res.lo,Upper=theta.res.hi)
}

p.summary <- function(RES){
	theta.res <- as.matrix(aggregate(RES,list(OD=OD),sum))[,-1]
	p.res <- RES/theta.res[OD,]
	p.res.hi <- apply(p.res,1,q975)
	p.res.lo <- apply(p.res,1,q025)
	p.res.mean <- apply(p.res,1,mean)
	list(Mean=p.res.mean,Lower=p.res.lo,Upper=p.res.hi)
}

psi.summary <- function(RES,omega,OD,Afull,Alpha,Beta,routes.per.OD){
	psi <- matrix(NA,nrow=ncol(Afull),ncol=length(omega))
	theta.res <- as.matrix(aggregate(RES,list(OD=OD),sum))[,-1]
	x.ini <- NULL
	for (i in 1:length(omega)){
		psi[,i] <- TrafficAssignment(ODdemand=theta.res[,i],ODpair=OD,A=Afull,Alpha=Alpha,Beta=Beta,logit.par=omega[i],x.ini=x.ini,routes.per.OD=routes.per.OD,FIX=F)
		if (i < length(omega)) x.ini <- psi[,i]*rep(theta.res[,i+1],routes.per.OD)
	}
	psi.res.hi <- apply(psi,1,q975)
	psi.res.lo <- apply(psi,1,q025)
	psi.res.mean <- apply(psi,1,mean)
	list(Mean=psi.res.mean,Lower=psi.res.lo,Upper=psi.res.hi)
}


panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.25) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = grey(0.8), ...)
}

panel.16.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.25) )
    h <- hist(x, plot = FALSE, breaks=16)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = grey(0.8), ...)
}

panel.log.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.25) )
    h <- hist(log(x), plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(exp(breaks[-nB]), 0, exp(breaks[-1]), y, col = grey(0.8), ...)
}

panel.dens <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, .01) )
  hist(x, freq = FALSE, col="cyan", add=TRUE) 
  lines(density(x))
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

