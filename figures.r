

links <- read.table(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/links.txt",header=T) 
monitored.links <- scan(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/monitored-links.txt")
D <- scan(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/D.txt") 
nodes <- read.table(file="https://raw.githubusercontent.com/MartinLHazelton/Leicester/main/nodes.txt",header=T)


library(transportation)

# FIGURE 1

node.x <- 1.25*c(2,4,6,8,10,12,14,16,18,16,2,0,2,4,6,10,16,2,6,10,2)
node.y <- 1.25*(10-c(2,2,2,2,2,2,2,2,2,0,4,6,6,6,6,6,6,8,8,8,10))
width.to.height <- 1.5
COL <- rep(grey(0.6),nrow(links))
COL[monitored.links] <- "black"
LWD <- rep(1,nrow(nodes))
LWD[unique(D)] <- 3
pdf(file="../tex/art/fig1.pdf",width=8,height=8/width.to.height)
par(mar=c(0.2,0.2,0.2,0.2))
plotnet(node.x,node.y,as.matrix(links[,2:3]),ARROW.LENGTH=0.2,link.sep=1/6,LABEL.CEX=0.7,link.lty=1,node.lwd=LWD,link.col=COL)
text(11,11.1,"London Road",cex=0.8)
text(7.5,1.6,"Lancaster Road",cex=0.8)
text(1,7.5,"Waterloo Way",srt=90,cex=0.8)
text(21.4,7.5,"Granville Road",srt=90,cex=0.8)
dev.off()


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

# FIGURE 2

pdf(file="../tex/art/fig2.pdf",width=6,height=6)
par(mar=c(0.2,0.2,0.2,0.2))
CEX <- 1
pairs(mean.theta.df, cex = CEX, pch = 1, horOdd=TRUE, diag.panel = panel.16.hist, cex.labels = 1.5, font.labels = 2,  labels=c(expression(gamma==0),expression(gamma==4),expression(gamma==16),expression(gamma==100))) #,log="xy")
dev.off()


# FIGURE 3

INDX <- which(priormean$Prior.Mean>120)
pdf(file="../tex/art/fig3.pdf",width=8,height=4)
par(mar=c(4.2,4.2,0.2,0.2))
plot(c(0.5,10.5),c(0,800),type="n",xlab="OD pair",ylab=expression(paste("Mean traffic volume,","  ",theta)),axes=F)
axis(1,at=1:10,labels=as.character(INDX),tick=FALSE)
axis(2)
box()
for (i in 1:length(INDX)){
	segments(i-0.2,qgamma(0.025,priormean$Prior.Mean[INDX[i]],1),i-0.2,qgamma(0.975,priormean$Prior.Mean[INDX[i]],1),lwd=2)
	segments(i-0.1,lo.theta.df[INDX[i],1],i-0.1,hi.theta.df[INDX[i],1],lwd=2)
	segments(i,lo.theta.df[INDX[i],2],i,hi.theta.df[INDX[i],2],lwd=2)	
	segments(i+0.1,lo.theta.df[INDX[i],3],i+0.1,hi.theta.df[INDX[i],3],lwd=2)	
	segments(i+0.2,lo.theta.df[INDX[i],4],i+0.2,hi.theta.df[INDX[i],4],lwd=2)
#	points(i+c(-2,-1,0,1,2)/10,c(priormean$Prior.Mean[INDX[i]],mean.theta.df[INDX[i],]),pch=1:5)
#	points(i+c(-2,-1,0,1,2)/10,c(priormean$Prior.Mean[INDX[i]],mean.theta.df[INDX[i],]),pch=19,col=c(grey(0.1),grey(0.3),grey(0.5),grey(0.7),grey(0.9)))
points(i+c(-2,-1,0,1,2)/10,c(priormean$Prior.Mean[INDX[i]],mean.theta.df[INDX[i],]),pch=19,col=1:5)
}
# legend("topright",pch=1:5,legend=c("Prior",expression(gamma==0),expression(gamma==4),expression(gamma==16),expression(gamma==100)))
# legend("topright",pch=19,legend=c("Prior",expression(gamma==0),expression(gamma==4),expression(gamma==16),expression(gamma==100)),col=c(grey(0.1),grey(0.3),grey(0.5),grey(0.7),grey(0.9)))
legend("topright",pch=19,legend=c("Prior",expression(gamma==0),expression(gamma==4),expression(gamma==16),expression(gamma==100)),col=1:5)
dev.off()



# FIGURE 4

pdf(file="../tex/art/fig4.pdf",width=6,height=6)
par(mar=c(0.2,0.2,0.2,0.2))
theta <- priormean$Prior.Mean
CEX <- 2*(rep(theta,table(OD))/max(theta))^0.5
pairs(mean.p.df, cex = CEX, pch = 19, horOdd=TRUE, diag.panel = panel.hist, cex.labels = 1.5, font.labels = 2,  labels=c(expression(gamma==0),expression(gamma==4),expression(gamma==16),expression(gamma==100)))
dev.off()


# FIGURE 5

xx <- seq(0.003,0.035,length=200)
fxx <- dnorm(xx,mean=0.02,sd=0.02/6)
pdf(file="../tex/art/fig5.pdf",width=6,height=6)
#par(mar=c(0.2,0.2,0.2,0.2))
par(mfrow=c(2,2))
plot(xx,fxx,main="Prior",xlab=expression(omega),ylab="density",type="l",xlim=range(xx))
abline(v=0.02,lty=2)
plot(density(leic.4$OTHER.PARS),main=expression(paste("Posterior with"," ",gamma==4)),xlab=expression(omega),ylab="density",xlim=range(xx))
abline(v=0.02,lty=2)
plot(density(leic.16$OTHER.PARS),main=expression(paste("Posterior with"," ",gamma==16)),xlab=expression(omega),ylab="density",xlim=range(xx))
abline(v=0.02,lty=2)
plot(density(leic.100$OTHER.PARS),main=expression(paste("Posterior with"," ",gamma==100)),xlab=expression(omega),ylab="density",xlim=range(xx))
abline(v=0.02,lty=2)
dev.off()




# FIGURE 6

pdf(file="../tex/art/fig6.pdf",width=7,height=7)
par(mfrow=c(2,2))
par(mar=c(4,4,2,1))
theta <- priormean$Prior.Mean
CEX <- 2*(rep(theta,table(OD))/max(theta))^0.5
plot(mean.psi.df$gamma0,mean.p.df$gamma0,cex=CEX,xlab=expression(psi),ylab=expression(p),pch=19)
title(expression(gamma==0))
plot(mean.psi.df$gamma4,mean.p.df$gamma4,cex=CEX,xlab=expression(psi),ylab=expression(p),pch=19)
title(expression(gamma==4))
plot(mean.psi.df$gamma16,mean.p.df$gamma16,cex=CEX,xlab=expression(psi),ylab=expression(p),pch=19)
title(expression(gamma==16))
plot(mean.psi.df$gamma100,mean.p.df$gamma100,cex=CEX,xlab=expression(psi),ylab=expression(p),pch=19)
title(expression(gamma==100))
dev.off()



par(mfrow=c(2,2))
par(mar=c(4,4,2,1))
theta <- priormean$Prior.Mean
theta <- rep(theta,routes.per.OD)
fitted.0 <- theta*mean.psi.df$gamma0
indx <- rep(routes.per.OD,routes.per.OD) > 1
pearson.resid.0 <- theta*(mean.p.df$gamma0-mean.psi.df$gamma0)/sqrt(theta*mean.psi.df$gamma0)
plot(fitted.0[indx],pearson.resid.0[indx],xlab="SUE route flow",ylab="Discrepancy",pch=19,log="x")
title(expression(gamma==0))
fitted.4 <- theta*mean.psi.df$gamma4
pearson.resid.4 <- theta*(mean.p.df$gamma4-mean.psi.df$gamma4)/sqrt(theta*mean.psi.df$gamma4)
plot(fitted.4[indx],pearson.resid.4[indx],xlab="SUE route flow",ylab="Discrepancy",pch=19,log="x")
title(expression(gamma==4))
fitted.16 <- theta*mean.psi.df$gamma16
pearson.resid.16 <- theta*(mean.p.df$gamma16-mean.psi.df$gamma16)/sqrt(theta*mean.psi.df$gamma16)
plot(fitted.16[indx],pearson.resid.16[indx],xlab="SUE route flow",ylab="Discrepancy",pch=19,log="x")
title(expression(gamma==16))
fitted.100 <- theta*mean.psi.df$gamma100
pearson.resid.100 <- theta*(mean.p.df$gamma100-mean.psi.df$gamma100)/sqrt(theta*mean.psi.df$gamma100)
plot(fitted.100[indx],pearson.resid.100[indx],xlab="SUE route flow",ylab="Discrepancy",pch=19,log="x")
title(expression(gamma==100))




