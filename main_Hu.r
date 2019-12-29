### Package ###############################################################################
library(xts) # xts
library(car) # qqPlot
library(QRM) # ESnorm
library(copula)
library(PerformanceAnalytics) # plot for p
library(EnvStats)




### (a) Data ###############################################################################
setwd("C:/Users/Eli/Desktop/429")
tmp = read.csv("429stocks.csv", stringsAsFactors=FALSE)
tmp$TICKER[tmp$TICKER=="SPWRA"] = "SPWR"
# obtain tickers and dates
tickers = unique(tmp$TICKER)
dates = as.Date(as.character( unique(tmp$date) ), format="%Y%m%d")
# returns for stocks
RET = c()
for (t in tickers){
  RET = cbind(RET, tmp$RET[tmp$TICKER==t])
}
colnames(RET) = tickers
MRET=RET
RET = log(RET+1) # RET_{t}:=S_{t+1}/S_{t}-1, X_{t}=log(S_{t+1}/S_{t})=log(RET_{t}+1)
RET = xts(RET, order.by=dates)
# obtain marketcap
tmp = tmp[tmp$date=="20181231",]
marketcap = matrix(tmp$PRC*tmp$SHROUT, nrow=1)/1000000
colnames(marketcap) = tickers
marketcap
# obtain statistics for RET
apply(RET,2,mean)*10000
apply(RET,2,sd)*100
apply(RET,2,skewness)
apply(RET,2,kurtosis)

# SP500 return
tmp = read.csv("429sp.csv", stringsAsFactors=FALSE)
SP500 = xts(tmp$vwretd, order.by=dates)
MSP500=SP500
SP500 = log(SP500+1)



### (b) Log Return ###############################################################################
par(mfrow=c(1,1))
plot.zoo(cbind(SP500,RET[,1:5]), screens=1:10, col=c("royalblue",rep(1,5)), las=1, main="Log Returns", xlab="",cex.main=2)
plot.zoo(cbind(SP500,RET[,6:10]), screens=1:10, col=c("royalblue",rep(1,5)), las=1, main="Log Returns", xlab="",cex.main=2)



### (c) Portfolio ###############################################################################

#' @title function for computing total loss of a portfolio
#' @param RET vector or matrix, returns
#' @param weights vector, weights of each asset
#' @param value numeric, total portfolio value
#' @param linear logical, indicator of whether the RET is linear-return or log-return
loss.portfolio = function(RET, weights, value, linear=FALSE){
  if (linear) { RET = RET } else { RET = (exp(RET)-1) }
  Profit = RET * value * weights
  return( -rowSums(Profit) )
}
# compute loss of portfolio
pf.value = 10000
pf.weights = rep(0.1,10)
alpha = c(0.95,0.975,0.99,0.995,0.999,0.9999)
loss = loss.portfolio(RET, pf.weights, pf.value)
# check normality
qqPlot(loss)
# normal VaR and ES
mu.hat = colMeans(RET)
sigma.hat = var(RET) #*(nrow(RET)-1)/nrow(RET)
meanloss = -sum(pf.weights*mu.hat) * pf.value
varloss = pf.value^2 * as.numeric(t(pf.weights) %*% sigma.hat %*% pf.weights)
VaR.normal = meanloss + sqrt(varloss) * qnorm(alpha)
ES.normal = meanloss + sqrt(varloss) * dnorm(qnorm(alpha))/(1-alpha) # ESnorm(0.95, mu=meanloss, sd=sqrt(varloss))
# historical VaR and ES
VaR.hs = quantile(loss,alpha)
ES.hs = apply(as.array(VaR.hs), 1, function(VaR){ mean(loss[loss > VaR]) } )
# histogram and comparing
par(mfrow=c(2,1))
par(mar = c(4.6, 4.1, 3.1, 2.1))
hist(loss,nclass=100, prob=TRUE, xlab="Portfolio Loss", main=paste("VaR at", paste(names(VaR.hs),collapse=", ")))
abline(v=VaR.hs,col=1,lty=2)
abline(v=VaR.normal,col=2,lty=5)
legend("topleft", legend=c("normal","HS"), col=1:2, lty=1, bty="n")
hist(loss,nclass=100, prob=TRUE, xlab="Portfolio Loss", main=paste("ES at", paste(names(VaR.hs),collapse=", ")))
abline(v=ES.hs,col=1,lty=2)
abline(v=ES.normal,col=2,lty=5)
legend("topleft", legend=c("normal","HS"), col=1:2, lty=1, bty="n")



### (d) Copulas ###############################################################################

## (d.1) Pseudo-observations ------------------------------------
RET = RET[,c("GE","ITT","PCG","IDA","SPWR","FSLR","TSLA")]
Uret = apply(RET, 2, edf, adjust=1) # pobs(RET, ties.method = "max")
# scatterplot
pairs2(RET, cex=0.3, col=adjustcolor("black",alpha.f=0.3))
pairs2(Uret, cex=0.3, col=adjustcolor("black",alpha.f=0.3)) # plot(as.matrix(Uret[,4:5]), pch=1) # CVA


## (d.2) Fit Copula --------------------------------------------
#' @title try to fitCopula or return NA
tryFitCopula = function(...){
  tryCatch(fitCopula(...), error=function(err) NA)
}
#' @title fit bivarite copulas including gauss,t,gumbel,clayton and frank
#' @param U matrix, bivarite pseudo observations
#' @return list, consists of different fitCopula class objects
fit.copulas = function(U){
  if (ncol(U)!=2) stop("ncol(U) is not 2")
  # gauss copula
  fit.gauss = tryFitCopula(normalCopula(dim=2,dispstr="un"), data=U, method="mpl")
  # gauss copula - method of moments: P^ = 2*sin(pi * Spearman/6))
  fit.gauss.spearman = tryFitCopula(normalCopula(dim=2,dispstr = "un"), data = U, method = "irho")
  fit.gauss.spearman@loglik = sum( dCopula(U, normalCopula(param=fit.gauss.spearman@estimate,dim=2,dispstr="un"), log=T) )
  # t copula
  fit.t = tryFitCopula(tCopula(dim=2,dispstr="un"), data=U, method="mpl")
  # t copula - method of moments: P^ = sin(pi * Kendall/2)
  fit.t.tau = tryFitCopula(tCopula(dim=2,dispstr="un"), data=U, method="itau.mpl")
  # AC gumbel copula
  fit.gumbel = tryFitCopula(gumbelCopula(dim=2), data=U, method="mpl")
  # AC gumbel copula - method of moments: tau
  fit.gumbel.tau = tryFitCopula(gumbelCopula(dim=2), data=U, method="itau")
  fit.gumbel.tau@loglik = sum( dCopula(U, gumbelCopula(param=fit.gumbel.tau@estimate,dim=2), log=T) )
  # AC clayton copula
  fit.clayton = tryFitCopula(claytonCopula(dim=2), data=U, method="mpl")
  # AC clayton copula - method of moments: tau
  fit.clayton.tau = tryFitCopula(claytonCopula(dim=2), data=U, method="itau")
  fit.clayton.tau@loglik = sum( dCopula(U, claytonCopula(param=fit.clayton.tau@estimate,dim=2), log=T) )
  # AC frank copula
  fit.frank = tryFitCopula(frankCopula(dim=2), data=U, method="mpl")
  # AC frank copula - method of moments: tau
  fit.frank.tau = tryFitCopula(frankCopula(dim=2), data=U, method="itau")
  fit.frank.tau@loglik = sum( dCopula(U, frankCopula(param=fit.frank.tau@estimate,dim=2), log=T) )
  res = list(gauss=fit.gauss, gauss.rho=fit.gauss.spearman
             , t=fit.t, t.tau=fit.t.tau
             , frank=fit.frank, frank.tau=fit.frank.tau
             , gumbel=fit.gumbel, gumbel.tau=fit.gumbel.tau
             , clayton=fit.clayton, clayton.tau=fit.clayton.tau)
}
# function to obtain p.value of gofCopula of fit using data U
gof.apply = function(fit, U){
  copula = tryCatch(fit@copula, error=function(err) NA)
  # for t-copula, df is fixed and round to the nearest integer
  if (class(copula)=="tCopula") {
    copula@df.fixed=TRUE
    copula@parameters[which(copula@param.names=="df")] = as.integer(round( copula@parameters[which(copula@param.names=="df")] ))
    attr(copula@parameters,"fixed")[which(copula@param.names=="df")] = TRUE
  }
  # set different estim.method
  method = tryCatch(fit@method, error=function(err) NA)
  if ( grepl("rho",method) ) { method="irho" }
  else if ( grepl("tau",method) ) { method="itau" }
  else { method="mpl" }
  tryCatch(gofCopula(copula,U,200,estim.method=method,simulation="mult",verbose=F)$p.value, error=function(err) NA)
}
# function to obtain loglik of fit
loglik.apply = function(fit){ tryCatch(fit@loglik, error=function(err) NA) }
# function to obtain lambda of fit
lambda.apply = function(fit){ tryCatch(lambda(fit@copula), error=function(err) rep(NA,2)) }
# function to obtain rho of fit
rho.apply = function(fit){ tryCatch(rho(fit@copula), error=function(err) NA) }

# check scatterplots
par(mfrow=c(2,2))
par(mar=c(4.1, 4.1, 1.1, 2.1))
plot(Uret[,c("GE","ITT")], cex=0.8, col=adjustcolor("black",alpha.f=0.3))
plot(Uret[,c("SPWR","FSLR")], cex=0.8, col=adjustcolor("black",alpha.f=0.3))
plot(Uret[,c("FSLR","TSLA")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot(Uret[,c("PCG","SPWR")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))

# GE:ITT
Fit.GE.ITT = fit.copulas( Uret[,c("GE","ITT")] )
( gof.GE.ITT = sapply(Fit.GE.ITT, gof.apply, Uret[,c("GE","ITT")]) )
sapply(Fit.GE.ITT, loglik.apply)
sapply(Fit.GE.ITT, lambda.apply)
sapply(Fit.GE.ITT, rho.apply)
# SPWR:FSLR
Fit.SPWR.FSLR = fit.copulas( Uret[,c("SPWR","FSLR")] )
( gof.SPWR.FSLR = sapply(Fit.SPWR.FSLR, gof.apply, Uret[,c("SPWR","FSLR")]) )
sapply(Fit.SPWR.FSLR, loglik.apply)
sapply(Fit.SPWR.FSLR, lambda.apply)
sapply(Fit.SPWR.FSLR, rho.apply)
# FSLR:TSLA
Fit.FSLR.TSLA = fit.copulas( Uret[,c("FSLR","TSLA")] )
( gof.FSLR.TSLA = sapply(Fit.FSLR.TSLA, gof.apply, Uret[,c("FSLR","TSLA")]) )
sapply(Fit.FSLR.TSLA, loglik.apply)
sapply(Fit.FSLR.TSLA, lambda.apply)
sapply(Fit.FSLR.TSLA, rho.apply)
# PCG:SPWR
Fit.PCG.SPWR = fit.copulas( Uret[,c("PCG","SPWR")] )
( gof.PCG.SPWR = sapply(Fit.PCG.SPWR, gof.apply, Uret[,c("PCG","SPWR")]) )
sapply(Fit.PCG.SPWR, loglik.apply)
sapply(Fit.PCG.SPWR, lambda.apply)
sapply(Fit.PCG.SPWR, rho.apply)



### (e) PCA #############################################################################

## (e.1) PCA and Construct New Index -----------------------------------
pca =  princomp(RET, cor = TRUE)
#reweight and plot
new_weight=c(rep(NA,10))
for(i in 1:10){
  new_weight[i]=pca$loadings[i,1]/sum(pca$loadings[,1])
}
sum(new_weight*RET[1,])
MRET=matrix(,nrow=2012,ncol=1)
for(i in 1:2012){
  MRET[i,1]=sum(new_weight*RET[i,])
}
MRET = xts(MRET, order.by=dates)
new_weight2=c(rep(NA,10))
for(i in 1:10){
  new_weight2[i]=pca$loadings[i,2]/sum(pca$loadings[,2])
}

MRET2=matrix(,nrow=2012,ncol=1)
for(i in 1:2012){
  MRET2[i,1]=sum(new_weight2*RET[i,])
}
MRET2 = xts(MRET2, order.by=dates)
Mdata = matrix(,nrow=2012,ncol=3)
for(i in 1:2012){
  Mdata[i,1]=Return.cumulative(MRET[1:i,], geometric = TRUE)[1]
  Mdata[i,2]=Return.cumulative(MSP500[1:i,], geometric = TRUE)[1]
  Mdata[i,3]=Return.cumulative(MRET2[1:i,], geometric = TRUE)[1]
}
Mdata=as.data.frame(Mdata)
timeline=1:2012
ggplot(Mdata, aes(x = timeline)) +
  geom_line(aes(y=Mdata$V1,color='First component'))+geom_line(aes(y=Mdata$V2,color='SP500'))+ylab('Cumulative Return')+theme_bw()+
  geom_line(aes(y=Mdata$V3,color='Second component'))+theme(legend.position = c(0.75,0.25))

# screeplot and biplot
fviz_eig(pca, barfill="grey", barcolor="black")
fviz_pca_biplot(pca, repel = TRUE, geom = "point", col.var = "black", col.ind = "grey" )

# construct new index and pseudo observations
NewIndex = RET %*% pca$loadings[,1:2]
NewIndex = cbind(NewIndex, SP500)
Uind = apply(NewIndex, 2, edf, adjust=1)
# check sample spearman's rho
cor(Uind, method="spearman")
# scatterplot of X and U
pairs2(NewIndex, cex=0.1, col=adjustcolor("black",alpha.f=0.3))
pairs2(Uind, cex=0.1, col=adjustcolor("black",alpha.f=0.3))


## (e.2) Analysis --------------------------------------------------------
# Comp.1: market factor

# Comp.2 verion 1: sentiment factor for Solar sector
RET[index( NewIndex[NewIndex[,"Comp.2"]< -0.1,] ), c("PCG","IDA","SPWR","FSLR")]
RET[index( NewIndex[NewIndex[,"Comp.2"]> 0.1,] ), c("PCG","IDA","SPWR","FSLR")]
# when Comp.2>>0, solar stocks plunged
# when Comp.2<<0, solar stocks surgerd
# 2013-02-14: https://www.fool.com/investing/general/2013/02/14/3-reasons-sunpower-is-leading-the-charge-in-solar.aspx
# 2013-04-09: https://www.marketwatch.com/story/first-solar-guidance-shocks-to-the-high-side-2013-04-09
# 2016-08-10: https://www.fool.com/investing/2016/08/10/why-shares-of-canadian-solar-inc-plunged-10-today.aspx

# Comp.2 version 2: risk factor related to the stock
apply(RET,2,sd)
apply(RET,2,skewness)
# For SPWR, FSLR, TSLA, sd > 0.03
# For SPWR, FSLR, TSLA, skewness >= 0.3
# Therefore combine volatility and skewness as one risk factor


## (e.3) Fit Copulas --------------------------------------------------------
# Copula for Comp.1 and sp500
Fit.c1 = fit.copulas( Uind[,c(1,3)] )
# check gof, loglik and lambda
( gof.c1 = sapply(Fit.c1, gof.apply, U=Uind[,c(1,3)]) )
sapply(Fit.c1, loglik.apply)
sapply(Fit.c1, lambda.apply)

# Copula for Comp.2 and sp500
Uc2rot = Uind[,c(2,3)]
plot(Uc2rot) # negatively correlated
Uc2rot[,1] = 1-Uc2rot[,1] # rotate s.t. able to fit Gumbel.etc
plot(Uc2rot) # now positively correlated
Fit.c2 = fit.copulas( Uc2rot )
# check gof, loglik and lambda
( gof.c2 = sapply(Fit.c2, gof.apply, U=Uc2rot) )
sapply(Fit.c2, loglik.apply)
sapply(Fit.c2, lambda.apply)


### (f) Marshall-Olkin Copula ---------------------------------------------------------------------------
test_data=as.matrix(Uret[,1:2])
eexp(test_data[,1], method = "mle", ci = FALSE, ci.type = "two-sided",
     ci.method = "exact", conf.level = 0.95)
eexp(test_data[,2], method = "mle", ci = FALSE, ci.type = "two-sided",
     ci.method = "exact", conf.level = 0.95)
#test_data2=matrix(c(0.2,0.7),nrow = 1,ncol=2)
alpha <- c(0.2, 0.8)#do no know how to find alpha
#wait for computing:http://home.iitk.ac.in/~kundu/paper141.pdf
#some explain:http://www.business.uwm.edu/gdrive/Soofi_E/995_Info/Dobrowolski-kumar%20MO.pdf
# Define the MO copula
C = function(u, alpha) {
  if(!is.matrix(u)) u <- rbind(u)
  stopifnot(0 <= alpha, alpha <= 1, length(alpha) == 2)
  pmin(u[,1]*u[,2]^(1-alpha[2]), u[,1]^(1-alpha[1])*u[,2])
}
# Define the singular component
S.C <- function(u, alpha) {
  stopifnot(0 <= u, u <= 1,
            0 <= alpha, alpha <= 1, length(alpha) == 2)
  tau <- alpha[1] * alpha[2] / (alpha[1] + alpha[2] - alpha[1] * alpha[2])
  tau * pmin(u[,1]^alpha[1], u[,2]^alpha[2])^(1/tau)
}

# Define the absolutely continuous component
A.C <- function(u, alpha) C(u, alpha) - S.C(u, alpha)
up <- test_data[,1]^(alpha[1]*(1/alpha[2]-1))
low <- test_data(1-alpha[1])*up

## Check the margins
plot(test_data[,1], ylab = expression(test_data[1]))
plot(test_data[,2], ylab = expression(test_data[2]))

## Evaluate the copula (and its singular and absolutely continuous components)
grid <- expand.grid("u[1]" = test_data, "u[2]" = test_data) # build a grid
val.C <- cbind(grid, "C(u[1],u[2])" = C(grid, alpha = alpha)) # append C
val.S <- cbind(grid, "S[C](u[1],u[2])" = S.C(grid, alpha = alpha)) # append S.C
val.A <- cbind(grid, "A[C](u[1],u[2])" = A.C(grid, alpha = alpha)) # append S.C
## Copula (wire frame and level curves)
#CARSH WARNING 30min not working
#wireframe2(val.C) # wire frame plot
#wireframe2(val.S) # singular component
#wireframe2(val.A) # absolutely continuous component