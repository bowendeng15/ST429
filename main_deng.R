### Package ###############################################################################
library(xts) # xts
library(car) # qqPlot
library(QRM) # ESnorm
library(copula)
library(factoextra) # plot for pca



### (a) Data ###############################################################################
setwd("/Users/Bowen.Deng/Desktop/LSE/ST429/project_025")
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
RET = log(RET+1) # RET_{t}:=S_{t+1}/S_{t}-1, X_{t}=log(S_{t+1}/S_{t})=log(RET_{t}+1)
RET = xts(RET, order.by=dates)
head(RET)
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
tickers = c("GE","ITT","PCG","IDA","SPWR","FSLR","TSLA")
Uret = apply(RET, 2, edf, adjust=1) # pobs(RET, ties.method = "max")
# scatterplot
pairs2(RET[,tickers], cex=0.3, col=adjustcolor("black",alpha.f=0.3))
pairs2(Uret[,tickers], cex=0.3, col=adjustcolor("black",alpha.f=0.3)) # plot(as.matrix(Uret[,4:5]), pch=1) # CVA


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
# function to obtain bic of fit
bic = function(fit){
  k = length(fit@estimate) # number of paramters
  if ( grepl("rho", fit@method) | grepl("tau", fit@method) ) {k=k-1}
  return( -2*fit@loglik+log(fit@nsample)*k )
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
# function to obtain bic of fit
bic.apply = function(fit){ tryCatch(bic(fit), error=function(err) NA) }
# function to obtain lambda of fit
lambda.apply = function(fit){ tryCatch(lambda(fit@copula), error=function(err) rep(NA,2)) }
# function to obtain rho of fit
rho.apply = function(fit){ tryCatch(rho(fit@copula), error=function(err) NA) }

# check scatterplots
par(mfrow=c(2,2))
par(mar=c(4.1, 4.1, 1.1, 2.1))
plot.default(RET[,c("GE","ITT")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot.default(RET[,c("SPWR","FSLR")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot.default(RET[,c("FSLR","TSLA")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot.default(RET[,c("PCG","SPWR")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))

plot(Uret[,c("GE","ITT")], cex=0.8, col=adjustcolor("black",alpha.f=0.3))
plot(Uret[,c("SPWR","FSLR")], cex=0.8, col=adjustcolor("black",alpha.f=0.3))
plot(Uret[,c("FSLR","TSLA")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot(Uret[,c("PCG","SPWR")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))

# GE:ITT
Fit.GE.ITT = fit.copulas( Uret[,c("GE","ITT")] )
( gof.GE.ITT = sapply(Fit.GE.ITT, gof.apply, Uret[,c("GE","ITT")]) )
sapply(Fit.GE.ITT, loglik.apply)
sapply(Fit.GE.ITT, bic.apply)
sapply(Fit.GE.ITT, lambda.apply)
sapply(Fit.GE.ITT, rho.apply)
# SPWR:FSLR
Fit.SPWR.FSLR = fit.copulas( Uret[,c("SPWR","FSLR")] )
( gof.SPWR.FSLR = sapply(Fit.SPWR.FSLR, gof.apply, Uret[,c("SPWR","FSLR")]) )
sapply(Fit.SPWR.FSLR, loglik.apply)
sapply(Fit.SPWR.FSLR, bic.apply)
sapply(Fit.SPWR.FSLR, lambda.apply)
sapply(Fit.SPWR.FSLR, rho.apply)
# FSLR:TSLA
Fit.FSLR.TSLA = fit.copulas( Uret[,c("FSLR","TSLA")] )
( gof.FSLR.TSLA = sapply(Fit.FSLR.TSLA, gof.apply, Uret[,c("FSLR","TSLA")]) )
sapply(Fit.FSLR.TSLA, loglik.apply)
sapply(Fit.FSLR.TSLA, bic.apply)
sapply(Fit.FSLR.TSLA, lambda.apply)
sapply(Fit.FSLR.TSLA, rho.apply)
# PCG:SPWR
Fit.PCG.SPWR = fit.copulas( Uret[,c("PCG","SPWR")] )
( gof.PCG.SPWR = sapply(Fit.PCG.SPWR, gof.apply, Uret[,c("PCG","SPWR")]) )
sapply(Fit.PCG.SPWR, loglik.apply)
sapply(Fit.PCG.SPWR, bic.apply)
sapply(Fit.PCG.SPWR, lambda.apply)
sapply(Fit.PCG.SPWR, rho.apply)


### (e) PCA #############################################################################

## (e.1) PCA and Construct New Index -----------------------------------
pca =  princomp(RET, cor = TRUE)

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
# plot of cumulative returns to compare trends
par(mfrow=c(1,1))
cumret = cumprod(exp(NewIndex))-1
cumret = scale(cumret, center=c(0,0,0), scale=c(1,1,1/8))
plot(cumret, legend.loc="topleft", col=c("red","blue","black"), main="Scaled Cumulative Returns")


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
# For SPWR, FSLR, TSLA, sd > 0.03
# Therefore take volatility as a risk factor


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
moCopula() # Marshall-Olkin copulas do not have densities (in "copula" package)

dcopula.MO = function (u, par, l13, l23, log = TRUE) {
  if ( length(par)!=1 ) stop("need 3 params")
  if ( ncol(u)!=2 ) stop("Marshall-Olkin copulas only available in the bivariate case")
  a = par/l13; b = par/l23
  res = apply(u, 1, function(row){
    u = row[1]; v = row[2]
    if (u^a > v^b) { res = log(1-a)-a*log(u) }
    else { res = log(1-b)-b*log(v) }
  })
  res = res + log(dexp(qexp(u,l13),l13)) + log(dexp(qexp(u,l23),l23))
  if (!log) { res = exp(res) }
  return(res)
}

fit.MO = function (Udata, initial = 0.5, ...) {
  l13 = 1; l23=1
  cat("l13", l13)
  cat("l23", l23)
  fn = function(par) { -sum( dcopula.MO(Udata, par, l13, l23, log = TRUE) ) }
  fit = optim(initial, fn=fn, lower=0, upper=min(l13,l23)
              , method="L-BFGS-B"
              , ...)
  # param = c(fit$par[3]/(fit$par[1]+fit$par[3]), fit$par[3]/(fit$par[2]+fit$par[3]))
  # res = list(param=param, loglik=-fit$value, fit=fit)
  res = fit
  return(res)
}

u = rCopula(1000, moCopula(c(0.3,0.1)))
sum(dcopula.MO(u, 0, 2,2))
fit.MO(Uret[,c("GE","ITT")])

library(MASS)
t = fitdistr(rexp(1000,10), "exponential")
t$estimate
