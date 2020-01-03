### Package ###############################################################################
library(xts) # xts
library(car) # qqPlot
library(QRM) # ESnorm
library(copula)
library(factoextra) # plot for pca
library(tseries) # jarque.bera.test


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
RET = xts(RET, order.by=dates)
head(RET)
# obtain marketcap ( price * outstanding shares )
tmp = tmp[tmp$date=="20181231",]
marketcap = matrix(tmp$PRC*tmp$SHROUT, nrow=1)/1000000
colnames(marketcap) = tickers
marketcap
# obtain beta ( cov(r_i, r_market) / var(r_market) )
apply(RET,2,cov,y=SP500)/var(as.vector(SP500))

# SP500 return
tmp = read.csv("429sp.csv", stringsAsFactors=FALSE)
SP500 = xts(tmp$vwretd, order.by=dates)

RET = log(RET+1) # RET_{t}:=S_{t+1}/S_{t}-1, X_{t}=log(S_{t+1}/S_{t})=log(RET_{t}+1)
SP500 = log(SP500+1)

for (i in 1:ncol(RET)){
  print((exp(RET[,colnames(RET)[i]])-1)[which.min(RET[,colnames(RET)[i]])])
  print((exp(RET[,colnames(RET)[i]])-1)[which.max(RET[,colnames(RET)[i]])])
  cat("\n")
}
(exp(SP500)-1)[which.min(SP500)]
(exp(SP500)-1)[which.max(SP500)]

plot(SP500["2015-01-01/2017-01-01"])


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
par(mfrow=c(1,4), cex=.9)
qqPlot(loss, id=F)
qqPlot(as.vector(RET$PCG), ylab="PCG", id=F)
qqPlot(as.vector(RET$UNFI), ylab="UNFI", id=F)
qqPlot(as.vector(RET$FSLR), ylab="FSLR", id=F)
jarque.bera.test(loss)
skewness(loss)
kurtosis(loss)
MardiaTest(as.matrix(RET))
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
legend("topleft", legend=c("HS","normal"), col=1:2, lty=1, bty="n")
hist(loss,nclass=100, prob=TRUE, xlab="Portfolio Loss", main=paste("ES at", paste(names(VaR.hs),collapse=", ")))
abline(v=ES.hs,col=1,lty=2)
abline(v=ES.normal,col=2,lty=5)
legend("topleft", legend=c("HS","normal"), col=1:2, lty=1, bty="n")



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
bic.apply = function(fit){ tryCatch(BIC(fit), error=function(err) NA) }
# function to obtain lambda of fit
lambda.apply = function(fit){ tryCatch(lambda(fit@copula), error=function(err) rep(NA,2)) }
# function to obtain rho of fit
rho.apply = function(fit){ tryCatch(rho(fit@copula), error=function(err) NA) }
# function to obtain tau of fit
tau.apply = function(fit){ tryCatch(tau(fit@copula), error=function(err) NA) }

# check scatterplots
par(mfrow=c(3,4))
par(mar=c(4.1, 4.1, 1.1, 2.1))
par(cex=0.7)
plot.default(RET[,c("GE","ITT")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot.default(RET[,c("SPWR","FSLR")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot.default(RET[,c("FSLR","TSLA")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot.default(RET[,c("PCG","SPWR")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))

plot(Uret[,c("GE","ITT")], cex=0.8, col=adjustcolor("black",alpha.f=0.3))
plot(Uret[,c("SPWR","FSLR")], cex=0.8, col=adjustcolor("black",alpha.f=0.3))
plot(Uret[,c("FSLR","TSLA")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot(Uret[,c("PCG","SPWR")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
set.seed(123)
plot(Fit.GE.ITT$gauss.rho@copula, n=2012, cex=0.8, col=adjustcolor("black",alpha.f=0.3),main="", xlab="gauss.rho U1",ylab="gauss.rho U2")
plot(Fit.SPWR.FSLR$t@copula, n=2012, cex=0.8, col=adjustcolor("black",alpha.f=0.3),main="", xlab="t U1",ylab="t U2")
plot(Fit.FSLR.TSLA$frank@copula, n=2012, cex=0.8, col=adjustcolor("black",alpha.f=0.5),main="", xlab="frank U1",ylab="frank U2")
plot(Fit.GE.ITT$gauss@copula, n=2012, cex=0.8, col=adjustcolor("black",alpha.f=0.5),main="", xlab="gauss U1",ylab="gauss U2")



options(digits = 4)
# GE:ITT
Fit.GE.ITT = fit.copulas( Uret[,c("GE","ITT")] )
( gof.GE.ITT = sapply(Fit.GE.ITT, gof.apply, Uret[,c("GE","ITT")]) )
sapply(Fit.GE.ITT, bic.apply)
sapply(Fit.GE.ITT, lambda.apply)
sapply(Fit.GE.ITT, rho.apply)
# SPWR:FSLR
Fit.SPWR.FSLR = fit.copulas( Uret[,c("SPWR","FSLR")] )
( gof.SPWR.FSLR = sapply(Fit.SPWR.FSLR, gof.apply, Uret[,c("SPWR","FSLR")]) )
sapply(Fit.SPWR.FSLR, bic.apply)
sapply(Fit.SPWR.FSLR, lambda.apply)
sapply(Fit.SPWR.FSLR, rho.apply)
# FSLR:TSLA
Fit.FSLR.TSLA = fit.copulas( Uret[,c("FSLR","TSLA")] )
( gof.FSLR.TSLA = sapply(Fit.FSLR.TSLA, gof.apply, Uret[,c("FSLR","TSLA")]) )
sapply(Fit.FSLR.TSLA, bic.apply)
sapply(Fit.FSLR.TSLA, lambda.apply)
sapply(Fit.FSLR.TSLA, rho.apply)
# PCG:SPWR
Fit.PCG.SPWR = fit.copulas( Uret[,c("PCG","SPWR")] )
( gof.PCG.SPWR = sapply(Fit.PCG.SPWR, gof.apply, Uret[,c("PCG","SPWR")]) )
sapply(Fit.PCG.SPWR, bic.apply)
sapply(Fit.PCG.SPWR, lambda.apply)
sapply(Fit.PCG.SPWR, rho.apply)


### (e) PCA #############################################################################

## (e.1) PCA and Construct New Index -----------------------------------
pca =  princomp(RET, cor = TRUE)
# percentage of variance explained
pca$sdev^2/sum(pca$sdev^2)

# screeplot and biplot
fviz_eig(pca, barfill="grey", barcolor="black")
fviz_pca_biplot(pca, repel = TRUE, geom = "point", col.var = "black", col.ind = "grey" )

# construct new index and pseudo observations
NewIndex = RET %*% pca$loadings[,1:2]
NewIndex = cbind(NewIndex, SP500)
Uind = apply(NewIndex, 2, edf, adjust=1)
# check sample spearman's rho
cor(Uind, method="spearman")
cor(Uind, method="kendall")
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
# check gof, bic
( gof.c1 = sapply(Fit.c1, gof.apply, U=Uind[,c(1,3)]) )
sapply(Fit.c1, bic.apply)

# Copula for Comp.2 and sp500
Uc2rot = Uind[,c(2,3)]
plot(Uc2rot) # negatively correlated
Uc2rot[,1] = 1-Uc2rot[,1] # rotate s.t. able to fit Gumbel.etc
plot(Uc2rot) # now positively correlated
Fit.c2 = fit.copulas( Uc2rot )
# check gof, bic
( gof.c2 = sapply(Fit.c2, gof.apply, U=Uc2rot) )
sapply(Fit.c2, bic.apply)

# check dependence measures
sapply(Fit.c1, lambda.apply)
sapply(Fit.c2, lambda.apply)
cor(Uind, method="spearman")
sapply(Fit.c1, rho.apply)
sapply(Fit.c2, rho.apply)
cor(Uind, method="kendall")
sapply(Fit.c1, tau.apply)
sapply(Fit.c2, tau.apply)

par(mfrow=c(2,3))
par(cex=0.9)
plot(as.matrix(NewIndex[,c(1,3)]), col=adjustcolor("black",alpha.f=0.2))
plot(Uind[,c(1,3)], col=adjustcolor("black",alpha.f=0.2))
set.seed(2)
plot(rCopula(Fit.c1$gauss.rho@copula, n=2012), col=adjustcolor("black",alpha.f=0.2), xlab="gauss.rho U1", ylab="gauss.rho U2")
set.seed(2)
tmp = rCopula(Fit.c2$frank.tau@copula, n=2012)
tmp[,1] = 1-tmp[,1]
plot(as.matrix(NewIndex[,c(2,3)]), col=adjustcolor("black",alpha.f=0.5))
plot(Uind[,c(2,3)], col=adjustcolor("black",alpha.f=0.5))
plot(tmp, col=adjustcolor("black",alpha.f=0.5), xlab="frank.tau U1", ylab="frank.tau U2")



### (f) Marshall-Olkin Copula ---------------------------------------------------------------------------

par(mfrow=c(1,3))
par(cex=.9)
plot(moCopula(c(.0, .0)), n = 10000, xaxs="i", yaxs="i", pch = 16, col = adjustcolor("black", 0.2))
plot(moCopula(c(.9, .2)), n = 10000, xaxs="i", yaxs="i", pch = 16, col = adjustcolor("black", 0.2))
plot(moCopula(c(1, 1)), n = 10000, xaxs="i", yaxs="i", pch = 16, col = adjustcolor("black", 0.2))

par(mfrow=c(1,4))
par(cex=.9)
tmp = mvrnorm(1000, mu=c(0,0), Sigma=matrix(c(1,0.6,0.6,1),nrow=2))
plot(tmp, xlab="X", ylab="Y")
plot(exp(tmp), xlab="exp(X)", ylab="exp(Y)")
plot(pobs(tmp,ties.method = "max"), xlab="F(X)", ylab="F(Y)")
plot(pobs(exp(tmp),ties.method = "max"), xlab="F(exp(X))", ylab="F(exp(Y))")


dcopula.MO = function (u, par, log = TRUE) {
  if ( length(par)!=2 ) stop("need 2 params")
  if ( ncol(u)!=2 ) stop("Marshall-Olkin copulas only available in the bivariate case")
  a = par[1]; b = par[2]
  res = apply(u, 1, function(row){
    u = row[1]; v = row[2]
    if (u^a > v^b) { log(1-a)-a*log(u) } else { log(1-b)-b*log(v) }
  })
  if (!log) { res = exp(res) }
  return(res)
}

fit.MO = function (data, initial = c(.5,.5), ...) {
  Udata = pobs(data, ties.method = "max")
  fn = function(par) { -sum( dcopula.MO(Udata, par, log = TRUE) ) }
  fit = optim(initial, fn=fn, lower=c(0,0), upper=c(1,1), method="L-BFGS-B",...)
  return(fit)
}

t = rCopula(1000, moCopula(0.4, 0.2))

fit.MO(exp(as.matrix(RET[,c("GE","ITT")])))
fit.MO(exp(as.matrix(RET[,c("SPWR","FSLR")])))
fit.MO(exp(as.matrix(RET[,c("FSLR","TSLA")])))
fit.MO(exp(as.matrix(RET[,c("PCG","SPWR")])))
fit.MO(exp(as.matrix(NewIndex[,c("Comp.1","SP500")])))
fit.MO(exp(as.matrix(NewIndex[,c("Comp.2","SP500")])))


# dcopula.MO.f = function (Udata, par, l13, l23, log = TRUE) {
#   if ( length(par)!=1 ) stop("need 1 param")
#   if ( ncol(Udata)!=2 ) stop("Marshall-Olkin copulas only available in the bivariate case")
#   a = par/l13; b = par/l23
#   res = apply(Udata, 1, function(row){
#     u = row[1]; v = row[2]
#     if (u^a > v^b) { c = log(1-a)-a*log(u) }
#     else { c = log(1-b)-b*log(v) }
#     return(c)
#   })
#   fx = dexp(qexp(Udata[,1],l13), l13, log=T)
#   fy = dexp(qexp(Udata[,2],l23), l23, log=T)
#   res = res + fx + fy
#   if (!log) { res = exp(res) }
#   return(res)
# }
# fit.MO = function (data, ...) {
#   l13 = fitdistr(data[,1], "exponential")$estimate
#   l23 = fitdistr(data[,2], "exponential")$estimate
#   Udata = pobs(data)
#   fn = function(par) { -sum( dcopula.MO.f(Udata, par, l13, l23, log = TRUE) ) }
#   fit = optim(par=min(l13,l23)/4, fn=fn, lower=0, upper=min(l13,l23)
#               , method="L-BFGS-B",...)
#   # param = c(fit$par[3]/(fit$par[1]+fit$par[3]), fit$par[3]/(fit$par[2]+fit$par[3]))
#   # res = list(param=param, loglik=-fit$value, fit=fit)
#   res = fit
#   return(res)
# }




# u = rCopula(1000, moCopula(c(0.9,0.3)))
# par(mfrow=c(1,1))
# plot(u)
# fit.MO(u)
#
# res = apply(u, 1, function(row){
#   u = row[1]; v = row[2]
#   if (u^1e-12 > v^1e-12) { c = 1 }
#   else { c = 0 }
#   return(c)
# })
#
# sum( dcopula.MO(u, par=c(1e-12,1e-12), log = TRUE) )
# sum( dcopula.MO(u, par=c(.5,.03), log = TRUE) )
# fit.MO(Uret[,c("PCG","SPWR")])
# fit.MO(Uind[,c("Comp.1","SP500")])
#
# sum(pobs(pobs(u))==1)
#
# library(MASS)
# t = fitdistr(NewIndex[,c(1,3)], "exponential")$estimate
# t$estimate
#
# rCopula()
# a = seq(0,1,length.out = 50)
# f = (1-a)/exp(log(0.2)*a)
# plot(a,f)
#
# u=v=0.1
# a = 0.3; b=0.1
# ch1 = 0.9^1.7
# ch2 = 0.9^1.9
#
# Y=1:10
# ( U = pobs(Y) )
# pobs(U)
#
# ch1-0.9-0.9+1
# ch2-0.9-0.9+1
#
# 0.1^1.7
# 0.1^1.9
#
# fit.gausscopula( rcopula.gauss(1000, cbind(c(1,0.4),c(0.4,1))) )
#
# cop <- gumbelCopula(2)
# contourplot2(cop, dCopula)
