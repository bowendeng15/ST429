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
tmp$TICKER[tmp$TICKER=="SPWRA"] = "SPWR" # same company
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
# obtain marketcap ( price * outstanding shares )
tmp = tmp[tmp$date=="20181231",]
matrix(tmp$PRC*tmp$SHROUT, nrow=1, dimnames=list("marketcap", tickers))/1000000

# SP500 return
tmp = read.csv("429sp.csv", stringsAsFactors=FALSE)
SP500 = xts(tmp$vwretd, order.by=dates)
# obtain beta ( cov(r_i, r_market) / var(r_market) )
apply(RET,2,cov,y=SP500)/var(as.vector(SP500))

# log returns # RET_{t}:=S_{t+1}/S_{t}-1, X_{t}=log(S_{t+1}/S_{t})=log(RET_{t}+1)
RET = log(RET+1)
SP500 = log(SP500+1)



### (b) Log Returns ###############################################################################
par(mfrow=c(1,1))
ylab = c("SP500", tickers[1:5], "SP500", tickers[6:10])
plot.zoo(cbind(SP500,RET[,1:5], SP500, RET[,6:10]), screens=1:12, col=c("royalblue",rep(1,5))
         , las=1, main="", xlab="", ylab=ylab, cex.main=1
         , panel=function(x,y,...){abline(v=as.Date(c("2011-08-08","2015-08-21","2018-02-05")), col="salmon",lty=2);lines(x, y,...)})



### (c) Portfolio ###############################################################################

#' @title function for computing loss, normal VaR and ES, HS VaR and ES of a portfolio
#' @param RET vector or matrix, returns
#' @param weights vector, weights of each asset
#' @param value numeric, total portfolio value
#' @param alpha vector, significance level of VaR and ES
#' @param linear logical, indicator of whether the RET is linear-return or log-return
#' @return a list of loss, normal VaR and ES, HS VaR and ES
summary.portfolio = function(RET, pf.weights, pf.value, alpha, linear=FALSE){
  # loss
  loss = loss.portfolio(RET, pf.weights, pf.value, linear=linear)
  # normal VaR and ES
  mu.hat = colMeans(RET)
  sigma.hat = var(RET)
  meanloss = -sum(pf.weights*mu.hat) * pf.value
  varloss = pf.value^2 * as.numeric(t(pf.weights) %*% sigma.hat %*% pf.weights)
  VaR.normal = meanloss + sqrt(varloss) * qnorm(alpha)
  ES.normal = meanloss + sqrt(varloss) * dnorm(qnorm(alpha))/(1-alpha)
  # HS VaR and ES
  VaR.hs = quantile(loss,alpha)
  ES.hs = apply(as.array(VaR.hs), 1, function(VaR){ mean(loss[loss > VaR]) } )
  list(loss=loss, VaR.normal=VaR.normal, ES.normal=ES.normal, VaR.hs=VaR.hs, ES.hs=ES.hs)
}
#' function for computing total loss of a portfolio
loss.portfolio = function(RET, weights, value, linear=FALSE){
  if (linear) { RET = RET } else { RET = (exp(RET)-1) }
  Profit = RET * value * weights
  -rowSums(Profit)
}
#' function for plotting the histogram of loss with VaR or ES
plot.portfolio = function(portfolio, type=c("VaR","ES"),...){
  hist(portfolio$loss, nclass=50, prob=TRUE, xlab="Loss", ...)
  measure.hs = switch(type, VaR=portfolio$VaR.hs, ES=portfolio$ES.hs)
  measure.normal = switch(type, VaR=portfolio$VaR.normal, ES=portfolio$ES.normal)
  abline(v=measure.hs, col=1, lty=2)
  abline(v=measure.normal, col=2, lty=5)
  legend("topleft", legend=c("HS","normal"), col=1:2, lty=1, bty="n")
}

# compute loss, normal VaR and ES, HS VaR and ES
pf1 = summary.portfolio(RET, pf.weights=rep(0.1,10), pf.value=10000, alpha=alpha)
pf2 = summary.portfolio(RET[,1:5], pf.weights=rep(0.2,5), pf.value=5000, alpha=alpha)
pf3 = summary.portfolio(RET[,1:3], pf.weights=rep(1/3,3), pf.value=3000, alpha=alpha)
# check normality
par(mfrow=c(1,3), cex=.9)
qqPlot(pf1$loss, id=F, ylab="loss", main="Portfolio 1")
qqPlot(pf2$loss, id=F, ylab="loss", main="Portfolio 2")
qqPlot(pf3$loss, id=F, ylab="loss", main="Portfolio 3")
jarque.bera.test(pf1$loss)$p.value
jarque.bera.test(pf2$loss)$p.value
jarque.bera.test(pf3$loss)$p.value
MardiaTest(as.matrix(RET))
MardiaTest(as.matrix(RET[,1:5]))
MardiaTest(as.matrix(RET[,1:3]))
# histogram and comparing
par(mfrow=c(2,3), cex=.9)
par(mar = c(3.1, 1.1, 3.1, 0.1))
plot.portfolio(pf1, type="VaR", main="VaR of Portfolio 1", yaxt="n")
plot.portfolio(pf2, type="VaR", main="VaR of Portfolio 2", yaxt="n")
plot.portfolio(pf3, type="VaR", main="VaR of Portfolio 3", yaxt="n")
plot.portfolio(pf1, type="ES", main="ES of Portfolio 1", yaxt="n")
plot.portfolio(pf2, type="ES", main="ES of Portfolio 2", yaxt="n")
plot.portfolio(pf3, type="ES", main="ES of Portfolio 3", yaxt="n")



### (d) Copulas ###############################################################################

## (d.1) Pseudo-observations ------------------------------------
Uret = pobs(RET, ties.method = "max") # apply(RET, 2, edf, adjust=1)
# scatterplot
tickers = c("GE","ITT","PCG","SPWR","FSLR","TSLA")
pairs2(RET[,tickers], cex=0.3, col=adjustcolor("black",alpha.f=0.3))
pairs2(Uret[,tickers], cex=0.3, col=adjustcolor("black",alpha.f=0.3))


## (d.2) Sample Correlations ------------------------------------
cor(Uret[,tickers], method="spearman")
cor(Uret[,tickers], method="kendall")


## (d.3) Fit Copula --------------------------------------------

#' @title fit bivarite copulas including gauss,t,gumbel,clayton and frank copulas
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
# try to fitCopula or return NA
tryFitCopula = function(...){
  tryCatch(fitCopula(...), error=function(err) NA)
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
  # 200 bootstrap are used
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


# Fitting Copulas
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


# check scatterplots
par(mfrow=c(3,4), cex=0.7)
par(mar=c(4.1, 4.1, 1.1, 2.1))
# log returns
plot.default(RET[,c("GE","ITT")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot.default(RET[,c("SPWR","FSLR")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot.default(RET[,c("FSLR","TSLA")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot.default(RET[,c("PCG","SPWR")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
# pseudo observations
plot.default(Uret[,c("GE","ITT")], cex=0.8, col=adjustcolor("black",alpha.f=0.3))
plot.default(Uret[,c("SPWR","FSLR")], cex=0.8, col=adjustcolor("black",alpha.f=0.3))
plot.default(Uret[,c("FSLR","TSLA")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
plot.default(Uret[,c("PCG","SPWR")], cex=0.8, col=adjustcolor("black",alpha.f=0.5))
# simulations
set.seed(123)
plot(Fit.GE.ITT$gauss.rho@copula, n=2012
     , cex=0.8, col=adjustcolor("black",alpha.f=0.3),main="", xlab="gauss.rho U1",ylab="gauss.rho U2")
plot(Fit.SPWR.FSLR$t@copula, n=2012
     , cex=0.8, col=adjustcolor("black",alpha.f=0.3),main="", xlab="t U1",ylab="t U2")
plot(Fit.FSLR.TSLA$frank@copula, n=2012
     , cex=0.8, col=adjustcolor("black",alpha.f=0.5),main="", xlab="frank U1",ylab="frank U2")
plot(Fit.GE.ITT$gauss@copula, n=2012
     , cex=0.8, col=adjustcolor("black",alpha.f=0.5),main="", xlab="gauss U1",ylab="gauss U2")



### (e) PCA #############################################################################

## (e.1) PCA and Construct New Index -----------------------------------
pca = princomp(RET, cor = TRUE)
# percentage of variance explained
pca$sdev^2/sum(pca$sdev^2)
# biplot
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


## (e.2) Fit Copulas --------------------------------------------------------
# Copula for Comp.1 and sp500
Fit.c1 = fit.copulas( Uind[,c(1,3)] )
# check gof, bic, dependence measures
( gof.c1 = sapply(Fit.c1, gof.apply, U=Uind[,c(1,3)]) )
sapply(Fit.c1, bic.apply)
sapply(Fit.c1, lambda.apply)
sapply(Fit.c1, rho.apply)
sapply(Fit.c1, tau.apply)

# Copula for Comp.2 and sp500
Uc2rot = Uind[,c(2,3)]
plot(Uc2rot) # negatively correlated
Uc2rot[,1] = 1-Uc2rot[,1] # rotate s.t. able to fit Gumbel.etc
plot(Uc2rot) # now positively correlated
Fit.c2 = fit.copulas( Uc2rot )
# check gof, bic, dependence measures
( gof.c2 = sapply(Fit.c2, gof.apply, U=Uc2rot) )
sapply(Fit.c2, bic.apply)
sapply(Fit.c2, lambda.apply)
sapply(Fit.c2, rho.apply)
sapply(Fit.c2, tau.apply)

# plot log returns, pseudo observations and simulations
par(mfrow=c(2,3), cex=0.9)
set.seed(2)
plot(as.matrix(NewIndex[,c(1,3)]), col=adjustcolor("black",alpha.f=0.2))
plot(Uind[,c(1,3)], col=adjustcolor("black",alpha.f=0.2))
plot(rCopula(Fit.c1$gauss.rho@copula, n=2012), col=adjustcolor("black",alpha.f=0.2), xlab="gauss.rho U1", ylab="gauss.rho U2")
set.seed(2)
tmp = rCopula(Fit.c2$frank.tau@copula, n=2012)
tmp[,1] = 1-tmp[,1]
plot(as.matrix(NewIndex[,c(2,3)]), col=adjustcolor("black",alpha.f=0.5))
plot(Uind[,c(2,3)], col=adjustcolor("black",alpha.f=0.5))
plot(tmp, col=adjustcolor("black",alpha.f=0.5), xlab="frank.tau U1", ylab="frank.tau U2")



### (f) Marshall-Olkin Copula ---------------------------------------------------------------------------

# scatter plot of MO copula simulations
par(mfrow=c(1,3))
plot(moCopula(c(.0, .0)), n = 10000, xaxs="i", yaxs="i", pch = 16, col = adjustcolor("black", 0.2))
plot(moCopula(c(.05, .05)), n = 10000, xaxs="i", yaxs="i", pch = 16, col = adjustcolor("black", 0.2))
plot(moCopula(c(1, 1)), n = 10000, xaxs="i", yaxs="i", pch = 16, col = adjustcolor("black", 0.2))

# scatter plot of exponential transformation
par(mfrow=c(1,4))
tmp = mvrnorm(1000, mu=c(0,0), Sigma=matrix(c(1,0.6,0.6,1),nrow=2))
plot(tmp, xlab="X", ylab="Y")
plot(exp(tmp), xlab="exp(X)", ylab="exp(Y)")
plot(pobs(tmp,ties.method = "max"), xlab="F(X)", ylab="F(Y)")
plot(pobs(exp(tmp),ties.method = "max"), xlab="F(exp(X))", ylab="F(exp(Y))")

#' @title function for computing copula density for MO
#' @param u matrix, pseudo observations
#' @param par vector of length 2, alpha and beta
#' @param log logical, whether to compute log-likelihood
#' @return a vector of copula density
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
#' @title function for fitting MO copula using pseudo MLE
#' @param data matrix, original data
#' @param initial vector of length 2, intial values for alpha and beta
#' @return fit information
fit.MO = function (data, initial = c(.5,.5), ...) {
  Udata = pobs(data, ties.method = "max")
  fn = function(par) { -sum( dcopula.MO(Udata, par, log = TRUE) ) }
  fit = optim(initial, fn=fn, upper=c(1-1e-12,1-1e-12), method="L-BFGS-B",...) #, lower=c(0,0), upper=c(1,1)
  return(fit)
}

fit.MO( exp(as.matrix(RET[,c("GE","ITT")])) )
fit.MO( exp(as.matrix(RET[,c("SPWR","FSLR")])) )
fit.MO( exp(as.matrix(RET[,c("FSLR","TSLA")])) )
fit.MO( exp(as.matrix(RET[,c("PCG","SPWR")])) )
fit.MO( exp(as.matrix(NewIndex[,c("Comp.1","SP500")])) )
fit.MO( exp(as.matrix(NewIndex[,c("Comp.2","SP500")])) )
