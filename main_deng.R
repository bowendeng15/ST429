### Package ----------------------------------------------------------------------------
library(xts) # xts
library(car) # qqPlot
library(QRM) # ESnorm
library(copula)



### Data ----------------------------------------------------------------------------
setwd("/Users/Bowen.Deng/Desktop/LSE/ST429/project_025")
tmp = read.csv("429stocks.csv", stringsAsFactors=FALSE)
tmp$TICKER[tmp$TICKER=="SPWRA"] = "SPWR"
# extract tickers and dates
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
# SP500 return
tmp = read.csv("429sp.csv", stringsAsFactors=FALSE)
SP500 = xts(tmp$vwretd, order.by=dates)
SP500 = log(SP500+1)



### Log Return ----------------------------------------------------------------------------
par(mfrow=c(1,1))
plot.zoo(cbind(SP500,RET[,1:5]), screens=1:10, col=c("royalblue",rep(1,5)), las=1, main="Log Returns", xlab="",cex.main=2)
plot.zoo(cbind(SP500,RET[,6:10]), screens=1:10, col=c("royalblue",rep(1,5)), las=1, main="Log Returns", xlab="",cex.main=2)



### Portfolio ---------------------------------------------------------------------------

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



### Copulas ---------------------------------------------------------------------------
Uret = apply(RET, 2, edf, adjust=1) # pobs(RET, ties.method = "max")
# scatterplot
pairs2(RET, cex=0.1, col=adjustcolor("black",alpha.f=0.3))
pairs2(Uret, cex=0.1, col=adjustcolor("black",alpha.f=0.3)) # plot(as.matrix(Uret[,4:5]), pch=1) # CVA


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
  res = list(fit.gauss=fit.gauss, fit.gauss.spearman=fit.gauss.spearman
             , fit.t=fit.t, fit.t.tau=fit.t.tau
             , fit.frank=fit.frank, fit.frank.tau=fit.frank.tau
             , fit.gumbel=fit.gumbel, fit.gumbel.tau=fit.gumbel.tau
             , fit.clayton=fit.clayton, fit.clayton.tau=fit.clayton.tau)
}

#' @title fit multivariate copulas for each bivariate pair of possible combinations
#' @param U matrix, multivariate pseudo observations
#' @return list, have choose(10,2) slots of object returned by fit.copulas
pairwise.fit.copulas = function(U){
  d = ncol(U)
  comb = combn(1:d, 2)
  res = list()
  for (i in 1:ncol(comb)){
    cat("fitting pairs", comb[,i], ", ")
    res[[i]] = fit.copulas(U[,comb[,i]])
  }
  return(res)
}


# fit pairwisely - take 5 minutes
FIT = pairwise.fit.copulas(Uret)
load("FIT.RData") # load the result directly

tickerpairs = apply(combn(tickers,2), 2, paste, collapse=":")
copulanames = c("gauss","gauss.irho","t","t.tau","frank","frank.tau","gumbel","gumbel.tau","clayton","clayton.tau")


# loglikelihood
LOGLIK = c()
for (i in 1:length(FIT)){
  tickerpair = tickerpairs[i]
  tmp = sapply(FIT[[i]], function(fit){tryCatch(fit@loglik, error=function(err) NA)} )
  tmp = matrix(tmp, nrow=1)
  rownames(tmp) = tickerpair
  LOGLIK = rbind(LOGLIK, tmp)
}
colnames(LOGLIK) = copulanames
tmp = LOGLIK
tmp[is.na(tmp)]=0
tmp = apply(tmp, 1, function(row){copulanames[which.max(row)]})
# print
options(digits=3)
cbind(data.frame(LOGLIK), winner=tmp)


# tail dependence(lambda)
LAMBDA = c()
for (i in 1:length(FIT)){
  tickerpair = tickerpairs[i]
  tmp = sapply(FIT[[i]], function(fit){tryCatch(lambda(fit@copula), error=function(err) rep(NA,2))} )
  rownames(tmp) = paste(tickerpair,c("l","u"),sep="-")
  LAMBDA = rbind(LAMBDA, tmp)
}
colnames(LAMBDA) = copulanames
LAMBDA # print


# Spearman's rho
RHO = c()
for (i in 1:length(FIT)){
  tickerpair = tickerpairs[i]
  tmp = sapply(FIT[[i]], function(fit){tryCatch(rho(fit@copula), error=function(err) NA)} )
  tmp = matrix(tmp, nrow=1)
  rownames(tmp) = tickerpair
  # rownames(tmp) = tickerpair
  RHO = rbind(RHO, tmp)
}
colnames(RHO) = copulanames
RHO # print



### Marshall-Olkin Copula ---------------------------------------------------------------------------
moCopula() # Marshall-Olkin copulas do not have densities (in "copula" package)
