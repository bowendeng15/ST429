### Package ----------------------------------------------------------------------------
library(xts) # xts
library(QRM) # ESnorm
library(car) # qqPlot

#' @title function for computing total loss of a portfolio
#' @param RET vector or matrix, returns
#' @param weights vector, weights of each asset
#' @param value numeric, total portfolio value
#' @param linear logical, indicator of whether the RET is linear-return or log-return
loss.portfolio = function(RET, weights, value, linear=FALSE){
  if (is.matrix(RET)) { Th = nrow(RET) } else { Th = 1 } # time horizon
  weight.mat = matrix(weights, nrow=Th, ncol=length(weights), byrow=TRUE) # weights in a matrix
  Vw.mat = value * weight.mat # values of each assets in a matriix
  if (linear) { summand = RET*Vw.mat } else { summand = (exp(RET)-1)*Vw.mat }
  loss = -rowSums(summand)
  return(loss)
}



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
RET = log(RET+1)
RET = xts(RET, order.by=dates)
head(RET)
# SP500 return
tmp = read.csv("429sp.csv", stringsAsFactors=FALSE)
SP500 = xts(tmp$vwretd, order.by=dates)
SP500 = log(SP500+1)



### Log Return ----------------------------------------------------------------------------
# library(viridis) # viridis, magma, inferno, plasma
# COLORS = c("black","darkblue","royalblue","maroon3","pink4","darkred","darkorange2")
plot.zoo(cbind(SP500,RET[,1:5]), screens=1:10, col=c("royalblue",rep(1,5)), las=1, main="Log Returns", xlab="")
plot.zoo(cbind(SP500,RET[,6:10]), screens=1:10, col=c("royalblue",rep(1,5)), las=1, main="Log Returns", xlab="")



### Portfolio ---------------------------------------------------------------------------
pf.value = 10000
pf.weights = rep(0.1,10)
alpha = c(0.95,0.975,0.99,0.995,0.999,0.9999)
# compute loss of portfolio
loss = loss.portfolio(RET, pf.weights, pf.value)
length(loss)
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
CopulaRET = apply(RET, 2, edf, adjust=1)
# scatterplot
pairs(as.matrix(RET), pch=".")
pairs(as.matrix(CopulaRET), pch=".") # plot(as.matrix(CopulaRET[,4:5]), pch=1) # CVA

## fit -----------------------------------------------
# gauss copula
fit.gauss = fit.gausscopula(CopulaRET)
fit.gauss$ll.max
# gauss copula - method of moments: P^ = Spearman
( ll.gauss.spearman = sum(dcopula.gauss(CopulaRET, Spearman(CopulaRET), log=TRUE)) )
# gauss copula - method of moments: P^ = sin(pi * Kendall/2)
( ll.gauss.kendall = sum(dcopula.gauss(CopulaRET, sin(pi*Kendall(CopulaRET)/2), log=TRUE)) )
# t copula
fit.t = fit.tcopula(CopulaRET)
fit.t$ll.max # sum(dcopula.t(CopulaRET, df=fit.t$nu, Sigma=fit.t$P,log=TRUE))
# t copula - method of moments: P^ = Spearman
fit.t.spearman = fit.tcopula(CopulaRET, method="Spearman")
fit.t.spearman$ll.max
# t copula - method of moments: P^ = sin(pi * Kendall/2)
fit.t.kendall = fit.tcopula(CopulaRET, method="Kendall")
fit.t.kendall$ll.max

