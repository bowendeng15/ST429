### Package ----------------------------------------------------------------------------
library(xts) # xts
library(car) # qqPlot
library(QRM) # ESnorm
library(copula) 

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
RET = log(RET+1) # RET_{t}:=S_{t+1}/S_{t}-1, X_{t}=log(S_{t+1}/S_{t})=log(RET_{t}+1)
RET = xts(RET, order.by=dates)
head(RET)
# SP500 return
tmp = read.csv("429sp.csv", stringsAsFactors=FALSE)
SP500 = xts(tmp$vwretd, order.by=dates)
SP500 = log(SP500+1) 



### Log Return ----------------------------------------------------------------------------
plot.zoo(cbind(SP500,RET[,1:5]), screens=1:10, col=c("royalblue",rep(1,5)), las=1, main="Log Returns", xlab="",cex.main=2)
plot.zoo(cbind(SP500,RET[,6:10]), screens=1:10, col=c("royalblue",rep(1,5)), las=1, main="Log Returns", xlab="",cex.main=2)



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
Uret = apply(RET, 2, edf, adjust=1)
# scatterplot
pairs2(RET, cex=0.1, col=adjustcolor("black",alpha.f=0.3))
pairs2(Uret, cex=0.1, col=adjustcolor("black",alpha.f=0.3)) # plot(as.matrix(Uret[,4:5]), pch=1) # CVA



## "QRM" version -----------------------------------------------
# # gauss copula
# fit.gauss = fit.gausscopula(Uret)
# fit.gauss$ll.max
# # gauss copula - method of moments: P^ = 2*sin(pi * Spearman/6))
# sum(dcopula.gauss(Uret, 2*sin(pi*Spearman(Uret)/6), log=TRUE))
# # t copula
# fit.t = fit.tcopula(Uret)
# fit.t$ll.max # sum(dcopula.t(Uret, df=fit.t$nu, Sigma=fit.t$P,log=TRUE))
# # t copula - method of moments: P^ = Spearman
# fit.t.spearman = fit.tcopula(Uret, method="Spearman")
# fit.t.spearman$ll.max
# # t copula - method of moments: P^ = sin(pi * Kendall/2)
# fit.t.kendall = fit.tcopula(Uret, method="Kendall")
# fit.t.kendall$ll.max
# # AC clayton
# fit.clayton = fit.AC(Uret, "clayton")
# fit.clayton$ll.max


## "copula" version -----------------------------------------------
# example: find dependenc structure of GE and ITT
U = Uret[,c(8,9)] 
# gauss copula
fit.c.gauss = fitCopula(normalCopula(dim=2,dispstr="un"), data=U, method="mpl")
fit.c.gauss@loglik
# gauss copula - method of moments: P^ = 2*sin(pi * Spearman/6))
fit.c.gauss.spearman = fitCopula(normalCopula(dim=2,dispstr = "un"), data = U, method = "irho")
fit.c.gauss.spearman@loglik = sum( dCopula(U, normalCopula(param=fit.c.gauss.spearman@estimate,dim=2,dispstr="un"), log=T) )
fit.c.gauss.spearman@loglik
# t copula
fit.c.t = fitCopula(tCopula(dim=2,dispstr="un"), data=U, method="mpl")
fit.c.t@loglik
# t copula - method of moments: P^ = sin(pi * Kendall/2)
fit.c.t.tau = fitCopula(tCopula(dim=2,dispstr="un"), data=U, method="itau.mpl")
fit.c.t.tau@loglik
# AC gumbel copula
fit.c.gumbel = fitCopula(gumbelCopula(dim=2), data=U, method="mpl")
fit.c.gumbel@loglik
# AC gumbel copula - method of moments: tau
fit.c.gumbel.tau = fitCopula(gumbelCopula(dim=2), data=U, method="itau")
fit.c.gumbel.tau@loglik = sum( dCopula(U, gumbelCopula(param=fit.c.gumbel.tau@estimate,dim=2), log=T) )
fit.c.gumbel.tau@loglik
# AC clayton copula
fit.c.clayton = fitCopula(claytonCopula(dim=2), data=U, method="mpl")
fit.c.clayton@loglik
# AC clayton copula - method of moments: tau
fit.c.clayton.tau = fitCopula(claytonCopula(dim=2), data=U, method="itau")
fit.c.clayton.tau@loglik = sum( dCopula(U, claytonCopula(param=fit.c.clayton.tau@estimate,dim=2), log=T) )
fit.c.clayton.tau@loglik
# AC frank copula
fit.c.frank = fitCopula(frankCopula(dim=2), data=U, method="mpl")
fit.c.frank@loglik
# AC frank copula - method of moments: tau
fit.c.frank.tau = fitCopula(frankCopula(dim=2), data=U, method="itau")
fit.c.frank.tau@loglik = sum( dCopula(U, frankCopula(param=fit.c.frank.tau@estimate,dim=2), log=T) )
fit.c.frank.tau@loglik

res = list(fit.c.gauss, fit.c.gauss.spearman
           , fit.c.t, fit.c.t.tau
           , fit.c.frank, fit.c.frank.tau
           , fit.c.gumbel, fit.c.gumbel.tau
           , fit.c.clayton, fit.c.clayton.tau)
funcloglik = function(fit){
  fit@loglik
}
funclambda = function(fit){
  lambda(fit@copula)
}

# loglikelihood of all copula fit
tmp = matrix( sapply(res, funcloglik), nrow=1 )
colnames(tmp) = c("gauss","gauss.irho","t","t.tau","frank","frank.tau","gumbel","gumbel.tau","clayton","clayton.tau")
rownames(tmp) = "loglik"
tmp

# tail dependence(lambda) of all copula fit
tmp = sapply(res, funclambda)
colnames(tmp) = c("gauss","gauss.irho","t","t.tau","frank","frank.tau","gumbel","gumbel.tau","clayton","clayton.tau")
tmp


### Marshall-Olkin Copula ---------------------------------------------------------------------------
moCopula() # Marshall-Olkin copulas do not have densities (in "copula" package)
