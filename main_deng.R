### Package ----------------------------------------------------------------------------
library(xts) # xts
library(QRM) # ESnorm
library(car) # qqPlot

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

loss = lo.fn(RET, pf.weights, pf.value)

mu.hat = colMeans(RET)
sigma.hat = var(RET) #*(nrow(RET)-1)/nrow(RET)

meanloss = -sum(pf.weights*mu.hat) * pf.value
varloss = pf.value^2 * as.numeric(t(pf.weights) %*% sigma.hat %*% pf.weights)

VaR.normal = meanloss + sqrt(varloss) * qnorm(0.95)
ES.normal = meanloss + sqrt(varloss) * dnorm(qnorm(0.95))/(1-0.95) # ESnorm(0.95, mu=meanloss, sd=sqrt(varloss))

qqPlot(loss)

VaR.hs = quantile(loss,0.95)
ES.hs <- mean(loss[loss>VaR.hs])

hist(loss,nclass=100, prob=TRUE, xlab="Loss Distribution", main="Historical simulation")
abline(v=c(VaR.normal,ES.normal),col=1,lty=2)
abline(v=c(VaR.hs,ES.hs),col=2,lty=5)
legend("topleft", legend = c("normal","HS"), col=1:2, lty=1, bty="n") 
