### Package ###############################################################################
library(xts) # xts
library(car) # qqPlot
library(QRM) # ESnorm
library(copula)
library(factoextra) # plot for pca


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

### (d) Copulas ###############################################################################
## (d.1) Pseudo-observations ------------------------------------
RET = RET[,c("PCG","IDA","SPWR","FSLR","TSLA")]
Uret = apply(RET, 2, edf, adjust=1) # pobs(RET, ties.method = "max")

### (f) Marshall-Olkin Copula ---------------------------------------------------------------------------
test_data=as.matrix(Uret[,1:2])
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

