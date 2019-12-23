### Package ----------------------------------------------------------------------------
library("xts")


### Data ----------------------------------------------------------------------------
setwd("/Users/Bowen.Deng/Desktop/LSE/ST429/project_025")
tmp = read.csv("429stocks.csv", stringsAsFactors=FALSE)
tmp$TICKER[tmp$TICKER=="SPWRA"] = "SPWR"

data = c()
for (t in unique(tmp$TICKER)){
  data = cbind(data, tmp$RET[tmp$TICKER==t])
}
colnames(data) = unique(tmp$TICKER)

tmp = read.csv("429sp.csv", stringsAsFactors=FALSE)
data = cbind(data, SP=tmp$vwretd)

data = xts(data, order.by=as.Date(as.character( unique(tmp$caldt) ), format="%Y%m%d") )
head(data)

rm(tmp)



### Log Return ----------------------------------------------------------------------------
plot(log(data[,c(11,1:5)]+1),main="",yaxis.right=F,legend.loc="topleft") # SP500 and 1:5
plot(log(data[,c(11,6:10)]+1),main="",yaxis.right=F,legend.loc="topleft") # SP500 and 6:10

