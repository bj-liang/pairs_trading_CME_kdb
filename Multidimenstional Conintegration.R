# \c 76 178
# h: hopen `:orf474:6000
# h"instinfo"
closeAllConnections()
rm(list=ls())
library(rkdb)

# set your directory to determine where things will be outputed
# no input files necesary
setwd("C:\\Users\\benny\\OneDrive\\Desktop\\github\\Pairs PCA CME")

midpoint_price <- function(
  db, 
  sym, 
  date,
  tmin, # where you start
  tmax, # this is the shift of when you end
  dt = 60) #increments in seconds
{
  date = as.Date(date, format="%Y-%m-%d")
  startDate= as.character(date, format = '%Y.%m.%d')
  q <- paste0('{[d;s;dt]
    t:([] time:`timespan$1000000000*((`int$-', tmin*60*60, ')+dt*til 1+`int$', tmax*60*60, '%dt));
    m: (select mid:0.5*(last bid+last ask) by time, sym from quote where date=d, sym=s);
    select from aj[`time; t; m]
                }[', startDate, ';`', sym, ';', dt,']')
  
  execute(db,q)
}

pred_price <- function(corMatrix, mu, sigma, realPrice, 
                       turnoff = append(c(1,1),numeric(nrow(corMatrix)-2)))
{
  V = svd(corMatrix)$u
  newCoor = t(V)%*%((realPrice-mu)/sigma)
  settoZero = turnoff * newCoor
  diag(sigma)%*%V%*%settoZero+mu
}

#download all the data beforehand
get_dataset <- function(symbol_list, db, date)
{
dataset = data.frame()
size = c()
for (s in symbol_list) {
  tmp <- midpoint_price(db, s, date, 7, 23)
  dataset = rbind(dataset, tmp)
  size = append(size, nrow(tmp))
}
if (anyNA(dataset)) {print('missing values in query')}
else if (any(size != median(size))) {print("size error in days")}
else {dataset}
}

EuroDollars = c('GEH8',
         'GEM8',
         'GEU8',
         'GEZ8',
         'GEH9',
         'GEM9',
         'GEU9',
         'GEZ9')
# what data do you want? and how are you going to source it
h1 <- open_connection('orf474.princeton.edu', 6007) # Interest rate futures
h2 <- open_connection('orf474.princeton.edu', 7007) # Non-interest rate futures
database = h1
d = get_dataset(EuroDollars, database, '2017-10-10')
#-------------------------------------------------------------------------------
#this is where the sausage is made 
wholeDaysSignal <- function(dataset, productorder,
                            turnoff = append(c(1,1),numeric(length(productorder)-2)),
                            m = 3600)
{
  eight = length(productorder) # this occurs a lot
  x = matrix(nrow = eight)
  meanA = as.matrix(numeric(eight))
  meanB = 0
  sigma = c()
  covA = 0
  predPrice_stor = NA
  open = 840
  close = 1320
  
  for (t in seq(nrow(dataset)/eight)) {
    ticker = 0 
    for (ge in productorder) {
      ticker = ticker + 1 
      x[ticker,1] = dataset[dataset$sym == ge,]$mid[t]
      meanA[ticker,1] = exp(-1/m)*meanA[ticker,1] + x[ticker,1]
    }
    meanB = exp(-1/m)*meanB + 1 
    
    if (t == 1){realPrice_stor = t(x)}else{realPrice_stor = rbind(realPrice_stor, t(x))}
    
    # stuff to do before feed into function
    meanie = meanA/meanB
    covA = exp(-1/m)*covA + (x%*%t(x))
    covM = covA/meanB - (meanie %*% t(meanie))
    if (t > open && t < close){ # main trading hours
      siggy = sqrt(diag(covM))
      corM = matrix(nrow = eight, ncol = eight)
      for (r in seq(eight)) {
        for (c in seq(eight)) {
          corM[r,c] = covM[r,c]/(siggy[r]*siggy[c])
        }
      }
      x_p = pred_price(corM, meanie, siggy, x, turnoff)
      # once again gonna want to store this data to test
      # remember this data only has the pred price of trading hours
      if (length(predPrice_stor) == 1){predPrice_stor = t(x_p)}else{predPrice_stor = rbind(predPrice_stor, t(x_p))}
    }
  }
  list(rp = realPrice_stor, pp = predPrice_stor, pointer = covM)
}

open = 840
close = 1320
se_ratio = data.frame()
correlationTable = data.frame()
slope_table = data.frame()
error_table = data.frame()
eight = length(EuroDollars)
tf = numeric(eight)
 for (i in seq(eight)) {
   tf[i] = 1
   eta = wholeDaysSignal(d, EuroDollars, tf)
   for (r in seq(eight)) {
     sym = EuroDollars[r]
     minpxincr = as.numeric(execute(database, paste0("select minpxincr from instinfo where inst=`",substr(sym,1,2))))
     rp = eta$rp[c((open+1):close),] # keep one extra for the forward price
     pp = eta$pp
     signal = ((pp - head(rp,-1))/minpxincr)[,r]
     fp = diff(rp)[,r]/minpxincr
     linear_model = lm(fp ~ signal -1)
     slope = summary(linear_model)$coefficients[1]
     error = summary(linear_model)$coefficients[2]
     ratio = error/slope
     se_ratio[i,r] = abs(ratio)
     slope_table[i,r] = slope
     error_table[i,r] = error 
     correlationTable[i,r] = cor(fp, signal)
   }
 }

sumMeasure_ratio = c()
sumMeasure_cor = c()
for (i in seq(eight)){
  sumMeasure_ratio[i] = sum(se_ratio[i,])
  sumMeasure_cor[i] = sum(correlationTable[i,])
}
colnames(correlationTable) = EuroDollars
colnames(se_ratio) = EuroDollars
etalabeler =c("eta1","eta2","eta3","eta4","eta5","eta6","eta7","eta8")
rownames(correlationTable) = etalabeler
rownames(se_ratio) = etalabeler

write.csv(correlationTable, 'correlationtble.csv')
write.csv(se_ratio, "slopeandseratio.csv")
write.csv(sumMeasure_cor, "correlation sum.csv")
write.csv(sumMeasure_ratio, "ratio sum.csv")

eta = wholeDaysSignal(d, EuroDollars, c(1,0,0,0,0,0,0,0))
r = 4
sym = EuroDollars[r]
minpxincr = as.numeric(execute(database, paste0("select minpxincr from instinfo where inst=`",substr(sym,1,2))))
rp = eta$rp[c((open+1):close),] 
pp = eta$pp
signal = ((pp - head(rp,-1))/minpxincr)[,r]
fp = diff(rp)[,r]/minpxincr

pdf(paste0("GEZ8 2017-10-10.pdf"))
plot.new()
par(las=1)
plot.window( xlim=c(min(signal)*1.05,max(signal)*1.05), ylim=c(min(fp)*1.05,max(fp)*1.05))
axis(1)
axis(2)
title( xlab="signals in ticks", ylab="foward price in ticks", main = paste("GEZ8 2017-10-10") )
for(i in seq(length(signal))){
  points(as.numeric(signal[i]), fp[i], col = "orange", pch = 25)
}
abline(b = slope_table[1,4],a=0, col = "green")
legend(0.4,-0.75,
       legend = c(paste0(round(slope_table[1,4],4), " \u00B1 ",round(error_table[1,4],4))),
       col = c("orange"),
       pch = c(25))
dev.off()