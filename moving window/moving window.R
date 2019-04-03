rm(list = ls())
setwd("C:/Users/wangbing.hub/Desktop/wangbing")
library(psych)
library(splines)
library(nlme)
library(mgcv)
library(lattice)
library(Matrix)
library(MASS)
library(tseries)
library(ggplot2)
library(plm)
library(survival)
library(stats)
source("helpfun.r")
#load data
load("fdata.RData")
MarketCap1 = read.csv("MarketCap_f.csv")
Price_to_book1 = read.csv("PB_f.csv")
Price_to_earnings1 = read.csv("PE_f.csv")
return1 = read.csv("exc_return_f.csv")
return25 = read.csv("exc_return_n(25).csv")
momentum1 = read.csv("momentum_diff_f.csv")
volatility1 = read.csv("volatility_diff_f.csv")
primium1 =  read.csv("mk_primium_f.csv")



#dimension
MarketCap2 = t(apply(MarketCap_n[,-1],1,as.numeric))
Price_to_book2= t(apply(PB_n[,-1],1,as.numeric))
Book_to_price2= 1/Price_to_book2
Price_to_earnings2 = t(apply(PE_n[,-1],1,as.numeric))
return2 = t(apply(return25[,-c(1:2)],1,as.numeric))
momentum2  = t(apply(momentum_n[,-c(1:2)],1,as.numeric))
volatility2 = t(apply(volatility_n[,-1],1,as.numeric))
primium2 = t(apply(mk_primium_n[,-1],1,as.numeric))


#plot

PRICE = read.csv("marketIn.csv")
# Scatterplot Matrices
vec_tr  <- function(dt)
{
  nr = nrow(dt)
  nl = ncol(dt)
  eps_vec = matrix(0, nr*nl, 1)
  for(i in 1:nr){
    eps_vec[(nl*(i-1)+1):(nl*i)] = dt[i,]
  }
  return(eps_vec)
}


win = c(2560:2565)

scatter_dt = array(0,dim=c(447*length(win),7))
scatter_dt[,1]      =  vec_tr(primium2[,win])
scatter_dt[,2]      =  vec_tr(MarketCap2[,win]/10000)
scatter_dt[,3]      =  vec_tr(1/Price_to_book2[,win])
scatter_dt[,4]      =  vec_tr(1/Price_to_earnings2[,win])
scatter_dt[,5]      =  vec_tr(momentum2[,win])
scatter_dt[,6]      =  vec_tr(volatility2[,win])
scatter_dt[,7]      =  vec_tr(month_return2[,win+5])

dt = as.data.frame(scatter_dt[,2:7])

p <- plot_ly(data = dt, x = ~V2, y = ~V7)
p1 <- plot_ly(data = dt, x = ~V3, y = ~V7)
p2 <- plot_ly(data = dt, x = ~V4, y = ~V7)
p3 <- plot_ly(data = dt, x = ~V5, y = ~V7)
p4 <- plot_ly(data = dt, x = ~V6, y = ~V7)
subplot(p,p1,p2,p3,p4,nrows = 2, shareX = TRUE, shareY = TRUE,
        titleY = FALSE, titleX = FALSE)

par(mfrow = c(1,1))
summary(dt)
plot(density(scatter_dt[,4]))
boxplot(dt)
skewness(scatter_dt)
kurtosis(scatter_dt)

par(mfrow = c(2,3))
#plot(scatter_dt[,1],scatter_dt[,7],xlab = "market premium",
#ylab= "excess return")
plot(scatter_dt[,2],scatter_dt[,7],xlab = "market cap",
     ylab= "excess return")
plot(scatter_dt[,3],scatter_dt[,7],xlab = "book_to_market",
     ylab= "excess return")
plot(scatter_dt[,4],scatter_dt[,7],xlab = "earnings_to_price",
     ylab= "excess return")
plot(scatter_dt[,5],scatter_dt[,7],xlab = "momentum",
     ylab= "excess return")
plot(scatter_dt[,6],scatter_dt[,7],xlab = "volatility",
     ylab= "excess return")

win2 = c(1:2577)
scatter_dt1 = array(0,dim=c(447,length(win2),7))
scatter_dt1[,,1]      =  primium2[,win2]
scatter_dt1[,,2]      =  MarketCap2[,win2]/10000
scatter_dt1[,,3]      =  1/Price_to_book2[,win2]
scatter_dt1[,,4]      =  1/Price_to_earnings2[,win2]
scatter_dt1[,,5]      =  momentum2[,win2]
scatter_dt1[,,6]      =  volatility2[,win2]
scatter_dt1[,,7]      =  month_return2[,win2+5]
par(mfrow = c(2,3))
plot(scatter_dt1[,2562,1],scatter_dt1[,2567,7],xlab = "market premium",
     ylab= "excess return")
plot(scatter_dt1[,2562,2],scatter_dt1[,2567,7],xlab = "market cap",
     ylab= "excess return")
plot(scatter_dt1[,2562,3],scatter_dt1[,2567,7],xlab = "book_to_market",
     ylab= "excess return")
plot(scatter_dt1[,2562,4],scatter_dt1[,2567,7],xlab = "earnings_to_price",
     ylab= "excess return")
plot(scatter_dt1[,2562,5],scatter_dt1[,2567,7],xlab = "momentum",
     ylab= "excess return")
plot(scatter_dt1[,2562,6],scatter_dt1[,2567,7],xlab = "volatility",
     ylab= "excess return")


Ini_par   = list()
Ini_par$T = 50 
Ini_par$N = 5
Ini_par$K = 3
Ini_par$J = 4
Ini_par$T.list = c(10, 50, 300, 400)
Ini_par$p.list = c(50,100,200,300) 
Ini_par$SDGamma= 0.0027


#estimation function

roll_fun = function(width, gap){

beta_result = c()
test.stat = c()
#loop
for (u in seq(from=1,to=2582,by=gap)){
#covariates
if( (u+2*width-1)<2582 ){
dt_ra = c(u:(u+width-1))
p_s       = nrow(return2[,dt_ra])
T_s       = ncol(return2[,dt_ra])

X.array      <-   array(0,dim=c(p_s,T_s,Ini_par$N))

#X.array[,,1]     =  primium2[,dt_ra]
X.array[,,1]      =  MarketCap2[,dt_ra]/10000
X.array[,,2]      =  1/Price_to_book2[,dt_ra]
X.array[,,3]      =  1/Price_to_earnings2[,dt_ra]
X.array[,,4]      =  momentum2[,dt_ra]
X.array[,,5]      =  volatility2[,dt_ra]


#X.array[,,2]      =  apply(MarketCap2[,dt_ra],2,scale,center = T,scale=T)
#X.array[,,3]      =  apply(Price_to_book2[,dt_ra],2,scale,center = T,scale=T)
#X.array[,,4]      =  apply(Price_to_earnings2[,dt_ra],2,scale,center = T,scale=T)
#X.array[,,5]      =  apply(momentum2[,dt_ra],2,scale,center = T,scale=T)
#X.array[,,6]      =  apply(volatility2[,dt_ra],2,scale,center = T,scale=T)



sum.X  <-matrix(0, p_s,Ini_par$N)
for (t in 1: T_s)
{
  sum.X        = sum.X+ X.array[,t,]    
}


Xbar = sum.X/T_s

ra         = range(Xbar)#boundary knots
numk       = 3#number of factors
mid        = seq(ra[1], ra[2], length.out = numk +2)[2: (numk+1)]#inner knots(quantiles)
deg        = 3#degree
numb       = 3+ deg + 1
designX         = array(0,dim = c( p_s, numb, Ini_par$N))


for (i in 1: Ini_par$N)
{
  designX[,,i]     = bs(Xbar[,i], df = NULL, knots = mid, degree = deg, intercept = T,  Boundary.knots = range(Xbar))
}  

designXnew = matrix(0,p_s,numb*Ini_par$N)

for (i in 1:p_s)
{
  designXnew[i,]    = as.vector(designX[i,,])
}



Ymatrix    = month_return2[,dt_ra+width]


A = solve(t(designXnew) %*% designXnew+ 0.00001*diag(numb*Ini_par$N))
P = designXnew %*% A %*% t(designXnew)


PPY     = (diag(p_s) - P)%*% Ymatrix 
sumaa   = matrix(0, Ini_par$N, Ini_par$N)
sumbb   = matrix(0, Ini_par$N,1)

for (i in 1: T_s )
{
  
  aa =  (t(X.array[,i,])%*%(diag(p_s)- P)%*%X.array[,i,])
  bb =  t(X.array[,i,])%*%(diag(p_s)- P)%*%Ymatrix[,i]
  
  sumaa = aa+ sumaa
  sumbb = bb+ sumbb
} 

betahat = solve(sumaa)%*%sumbb

Iahat = matrix(0, p_s, T_s )

for (i in 1:p_s)
{
  
  Iahat[i,] = X.array[i,,]%*%betahat
  
}

YY = Ymatrix - Iahat


eigObj = eigen(t(YY) %*% P %*%YY)
F = eigObj$vectors[,1:Ini_par$K] * sqrt(T_s)
B = 1/T_s * A %*% t(designXnew) %*%YY %*% F


### this is G(bar{X})
J = numb

fitted = designXnew %*% B


Gamma = 1/T_s * YY %*% F - fitted
SDGamma = sd(Gamma)


###TSET STATISTIC
#C.array = array(0, dim = c(numb,Ini_par$N,Ini_par$N))
#for (i in 1:Ini_par$N){
#C = solve(A)%*%t(designXnew)%*% Xbar[,i]
#C.array[,,i] = C
#i = i+1
#}

#C.arraynew = array(0,dim = c(numb*Ini_par$N,Ini_par$N))
#for (i in 1:Ini_par$N){
# C.arraynew[,i]    = t(as.vector(C.array[,,i]))
#  i = i+1
#}

#appr = designXnew %*% as.matrix(C.arraynew)
#pi_t = array(0, dim = c(p_s,T_s,Ini_par$N))
#for (i in 1:T_s){
#  pi_t[,i,] = P %*% X.array[,i,] - appr
#i = i+1
#}

U = YY - fitted %*% t(F)
sigma_u = 0
for (i in 1:T_s){
sigma_u = sigma_u + tr(U[,1]%*%t(U[,1]))
i = i+1
}
V_b = 1/T_s * sigma_u
V_lu1=tr(Gamma %*% t(Gamma)) + V_b

V_pi = 1/(p_s * T_s) * sumaa

V_lu = 1/p_s * V_lu1 * V_pi

V_approx = solve(V_pi)%*% V_lu %*% solve(V_pi)

t.stat = c()
for (i in 1:Ini_par$N){
  t.stat[i] = abs(betahat[i]) * sqrt(p_s * T_s)/sqrt(V_approx[i,i])
i =i+1}


beta_result = cbind(beta_result,betahat)
test.stat = cbind(test.stat,t.stat)
}
}
return(list(beta_result, test.stat))
}


###################################################
load("10 data.RData")

result1 = roll_fun(10,10)
result2 = roll_fun(50,10)
result3 = roll_fun(126,10)
result4 = roll_fun(250,10)
result5 = roll_fun(400,10)

result6 = roll_fun(10,10)
result7 = roll_fun(50,10)
result8 = roll_fun(126,10)
result9 = roll_fun(250,10)
result10 = roll_fun(400,10)

par(mfrow = c(1,1))

lines(c(1:179),result1[[2]][4,-c(1:78)],col = "black",type="l",xlab = "",
     ylab= "beta value")
#lines(c(1:179),result2[[2]][2,-c(1:70)],col="black",type="l")
lines(c(1:179), result3[[2]][4,-c(1:54)],col="red",type="l",pch=18,xlab = "",
     ylab= "t statistic")
#lines(c(1:179), result4[[2]][2,-c(1:30)],col="grey",type="l")
plot(c(1:179), result5[[2]][4,],col="blue",type="l",pch=18,xlab = "",
     ylab= "t statistic")
abline(h=1.96,col="red")
legend(5,25, legend=c("T=10","T=126","T=400"),
       col=c("black","red", "blue"), lty=1:2, cex=0.8)

industry= read.csv("industry.csv")
names=cap_raw[c(pointer_a),2]
write.csv(names,"C:/Users/wangbing.hub/Desktop/wangbing/names1.csv")
ylim = c(0,20)
#################################################

#predict
GF=fitted%*% t(F)
Y_p = Iahat  + GF

#################################################
plot(c(1:100),return2[199,51:150],"l",col = "red", lwd = 1.5,main = "prediction result without factor", ylab = "excess return", xlab = "date")
lines(c(1:100),Iahat[199,51:150], "l",col = "blue", lwd = 1.5)

plot(c(1:100),return2[199,51:150], "l",col = "red", lwd = 1.5,main = "prediction result with factor", ylab = "excess return", xlab = "date")
lines(c(1:100),GF[199,51:150],"l",col = "blue", lwd = 1.5)


plot(c(1:10),return2[199,1:10], "l",col = "blue",lwd = 1.5,main = "prediction result without factor", ylab = "excess return", xlab = "date")
lines(c(1:10),Iahat2[199,],"l",col = "red")

plot(c(1:10),return2[199,1:10], "l",col = "blue", lwd = 1.5,main = "prediction result with factor", ylab = "excess return", xlab = "date")
lines(c(1:10),GF2[199,1:10],"l",col = "red")

#################################################################3






par(mfrow = c(1,1))

fitted.Size = designXnew[,1:J] %*% B[1:J,]
fitted.Size1 = fitted.Size[order(Xbar[,1]),]
S = sort(Xbar[,1])
plot(S, fitted.Size1[,1], "l", main = "Size Characteristics", ylab = "Size effect", xlab = "")
lines(S, fitted.Size1[,2], lty="dashed")
lines(S, fitted.Size1[,3], lty="dotted")


fitted.Size = designXnew[,2:J] %*% B[2:J,]
fitted.Size1 = fitted.Size[order(Xbar[,2]),]
S = sort(Xbar[,2])
plot(S, fitted.Size1[,1], "l", main = "Size Characteristics", ylab = "Size effect", xlab = "", ylim = c(-0.02,0.02))
lines(S, fitted.Size1[,2], lty="dashed")
lines(S, fitted.Size1[,3], lty="dotted")
#legend(1.5, 0.18, c("Facter 1", "Factor 2", "Factor 3"), lty = c("solid", "dashed", "dotted"))

fitted.Val = designXnew[,(2*J+1):(3*J)] %*% B[(2*J+1):(3*J),]
fitted.Val = fitted.Val[order(Xbar[,3]),]
Va = sort(Xbar[,3])
plot(Va, fitted.Val[,1], "l", 
     main = "Value Characteristics", ylab = "Value effect", xlab = "", ylim = c(-0.5,0.5))
lines(Va, fitted.Val[,2], lty="dashed")
lines(Va, fitted.Val[,3], lty="dotted")
#legend(1.5, 0.18, c("Facter 1", "Factor 2", "Factor 3"), lty = c("solid", "dashed", "dotted"))

fitted.Val = designXnew[,(3*J+1):(4*J)] %*% B[(3*J+1):(4*J),]
fitted.Val = fitted.Val[order(Xbar[,4]),]
Va = sort(Xbar[,4])
plot(Va, fitted.Val[,1], "l", 
     main = "Profitability Characteristics", ylab = "Value effect", xlab = "",ylim = c(-0.5,0.5))
lines(Va, fitted.Val[,2], lty="dashed")
lines(Va, fitted.Val[,3], lty="dotted")


fitted.Mom = designXnew[,(4*J+1):(5*J)] %*% B[(4*J+1):(5*J),]
fitted.Mom = fitted.Mom[order(Xbar[,5]),]
M = sort(Xbar[,5])
plot(M, fitted.Mom[,1], "l", 
     main = "Momentum Characteristics", ylab = "Momentum effect", xlab = "", ylim = c(-0.5,0.5))
lines(M, fitted.Mom[,2], lty="dashed")
lines(M, fitted.Mom[,3], lty="dotted")
#legend(1.5, -0.12, c("Facter 1", "Factor 2", "Factor 3"), lty = c("solid", "dashed", "dotted"))

fitted.Vol = designXnew[,(5*J+1):(6*J)] %*% B[(5*J+1):(6*J),]
fitted.Vol = fitted.Vol[order(Xbar[,6]),]
Vo = sort(Xbar[,6])
plot(Vo, fitted.Vol[,1], "l", 
     main = "Volatility Characteristics", ylab = "Volatility effect", xlab = "", ylim = c(-0.5,0.5))
lines(Vo, fitted.Vol[,2], lty="dashed")
lines(Vo, fitted.Vol[,3], lty="dotted")


#################################################################
purtest(volatility2[,1:500], test = "levinlin",lags = "SIC", pmax = 10)

names = read.csv("names.csv")
names1 = names[,2]
dimnames(X.array)<-list(firm = names1,dates = c(1:500),c("prim","size","value","prof", "mom", "vol"))
dimnames(Ymatrix)<-list(firm = names1,dates = c(1:500))

array.all = array(0,dim=c(p_s,T_s,Ini_par$N+1))
array.all[,,1]      =  return2
array.all[,,2]      =  primium2
array.all[,,3]      =  apply(MarketCap2,2,scale,center = T,scale=T)
array.all[,,4]      =  apply(Price_to_book2,2,scale,center = T,scale=T)
array.all[,,5]      =  apply(Price_to_earnings2,2,scale,center = T,scale=T)
array.all[,,6]      =  apply(momentum2,2,scale,center = T,scale=T)
array.all[,,7]      =  apply(volatility2,2,scale,center = T,scale=T)


dimnames(array.all)<-list(firm = names1,dates = c(1:T_s),c("exc_return","prim","size","value","prof", "mom", "vol"))


d = data.frame(return =Ymatrix, prim = X.array[,,"prim"], size =X.array[,,"size"], 
               value =  X.array[,,"value"], prof = X.array[,,"prof"], mom = X.array[,,"mom"])
write.csv(X.array)
pgr = pdata.frame(X.array,index = c(firm,dates))

form = as.formula(return~ prim +size +value +prof +mom +vol)
wi<-plm(form,data=d,model="within")
re<-plm(form,data=d,model="random")

