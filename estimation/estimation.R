load("fdata.RData")

Ini_par   = list()
Ini_par$T = 50 
Ini_par$N = 5
Ini_par$K = 3
Ini_par$J = 4
Ini_par$T.list = c(10, 50, 300, 400,126)
Ini_par$p.list = c(50,100,200,300,126) 
Ini_par$SDGamma= 0.0027


Returns   = month_return2
T_s        = Ini_par$T.list[5]
## number of observed covariates
p_s        = 447

X.array      <-   array(0,dim=c(p_s,T_s,Ini_par$N))
X.array[,,1]      =  apply(MarketCap2[,c(2001:2126)]/10000,2,scale,center = T,scale=T)
X.array[,,2]      =  apply(1/Price_to_book2[,c(2001:2126)],2,scale,center = T,scale=T)
X.array[,,3]      =  apply(1/Price_to_earnings2[,c(2001:2126)],2,scale,center = T,scale=T)
X.array[,,4]      =  apply(momentum2[,c(2001:2126)],2,scale,center = T,scale=T)
X.array[,,5]      =  apply(volatility2[,c(2001:2126)],2,scale,center = T,scale=T)



sum.X  <-matrix(0, p_s,Ini_par$N)
for (t in 1: T_s)
{
  
  sum.X        = sum.X+ X.array[,t,]    
  
}


Xbar = sum.X/T_s


ra         = range(Xbar)
numk       = 3
mid        = seq(ra[1], ra[2], length.out = numk +2)[2: (numk+1)]
deg        = 3
numb       = 3+ deg + 1
designX         = array(0,dim = c(  p_s, numb, Ini_par$N))

for (i in 1: Ini_par$N)
{
  designX[,,i]     = bs(Xbar[,i], df = NULL, knots = mid, degree = deg, intercept = T,  Boundary.knots = range(Xbar))
}  

designXnew = matrix(0,p_s,numb*Ini_par$N)

for (i in 1:p_s)
{
  designXnew[i,]    = as.vector(designX[i,,])
}



coef                     = seq(0.00001,0.01, length.out = numb*Ini_par$K*Ini_par$N)
coef.matrix              = matrix(coef,numb*Ini_par$N,Ini_par$K )



Ymatrix    = Returns[,c(1:T_s)]


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


### this is G(bar{X}
J = numb

fitted = designXnew %*% B


Gamma = 1/T_s * YY %*% F - fitted
SDGamma = sd(Gamma)


#predict
GF=fitted%*% t(F)
Y_p = Iahat  + GF

plot(c(1:126),Returns[199,c(1:T_s)],"l",col = "black", lwd = 1.5,main = "prediction result", ylab = "excess return", xlab = "date")
lines(c(1:126),Y_p[199,], "l",col = "orange", lwd = 1.5)

plot(c(1:126),Returns[199,c(1:T_s)], "l",col = "black", lwd = 1.5,main = "prediction result with factor", ylab = "excess return", xlab = "date")
lines(c(1:126),Iahat[199,],"l",col = "red", lwd = 1.5)
lines(c(1:126),GF[199,],"l",col = "blue", lwd = 1.5)
legend(5,- 0.4, legend=c("returns","predicted returns","linear part","non-linear part"),
       col=c("black","orange","red", "blue"),lty="solid", cex=0.8 )




par(mfrow = c(2,3))

fitted.Size = designXnew[,1:J] %*% B[1:J,]
fitted.Size1 = fitted.Size[order(Xbar[,1]),]
S = sort(Xbar[,1])
plot(S, fitted.Size1[,1], "l", main = "Size Characteristics", ylab = "Size effect", xlab = "",ylim = c(-0.8,0.8))
lines(S, fitted.Size1[,2], lty="dashed")
lines(S, fitted.Size1[,3], lty="dotted")


fitted.Val = designXnew[,(J+1):(2*J)] %*% B[(J+1):(2*J),]
fitted.Val = fitted.Val[order(Xbar[,2]),]
Va = sort(Xbar[,2])
plot(Va, fitted.Val[,1], "l", 
     main = "Value Characteristics", ylab = "Value effect", xlab = "")
lines(Va, fitted.Val[,2], lty="dashed")
lines(Va, fitted.Val[,3], lty="dotted")
#legend(1.5, 0.18, c("Facter 1", "Factor 2", "Factor 3"), lty = c("solid", "dashed", "dotted"))

fitted.prof = designXnew[,(2*J+1):(3*J)] %*% B[(2*J+1):(3*J),]
fitted.prof = fitted.prof[order(Xbar[,3]),]
Va = sort(Xbar[,3])
plot(Va, fitted.prof[,1], "l", 
     main = "Profitability Characteristics", ylab = "Value effect", xlab = "")
lines(Va, fitted.prof[,2], lty="dashed")
lines(Va, fitted.prof[,3], lty="dotted")


fitted.Mom = designXnew[,(3*J+1):(4*J)] %*% B[(3*J+1):(4*J),]
fitted.Mom = fitted.Mom[order(Xbar[,4]),]
M = sort(Xbar[,4])
plot(M, fitted.Mom[,1], "l", 
     main = "Momentum Characteristics", ylab = "Momentum effect", xlab = "",ylim = c(-0.5,0.4))
lines(M, fitted.Mom[,2], lty="dashed")
lines(M, fitted.Mom[,3], lty="dotted")
#legend(1.5, -0.12, c("Facter 1", "Factor 2", "Factor 3"), lty = c("solid", "dashed", "dotted"))

fitted.Vol = designXnew[,(4*J+1):(5*J)] %*% B[(4*J+1):(5*J),]
fitted.Vol = fitted.Vol[order(Xbar[,5]),]
Vo = sort(Xbar[,5])
plot(Vo, fitted.Vol[,1], "l", 
     main = "Volatility Characteristics", ylab = "Volatility effect", xlab = "",ylim = c(-1.0,1.0))
lines(Vo, fitted.Vol[,2], lty="dashed")
lines(Vo, fitted.Vol[,3], lty="dotted")

