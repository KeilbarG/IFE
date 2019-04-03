rm(list = ls())
setwd("C:/Users/wangbing.hub/Desktop/wangbing")
library(mvtnorm)
library(yaml)
library("splines")
library(nlme)
library(mgcv)
library(lattice)
library(Matrix)
library(MASS)
source("helpfun.r")
library("ggplot2")
install.packages("ggfortify")
load("C:/Users/wangbing.hub/Desktop/wangbing/simdata1,1.RData")
###simulate observations

Ini_par   = list()
Ini_par$q.list = seq(4, 10, by = 1)
Ini_par$K = 3
Ini_par$J.list = seq(4, 10, by = 1)
Ini_par$T.list = c(3, 5, 10, 50, 126, 400)
Ini_par$N.list = c(3, 5, 10, 50, 126, 400)
Ini_par$SDGamma= 0.0027



  ## Generate F
      ## number of time point 
    T_s        = Ini_par$T.list[5]
      ## number of observed covariates
    N_s        = Ini_par$N.list[5]
      ## number of individuals
    q_s        = Ini_par$q.list[3]
    Ini_par$J  = Ini_par$J.list[1]
    
    
    
  estimate_fun = function(sample, T_s, N_s, q_s){
    
    beta_re = 0
    f_re =  0
    y_re =  0
    
    for (p in 1:sample){
    ### from Fan's paper
    veccov_f    = c(0.9076,0.0049,0.0230,0.0049,0.8737,0.0403,0.0230,0.0403,0.9266)
    Sigma_eps   = matrix(veccov_f, 3, 3)
    vecphi      = c(-0.0371,-0.1226,-0.1130,-0.2339,0.1060,-0.2793,0.2803,-.0755,-0.0529)
    phi         = matrix(vecphi, 3, 3, byrow = TRUE)
    mu_f        = rep(0,Ini_par$K)
    mu          = rep(0,Ini_par$K)
    ## mu_f    = solve(diag(rep(1,Ini_par$K)) - phi) %*% mu
    
  ## K is the number of factors
    F_s     = matrix(0, T_s, Ini_par$K)
    

  ####calculate the variance of f
    Sigma_eps_vec = matrix(0, Ini_par$K^2, 1)
    Sigma_eps_vec = as.vector(Sigma_eps)
    phi_tensor    = tensor(phi)
    cov_f_vec     = solve(diag(rep(1,Ini_par$K^2)) - phi_tensor) %*% Sigma_eps_vec
    cov_f = matrix(0, Ini_par$K,Ini_par$K)
    for(i in 1:Ini_par$K){
     cov_f[,i] = cov_f_vec[(Ini_par$K*(i-1)+1):(Ini_par$K*i)]
    }

  
  ### simulate F
    F_s[1,] =mvrnorm(1, mu_f, cov_f)
 
    for(t in 2:T_s){
      eps_s   = mvrnorm(1, rep(0,Ini_par$K), Sigma_eps)
      F_s[t,] = mu + phi %*% F_s[t-1,] + eps_s
      }
    

  #### simulate design
    X.array<-array(0,dim=c(N_s,T_s,q_s))
    sum.X  <-matrix(0, N_s,q_s)
    for (t in 1: T_s){
      
      X.array[,t,] = mvrnorm(N_s,rep(0,q_s), diag(rep(0.1,q_s)))
      
      sum.X        = sum.X+ X.array[,t,]    
    }
  
    
    ## Xbar is X average over time
    Xbar = sum.X/T_s
    
    ra         = range(Xbar)
    numk       = 3
    mid        = seq(ra[1], ra[2], length.out = numk +2)[2: (numk+1)]
    deg        = 2
    numb       = 3+ deg + 1
    
    designX         = array(0,dim = c(  N_s, numb, q_s))
    
    for (i in 1: q_s)
    {
      designX[,,i]     = bs(Xbar[,i], df = NULL, knots = mid, degree = deg, intercept = T,  Boundary.knots = range(Xbar))
    }  
    

    designXnew = matrix(0,N_s,numb*q_s)
    
    for (i in 1:N_s)
    {
      designXnew[i,]   = as.vector(designX[i,,])
    }
    
    
    coef                     = seq(0.00001,0.01, length.out = numb*Ini_par$K*q_s)
    coef.matrix              = matrix(coef,numb*q_s ,Ini_par$K )
    
    
    R_xbar                   = t(mvrnorm(Ini_par$K, rep(0,N_s), diag(rep(0.05,N_s))))
    G.array                  = designXnew%*%coef.matrix + R_xbar
    
  
  #### simulate g matrix
     ###G.array    = array(0,dim=c(p_s,Ini_par$N,Ini_par$K))
     ###G.array[,,1] = Xbar
     ###G.array[,,2] = Xbar^2
     ###G.array[,,3] = Xbar^3-2*Xbar
    

    
    ### simulate error 
    Sigma_u = diag(rep(0.5,N_s))
    U_s = t(mvrnorm(T_s, rep(0,N_s), Sigma_u))
    #U_s = t(rmvt(T_s, Sigma_u, rep(0,p_s) ))
    Gamma_s = t(mvrnorm(Ini_par$K, rep(0,N_s), diag(rep(Ini_par$SDGamma^2,N_s))))
    #Gamma_s = t(rmvt(Ini_par$K, diag(rep(Ini_par$SDGamma^2,p_s)), rep(0,p_s)))
    
    
    Y_s = (G.array + Gamma_s) %*% t(F_s) + U_s

    ### interactive fixed effect parameter


    #beta = matrix(seq(-0.5,0.5, length.out = Ini_par$N ), Ini_par$N, 1)
    set.seed(123)
    beta = runif(q_s,-1,1)
    
    Ia   = matrix(0, N_s, T_s)
    
    for (i in 1:N_s)
    {
      Ia[i,] = X.array[i,,]%*%beta
    }
    
    
    Ymatrix = Y_s + Ia
       

   ##### estimation

    A = solve(t(designXnew) %*% designXnew+ 0.00001*diag(numb*q_s))
    P = designXnew %*% A %*% t(designXnew)
    
    
    PPY     = (diag(N_s) - P)%*% Ymatrix 
    sumaa   = matrix(0, q_s, q_s)
    sumbb   = matrix(0, q_s,1)
    
    for (i in 1: T_s )
    {
      
      aa =  (t(X.array[,i,])%*%(diag(N_s)- P)%*%X.array[,i,])
      bb =  t(X.array[,i,])%*%(diag(N_s)- P)%*%Ymatrix[,i]
      
      sumaa = aa+ sumaa
      sumbb = bb+ sumbb
    } 
    
    betahat = solve(sumaa)%*%sumbb
    
    mse_beta     = sum((beta- betahat)^2)/4
    
    rmse_beta    = mse_beta^(1/2)
    
    ##### calculate the factors and loadings
    Iahat = matrix(0, N_s, T_s )
    
    for (i in 1:N_s)
    {
      Iahat[i,] = X.array[i,,]%*%betahat
    }
    
    YY = Ymatrix - Iahat 
    
    
    eigObj = eigen(t(YY) %*% P %*%YY)
    F = eigObj$vectors[,1:Ini_par$K] * sqrt(T_s)
    B = 1/T_s * A %*% t(designXnew) %*%YY %*% F
    
    #predict
    fitted = designXnew %*% B
    GF = fitted%*% t(F)
    Y_p = Iahat  + GF
    mse_y = sum(t(Ymatrix - Y_p)%*%(Ymatrix - Y_p))/(T_s*N_s)
    rmse_y = mse_y^(1/2)
    
    mse_f = sum(t(F - F_s)%*%(F - F_s))/(T_s*Ini_par$K)
    rmse_f = mse_f^(1/2)

    beta_re = beta_re + rmse_beta
    f_re = f_re + rmse_f
    y_re = y_re + rmse_y
    
    p= p+1
    }
    return(list(beta_re/sample, f_re/sample, y_re/sample))
  }
  
#monte carlo simulation 
  re34 = estimate_fun(1000,3,400,6)
  re35 = estimate_fun(1000,5,400,6)
  re36 = estimate_fun(1000,10,400,6)
  re37 = estimate_fun(1000,50,400,6)
  re38 = estimate_fun(1000,126,400,6)
  re39 = estimate_fun(1000,400,400,6)
  re40 = estimate_fun(1000,400,3,6)
  re41 = estimate_fun(1000,400,5,6)
  re42 = estimate_fun(1000,400,10,6)
  re43 = estimate_fun(1000,400,50,6)
  re44 = estimate_fun(1000,400,126,6)
 
  simu = rbind(t(re34),t(re35),t(re36),t(re37),t(re38),t(re39),
               t(re40),t(re41),t(re42),t(re43),t(re44))
  sim3 = as.data.frame(simu)
  simre3 = data.frame(lapply(sim3, as.numeric))
  write.csv(simre3, "C:/Users/wangbing.hub/Desktop/wangbing/simre2.csv")
  #plot F
    library(ggfortify)
    df <-t(YY) %*% P %*%YY
    autoplot(prcomp(df))
    df1=G.array[,1]
    df2=G.array[,2]
    plot(df1,df2)
    ### this is G(bar{X})
    
    fitted = designXnew %*% B

    ### calibrating Gamma

    Gamma = 1/T_s * YY %*% F - fitted
    SDGamma = sd(Gamma)
    Gammareal = Re(Gamma)
    hist(Gammareal, freq=FALSE, breaks = 20)

    ### qqplot(Gammareal)
    qqline(Gammareal)



     
    

