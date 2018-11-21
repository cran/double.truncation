PMLE.Weibull <-
function(u.trunc,y.trunc,v.trunc,epsilon = 1e-5,D1=2,D2=2,d1=2,d2=2){
  
  n=length(y.trunc)
  p=2 #number of parameter
  
  # define functions for fi
  fi_func = function(w){exp( w - exp(w) )}
  Fi_func = function(w){1 - exp( -exp(w) )}
  fi1_func = function(w){(1 - exp(w)) * exp(w - exp(w))}
  fi2_func = function(w){( -exp(w) + (1 - exp(w))^2) * exp(w - exp(w))}
  
  # define functions for Hui, Hvi, Kui, Kvi
  Hvi = function(para){
    m = para[1]
    s = exp( para[2] )
    Zy = ( y.trunc - m ) / s
    Zu = ( u.trunc - m ) / s
    Zv = ( v.trunc - m ) / s
    A = fi_func(Zv)
    B = Fi_func(Zv) - Fi_func(Zu)
    return(A/B)
  }
  
  Hui = function(para){
    m = para[1]
    s = exp( para[2] )
    Zy = ( y.trunc - m ) / s
    Zu = ( u.trunc - m ) / s
    Zv = ( v.trunc - m ) / s
    A = fi_func(Zu)
    B = Fi_func(Zv) - Fi_func(Zu)
    return(A/B)
  }
  
  Kvi = function(para){
    m = para[1]
    s = exp( para[2] )
    Zy = ( y.trunc - m ) / s
    Zu = ( u.trunc - m ) / s
    Zv = ( v.trunc - m ) / s
    A = fi1_func(Zv)
    B = Fi_func(Zv) - Fi_func(Zu)
    return(A/B)
  }
  
  Kui = function(para){
    m = para[1]
    s = exp( para[2] )
    Zy = ( y.trunc - m ) / s
    Zu = ( u.trunc - m ) / s
    Zv = ( v.trunc - m ) / s
    A = fi1_func(Zu)
    B = Fi_func(Zv) - Fi_func(Zu)
    return(A/B)
  }
  
  # likelihood function
  lik = function(para){
    m = para[1]
    s = exp( para[2] )
    Zy = ( y.trunc - m ) / s
    Zu = ( u.trunc - m ) / s
    Zv = ( v.trunc - m ) / s 
    A = length(y.trunc) * log(s)
    B = sum( log(fi_func(Zy)) )
    C = sum( log(Fi_func(Zv)  - Fi_func(Zu)) )
    return( - A + B - C )
  }
  
  # score function 
  first_mu = function(para){
    m = para[1]
    s = exp( para[2] )
    Zy = ( y.trunc - m ) / s
    Zu = ( u.trunc - m ) / s
    Zv = ( v.trunc - m ) / s
    A = sum( fi1_func(Zy) / fi_func(Zy) )
    B = sum( Hvi(para) - Hui(para) ) 
    total = (1/s) * (-A + B)
    return(total)	
  }
  
  first_sig = function(para){
    m = para[1]
    s = exp( para[2] )
    Zy = ( y.trunc - m ) / s
    Zu = ( u.trunc - m ) / s
    Zv = ( v.trunc - m ) / s
    A = length(y.trunc) / s
    B = sum( Zy * fi1_func(Zy) / fi_func(Zy) ) / s
    C = sum( Zv * Hvi(para) - Zu*Hui(para) ) / s
    total = -A - B + C
    return(total)
  }
  
  score_fun = function(para){
    SF1 = first_mu(para)
    SF2 = first_sig(para) * exp( para[2] )
    A = matrix(c(SF1, SF2), nrow = 2, ncol = 1)
    return(A)	
  }
  
  # Hessian matrix
  second_mu = function(para){
    m = para[1]
    s = exp( para[2] )
    Zy = ( y.trunc - m ) / s
    Zu = ( u.trunc - m ) / s
    Zv = ( v.trunc - m ) / s
    A = sum( fi2_func(Zy) / fi_func(Zy) )
    B = sum( (fi1_func(Zy) / fi_func(Zy))^2 )
    C = sum( Kvi(para) - Kui(para) )
    D = sum( (Hvi(para) - Hui(para))^2 )
    total = (1/s^2) * ( A - B - C + D )
    return(total)
  }
  
  second_sig = function(para){
    m = para[1]
    s = exp( para[2] )
    Zy = ( y.trunc - m ) / s
    Zu = ( u.trunc - m ) / s
    Zv = ( v.trunc - m ) / s
    A = length(y.trunc)
    B = sum( (2 * Zy * fi1_func(Zy) + Zy^2 * fi2_func(Zy)) / fi_func(Zy) )
    C = sum( (Zy * fi1_func(Zy) / fi_func(Zy))^2 )
    D = 2 * sum( Zv * Hvi(para) - Zu * Hui(para) )
    E = sum( Zv^2 * Kvi(para) - Zu^2 * Kui(para) )
    F = sum( (Zv * Hvi(para) - Zu * Hui(para))^2 )
    total = (1/s^2) * ( A + B - C - D - E + F )
    return(total)	
  }
  
  second_musig = function(para){
    m = para[1]
    s = exp( para[2] )
    Zy = ( y.trunc - m ) / s
    Zu = ( u.trunc - m ) / s
    Zv = ( v.trunc - m ) / s
    A = sum( (fi1_func(Zy) + Zy * fi2_func(Zy)) / fi_func(Zy) )
    B = sum( Zy * (fi1_func(Zy) / fi_func(Zy))^2 )
    C = sum( Zv * Kvi(para) - Zu * Kui(para) )
    D = sum( (Hvi(para) - Hui(para)) * (Zv * Hvi(para) - Zu * Hui(para) - 1) )
    total = (1/s^2) * ( A - B - C + D )
    return(total)
  }
  
  hessian_fun = function(para){	
    H11 = second_mu(para)
    H22 = second_sig(para) * exp(2*para[2]) + score_fun(para)[2]
    H12 = second_musig(para) * exp( para[2] )
    H21 = H12
    A = matrix( c(H11, H12, H21, H22), nrow = 2, ncol = 2 )
    return(A)
  }

  # RNR-algorithm
  random = 0
  theta = matrix(,1, 2)
  theta_old = matrix(, 1, 2)
  theta_new = matrix(, 1, 2)
  
  # initial values 
  alpha_s = pi / ( sqrt(6) * sd(y.trunc))
  lambda_s = 1 / mean( exp(y.trunc*alpha_s) )
  theta_0 = c( -log(lambda_s)/alpha_s, log(1/alpha_s) )
  #theta_0 = c( median(y.trunc), log(IQR(y.trunc)/2) )
  
  # if we change initial value, the initial value in repeat function shoud also change
  theta[1, ] = theta_0
  
  h = 1
  repeat{	
    theta_old = theta[1:h, ]
    theta_new = theta[h, ]-solve(hessian_fun(theta[h, ]))%*%score_fun(theta[h, ])
    theta = matrix(,h+1,2)
    theta[1:h, ] = theta_old
    theta[h+1, ] = c(theta_new[1, ], theta_new[2, ])
    if (1*is.nan(theta)[h+1, 1]==1){
      theta[1, ] = theta_0+c(runif(1,-d1,d1), runif(1,-d2,d2))
      h = 0
      random = random + 1
    }else if (abs(theta[h+1,1]-theta[h,1])>D1 | abs(theta[h+1,2]-theta[h,2])>D2){
      theta[1, ] = theta_0 + c(runif(1, -d1, d1), runif(1, -d2, d2))
      h = 0
      random = random + 1
    }else if (max(abs(theta[h+1, ]-theta[h, ]))<epsilon){	
      break
    }
    h = h + 1
  }
  theta
  muhat = theta[h+1, 1]
  sighat = exp(theta[h+1, 2])
  
  # likelihood function (without transformation)
  likn = function(para){
    m = para[1]
    s = para[2]
    Zy = ( y.trunc - m ) / s
    Zu = ( u.trunc - m ) / s
    Zv = ( v.trunc - m ) / s
    A = length(y.trunc) * log(s)
    B = sum( log( fi_func(Zy) ) )
    C = sum( log( Fi_func(Zv)  - Fi_func(Zu) ) )
    return( - A + B - C )
    
  }
  logL = likn(c(muhat, sighat))
 
  # transform the hessian to (mu, sigma)
  hessian_ms = function(para){
    m = para[1]
    s = para[2]
    trans = matrix(c(1,0,0,s^(-1)), ncol = 2, nrow = 2)
    A = matrix(, 2, 2)
    A[1, 1] = hessian_fun(c(m, log(s)))[1, 1]
    A[2, 2] = (hessian_fun(c(m, log(s)))[2,2] - score_fun(c(m, log(s)))[2])/s^2
    A[2, 1] = hessian_fun(c(m, log(s)))[2, 1]/s
    A[1, 2] = A[2, 1]
    B = trans%*%A%*%trans
    return(A)
  }
  
  ## Interval estimates ##
  info = - hessian_ms(c(muhat, sighat))
  alpha=0.05
  SE_mu=sqrt(solve(info)[1, 1])
  Low_mu=muhat+qnorm(alpha/2,0,1)*SE_mu
  Up_mu=muhat+qnorm(1-alpha/2,0,1)*SE_mu
  SE_sig=sqrt(solve(info)[2, 2])
  Low_sig=sighat*exp(qnorm(alpha/2, 0, 1)*SE_sig/sighat)
  Up_sig=sighat*exp(qnorm(1-alpha/2, 0, 1)*SE_sig/sighat)
  
  mu_res=c(estimate=muhat,SE=SE_mu,Low=Low_mu,Up=Up_mu)
  sig_res=c(estimate=sighat,SE=SE_sig,Low=Low_sig,Up=Up_sig)
  
  p = 2
  AIC = -2 * logL + 2 * p
  convergence_res = c(logL = logL, DF = p, AIC = AIC, No.of.iterations = h)
  list(estimate = c(mu = muhat, sigma =sighat), 
       SE = c(mu = SE_mu, sigma = SE_sig), 
       convergence = convergence_res, Score = as.vector(score_fun(theta[h, ])), 
       Hessian = -info)
  #round(cbind(theta[,1],exp(theta[,2])),2)
}
