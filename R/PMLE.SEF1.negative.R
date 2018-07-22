PMLE.SEF1.negative <-
function(u.trunc,y.trunc,v.trunc, 
                            tau1=min(y.trunc),epsilon=0.0001){
  
  n=length(y.trunc)
  u.max=pmax(u.trunc,tau1)
  
  score_func=function(eta){
    S0=-exp(eta*v.trunc)+exp(eta*u.max)
    S1=-v.trunc*exp(eta*v.trunc)+u.max*exp(eta*u.max)
    n/eta+sum(y.trunc)-sum(S1/S0)
  }
  
  Hessian_func=function(eta){
    S0=-exp(eta*v.trunc)+exp(eta*u.max)
    S1=-v.trunc*exp(eta*v.trunc)+u.max*exp(eta*u.max)
    S2=-(v.trunc^2)*exp(eta*v.trunc)+u.max^2*exp(eta*u.max)
    -n/eta^2-sum(S2/S0)+sum( (S1/S0)^2 )
  }
  
  #Newton-Raphson
  eta=c()
  eta[1]=1/(tau1-mean(y.trunc))
  k=1
  repeat{
    eta[k+1]=eta[k]-score_func(eta[k])/Hessian_func(eta[k])
    if(abs(eta[k+1]-eta[k])<epsilon) break
    k=k+1
  }
  
  eta_hat=eta[k+1] # estimate
  V=-1/Hessian_func(eta[k+1])
  eta_SE=sqrt( V )
  
  #log-likelihood function
  logL_func=function(eta){
    S=sum( log( -exp(eta*v.trunc)+exp(eta*u.max) ) )
    n*log(-eta)+eta*sum(y.trunc)-S
  }
  logL=logL_func(eta_hat)
  p=1 #number of parameter
  AIC=-2*logL+2*p
  
  convergence_res=c(logL=logL,DF=p,AIC=AIC,No.of.iterations=k)
  
  list(eta1=eta_hat,SE1=eta_SE,
       convergence=convergence_res,
       Score=score_func(eta_hat),
       Hessian=Hessian_func(eta_hat)
  )

}
