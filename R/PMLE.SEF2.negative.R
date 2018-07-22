PMLE.SEF2.negative <-
function(u.trunc,y.trunc,v.trunc,epsilon=0.0001){
  
  n=length(y.trunc)
  p=2 #number of parameter
  
  #Define function
  H1=function(eta1,eta2){
    den1=pnorm((v.trunc+eta1/(2*eta2))/sqrt(-1/(2*eta2)))
    den2=pnorm((u.trunc+eta1/(2*eta2))/sqrt(-1/(2*eta2)))
    num=dnorm((v.trunc+eta1/(2*eta2))/sqrt(-1/(2*eta2)))
    num/(den1-den2)
  }
  
  H2=function(eta1,eta2){
    den1=pnorm((v.trunc+eta1/(2*eta2))/sqrt(-1/(2*eta2)))
    den2=pnorm((u.trunc+eta1/(2*eta2))/sqrt(-1/(2*eta2)))
    num=dnorm((u.trunc+eta1/(2*eta2))/sqrt(-1/(2*eta2)))
    num/(den1-den2)
  }
  
  CR1=function(eta1,eta2){  
    (-v.trunc/sqrt(-2*eta2))-sqrt(-2*eta2)*eta1/(4*eta2^2)
  }
  CR2=function(eta1,eta2){
    (-u.trunc/sqrt(-2*eta2))-sqrt(-2*eta2)*eta1/(4*eta2^2)
  }
  
  CR3=function(eta1,eta2){
    (v.trunc+eta1/(2*eta2))/sqrt(-1/(2*eta2))
  }
  
  CR4=function(eta1,eta2){
    (u.trunc+eta1/(2*eta2))/sqrt(-1/(2*eta2))
  }
  
  CR5=function(eta1,eta2){
    -v.trunc/(sqrt((-2*eta2)^3))-3*eta1/(sqrt((-2*eta2)^5))
  }
  
  CR6=function(eta1,eta2){
    -u.trunc/(sqrt((-2*eta2)^3))-3*eta1/(sqrt((-2*eta2)^5))
  }
  
  #score function
  SF=matrix(,2,1)
  score_func=function(ETA){
    eta1=ETA[1]
    eta2=ETA[2]
    h1=H1(eta1,eta2)
    h2=H2(eta1,eta2)
    cr1=CR1(eta1,eta2)
    cr2=CR2(eta1,eta2)
    
    SF[1,]=sum(y.trunc)+n*eta1/(2*eta2)+(1/sqrt(-2*eta2))*sum(h1-h2)
    SF[2,]=sum(y.trunc^2)-n*eta1^2/(4*eta2^2)+n/(2*eta2)-sum(h1*cr1-h2*cr2)
    return(SF)
  }
  
  #Hessian matrix
  J=matrix(,2,2)  
  Hessian_func=function(ETA){
    eta1=ETA[1];eta2=ETA[2]
    h1=H1(eta1,eta2); h2=H2(eta1,eta2)
    cr1=CR1(eta1,eta2); cr2=CR2(eta1,eta2)
    cr3=CR3(eta1,eta2); cr4=CR4(eta1,eta2)
    cr5=CR5(eta1,eta2); cr6=CR6(eta1,eta2)
    v1=cr3*(cr1^2)-cr5; v2=cr4*(cr2^2)-cr6
    v3=h1*cr1-h2*cr2; v4=-h1*cr3*cr1+h2*cr4*cr2
    J[1,1]=n/(2*eta2)-(1/(2*eta2))*sum(h1*cr3-h2*cr4+(h1-h2)^2)
    J[2,2]=n*eta1^2/(2*eta2^3)-n/(2*eta2^2)+sum(h1*v1-h2*v2+v3^2)
    J[1,2]=-n*eta1/(2*eta2^2)+((-2*eta2)^(-3/2))*sum(h1-h2)+(1/sqrt(-2*eta2))*sum(v4-v3*(h1-h2))
    J[2,1]=J[1,2]
    return(J)
  }
  
  #Newton-Raphson method
  Eta=matrix(,1,2)
  Eta_old=Eta_new=matrix(,1,2)
  Eta[1,]=c(mean(y.trunc)/var(y.trunc),-1/(2*var(y.trunc)))
  
  k=1
  repeat{
    Eta_old=Eta[1:k,]
    Eta_new=Eta[k,]-solve(Hessian_func(Eta[k,]))%*%score_func(Eta[k,])
    Eta=matrix(,k+1,2)
    Eta[1:k,]=Eta_old
    Eta[k+1,]=c(Eta_new[1,],Eta_new[2,])
    Err1=abs(Eta[k+1,1]-Eta[k,1])
    Err2=abs(Eta[k+1,2]-Eta[k,2])
    if( (Err1<epsilon)&(Err2<epsilon) ){break}
    k=k+1
  }

  #estimator eta1,eta2,mu and sigma
  Eta_hat=Eta[k+1,]
  V=-solve(Hessian_func(Eta[k,]))
  SE1=sqrt(V[1,1])
  SE2=sqrt(V[2,2])
  mu_hat=-Eta_hat[1]/(2*Eta_hat[2])
  sigma_hat=sqrt(-1/(2*Eta_hat[2]))
  
  #log-likelihood function
  logL_func=function(eta1,eta2){
    a=pnorm((v.trunc+eta1/(2*eta2))/sqrt(-1/(2*eta2)))
    b=pnorm((u.trunc+eta1/(2*eta2))/sqrt(-1/(2*eta2)))
    c=sum(eta1*y.trunc+eta2*y.trunc^2)+(n*eta1^2)/(4*eta2)+(n/2)*log(-eta2)-(n/2)*log(pi)-sum(log(a-b))
    return(c)  
  }
  
  #value of log-likelihood
  eta1_hat=Eta_hat[1]
  eta2_hat=Eta_hat[2]
  logL=logL_func(eta1_hat,eta2_hat)
  
  #AIC
  p=2
  AIC=-2*logL+2*p
  convergence_res=c(logL=logL,DF=p,AIC=AIC,No.of.iterations=k)
  
  list(eta=c(eta1=eta1_hat,eta2=eta2_hat),
       SE=c(eta1=SE1,eta2=SE2),
       normal=c(mean=mu_hat,SD=sigma_hat),
       convergence=convergence_res,
       Score=score_func(Eta_hat),
       Hessian=Hessian_func(Eta_hat)
  )

}
