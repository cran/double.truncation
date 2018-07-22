PMLE.SEF3.positive <-
function(u.trunc,y.trunc,v.trunc,
                        tau2=max(y.trunc),epsilon=0.0001){
  
  n=length(y.trunc)
  u.max=u.trunc
  v.min=pmin(v.trunc,tau2)

  myfun0=function(eta1,eta2,eta3){
    int=c()
    for(i in 1:n){
      myfun00=function(y){ exp(eta1*y+eta2*y^2+eta3*y^3) }
      int[i]=integrate(myfun00,u.max[i],v.min[i])$value
    }
    return(int)   
  }
  
  myfun1=function(eta1,eta2,eta3){
    int1=c()
    for(i in 1:n){
      myfun11=function(y){ y*exp(eta1*y+eta2*y^2+eta3*y^3) }
      int1[i]=integrate(myfun11,u.max[i],v.min[i])$value
    }
    return(int1)   
  }
  
  myfun2=function(eta1,eta2,eta3){
    int2=c()
    for(i in 1:n){
      myfun22=function(y){  (y^2)*exp(eta1*y+eta2*y^2+eta3*y^3) }
      int2[i]=integrate(myfun22,u.max[i],v.min[i])$value
    }
    return(int2)   
  }
  
  myfun3=function(eta1,eta2,eta3){
    int3=c()
    for(i in 1:n){
      myfun33=function(y){ (y^3)*exp(eta1*y+eta2*y^2+eta3*y^3) }
      int3[i]=integrate(myfun33,u.max[i],v.min[i])$value
    }
    return(int3)   
  }
  
  myfun4=function(eta1,eta2,eta3){
    int4=c()
    for(i in 1:n){
      myfun44=function(y){ (y^4)*exp(eta1*y+eta2*y^2+eta3*y^3)  }
      int4[i]=integrate(myfun44,u.max[i],v.min[i])$value
    }
    return(int4)   
  } 
  
  myfun5=function(eta1,eta2,eta3){
    int5=c()
    for(i in 1:n){
      myfun55=function(y){ (y^5)*exp(eta1*y+eta2*y^2+eta3*y^3) }
      int5[i]=integrate(myfun55,u.max[i],v.min[i])$value
    }
    return(int5)   
  }
  
  myfun6=function(eta1,eta2,eta3){
    int6=c()
    for(i in 1:n){
      myfun66=function(y){ (y^6)*exp(eta1*y+eta2*y^2+eta3*y^3) }
      int6[i]=integrate(myfun66,u.max[i],v.min[i])$value
    }
    return(int6)   
  }
  
  score_func=function(ETA){
SF=matrix(,3,1)
    sf0=myfun0(ETA[1],ETA[2],ETA[3])
    sf1=myfun1(ETA[1],ETA[2],ETA[3])
    sf2=myfun2(ETA[1],ETA[2],ETA[3])
    sf3=myfun3(ETA[1],ETA[2],ETA[3])
    SF[1,]=sum(y.trunc)-sum(sf1/sf0)
    SF[2,]=sum(y.trunc^2)-sum(sf2/sf0)
    SF[3,]=sum(y.trunc^3)-sum(sf3/sf0)
    return(SF)
  }
  
  
  Hessian_func=function(ETA){
    J=matrix(,3,3)  
    sf0=myfun0(ETA[1],ETA[2],ETA[3])
    sf1=myfun1(ETA[1],ETA[2],ETA[3])
    sf2=myfun2(ETA[1],ETA[2],ETA[3])
    sf3=myfun3(ETA[1],ETA[2],ETA[3])
    sf4=myfun4(ETA[1],ETA[2],ETA[3])
    sf5=myfun5(ETA[1],ETA[2],ETA[3])
    sf6=myfun6(ETA[1],ETA[2],ETA[3])
    J[1,1]=-sum(sf2/sf0)+sum((sf1/sf0)^2)
    J[2,2]=-sum(sf4/sf0)+sum((sf2/sf0)^2)
    J[3,3]=-sum(sf6/sf0)+sum((sf3/sf0)^2)
    J[1,2]=-sum(sf3/sf0)+sum((sf2/sf0)*(sf1/sf0))
    J[2,1]=J[1,2]
    J[1,3]=-sum(sf4/sf0)+sum((sf3/sf0)*(sf1/sf0))
    J[3,1]=J[1,3]
    J[2,3]=-sum(sf5/sf0)+sum((sf3/sf0)*(sf2/sf0))
    J[3,2]=J[2,3]
    return(J)
  }
  
  #Newton-Raphson
  Eta=matrix(,1,3)
  Eta_old=matrix(,1,3)
  Eta_new=matrix(,1,3)
  Eta[1,]=c(mean(y.trunc)/var(y.trunc),-1/(2*var(y.trunc)),0)

  k=1
  repeat{
    Eta_old=Eta[1:k,]
    Eta_new=Eta[k,]-solve(Hessian_func(Eta[k,]),tol=10^(-18))%*%score_func(Eta[k,])
    Eta=matrix(,k+1,3)
    Eta[1:k,]=Eta_old
    Eta[k+1,]=c(Eta_new[1,],Eta_new[2,],Eta_new[3,])
    Error1=abs(Eta[k+1,1]-Eta[k,1])
    Error2=abs(Eta[k+1,2]-Eta[k,2])
    Error3=abs(Eta[k+1,3]-Eta[k,3])
    if( (Error1<0.0001)&(Error2<0.0001)&(Error3<10^(-7)) ){break
    }else if( (Error1>20)|(Error2>10)|(Error3>1) ){
      Eta[1,]=c(mean(y.trunc)/var(y.trunc)+runif(1,-6,6),-1/(2*var(y.trunc))+runif(1,-0.5,0.5),0)
      k=0
    }
    k=k+1
  }
  
  Eta_hat=Eta[k+1,]
  
  se_eta1=sqrt(solve( -Hessian_func(Eta_hat),tol=10^(-18) )[1,1])
  se_eta2=sqrt(solve( -Hessian_func(Eta_hat),tol=10^(-18) )[2,2])
  se_eta3=sqrt( solve(-Hessian_func(Eta_hat),tol=10^(-18) )[3,3])
  
  #log-likelihood function
  logL_func=function(eta1,eta2,eta3){
    int=c()
    for(i in 1:n){
      myfun=function(y){
        exp(eta1*y+eta2*y^2+eta3*y^3)
      }
      int[i]=integrate(myfun,u.max[i],v.min[i])$value
    }
    sum(eta1*y.trunc+eta2*y.trunc^2+eta3*y.trunc^3)-sum(log(int))
  }
    
  logL=logL_func(Eta_hat[1],Eta_hat[2],Eta_hat[3])
  #AIC
  p=3
  AIC=-2*logL+2*p
  convergence_res=c(logL=logL,DF=p,AIC=AIC,No.of.iterations=k)
  
  list(eta=c(eta1=Eta_hat[1],eta2=Eta_hat[2],eta3=Eta_hat[3]),
       SE=c(eta1=se_eta1,eta2=se_eta2,eta3=se_eta3),
       convergence=convergence_res,
       Score=as.vector(score_func(Eta_hat)),
       Hessian=Hessian_func(Eta_hat)
  )
  
}
