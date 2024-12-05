GoF <- function(u.trunc,y.trunc,v.trunc,epsilon=1e-08,F0,B=500,F.plot=TRUE){
  n=length(y.trunc)
  res=NPMLE(u.trunc,y.trunc,v.trunc,epsilon=epsilon,detail=TRUE)
  F.est=res$F
  f.est=res$f
  V.mat=res$V
  T.mat=matrix(rep(sort(y.trunc)[-n],n-1),n-1,n-1,byrow=TRUE)
  W.mat=T.mat<=t(T.mat)
  C.est=sum(  (F.est[-n]-F0(sort(y.trunc))[-n])^2  )
  K.est=max(abs(F.est-F0(sort(y.trunc))))

  #Tau1=0.2
  #Tau2=0.8
  #temp=F.est[-n]
  #phi=sqrt( pmax(F.est[-n],Tau1)*(1-pmin(F.est[-n],Tau2)) )

  C.boot=K.boot=Z.boot=rep(0,B)
  for(i in 1:B){
    G=mvrnorm(n=1,mu=rep(0,n-1),Sigma=V.mat)
    #Z=G/phi
    C.boot[i]=sum(G^2)
    K.boot[i]=max(abs(G))
    #Z.boot[i]=max(abs(Z))
  }

  if(F.plot==TRUE){
    plot(F0(sort(y.trunc)),F.est,xlab="F0",ylab="NPMLE")
    abline(a=0,b=1,col="red")
  }

  C.res=c(C=C.est,P=mean(C.boot>C.est),
          boot90=sort(C.boot)[B*0.90],
          boot95=sort(C.boot)[B*0.95],
          boot99=sort(C.boot)[B*0.99])
  K.res=c(K=K.est,P=mean(K.boot>K.est),
          boot90=sort(K.boot)[B*0.90],
          boot95=sort(K.boot)[B*0.95],
          boot99=sort(K.boot)[B*0.99])
  list(CvM=C.res,KS=K.res)

}
