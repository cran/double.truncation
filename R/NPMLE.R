NPMLE <- function(u.trunc,y.trunc,v.trunc,epsilon=1e-08,detail=FALSE){
  n=length(y.trunc)
  u.trunc=u.trunc[order(y.trunc)]
  v.trunc=v.trunc[order(y.trunc)]
  y.trunc=y.trunc[order(y.trunc)]

  atr_U=matrix(y.trunc,n,n,byrow=TRUE)>=matrix(u.trunc,n,n)
  atr_V=matrix(y.trunc,n,n,byrow=TRUE)<=matrix(v.trunc,n,n)
  J_mat=atr_U&atr_V

  f_old=rep(1/n,n)
  k=0 ## the number of iterations
  repeat{
    k=k+1
    F_old=J_mat%*%f_old
    f_new=1/(  t(J_mat)%*%(1/F_old)  )
    f_new=f_new/sum(f_new)
    f_new=as.vector(f_new)
    Error=sum( abs(f_new-f_old) )
    if(Error<epsilon){break}
    f_old=f_new
  }
  f_est=f_new
  F_est=as.vector(J_mat%*%f_est) ## This is not CDF ##

  ##### Score & Fisher infomation #####
  D=cbind(diag(n-1),-rep(1,n-1))
  Score=D%*%(  1/(f_est)-t(J_mat)%*%(1/(J_mat%*%f_est))  )
  Info=D%*%(diag(1/f_est^2)-t(J_mat)%*%diag(1/F_est^2)%*%J_mat)%*%t(D)
  W=upper.tri(matrix(1,n-1,n-1), diag=TRUE)
  V.mat=t(W)%*%solve(Info)%*%W
  SE=c(sqrt( diag(V.mat) ),0)

  logL=sum(log(f_est)-log(F_est))
  convergence_res=c(logL=logL, No.of.iterations=k)

  if(detail==FALSE){
    list(f=f_new,F=cumsum(f_new),SE=SE,convergence=convergence_res)
  }else{
    list(f=f_new,F=cumsum(f_new),SE=SE,convergence=convergence_res,
         V=V.mat)
  }

}
