simu.Weibull <- function(n,mu,sigma,delta){

u=y=v=a=c()

i=1
repeat{
  y[i]=mu+sigma*log(-log(1-runif(1)))
  u[i]=mu-delta+sigma*log(-log(1-runif(1)))
  v[i]=mu+delta+sigma*log(-log(1-runif(1)))	
  if(u[i]<=y[i]&&y[i]<=v[i]){		
    a[i]=1
  }else{	
    a[i]=0
  }
  if(sum(a)==n){break}
  i=i+1
}

Inclu=mean(a)
print(c(Inclusion_probability=Inclu))
# truncated data 
u.trunc=u[a==1]
y.trunc=y[a==1]
v.trunc=v[a==1]
Dat=data.frame(u=u.trunc,y=y.trunc,v=v.trunc)
Dat

}




