\name{PMLE.SEF1.free}
\alias{PMLE.SEF1.free}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric inference for the one-parameter SEF model (free parameter space)}
\description{
Maximum likelihood estimates and their standard errors (SEs) are computed.
 Also computed are the likelihood value, AIC, and other qnantities.
}
\usage{
PMLE.SEF1.free(u.trunc, y.trunc, v.trunc,
 tau1 = min(y.trunc), tau2 = max(y.trunc), epsilon = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u.trunc}{lower truncation limit}
  \item{y.trunc}{variable of interest}
  \item{v.trunc}{upper truncation limit}
  \item{tau1}{lower support}
  \item{tau2}{upper support}
  \item{epsilon}{error tolerance for Newton-Raphson}
}
\details{
Details are seen from the references.
}

\value{
  \item{eta}{estimates}
  \item{SE}{standard errors}
  \item{convergence}{Log-likelihood, degree of freedom, AIC, the number of iterations}
  \item{Score}{score at the converged value}
  \item{Hessian}{Hessian at the converged value}
}
\references{
Hu YH, Emura T (2015) Maximum likelihood estimation for a special exponential family
 under random double-truncation, Computation Stat 30 (4): 1199-229

Emura T, Hu YH, Konno Y (2017) Asymptotic inference for maximum likelihood estimators
 under the special exponential family with double-truncation, Stat Pap 58 (3): 877-909

Dorre A, Emura T (2019) Analysis of Doubly Truncated Data, An Introduction,
 JSS Research Series in Statistics, Springer
}
\author{Takeshi Emura}

\examples{
### Data generation: see Appendix of Hu and Emura (2015) ###
eta_true=-3
eta_u=-9
eta_v=-1
tau=10
n=300

a=u=v=y=c()

j=1
repeat{
  u1=runif(1,0,1)
  u[j]=tau+(1/eta_u)*log(1-u1)
  u2=runif(1,0,1)
  v[j]=tau+(1/eta_v)*log(1-u2)
  u3=runif(1,0,1)
  y[j]=tau+(1/eta_true)*log(1-u3)
  if(u[j]<=y[j]&&y[j]<=v[j]) a[j]=1 else a[j]=0
  if(sum(a)==n) break
  j=j+1
}
mean(a) ## inclusion probability around 0.5
  
v.trunc=v[a==1]
u.trunc=u[a==1]
y.trunc=y[a==1]

PMLE.SEF1.free(u.trunc,y.trunc,v.trunc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Special exponential family }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
