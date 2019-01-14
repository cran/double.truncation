\name{PMLE.SEF2.negative}
\alias{PMLE.SEF2.negative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Inference for the two-parameter SEF model
 (negative parameter space for eta_2)
}
\description{
Maximum likelihood estimates and their standard errors (SEs) are computed.
 Also computed are the likelihood value, AIC, and other qnantities.
 Since this is the model, estimates for the mean and SD are also computed.
}
\usage{
PMLE.SEF2.negative(u.trunc, y.trunc, v.trunc, epsilon = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u.trunc}{lower truncation limit}
  \item{y.trunc}{variable of interest}
  \item{v.trunc}{upper truncation limit}
  \item{epsilon}{error tolerance for Newton-Raphson}
}
\details{
Details are seen from the references.
}
\value{
  \item{eta}{estimates}
  \item{SE}{standard errors}
  \item{convergence}{Log-likelihood, degree of freedom, AIC, the number of iterations}
  \item{Score}{score vector at the converged value}
  \item{Hessian}{Hessian matrix at the converged value}
}
\references{
Hu YH, Emura T (2015) Maximum likelihood estimation for a special exponential family
 under random double-truncation, Computation Stat 30 (4): 1199-229

Emura T, Hu YH, Konno Y (2017) Asymptotic inference for maximum likelihood estimators
 under the special exponential family with double-truncation, Stat Pap 58 (3): 877-909

Dorre A, Emura T (2019) Analysis of Doubly Truncated Data, An Introduction,
 JSS Research Series in Statistics, Springer
}
\author{Takeshi Emura, Ya-Hsuan Hu}

\examples{
### Data generation: see Appendix of Hu and Emura (2015)
n=300
eta1_true=30
eta2_true=-0.5
mu_true=30
mu_u=29.09
mu_v=30.91

a=u=v=y=c()

###generate n samples of (ui,yi,vi) subject to ui<=yi<=vi###
j=1
repeat{  
  u[j]=rnorm(1,mu_u,1)
  v[j]=rnorm(1,mu_v,1)
  y[j]=rnorm(1,mu_true,1)
  if(u[j]<=y[j]&&y[j]<=v[j]) a[j]=1 else a[j]=0
  if(sum(a)==n) break ###we need n data set###
  j=j+1
}
mean(a) ### inclusion probability around 0.5 ###

v.trunc=v[a==1]
y.trunc=y[a==1]
u.trunc=u[a==1]
  
PMLE.SEF2.negative(u.trunc,y.trunc,v.trunc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Special exponential family }% use one of  RShowDoc("KEYWORDS")
\keyword{ Normal distribution }% __ONLY ONE__ keyword per line
