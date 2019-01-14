\name{PMLE.SEF3.negative}
\alias{PMLE.SEF3.negative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Inference for the three-parameter SEF model
 (negative parameter space for eta_3)
}
\description{
Maximum likelihood estimates and their standard errors (SEs) are computed.
 Also computed are the likelihood value, AIC, and other qnantities.
}
\usage{
PMLE.SEF3.negative(u.trunc, y.trunc, v.trunc, tau1 = min(y.trunc),
 epsilon = 1e-04, D1=20, D2=10, D3=1, d1=6, d2=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u.trunc}{lower truncation limit}
  \item{y.trunc}{variable of interest}
  \item{v.trunc}{upper truncation limit}
  \item{tau1}{lower support}
  \item{epsilon}{error tolerance for Newton-Raphson}
  \item{D1}{Divergence condition for eta_1}
  \item{D2}{Divergence condition of eta_2}
  \item{D3}{Divergence condition of eta_3}
  \item{d1}{Range of randomization for eta_1}
  \item{d2}{Range of randomization for eta_2}
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
## The first 10 samples of the childhood cancer data ##
y.trunc=c(6,7,15,43,85,92,96,104,108,123)
u.trunc=c(-1643,-24,-532,-1508,-691,-1235,-786,-261,-108,-120)
v.trunc=u.trunc+1825
PMLE.SEF3.negative(u.trunc,y.trunc,v.trunc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Special exponential family }% use one of  RShowDoc("KEYWORDS")
\keyword{ Skew normal distribution }% __ONLY ONE__ keyword per line
