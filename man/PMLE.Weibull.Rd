\name{PMLE.Weibull}
\alias{PMLE.Weibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Inference for the Weibull model
}
\description{
Maximum likelihood estimates and their standard errors (SEs) are computed.
 Also computed are the likelihood value, AIC, and other qnantities.
}
\usage{
PMLE.Weibull(u.trunc, y.trunc, v.trunc,epsilon = 1e-5,D1=2,D2=2,d1=2,d2=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u.trunc}{lower truncation limit}
  \item{y.trunc}{variable of interest}
  \item{v.trunc}{upper truncation limit}
  \item{epsilon}{error tolerance for Newton-Raphson}
  \item{D1}{Randomize the intial value if |mu_h-mu_{h+1}|>D1}
  \item{D2}{Randomize the intial value if |sigma_h-sigma_{h+1}|>D2}
  \item{d1}{U(-d1,d1) is added to the intial value of mu}
  \item{d2}{U(-d2,d2) is added to the intial value of sigma}
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
Dorre A, Huang CY, Tseng YK, Emura T (2020-) Likelihood-based analysis of
 doubly-truncated data under the location-scale and AFT model, Computation Stat 
}
\author{Takeshi Emura}

\examples{
## A data example from Efron and Petrosian (1999) ## 
y.trunc=c(0.75, 1.25, 1.50, 1.05, 2.40, 2.50, 2.25)
u.trunc=c(0.4, 0.8, 0.0, 0.3, 1.1, 2.3, 1.3)
v.trunc=c(2.0, 1.8, 2.3, 1.4, 3.0, 3.4, 2.6)
PMLE.Weibull(u.trunc,y.trunc,v.trunc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Location-scale family }% use one of  RShowDoc("KEYWORDS")
\keyword{ Weibull distribution }% __ONLY ONE__ keyword per line
