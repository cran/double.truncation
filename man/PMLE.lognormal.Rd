\name{PMLE.lognormal}
\alias{PMLE.lognormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Inference for the lognormal model
}
\description{
Maximum likelihood estimates (MLEs) and their standard errors (SEs) are computed for the lognormal model
 based on doubly-truncated data (Dorre et al. 2021). Also computed are the likelihood value, AIC, and other qnantities.
}
\usage{
PMLE.lognormal(u.trunc, y.trunc, v.trunc,epsilon = 1e-5,D1=2,D2=2,d1=2,d2=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u.trunc}{a vector of lower truncation limits}
  \item{y.trunc}{a vector of variables of interest}
  \item{v.trunc}{a vector of upper truncation limits}
  \item{epsilon}{a small positive number for the error tolerance for Newton-Raphson iterations}
  \item{D1}{a positive number: Randomize the intial value for a divergent iteration (the updated amount for mu is greater than D1)}
  \item{D2}{a positive number: Randomize the intial value for a divergent iteration (the updated amount for sigma is greater than D2)}
  \item{d1}{a positive number: For a divergent iteration, U(-d1,d1) is added to the intial value of mu}
  \item{d2}{a positive number: For a divergent iteration, U(-d2,d2) is added to the intial value of log(sigma)}
}
\details{
A randomized Newton–Raphson algorithm (Section 3.2 of Dorre et al.(2021)) was employed to compute the MLE.
}
\value{
  \item{eta}{estimates}
  \item{SE}{standard errors}
  \item{convergence}{Log-likelihood, degree of freedom, AIC, the number of iterations}
  \item{Score}{score vector at the converged value}
  \item{Hessian}{Hessian matrix at the converged value}
}
\references{
Dorre A, Huang CY, Tseng YK, Emura T (2021) Likelihood-based analysis of
 doubly-truncated data under the location-scale and AFT model, Computation Stat 36(1): 375-408 
}
\author{Takeshi Emura}

\examples{
## A data example from Efron and Petrosian (1999) ## 
y.trunc=c(0.75, 1.25, 1.50, 1.05, 2.40, 2.50, 2.25)
u.trunc=c(0.4, 0.8, 0.0, 0.3, 1.1, 2.3, 1.3)
v.trunc=c(2.0, 1.8, 2.3, 1.4, 3.0, 3.4, 2.6)
PMLE.lognormal(u.trunc,y.trunc,v.trunc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Location-scale family }% use one of  RShowDoc("KEYWORDS")
\keyword{ Weibull distribution }% __ONLY ONE__ keyword per line
