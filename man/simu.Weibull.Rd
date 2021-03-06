\name{simu.Weibull}
\alias{simu.Weibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating doubly-truncated data from the Weibull model
}
\description{
A data frame is generated by simulated data from the Weibull model.
}
\usage{
simu.Weibull(n,mu,sigma,delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size}
  \item{mu}{location parameter}
  \item{sigma}{scale parameter}
  \item{delta}{a positive parameter controlling the inclusion probability}
}
\details{
The data are generated from the random vector (U,Y,V) subject to the inclusion criterion U<=Y<=V.
The random vector are defined as U=mu-delta+sigma*W, Y=mu+sigma*W, and U=mu+delta+sigma*W, where P(W>w)=exp(-exp(w)).
See Section 5.1 of Dorre et al. (2020-) for details. The inclusion probability is P(U<=Y<=V).
}
\value{
  \item{u}{lower truncation limits}
  \item{y}{log-transformed lifetimes}
  \item{v}{upper truncation limits}
}
\references{
Dorre A, Huang CY, Tseng YK, Emura T (2020-) Likelihood-based analysis of
 doubly-truncated data under the location-scale and AFT model, Computation Stat, DOI:10.1007/s00180-020-01027-6
}
\author{Takeshi Emura}

\examples{
## A simulation from Dorre et al.(2020) ##
simu.Weibull(n=100,mu=5,sigma=2,delta=2.08)

Dat=simu.Weibull(n=100,mu=5,sigma=2,delta=2.08)
PMLE.Weibull(Dat$u,Dat$y,Dat$v)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Location-scale family }% use one of  RShowDoc("KEYWORDS")
\keyword{ Weibull distribution }% __ONLY ONE__ keyword per line
