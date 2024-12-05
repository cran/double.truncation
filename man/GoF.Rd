\name{GoF}
\alias{GoF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Goodness-of-fit test based on the CvM and KS statistics}
\description{
Goodness-of-fit test statistics are computed based on the Cramér–von Mises (CvM) and Kolmogorov–Smirnov (KS) test statistics proposed in Emura et al. (2015). P-value and critical values with significance levels of 0.01, 0.05 and 0.10 are also computed.
}
\usage{
GoF(u.trunc, y.trunc, v.trunc,epsilon=1e-08,F0,B=500,F.plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u.trunc}{lower truncation limit}
  \item{y.trunc}{variable of interest}
  \item{v.trunc}{upper truncation limit}
  \item{epsilon}{error tolerance for the self-consistency algorithm}
  \item{F0}{a function for the null distribution function}
  \item{B}{the number of bootstrap resamples (B=500 is the default)}
  \item{F.plot}{model diagnostic plot}

}
\details{
Details are seen from Emura et al.(2015).
}

\value{
  \item{CvM}{Test statistics, P-value, and critical values for the Cramér–von Mises (CvM) test}
  \item{KS}{Test statistics, P-value, and critical values for the Kolmogorov–Smirnov (KS) test}
}
\references{
Emura T, Konno Y, Michimae H (2015). Statistical inference based on the 
nonparametric maximum likelihood estimator under double-truncation. 
Lifetime Data Analysis 21: 397-418
}
\author{Takeshi Emura}

\examples{
## A data example from Efron and Petrosian (1999) ## 
y.trunc=c(0.75, 1.25, 1.50, 1.05, 2.40, 2.50, 2.25)
u.trunc=c(0.4, 0.8, 0.0, 0.3, 1.1, 2.3, 1.3)
v.trunc=c(2.0, 1.8, 2.3, 1.4, 3.0, 3.4, 2.6)
F0=function(x){x/3}
GoF(u.trunc,y.trunc,v.trunc,F0=F0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Goodness-of-fit }% use one of  RShowDoc("KEYWORDS")
\keyword{ Model diagnostic }% __ONLY ONE__ keyword per line
