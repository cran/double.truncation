\name{NPMLE}
\alias{NPMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Nonparametric inference based on the self-consistency method}
\description{
Nonparametric maximum likelihood estimates are computed
 based on the self-consistency method (Efron and Petrosian 1999). The SE is computed from the asymptotic variance
 derived in Emura et al. (2015).
}
\usage{
NPMLE(u.trunc, y.trunc, v.trunc,epsilon=1e-08,detail=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u.trunc}{lower truncation limit}
  \item{y.trunc}{variable of interest}
  \item{v.trunc}{upper truncation limit}
  \item{epsilon}{error tolerance for the self-consistency algorithm}
  \item{detail}{if TRUE, show the details including the covariate matrix}
}
\details{
Details are seen from the references.
}

\value{
  \item{f}{density}
  \item{F}{cumulative distribution}
  \item{SE}{standard error} 
  \item{convergence}{Log-likelihood, and the number of iterations}
  \item{V}{covariance matrix for the NPMLE} 
}
\references{
Efron B, Petrosian V (1999). Nonparametric methods for doubly truncated data. 
J Am Stat Assoc 94: 824-834

Emura T, Konno Y, Michimae H (2015). Statistical inference based on the 
nonparametric maximum likelihood estimator under double-truncation. 
Lifetime Data Analysis 21: 397-418

Dorre A, Emura T (2019) Analysis of Doubly Truncated Data, An Introduction,
 JSS Research Series in Statistics, Springer
}
\author{Takeshi Emura}

\examples{
## A data example from Efron and Petrosian (1999) ## 
y.trunc=c(0.75, 1.25, 1.50, 1.05, 2.40, 2.50, 2.25)
u.trunc=c(0.4, 0.8, 0.0, 0.3, 1.1, 2.3, 1.3)
v.trunc=c(2.0, 1.8, 2.3, 1.4, 3.0, 3.4, 2.6)
NPMLE(u.trunc,y.trunc,v.trunc)
NPMLE(u.trunc,y.trunc,v.trunc,detail=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ NPMLE }% use one of  RShowDoc("KEYWORDS")
\keyword{ Self-consistency method }% __ONLY ONE__ keyword per line
