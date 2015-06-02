\name{eval.greedy}
\alias{eval.greedy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns a greedy histogram }
\description{
Returns a greedy histogram.
A greedy histogram is grown stagewise by minimizing an empirical risk
functional.
}
\usage{
eval.greedy(dendat, leaf, method = "loglik", minobs = NULL, bound = 0, 
suppo = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendat}{ n*d data matrix }
  \item{leaf}{ the (maximal) number of rectangles in the partition of 
   the histogram }
  \item{method}{ "loglik" or "projec"; the empirical risk is either the 
   log-likelihood or the L2 empirical risk }
  \item{minobs}{ non-negative integer; splitting of a bin will be continued if
   the bin containes "minobs" or more observations  }
  \item{bound}{ internal }
  \item{suppo}{2*d vector of real numbers;
  the rectangle to be splitted;
  the rectangle has to contain the data}
}
%\details{}
\value{
 An evaluation tree
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemela }
%\note{ ~~further notes~~ }
\seealso{ 
\code{\link{lstseq.greedy}},
\code{\link{partition}}
}

\examples{
library(denpro)
dendat<-sim.data(n=200,seed=5,type="mulmodII")
eva<-eval.greedy(dendat,leaf=15)

dp<-draw.pcf(eva,pnum=c(60,60))         
persp(dp$x,dp$y,dp$z,theta=-20,phi=30)

}

\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
