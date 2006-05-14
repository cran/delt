\name{eval.cart}
\alias{eval.cart}
\title{ Calculates a CART histogram }
\description{
Calculates a CART histogram.
The estimate is represented as an evaluation tree.
An CART histogram is a multivariate adaptive histogram
which is obtained by pruning an evaluation tree of an overfitting
histogram. 
}
\usage{
eval.cart(dendat, leaf, minobs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendat}{ n*d data matrix }
  \item{leaf}{ positive integer; the cardinality of the partition of the 
   histogram }
  \item{minobs}{ non-negative integer; splitting of a bin of the overfitting
   histogram will be continued if
   the bin containes "minobs" or more observations }
}
\details{
The partition of the histogram may not contain exactly "leaf" 
rectangles: the cardinality of the partition is as close as possible
to "leaf"
 }
\value{
An evaluation tree
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemelä }
%\note{ ~~further notes~~ }
\seealso{ 
\code{\link{lstseq.cart}},
\code{\link{densplit}} 
}
\examples{
dendat<-sim.data(n=600,seed=5,type="mulmodII")
eva<-eval.cart(dendat,16)

dp<-draw.pcf(eva,pnum=c(60,60))
persp(dp$x,dp$y,dp$z,theta=-20,phi=30)

}

\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
