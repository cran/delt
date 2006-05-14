\name{eval.stage}
\alias{eval.stage}
\title{ Returns a stagewise minimization estimate }
\description{
Returns a stagewise minimization estimate.
A stagewise minimization estimator is a convex combination of
greedy histograms.
The convex combination is constructed by a stagewise minimization of
an empirical risk functional.  
}
\usage{
eval.stage(dendat, leaf, M, pis = NULL, mcn = dim(dendat)[1], 
minobs = NULL, seedi = 1, method = "projec", bound = 0)
}
\arguments{
  \item{dendat}{ n*d data matrix }
  \item{leaf}{ the (maximal) number of rectangles in the partition of
   the greedy histograms }
  \item{M}{ the number of histograms in the convex combination }
  \item{pis}{ the vector of weights of the convex combination }
  \item{mcn}{ the size of the Monte Carlo sample used in the numerical
   integration in calculating the empirical risk functional }
  \item{minobs}{ non-negative integer; splitting of a bin of a greedy
   histogram will be continued if the bin containes "minobs" 
   or more observations }
  \item{seedi}{ the seed for the generation of the Monte Carlo sample }
  \item{method}{ "loglik" or "projec"; the empirical risk is either the
   log-likelihood or the L2 empirical risk }
  \item{bound}{ internal }
}
%\details{}
\value{
An evaluation tree
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemelä }
%\note{ ~~further notes~~ }
\seealso{ 
\code{\link{eval.greedy}},
\code{\link{eval.stage.gauss}} 
}
\examples{
dendat<-sim.data(n=200,seed=5,type="mulmodII")
leaf<-13  # the number of leafs of the greedy histograms
M<-5      # the number of greedy histograms

pcf<-eval.stage(dendat,leaf=leaf,M=M)

dp<-draw.pcf(pcf,pnum=c(120,120))
persp(dp$x,dp$y,dp$z,ticktype="detailed",phi=25,theta=-120)

}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
