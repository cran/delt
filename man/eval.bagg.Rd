\name{eval.bagg}
\alias{eval.bagg}
\title{ Returns a bootstrap aggregation of adaptive histograms }
\description{
Returns a bootstrap aggregation of CART-histograms or greedy histograms.
}
\usage{
eval.bagg(dendat, B, leaf, minobs = NULL, seed = 1, sample = "bagg", 
prune = "off", splitscan = 0, seedf = 1, scatter = 0, src = "c", 
method = "loglik")
}
\arguments{
  \item{dendat}{ n*d data matrix }
  \item{B}{ positive integer; the number of aggregated histograms }
  \item{leaf}{ the cardinality of the partitions of the aggregated histograms }
  \item{minobs}{ non-negative integer; 
   a property of aggregated histograms;
   splitting of a bin will be continued if
   the bin containes "minobs" or more observations }
  \item{seed}{ the seed for the random number generation of the
   random selection of the bootstrap sample }
  \item{sample}{ "bagg" or "worpl";
   the bootstrapping method;
   "worpl" for the n/2-out-of-n without replacement;
   "bagg" for n-out-of-n with replacement }
  \item{prune}{ "on" or "off";
   if "on", then CART-histograms will be aggregated;
   if "off", then greedy histograms will be aggregated }
  \item{splitscan}{ internal 
   (how many splits will be used for random split selection) }
  \item{seedf}{ internal (seed for random split selection) }
  \item{scatter}{ internal (random perturbation of observations) }
  \item{src}{ internal ("c" or "R" code) }
  \item{method}{ "loglik" or "projec";
   the empirical risk is either the log-likelihood or the L2 empirical risk }
}
%\details{}
\value{
 An evaluation tree
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemela }
%\note{ ~~further notes~~ }
\seealso{ 
\code{\link{lstseq.bagg}},
\code{\link{eval.cart}},
\code{\link{eval.greedy}}
}
\examples{
library(denpro)
dendat<-sim.data(n=600,seed=5,type="mulmodII")

leaf<-7     # number of leaves in the histograms
seed<-1     # seed for choosing bootstrap samples
sample="worpl" # without-replacement bootstrap
prune="on"  # we use CART-histograms
B<-5        # the number of histograms in the average

eva<-eval.bagg(dendat,B,leaf,seed=seed,sample=sample,prune=prune)

dp<-draw.pcf(eva,pnum=c(60,60))
persp(dp$x,dp$y,dp$z,theta=-20,phi=30)

}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
