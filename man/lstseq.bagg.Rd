\name{lstseq.bagg}
\alias{lstseq.bagg}
\title{ Calculates a scale of bootstrap aggregated histograms }
\description{
Calculates a scale of bootstrap aggregated histograms.
The estimates in the sequence are calculated with function "eval.bagg".
}
\usage{
lstseq.bagg(dendat, B, lstree=NULL, level = NULL, 
maxleaf = NULL, leafseq = NULL, 
minobs = NULL, seed = 1, sample = "bagg", prune = "off", 
splitscan = 0, seedf = 1, scatter = 0, src = "c", method = "loglik")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{dendat}{ n*d data matrix }
  \item{B}{ positive integer; the number of aggregated histograms }
  \item{maxleaf}{ the maximal cardinality of the partitions of the histograms
   in the sequence } 
  \item{lstree }{ if NULL, then level set trees are not calculated }
  \item{level}{ if NULL, then shape trees are not calculated;
   if positive number, then it is
   the level of the level sets for which the shape trees are calculated }
  \item{leafseq}{ a vector giving the cardinalities of the partitions 
   of the aggregated histograms }
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
A list with components
\item{lstseq}{a list of level set trees}
  \item{pcfseq}{a list of piecewise constant functions}
  \item{stseq}{a list of shape trees}
  \item{hseq}{a vector of smoothing parameters corresponding to the estimates
  in the sequence;
  the smoothing parameter is the cardinality of the partitions
  of the aggregated histograms}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemela }
%\note{ ~~further notes~~ }
\seealso{ 
\code{\link{eval.bagg}}
}
\examples{
dendat<-sim.data(n=200,seed=1,type="mulmodII")

seed<-1        # seed for choosing bootstrap samples
sample="worpl" # without-replacement bootstrap
prune="on"     # we use CART-histograms
B<-2           # the number of histograms in the average

estiseq<-lstseq.bagg(dendat,B,maxleaf=10,lstree=TRUE,
         seed=seed,sample=sample,prune=prune)

mt<-modegraph(estiseq)

plotmodet(mt)

#scaletable(estiseq)

}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
