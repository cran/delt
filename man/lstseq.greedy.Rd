\name{lstseq.greedy}
\alias{lstseq.greedy}
\title{ Calculates a scale of greedy histograms }
\description{
Calculates a scale of greedy histograms.
The histograms in the scale have a partition of growing cardinality. 
Returns a sequence of estimates as piecewise constant functions and
optionally as level set trees.
Optionally a shape tree for a level set of each estimate is calculated.
}
\usage{
lstseq.greedy(dendat, maxleaf, lstree = NULL, level = NULL)
}
\arguments{
  \item{dendat}{ n*d data matrix }
  \item{maxleaf}{ integer>1; the scale consists of histograms whose
   partitions have cardinality 1,...,maxleaf  }
  \item{lstree }{ if NULL, then level set trees are not calculated }
  \item{level}{ if NULL, then shape trees are not calculated;
   if positive number, then it is
   the level of the level sets for which the shape trees are calculated }
}
%\details{}
\value{
 A list with components
  \item{lstseq}{a list of level set trees}
  \item{pcfseq}{a list of piecewise constant functions}
  \item{stseq}{a list of shape trees}
  \item{hseq}{a vector of smoothing parameters corresponding to the members
  in the sequences; 
  the smoothing parameter is the cardinality of the partition}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemela }
%\note{ ~~further notes~~ }
\seealso{ 
\code{\link{eval.greedy}}, 
}
\examples{
dendat<-sim.data(n=200,seed=1,type="mulmodII")
estiseq<-lstseq.greedy(dendat,maxleaf=20,lstree=TRUE)

mt<-modegraph(estiseq)

plotmodet(mt)

#scaletable(estiseq)

}

\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
