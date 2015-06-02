\name{lstseq.cart}
\alias{lstseq.cart}
\title{ Calculates a scale of CART histograms }
\description{
Calculates a scale of CART histograms.
The histograms in the scale have partitions of growing cardinality.
Returns a sequence of estimates as piecewise constant functions and
optionally as level set trees.
Optionally a shape tree for a level set of each estimate is calculated.
}
\usage{
lstseq.cart(treeseq, maxleaf=NULL, lstree=NULL, level = NULL,
indvec = NULL)
}
\arguments{
  \item{treeseq}{ output of function "prune" }
  \item{maxleaf}{ the maximal cardinality of the partitions of the histograms
   in the sequence } 
  \item{lstree }{ if NULL, then level set trees are not calculated }
  \item{level}{ if NULL, then shape trees are not calculated;
   if positive number, then it is
   the level of the level sets for which the shape trees are calculated }
  \item{indvec}{ a vector of indeces;
   chooses a subset of the complete sequence of subhistograms
   of the overfitting histogram }
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
\code{\link{densplit}}, 
\code{\link{prune}} 
}
\examples{
library(denpro)
dendat<-sim.data(n=100,seed=1,type="mulmodII")
et<-densplit(dendat)
treeseq<-prune(et)

estiseq<-lstseq.cart(treeseq,maxleaf=20,lstree=TRUE)

mt<-modegraph(estiseq)

plotmodet(mt)

#scaletable(estiseq)

}

\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
