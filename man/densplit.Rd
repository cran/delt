\name{densplit}
\alias{densplit}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Calculation of an overfitting histogram}
\description{
The function returns an overfitting histogram when a data matrix
is given as an input.
The output is an evaluation tree which is grown with greedy growing.
The evaluation tree defines a partition of the sample space.
The evaluation tree may be pruned to get a density estimate.
}
\usage{
densplit(dendat, minobs=NULL, leaf=0, method="loglik", 
splitscan=0, seedf=1, suppo=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{dendat}{n*d data matrix}
 \item{minobs}{ non-negative integer; splitting of a bin will be continued if
 the bin containes "minobs" or more observations}
 \item{leaf}{internal (maximal number of leafs in the evaluation tree)}
 \item{method}{ "loglik" or "projec"; the contrast function }
 \item{splitscan}{internal (random selection of splits)}
 \item{seedf}{internal}
 \item{suppo}{ 2*d vector of real numbers; 
  the rectangle to be splitted; 
  the rectangle has to contain the data}
}
%\details{
%}
\value{
Returns an evaluation tree as a list of vectors. 
  \item{direc }{integer in 1,...,d; variable which is splitted }
  \item{split }{real number; splitting point}
  \item{mean }{nonnegative number; value of the histogram on the
   rectangle corresponding to the node}
  \item{nelem }{nonnegative integer; number of observations in the
   rectangle corresponding to the node}
  \item{ssr }{real number; value of the likelihood criterion}
  \item{volume }{non-negative number; volume of the rectangle
   corresponding to the node }
  \item{left }{non-negative integer; link to the left child, 0 if
   terminal node }
  \item{right }{non-negative integer; link to the right child, 0 if
   terminal node}
  \item{low }{the lower vertice of the rectangles}
  \item{upp }{the upper vertice of the rectangles}
  \item{N }{the number of grid points at each direction}
  \item{support}{the support of the histogram}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemela }
%\note{ ~~further notes~~ }


\seealso{  
\code{\link{prune}},
\code{\link{eval.pick}}
}

\examples{
dendat<-sim.data(n=200,seed=5,type="mulmodII")
et<-densplit(dendat)

treeseq<-prune(et)
treeseq$leafs
len<-length(treeseq$leafs)

leaf<-treeseq$leafs[len-10]
leaf
etsub<-eval.pick(treeseq,leaf=leaf)

dp<-draw.pcf(etsub)
persp(dp$x,dp$y,dp$z,phi=25,theta=-120)

}

\keyword{smooth}
\keyword{multivariate}
