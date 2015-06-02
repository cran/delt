\name{scaspa}
\alias{scaspa}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Finds the number of modes of histograms which
are obtained by pruning an overfitting histogram}
\description{
Function "densplit" returns an overfitting histogram as
an evaluation tree.
Function "prune" finds the candidate nodes for pruning.
Function "scaspa" finds the number of modes (local maxima)
in the histograms which are obtained by using these candidate 
nodes for pruning.
}
\usage{
scaspa(treeseq, bind, eind)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{treeseq}{ a list returned by function "prune"}
 \item{bind}{ integer in 1:subnum, where subnum is the number of
  elements in field "leafs" of "treeseq". We need bind<eind.}
 \item{eind}{integer in 1:subnum, where subnum is the number of
  elements in field "leafs" of "treeseq". We need bind<eind. }
}
%\details{}
\value{
List with the following vectors whose length is (eind-bind+1)
  \item{moodilkm }{number of local maxima for each subtree}
  \item{alfas }{value of the smoothing parameter alpha for each subtree}
  \item{leafnums }{number of leaves for each subtree}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemela }
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ 
\code{\link{densplit}},
\code{\link{prune}}  
}

\examples{
set.seed(1)
dendat<-matrix(rnorm(20),10)  
minlkm<-2
et<-densplit(dendat,minlkm)
treeseq<-prune(et)
treeseq$leafs

scaspa(treeseq,1,5)

}

\keyword{ smooth }

