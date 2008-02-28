\name{eval.pick}
\alias{eval.pick}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Returns a subtree of an evaluation tree}
\description{
Returns a subtree of an evaluation tree.
The subtree has a specified number of leafs.
The evaluation tree is calculated by "densplit" function.
To find out the possible values for the number of leaves
we use function "prune". 
}
\usage{
eval.pick(treeseq, leaf)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{treeseq}{ an overfitting evaluation tree with information on the
  possible pruning nodes; output of function "prune"}
 \item{leaf}{ positive integer; number of leaves in the subtree to
  be returned }
}
%\details{}
\value{
Returns an evaluation tree, see the documentation of function "eval.cart"
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



