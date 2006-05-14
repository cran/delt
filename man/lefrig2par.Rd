\name{lefrig2par}
\alias{lefrig2par}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ 
Transforms an evaluation tree so that it can be plotted with
the "plottree" function of package "denpro" 
}
\description{
Evaluation trees are trees which are implemented with 
"left" and "right" pointers.
We transform this tree representation to the representation with 
"parent" pointers (level set tree), so that it can be plotted with
the "plottree" function of package "denpro".
For example, functions "densplit" and "eval.pick" return evaluation trees.
}
\usage{
lefrig2par(et)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{et}{evaluation tree; 
  result of "densplit", "eval.pick", "eval.cart", ...;
  see the documentation of "eval.cart"}
}
%\details{}
\value{
Returns a level set tree:
  \item{parent }{parent links }
  \item{level }{height of the node}
  \item{center }{determines the ordering of the nodes}
  \item{volume }{determines the horizontal positioning of the nodes}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemelä }
%\note{ ~~further notes~~ }


\seealso{  
\code{\link{densplit}}, 
\code{\link{eval.pick}}
}

\examples{
dendat<-sim.data(n=200,seed=5,type="mulmodII")
et<-densplit(dendat)

lst<-lefrig2par(et)

plottree(lst)

}

\keyword{hplot}


