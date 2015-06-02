\name{makebina}
\alias{makebina}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Tranforms and evaluation tree to the tree object of R 
}
\description{
Evaluation trees are such trees that are implemented 
with "left" and "right" pointers.
We transform this tree representation to the tree object
of the package "tree",
so that it can be plotted by 
"plot.tree" function from package "tree" or 
by "draw.tree" function from package "maptree".
}
\usage{
makebina(et)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{et}{evaluation tree; 
  result of "densplit", "eval.pick", "eval.cart", ...;
  see the documentation of "eval.cart"}
}
%\details{
%}
\value{
Returns an object of class tree. 
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemela }
%\note{ ~~further notes~~ }

\seealso{  
\code{\link{densplit}}, 
\code{\link{eval.pick}}   
}

\examples{
library(denpro)
dendat<-sim.data(n=100,seed=5,type="mulmodII")
et<-densplit(dendat)
mb<-makebina(et)

set.seed(1)
dendat<-matrix(rnorm(20),10)  
et<-densplit(dendat,minobs=2)
mb<-makebina(et)

#library(tree)
#plot.tree(mb)

#library(maptree)
#draw.tree(mb)

}

\keyword{hplot}














