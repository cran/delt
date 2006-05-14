\name{plotparti}
\alias{plotparti}
\title{ Draws a partition }
\description{
Draws the partition calculated with function "partition". 
}
\usage{
plotparti(pa, d1 = NULL, d2 = NULL, dendat = NULL, restri = NULL, 
pch = 21, support = pa$support, col = "black")
}
\arguments{
  \item{pa}{ partition; output of function "partition" }
  \item{d1}{ integer 1,...,d;
   for the case the partition is a partition of a higher than 2 dimensional
   rectangle, "d1" is the first direction of the partition }
  \item{d2}{integer 1,...,d;
   for the case the partition is a partition of a higher than 2 dimensional
   rectangle, "d2" is the second direction of the partition }
  \item{dendat}{ n*d data matrix;
   if given as an argument it will also be plotted }
  \item{restri}{ internal }
  \item{pch}{ symbol for plotting "dendat";
   see function "points" }
  \item{support}{ the bounds of the partition }
  \item{col}{ color of the lines of the partition }
}
%\details{}
\value{
A plot at the graphics window
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemelä }
%\note{ ~~further notes~~ }
\seealso{ 
\code{\link{partition}} 
}

\examples{
dendat<-sim.data(n=200,seed=5,type="mulmodII")
et<-densplit(dendat)

pa<-partition(et)

plotparti(pa)

}

\keyword{ hplot }% at least one, from doc/KEYWORDS



