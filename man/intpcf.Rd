\name{intpcf}
\alias{intpcf}
\title{ Calculates the integral of a piecewise constant function }
\description{
Calculates the integral of a piecewise constant function.
}
\usage{
intpcf(pcf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{ piecewise constant function }
}
%\details{}
\value{
Real number; the integral of the piecewise cosntant function
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemela}
%\note{ ~~further notes~~ }
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
dendat<-sim.data(n=100,seed=5,type="mulmodII")
eva<-eval.greedy(dendat,leaf=5)
intpcf(eva)
 
}

\keyword{ math }% at least one, from doc/KEYWORDS



