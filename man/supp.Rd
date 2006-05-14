\name{supp}
\alias{supp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Returns the bounding box of observations}
\description{
Returns the smallest rectangle containing the observations.
The sides of the rectangle are parallel to the coordinate axis.
}
\usage{
supp(dendat, epsi=0, blown=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{dendat}{n*d data matrix}
 \item{epsi}{positive number: option to return the smallest rectangle
  such that the epsi-shrinkage of the rectangle containes the observations}
 \item{blown}{internal}
}
%\details{}
\value{
2*d vector: (2*i-1)-element, i=1,...,d, is the start of the i:th
interval and (2*i)-element is the end of the i:th interval
}
%\references{ http://www.rni.helsinki.fi/~jsk/delt }
\author{Jussi Klemelä}
%\note{ ~~further notes~~ }

%\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

\examples{
set.seed(1)
dendat<-matrix(rnorm(20),10)
supp(dendat)

}

\keyword{ multivariate }%-- one or more ...







