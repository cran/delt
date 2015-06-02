\name{cluster.lst}
\alias{cluster.lst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assigns labels to data points
}
\description{
Assigns labels to data points according to cluster membership, 
when the clusters are defined as high density regions 
}
\usage{
cluster.lst(dendat, h, N = NULL, cut = NULL, lambda = NULL, complete = FALSE,
type = "grid", labels = "number", nodes = NULL, minobs = 1)
}
\arguments{
  \item{dendat}{n*d matrix of real numbers; the data matrix.}
  \item{h}{positive real number; smoothing parameter of a kernel density estimator}
  \item{N}{d vector of positive integers;  
   a kernel estimate is evaluated on a regular grid which is such that
   in direction i there are N[i] points;
   N is needed only when type="grid".}
  \item{cut}{real number between 0 and 1; this parameter is used to determine
   the level "lambda" of the level set whose disconnected components determine the clusters.}
  \item{lambda}{positive real number between; "lambda" is the level of the level set 
   whose disconnected components determine the clusters.}
  \item{complete}{TRUE or FALSE; if complete=FALSE, then partial clustering is 
   performed, otherwise complete clustering is performed.}
  \item{type}{either "grid" or "adaptive"; if type="grid", then the density is
   estimated using a discretized kernel estimator with a regular grid;
   otherwise the density is estimated using a discretized kernel estimator 
   with an adaptive grid.}
  \item{labels}{if labels="number", then the cluster labels are integers 1,2,...,
   otherwise the cluster labels are colors.}
  \item{nodes}{a vector of positive integers; contains pointers to the nodes of
   a level set tree; the nodes indicate which disconnected components of level sets
   define the clusters.}
 \item{minobs}{a positive integer; this is a parameter of function "pcf.greedy.kernel".}
}
%\details{}
\value{
a vector of cluster labels;
the vector has length equal to the number of rows of the data matrix "dendat".
The cluster labels are either numbers or names of colors.
}
%\references{http://jussi.klemela.com/delt/}
\author{Jussi Klemela}
%\note{ ~~further notes~~ }

 \seealso{
\code{\link{pcf.greedy.kernel}}
}

\examples{
library(denpro)
# generate data
seed<-1
n<-50
d<-2
l<-3; D<-4; c<-D/sqrt(2)
M<-matrix(0,l,d); M[2,]<-c; M[3,]<--c
sig<-matrix(1,l,d)
p<-rep(1/l,l)
dendat<-sim.data(type="mixt",n=n,M=M,sig=sig,p=p,seed=seed)

# partial clustering with a fixed level 
h<-(4/(d+2))^(1/(d+4))*n^(-1/(d+4))*apply(dendat,2,sd)
N<-rep(20,d)
cl<-cluster.lst(dendat,h,N=N,labels="colors",type="grid",lambda=0.02)
#plot(dendat,col=cl)

# complete clustering with a fixed level
cl<-cluster.lst(dendat,h,N=N,complete=TRUE,labels="colors",type="grid",lambda=0.02)
#plot(dendat,col=cl)

# complete clustering with locally changing levels
N<-rep(20,d)
pcf<-pcf.kern(dendat,h,N)
lst<-leafsfirst(pcf)
nodes<-findbnodes(lst,modenum=3)
cl<-cluster.lst(dendat,h,N,nodes=nodes,complete=TRUE,labels="colors")
#plot(dendat,col=cl)
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{smooth}% __ONLY ONE__ keyword per line

