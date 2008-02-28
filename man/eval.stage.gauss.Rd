\name{eval.stage.gauss}
\alias{eval.stage.gauss}
\title{ Returns a 1D Gaussian mixture density estimate }
\description{
 Estimates a 1D density with a mixture of Gaussians.
 The mixture is found by minimizing the L2 empirical risk 
 in a stagewise manner.
}
\usage{
eval.stage.gauss(dendat, M, mugrid, siggrid = 1, sigeka = TRUE, src = "c",
sampstart=FALSE, boost=FALSE, N=60)
}
\arguments{
  \item{dendat}{n-vector of 1D observations }
  \item{M}{integer >= 1; the number of mixture components in the estimate }
  \item{mugrid}{a vector of real numbers; the range for the means of mixture 
                components }
  \item{siggrid}{a vector of real numbers; the range of possible standard 
                 deviations in the mixture components  }
  \item{sigeka}{ TRUE or FALSE; if TRUE, then the standard deviation of the 
   first mixture component is equal to 1, otherwise the standard deviation 
   of the first mixture component is found by minimization  }
  \item{src}{ "R" or "c"; if "R", then the R-code is used, otherwise the 
   c-code is used}
  \item{sampstart}{internal}
  \item{boost}{internal}
  \item{N}{postive integer; the number of evaluation points}
}
\value{
  A piecewise constant function with the additional components:
  \item{muut }{vector of real numbers; the means of the mixture components}
  \item{sigit }{vector of positive real numbers; the standard deviations of
  the mixture components}
  \item{curmix }{a probability vector; the weights of the mixture components}
}
\references{ Jussi Klemela (2005). 
Density Estimation with Stagewise Optimization of the Empirical Risk }
\author{ Jussi Klemela }

\seealso{ 
\code{\link{eval.stage}}
}
\examples{
dendat<-sim.data(n=100,type="1d2modal",seed=1) 

mugrid<-seq(-1,5,0.3)    # grid of mu-values
siggrid<-seq(0.2,2,0.2)  # grid of sigma-values
M<-17                     # number of mixture components
pcf<-eval.stage.gauss(dendat,M,mugrid,siggrid)

dp<-draw.pcf(pcf)
plot(dp$x,dp$y,type="l")

# draw the estimate with the help of package "denpro"
#N<-100
#pcf2<-pcf.func("mixt",N,sig=pcf$sigit,M=pcf$muut,p=pcf$curmix)  
#pnum<-100
#dm<-draw.pcf(pcf2,pnum=pnum) 
#plot(dm$x,dm$y,type="l")

}

\keyword{smooth }



