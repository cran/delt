\name{prune}
\alias{prune}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Prepares for pruning an overfitting evaluation tree}
\description{Finds a sequence of nodes of an overfitting evaluation tree
 which are candidates to be the pruning nodes.
 Pruning a tree means removing a branch starting from a node. 
}
\usage{
prune(et)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{et}{ an evaluation tree;
  output of "eval.cart", "densplit", ...}
}
%\details{}
\value{
 A list containing the following components. 
  \item{tree}{the original tree which was given as the input}
  \item{delnodes }{vector giving a sequence of nodes in the order in
   which we should prune the branches starting from these nodes  }
  \item{delend}{vector whose length is the number of subtrees of the
   original tree. With the help of "delend" we define the subtrees.
   Elements of "delend" define a sequence of nodes from "delnodes" in the 
   following way: (1:delend[1]) is the first sequence,
   (delend[1]+1:delend[2]) is the second sequence, and so on.
   Then, i:th subtree is the result of pruning branches away whose
   roots are the nodes which are the first delend[i] elements of delnodes.}
  \item{leafs}{vector whose length is the number of subtrees of the
   original tree; number of leafs of the subtrees}
  \item{alfa}{vector whose length is the number of subtrees of the
    original tree; value of the corresponding alfa (complexity parameter)
    for every subtree}
  \item{loglik}{vector whose length is the number of subtrees of the
   original tree; the value of the likelihood criterion for the subtree}

}
%\references{ ~put references to the literature/website}
\author{ Jussi Klemelä}
%\note{ ~~further notes~~ }

\seealso{ 
\code{\link{densplit}},
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





