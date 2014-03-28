\name{ordcont}
\alias{ordcont}
\title{
computes the "intermediate" correlation matrix for the multivariate standard normal in order to achieve the "target" correlation matrix for the ordinal/discrete variables
}
\description{
The function computes the correlation matrix of the k-dimensional standard normal r.v. yielding the desired correlation matrix \code{Sigma} for the k-dimensional r.v. with desired marginal distributions \code{marginal}
}
\usage{
ordcont(marginal,Sigma,support=list(),Spearman=FALSE,epsilon=1e-06,maxit=100)
}
\arguments{
  \item{marginal}{
a list of \eqn{k} elements, where \eqn{k} is the number of variables.
The \eqn{i}-th element of marginal is the vector of the cumulative probabilities defining the marginal distribution of the \eqn{i}-th component of the  multivariate variable. If the \eqn{i}-th component has \eqn{k_i} categories, the \eqn{i}-th element of marginal will contain \eqn{k_i-1} probabilities (the \eqn{k_i}-th is  obviously 1).
}
  \item{Sigma}{
the target correlation matrix of the ordinal/discrete variables
}
  \item{support}{
a list of \eqn{k} elements, where \eqn{k} is the number of variables. The \eqn{i}-th element of support contains the ordered values of the support of the \eqn{i}-th variable. By default, the support of the \eqn{i}-th variable is \eqn{1,2,...,k_i}
}
  \item{Spearman}{
if TRUE, the function finds Spearman's correlations (and it is not necessary to prvide \code{support}), if FALSE (default) Pearson's correlations
}
  \item{epsilon}{
the maximum tolerated error among target and actual correlations
}
  \item{maxit}{
the maximum number of iterations allowed for the algorithm
}
}
\value{
a list of four elements
  \item{SigmaC}{
the correlation matrix of the multivariate standard normal variable
}
  \item{SigmaO}{
the actual correlation matrix of the discretized variables (it should approximately coincide with the target correlation matrix \code{Sigma})
}
  \item{Sigma}{
the target correlation matrix of the ordinal/discrete variables
}
  \item{niter}{
a matrix containing the number of iterations performed by the algorithm, one for each pair of variables
}
  \item{maxerr}{
the actual maximum error (the absolute maximum deviation between actual and target correlations of the ordinal/discrete variables)
}
}
\author{
Alessandro Barbiero, Pier Alda Ferrari
}
\note{
For some choices of \code{marginal} and \code{Sigma}, the algorithm may not provide a feasible correlation matrix \code{SigmaC}. In this case, a warning is issued and the function supplies as \code{SigmaC} - just for consistency - the nearest correlation matrix, according to the algorithm discussed in Higham (2002). Obviuously, the corresponding \code{SigmaO} may sensibly differ from \code{Sigma}.
The value of the maximum tolerated absolute error \code{epsilon} on the elements of the correlation matrix for ordinal r.v. can be set by the user: a value between 0.000001 and 0.01 seems to be an acceptable compromise assuring both the precision of the results and the convergence of the algorithm; moreover, a maximum number of iteration can be chosen (\code{maxit}), in order to avoid possible endless loops in case of non-convergence
}
\references{
N. Higham (2002) Computing the nearest correlation matrix: A problem from finance, IMA Journal of Numerical Analysis, 22, 329-343.
}
\seealso{
\code{\link{contord}}, \code{\link{ordsample}}, \code{\link{corrcheck}}
}
\examples{
# consider a 4-dimensional ordinal variable
k<-4
# with different number of categories
kj<-c(2,3,4,5)
# and uniform marginal distributions
marginal<-list(0.5,c(1/3,2/3),c(1/4,2/4,3/4),c(1/5,2/5,3/5,4/5))
corrcheck(marginal)
# and the following correlation matrix
Sigma<-matrix(c(1,0.5,0.4,0.3,0.5,1,0.5,0.4,0.4,0.5,1,0.5,0.3,0.4,0.5,1),4,4,byrow=TRUE)
Sigma
# the correlation matrix of the standard 4-dimensional standard normal
# ensuring Sigma is
res<-ordcont(marginal, Sigma)
res[[1]]
# change some marginal distributions
marginal<-list(0.3,c(1/3,2/3),c(1/5,2/5,3/5),c(0.1,0.2,0.4,0.6))
corrcheck(marginal)
# and notice how the correlation matrix of the multivariate normal changes...
res<-ordcont(marginal, Sigma)
res[[1]]
# change Sigma, adding a negative correlation
Sigma[1,2]<--0.2
Sigma[2,1]<-Sigma[1,2]
Sigma
# checking whether Sigma is still positive definite
eigen(Sigma)$values # all >0, OK
res<-ordcont(marginal, Sigma)
res[[1]]
}
\keyword{htest}
\keyword{models}
\keyword{distribution}
\keyword{datagen}

