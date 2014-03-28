\name{contord}
\alias{contord}
\title{
Correlations of discretized variables
}
\description{
The function computes the correlation matrix of the \eqn{k} variables, with given marginal distributions, derived discretizing a \eqn{k}-variate standard normal variable with given correlation matrix
}
\usage{
contord(marginal, Sigma, support = list(), Spearman = FALSE)
}
\arguments{
  \item{marginal}{
a list of \eqn{k} elements, where \eqn{k} is the number of variables.
The \eqn{i}-th element of marginal is the vector of the cumulative probabilities defining the marginal distribution of the \eqn{i}-th component of the  multivariate variable. If the \eqn{i}-th component has \eqn{k_i} categories, the \eqn{i}-th element of marginal will contain \eqn{k_i-1} probabilities (the \eqn{k_i}-th is  obviously 1).
}
  \item{Sigma}{
the correlation matrix of the standard multivariate normal variable
}
  \item{support}{
a list of \eqn{k} elements, where \eqn{k} is the number of variables. The \eqn{i}-th element of support contains the ordered values of the support of the \eqn{i}-th variable. By default, the support of the \eqn{i}-th variable is \eqn{1,2,...,k_i}
}
  \item{Spearman}{
if TRUE, the function finds Spearman's correlations (and it is not necessary to provide \code{support}), if FALSE (default) Pearson's correlations
}
}
\value{
the correlation matrix of the discretized variables
}
\author{
Alessandro Barbiero, Pier Alda Ferrari
}
\seealso{
\code{\link{ordcont}}, \code{\link{ordsample}}, \code{\link{corrcheck}}
}
\examples{
# consider 4 discrete variables
k<-4
# with these marginal distributions
marginal<-list(0.4,c(0.3,0.6),c(0.25,0.5,0.75),c(0.1,0.2,0.8,0.9))
# generated discretizing a multivariate standard normal variable
# with correlation matrix
Sigma<-matrix(0.5,4,4)
diag(Sigma)<-1
# the resulting correlation matrix for the discrete variables is
contord(marginal,Sigma)
# note all the correlations are smaller than the original 0.6
# change Sigma, adding a negative correlation
Sigma[1,2]<--0.15
Sigma[2,1]<-Sigma[1,2]
Sigma
# checking whether Sigma is still positive definite
eigen(Sigma)$values # all >0, OK
contord(marginal, Sigma)
}
\keyword{htest}
\keyword{models}
\keyword{distribution}
\keyword{datagen}
