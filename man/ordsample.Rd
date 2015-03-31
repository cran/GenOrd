\name{ordsample}
\alias{ordsample}
\title{
Drawing a sample of ordinal/discrete data
}
\description{
The function draws a sample from a multivariate ordinal/discrete variable with correlation matrix \code{Sigma} and pre-specified marginals \code{marginal}}
\usage{
ordsample(n, marginal, Sigma, support=list(), Spearman=FALSE, cormat="ordinal")
}
\arguments{
  \item{n}{
the sample size
}
  \item{marginal}{
a list of \eqn{k} elements, where \eqn{k} is the number of variables.
The \eqn{i}-th element of \code{marginal} is the vector of the cumulative probabilities defining the marginal distribution of the \eqn{i}-th component of the  multivariate variable. If the \eqn{i}-th component has \eqn{k_i} categories, the \eqn{i}-th element of marginal will contain \eqn{k_i-1} probabilities (the \eqn{k_i}-th is obviously 1 and shall not be included).
}
  \item{Sigma}{
the target correlation matrix of the ordinal/discrete variables
}
  \item{support}{
a list of \eqn{k} elements, where \eqn{k} is the number of variables. The \eqn{i}-th element of \code{support} contains the ordered values of the support of the \eqn{i}-th variable. By default, the support of the \eqn{i}-th variable is \eqn{1,2,...,k_i}
}
  \item{Spearman}{
if \code{TRUE}, the function finds Spearman's correlations (and it is not necessary to prvide \code{support}), if \code{FALSE} (default) Pearson's correlations 
}
  \item{cormat}{
\code{ordinal} if the \code{Sigma} in input is the target correlation matrix of ordinal/discrete variables; \code{continuous} if the \code{Sigma} in input is the intermediate correlation matrix of the multivariate standard normal
}
}
\value{
a \eqn{n\times k} matrix of discrete/ordinal data drawn from the \eqn{k}-variate discrete/ordinal r.v. with the desired marginal distributions and correlation matrix
}
\author{
Alessandro Barbiero, Pier Alda Ferrari
}
\seealso{
\code{\link{contord}}, \code{\link{ordcont}}, \code{\link{corrcheck}}
}
\examples{
# Example 1

# draw a sample from a bivariate ordinal variable
# with 4 of categories and asymmetrical marginal distributions
# and correlation coefficient 0.6 (to be checked)
k <- 2
marginal <- list(c(0.1,0.3,0.6), c(0.4,0.7,0.9))
corrcheck(marginal) # check ok
Sigma <- matrix(c(1,0.6,0.6,1),2,2)
# sample size 1000
n <- 1000
# generate a sample of size n
m <- ordsample(n, marginal, Sigma)
head(m)
# sample correlation matrix
cor(m) # compare it with Sigma
cumsum(table(m[,1]))/n
cumsum(table(m[,2]))/n # compare it with the two marginal distributions

# Example 1bis

# draw a sample from a bivariate ordinal variable
# with 4 of categories and asymmetrical marginal distributions
# and Spearman correlation coefficient 0.6 (to be checked)
k <- 2
marginal <- list(c(0.1,0.3,0.6), c(0.4,0.7,0.9))
corrcheck(marginal, Spearman=TRUE) # check ok
Sigma <- matrix(c(1,0.6,0.6,1),2,2)
# sample size 1000
n <- 1000
# generate a sample of size n
m <- ordsample(n, marginal, Sigma, Spearman=TRUE)
head(m)
# sample correlation matrix
cor(rank(m[,1]),rank(m[,2])) # compare it with Sigma
cumsum(table(m[,1]))/n
cumsum(table(m[,2]))/n # compare it with the two marginal distributions

# Example 2

# draw a sample from a 4-dimensional ordinal variable
# with different number of categories and uniform marginal distributions
# and different correlation coefficients
k <- 4
marginal <- list(0.5, c(1/3,2/3), c(1/4,2/4,3/4), c(1/5,2/5,3/5,4/5))
corrcheck(marginal)
# select a feasible correlation matrix
Sigma <- matrix(c(1,0.5,0.4,0.3,0.5,1,0.5,0.4,0.4,0.5,1,0.5,0.3,0.4,0.5,1),4,4,byrow=TRUE)
Sigma
# sample size 100
n <- 100
# generate a sample of size n
set.seed(1)
m <- ordsample(n, marginal, Sigma)
# sample correlation matrix
cor(m) # compare it with Sigma
cumsum(table(m[,4]))/n # compare it with the fourth marginal
head(m)
# or equivalently...
set.seed(1)
res <- ordcont(marginal, Sigma)
res[[1]] # the intermediate correlation matrix of the multivariate normal
m <- ordsample(n, marginal, res[[1]], cormat="continuous")
head(m)
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}