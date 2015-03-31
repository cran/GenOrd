\name{corrcheck}
\alias{corrcheck}
\title{
Checking correlations for feasibility
}
\description{
The function returns the lower and upper bounds of the correlation coefficients of each pair of ordinal/discrete variables given their marginal distributions, i.e., returns the range of feasible bivariate correlations.
}
\usage{
corrcheck(marginal, support = list(), Spearman = FALSE)
}
\arguments{
  \item{marginal}{
a list of \eqn{k} elements, where \eqn{k} is the number of variables.
The \eqn{i}-th element of \code{marginal} is the vector of the cumulative probabilities defining the marginal distribution of the \eqn{i}-th component of the  multivariate variable. If the \eqn{i}-th component has \eqn{k_i} categories, the \eqn{i}-th element of \code{marginal} will contain \eqn{k_i-1} probabilities (the \eqn{k_i}-th is obviously 1 and shall not be included).
}
  \item{support}{
a list of \eqn{k} elements, where \eqn{k} is the number of variables. The \eqn{i}-th element of support contains the ordered values of the support of the \eqn{i}-th variable. By default, the support of the \eqn{i}-th variable is \eqn{1,2,...,k_i}
}
  \item{Spearman}{
\code{TRUE} if we consider Spearman's correlation, \code{FALSE} (default) if we consider Pearson's correlation
}
}
\value{
The functions returns a list of two matrices: the former contains the lower bounds, the latter the upper bounds of the feasible pairwise correlations (on the extra-diagonal elements)
}
\author{
Alessandro Barbiero, Pier Alda Ferrari
}
\seealso{
\code{\link{contord}}, \code{\link{ordcont}}, \code{\link{ordsample}}
}
\examples{
# four variables
k <- 4
# with 2, 3, 4, and 5 categories (Likert scales, by default)
kj <- c(2,3,4,5)
# and these marginal distributions (set of cumulative probabilities)
marginal <- list(0.4, c(0.6,0.9), c(0.1,0.2,0.4), c(0.6,0.7,0.8,0.9))
corrcheck(marginal) # lower and upper bounds for Pearson's rho
corrcheck(marginal, Spearman=TRUE) # lower and upper bounds for Spearman's rho
# change the supports
support <- list(c(0,1), c(1,2,4), c(1,2,3,4), c(0,1,2,5,10))
corrcheck(marginal, support=support) # updated bounds
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
\keyword{distribution}
