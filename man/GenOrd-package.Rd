\name{GenOrd-package}
\alias{GenOrd-package}
\docType{package}
\title{
Simulation of Ordinal and Discrete Variables with Given Correlation Matrix and Marginal Distributions
}
\description{
The package implements a procedure for generating samples from ordinal/discrete random variables with pre-specified correlation matrix and marginal distributions. The marginal distribution are linked together through a gaussian copula.
The procedure is developed in two steps: the first step (function \code{\link{ordcont}}) sets up the gaussian copula in order to achieve the desired correlation matrix on the target ordinal/discrete variables; the second step (\code{\link{ordsample}}) generates samples from the target variables.
The procedure can handle both Pearson's and Spearman's correlations, and any finite support for the discrete variables.
The intermediate function \code{\link{contord}} computes the correlations of discrete/ordinal variables derived from correlated normal variables through discretization.
Function \code{\link{corrcheck}} returns the lower and upper bounds of the correlation coefficient of each pair of ordinal/discrete variables given their marginal distributions, i.e., returns the range of feasible bivariate correlations.

This version has fixed some drawbacks of \code{\link{ordcont}}, which now stops and returns an error if the resulting correlation matrix of the gaussian copula is not a valid correlation matrix. Validity checks on the marginal distributions and correlation matrix provided by the user in the four functions have been added.
}
\details{
\tabular{ll}{
Package: \tab GenOrd\cr
Type: \tab Package\cr
Version: \tab 1.2.0\cr
Date: \tab 2015-03-31\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Alessandro Barbiero, Pier Alda Ferrari
Maintainer: Alessandro Barbiero <alessandro.barbiero@unimi.it>
}
\references{
P.A. Ferrari, A. Barbiero (2012) Simulating ordinal data, Multivariate Behavioral Research, 47(4), 566-589
}

\keyword{ package }
\seealso{
\code{\link{contord}}, \code{\link{ordcont}}, \code{\link{corrcheck}}, \code{\link{ordsample}}
}

