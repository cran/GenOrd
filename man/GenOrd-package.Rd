\name{GenOrd-package}
\alias{GenOrd-package}
\docType{package}
\title{
Simulation of ordinal and discrete variables with given correlation matrix and marginal distributions
}
\description{
The package implements a procedure for generating samples from ordinal/discrete random variables with pre-specified correlation matrix and marginal distributions. 
It is developed in two steps: the first step (function \code{\link{ordcont}}) sets up the original continuous variables in order to achieve the final  discrete/ordinal variables meeting the experimental conditions; the second step (\code{\link{ordsample}}) generates samples from the adjusted original variables and discretizes them, thus simulating samples from the target variables. The procedure can handle both Pearson's correlation and Spearman's rho, and any finite support for the discrete variables.
The intermediate function \code{\link{contord}} computes the correlations of discrete/ordinal variables derived from correlated normal variables through discretization.
Function \code{\link{corrcheck}} returns the lower and upper bounds of the correlation coefficients of ordinal/discrete variables given their marginal  distributions, i.e. returns the range of feasible pairwise correlations.
This version has fixed some inconsistencies regarding Spearman correlation coefficient that affected the two previous versions.
}
\details{
\tabular{ll}{
Package: \tab GenOrd\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2012-05-16\cr
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

