\name{GenOrd-package}
\alias{GenOrd-package}
\docType{package}
\title{
Simulation of Discrete Random Variables with Given Correlation Matrix and Marginal Distributions
}
\description{
The package implements a procedure for generating samples from a multivariate discrete random variable with pre-specified correlation matrix and marginal distributions. The marginal distributions are linked together through a gaussian copula.
The procedure is developed in two steps: the first step (function \code{\link{ordcont}}) sets up the gaussian copula in order to achieve the desired correlation matrix on the target random discrete components; the second step (\code{\link{ordsample}}) generates samples from the target variables.
The procedure can handle both Pearson's and Spearman's correlations, and any finite support for the discrete variables.
The intermediate function \code{\link{contord}} computes the correlations of the multivariate discrete variable derived from correlated normal variables through discretization.
Function \code{\link{corrcheck}} returns the lower and upper bounds of the correlation coefficient of each pair of discrete variables given their marginal distributions, i.e., returns the range of feasible bivariate correlations.

This version has fixed some drawbacks in terminology in the previous version; the only actual change concerns the parameter \code{cormat} in the \code{\link{ordsample}} function. Further examples of implementation have been added.
}
\details{
\tabular{ll}{
Package: \tab GenOrd\cr
Type: \tab Package\cr
Version: \tab 1.4.0\cr
Date: \tab 2015-09-11\cr
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

A. Barbiero, P.A. Ferrari (2014) Simulation of correlated Poisson variables. Applied Stochastic Models in Business and Industry, doi 10.1002/asmb.2072
}
\keyword{ package }
\seealso{
\code{\link{contord}}, \code{\link{ordcont}}, \code{\link{corrcheck}}, \code{\link{ordsample}}
}

