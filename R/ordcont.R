ordcont <-
function(marginal, Sigma, support=list(), Spearman=FALSE, epsilon=0.000001, maxit=100)
{
# cormat=target correlation matrix for ordinal r.v. (upper triangle)
# epsilon=maximum absolute error between target and actual correlations for ordinal r.v.
# maxit=maximum number of iterations allowed for the iterative step

# k=number of variables
k<-length(marginal)
# kj=number of categories for the k variables (vector of k integer numbers)
kj<-numeric(k)
for(i in 1:k)
{
kj[i]<-length(marginal[[i]])+1
}

## check on correlation matrix for positive semidefiniteness

for(g in 2:k)
{
if (det(Sigma[1:g,1:g])<=0)
{
# exit from the procedure
stop("Main minor number ",g," is not positive!","\n")
}
}

## check for feasibility of single elements of Sigma given the marginal distributions

mcmin<-corrcheck(marginal,support,Spearman)[[1]]
mcmax<-corrcheck(marginal,support,Spearman)[[2]]
if(sum(mcmin<=Sigma & Sigma<=mcmax)!=k^2)
{
# exit from the procedure
stop("Some correlation coefficients are not feasible!","\n",
"Please use function corrcheck to get lower and upper bounds!","\n")
}

# Sigma0 is now the target correlation matrix for ordinal r.v.
Sigma0<-Sigma
Sigmaord<-Sigma
# find the correlation coefficients for ordinal r.v. corresponding to
# correlation coefficients for normal r.v. with correlation matrix equal to
# the target correlation matrix for ordinal r.v.

Sigmaord<-contord(marginal,Sigma,support,Spearman)

Sigmaold<-Sigma
Sigmaordold<-Sigmaord

# iteration step

it<-0
# the loop ends when the maximum absolute error is smaller than the prescribed epsilon
# or when the number of iterations reaches the prescribed maxit
while(max(abs(Sigmaord-Sigma0)>epsilon) & it<maxit)
{
for(q in 1:(k-1))
{
for(r in (q+1):k)
{
### update the elements of the correlation matrix for the continuous r.v.
# null correlation
if (Sigma0[q,r]==0)
{
Sigma[q,r]<-0
}
else
{
# in order to avoid the updated elements are greater than 1, use this update...
# that mildly increases the correlations of the normal variables...
if (Sigma0[q,r]*(Sigma0[q,r]/Sigmaordold[q,r])>=1)
{
Sigma[q,r]<-Sigmaold[q,r]*(1+0.1*(1-Sigmaold[q,r])*sign(Sigma0[q,r]-Sigmaord[q,r]))
}
# update: the new correlation coefficient for continuous r.v.
# is equal to the old one multiplied by the ratio between the target and the actual
# correlation coefficient for ordinal r.v.
else
{
Sigma[q,r]<-Sigmaold[q,r]*(Sigma0[q,r]/Sigmaord[q,r])
}
}
# symmetry
Sigma[r,q]<-Sigma[q,r]
}
}
###
# if the correlation matrix built with the updated correlation coefficients
# is no more definite positive, use the nearPD algorithm to fix it
Sigma<-as.matrix(nearPD(Sigma, corr=TRUE)$mat)
###
# discretization: compute the probability mass functions of the correlated ordinal variables...
###
Sigmaord<-contord(marginal,Sigma,support,Spearman)

Sigmaold<-Sigma
# next iteration
it<-it+1
}
# maximum absolute error between the current correlation matrix and the target correlation matrix
emax<-max(abs(Sigmaord-Sigma0))
# return a list comprising the correlation matrix for continuous r.v., the correlation matrix for ordinal r.v.
# and the target ocrrelation matrix for ordinal r.v., the number of iterations, the maximum absolute error
list(SigmaC=Sigma, SigmaO=Sigmaord, Sigma=Sigma0, niter=it, maxerr=emax)
}

