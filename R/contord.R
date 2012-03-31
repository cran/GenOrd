contord <-
function(marginal,Sigma,support=list(),Spearman=FALSE)
{
len<-length(support)
k<-length(marginal)
kj<-numeric(k)
Sigmaord<-diag(k)
for(i in 1:k)
{
kj[i]<-length(marginal[[i]])+1

if(len==0)
{
support[[i]]<-1:kj[i]
}

if(Spearman)
{
support[[i]]<-c(marginal[[i]],1)
}
}
# compute the quantiles corresponding to the cumulative marginal distributions
# of the ordinal r.v.
L<-vector("list",k)
for(i in 1:k)
{
L[[i]]<-qnorm(marginal[[i]])
L[[i]]<-c(-Inf,L[[i]],+Inf)
}
#
#
for(q in 1:(k-1))
{
for(r in (q+1):k)
{
# compute the joint probability mass function of the bivariate ordinal variable (components q and r)
pij<-matrix(0,kj[q],kj[r])
for(i in 1:kj[q])
{
for (j in 1:kj[r])
{
# set all the lower bound equal to -Inf and the upper bounds equal to +Inf
low<-rep(-Inf,k)
upp<-rep(Inf,k)
# set the desired lower and upper bounds for components q and r
low[q]<-L[[q]][i]
low[r]<-L[[r]][j]
upp[q]<-L[[q]][i+1]
upp[r]<-L[[r]][j+1]
# pvnorm computes the probability of a rectangle for a k-variate normal
pij[i,j]<-pmvnorm(low,upp,rep(0,k),corr=Sigma)
low<-rep(-Inf,k)
upp<-rep(Inf,k)
}
}
# the sum of pij must be one
sum(pij)
# mean and variances of the two ordinal r.v.
# NOTE: implicity here we assume a Likert scale for both r.v.:
# their ordered categories are coded as 1, 2, 3,..., k1 and 1, 2, 3,..., k2
my<-sum(apply(pij,2,sum)*support[[r]])
sigmay<-sqrt(sum(apply(pij,2,sum)*support[[r]]^2)-my^2)
mx<-sum(apply(pij,1,sum)*support[[q]])
sigmax<-sqrt(sum(apply(pij,1,sum)*support[[q]]^2)-mx^2)
mij<-support[[q]]%*%t(support[[r]])
muij<-sum(mij*pij)
# covariance between the ordinal r.v.
covxy<-muij-mx*my
corxy<-covxy/(sigmax*sigmay)
# return the correlation coefficient between the ordinal r.v.
Sigmaord[q,r]<-corxy
}
}
as.matrix(forceSymmetric(Sigmaord))
}

