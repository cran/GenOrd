corrcheck <-
function(marginal,support=list(),Spearman=FALSE)
{
k<-length(marginal)
mcmax<-diag(1,k)
mcmin<-diag(1,k)
len<-length(support)

if(len==0)
{
for(i in 1:k)
{
support[[i]]<-1:(length(marginal[[i]])+1)
}
}

if(Spearman)
{
for(i in 1:k)
{
support[[i]]<-c(marginal[[i]],1)
}
}

for(i in 1:(k-1))
{
for (j in (i+1):k)
{
# variables i and j
P1<-c(0,marginal[[i]],1)
P2<-c(0,marginal[[j]],1)
l1<-length(P1)-1
l2<-length(P2)-1
p1<-numeric(0)
p2<-numeric(0)
# p1 and p2 will contain the probability mass function of the two variables i and j
# obtained by the corresponding cumulative marginal distributions
for(g in 1:l1)
{
p1[g]<-P1[g+1]-P1[g]
}
for(g in 1:l2)
{
p2[g]<-P2[g+1]-P2[g]
}
# expected values
E1<-sum(p1*support[[i]])
E2<-sum(p2*support[[j]])
# variances
V1<-sum(p1*support[[i]]^2)-E1^2
V2<-sum(p2*support[[j]]^2)-E2^2
# the maximum and minimum values of correlation coefficients
# are computed by building the cograduation and countergraduation tables
#
# cograduation table
y1<-1
y2<-1
lim<-0
E12<-0
PP1<-P1
PP2<-P2
PP1<-PP1[-1]
PP2<-PP2[-1]
while(length(PP1)>0)
{
E12<-E12+support[[i]][y1]*support[[j]][y2]*(min(PP1[1],PP2[1])-lim)
lim<-min(PP1,PP2)
if(PP1[1]==lim)
{
PP1<-PP1[-1]
y1<-y1+1
}
if(PP2[1]==lim)
{
PP2<-PP2[-1]
y2<-y2+1
}
}
# maximum correlation
c12<-(E12-E1*E2)/sqrt(V1*V2)
# countergraduation table
y1<-1
y2<-l2
lim<-0
E21<-0
PP1<-P1
PP2<-cumsum(rev(p2))
PP1<-PP1[-1]
while(length(PP1)>0)
{
E21<-E21+support[[i]][y1]*support[[j]][y2]*(min(PP1[1],PP2[1])-lim)
lim<-min(PP1,PP2)
if(PP1[1]==lim)
{
PP1<-PP1[-1]
y1<-y1+1
}
if(PP2[1]==lim)
{
PP2<-PP2[-1]
y2<-y2-1
}
}
# minimum correlation
c21<-(E21-E1*E2)/sqrt(V1*V2)
# bounds for element (i,j) saved in a matrix
mcmax[i,j]<-c12
mcmin[i,j]<-c21
}
}
mcmax<-forceSymmetric(mcmax)
mcmin<-forceSymmetric(mcmin)
# returns max e min correlation coefficients
list(mcmin,mcmax)
}

