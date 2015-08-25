###     create.logicrule  (2007-04-02)
###
###	Copyright (C) 2006, 2007  Di Camillo Barbara
###
###	This program is free software, part of the package netsim; you can redistribute it and/or
###	modify it under the terms of the GNU General Public License
###	as published by the Free Software Foundation; either version 2
###	of the License, or (at your option) any later version.

###	This program is distributed in the hope that it will be useful,
###	but WITHOUT ANY WARRANTY; without even the implied warranty of
###	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###	GNU General Public License for more details.

###	You should have received a copy of the GNU General Public License
###	along with this program; if not, write to the Free Software
###	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


create.logicrule<-function(v, f.pr.and=NULL)

{
 cat(file=fp, "\n***create_logicRule***\ninput:\n\tv = ", v, "\n\tf_pr_and = ...\n")
#create a binary tree that codify for regulatory interactions amnong the regulators v
# the output r is a vector of positive and negative integers corresponding to the ninary tree representation
#In particular, the root of the tree has position 1 in the vector and for every node of the tree in position i, the left son has position 2*i and the right son has position 2*i+1
#If a node corresponds to a regulator, it becomes a leaf; therefore, the maximum number of levels that the tree can have corresponds to the number of regulators.
#Functions and regulators can be positioned anywhere on the tree, with the only constraint that regulators must have a function as parent node, so, for example,
#the root of the tree is always a function, unless there is a single regulator; in this case the tree collapses in a single node corresponding to the regulator.
#Only fCOOPERATIVE and fSYNERGIC are used as function nodes

#-The COOPERATIVE rule is codified by -2
#-The SYNERGIC rule is codified by -3
#-The regulatory genes are codified by their index (i.e. an integer ranging from 1 to N).
#-The empty node by a -1
#For example,
#the rule fCOOPERATIVE(fSYNERGIC(x1,x2), fCOOPERATIVE(x3,x4)) is codified by the vector v
#v=c(-2,-3,-2,1,2,3,4)
#the rule fCOOPERATIVE(fSYNERGIC(x1,x2), fCOOPERATIVE(fCOOPERATIVE(x3,x5),x4)) is codified by the vector v
#v=c(-2,-3,-2,1,2,-2,4,-1,-1,-1,-1,3,5,-1,-1)



#\arguments{
#  \item{v}{(vector of positive integers) vector of regulators}
#  \item{f.pr.and}{Function with domain and codomain in [0-1] that expresses the probability to obtain a cooperative rather than a synergic rule as a function of the level of the binary tree used to implement the interactions among the regulators. Given a binary tree of maximum L levels, the probability to have a cooperative rule at level i is f.pr.and(i/L), whereas the probability to have a synergic rule is 1-f.pr.and(i/L). See details}
#  }

#In particular, if f.pr.and=NULL, the probability to have a cooperative or a synergic function is the same (50%).
#Alternatively, it is possible to assign f.pr.and a function with domain and codomain in [0-1], which expresses the probability
#to obtain a cooperative rather than a synergic rule as a function of the level of the binary tree used to implement the interactions among the regulators.
#Given a binary tree of maximum L levels, the probability to have a cooperative rule at level i is f.pr.and(i/L),
#whereas the probability to have a synergic rule is 1-f.pr.and(i/L).
#Examples of suitable functions are:
#f1<-function(norm.lev) {y<-1/(1+exp(-10*(norm.lev-0.5))); return(y)} #norm.lev=i/L (probability to have a cooperative function increases with the tree level)
#or
#f1<-function(norm.lev) {y<-0.8; return(y)} #probability to have a cooperative function is constant




l<-length(v)

x<-(1:l)/(l)
if (is.null(f.pr.and)) pr.and<-rep(0.5,length(x))
else pr.and<-f.pr.and(x)
pr.or<-1-pr.and

r<-c(sample_s(c(-2,-3),1, prob=c(pr.and[1],pr.or[1]), chi="create_logicRule"),0,0)
o<-c(2,3)
if (r[1]==(-2)) {pr.or<-rep(0,length(pr.or)); pr.and<-rep(1,length(pr.and)) }
black.p<-blacklist<-vector()

if (l>2){
	for(i in 1:(l-2)){
		lr<-length(r)
		lo<-length(o)
		e<-sample_s(seq(1,lo,1),1, chi="create_logicRule")
		p<-o[e]
		if (2*p>lr){
			nvect<-rep(-1,lr)
			r<-c(r,nvect)
		}
		o[e]<-2*p
		o<-c(o,(2*p+1))

		L.p<-length(black.p)
		if (L.p>0) {
			for (j in (1:L.p)) {
				blacklist<-vector();
				x<-ceiling(log(max(o)/black.p[j],2)) ;
				for (i in (1:x))
					blacklist<-c(blacklist,black.p[j]*(2^i)+(0:((2^i)-1)))
			}
		}
		livello<-floor(log(p,2))+1
		if (p%in%blacklist)    r[p]<-(-2)
    		else r[p]<-sample_s(c(-2,-3),1,prob=c(pr.and[livello],pr.or[livello]), chi="create_logicRule")
		if (r[p]==(-2)) black.p<-c(black.p,p)

    r[2*p]<-0
		r[2*p+1]<-0
	}
}
for(i in 1:length(o)){
	r[o[i]]<-v[i]
}
k<-length(r)
repeat{
	if (r[k]!=-1) break
	r<-r[-k]
	k<-k-1
}
cat(file=fp, "create_logicRule output:\n\tr = ")
print_list(r)
return(r)
}
