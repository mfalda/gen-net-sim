###     createRules  (2007-04-02)
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

createRules<-function(M,f.pr.and)
{
 cat(file=fp, "\n***createRules***\ninput:\n\tM = ", M, "\n\tf.pr.and = ...\n")
#This function creates for all the gene in the network a set of rules given as output in the list L
#In particular, for each gene i, the interactions among its regulators are codified in position  [[i]] of the list "L" as a vector representing a binary tree.
#see create.logicrule and the help of simulateprofiles for further details

# the function returns the list L of regulatory rules and the maximum number of nodes in the binary trees created to codify the rules

#\arguments{
#  \item{M}{(matrix of real numbers) : the connectivity matrix containing information on network topology and regulatory
#	     efficiency and sign: a non zero value in position (i,j) means that gene j regulates
#	     gene i with efficiency and sign indicated by the value itself;}
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



n<-dim(M)[1]
L<-list()
max.lengthL<-0

for(i in 1:n){
	op<-which(M[i,]!=0)
	l_op<-length(op)
	if (l_op>1)  logic_rule<-create.logicrule(op,f.pr.and)
	 else {if (l_op==1) logic_rule<-op else logic_rule<-NaN }
	lengthL<-length(logic_rule)
        max.lengthL<-max(max.lengthL,lengthL)

	L[[i]]<-logic_rule
}
cat(file=fp, "createRules output:\n\tlist(L, max.lengthL) = ")
print_list(L)
cat(file=fp, max.lengthL, "\n")
return(list(L, max.lengthL))
}


