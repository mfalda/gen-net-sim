###     next.op  (2007-04-02)
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

next.op<-function(r){
 cat(file=fp, "\n***next_op***\ninput:\n\tr = ")
print_list(r)
# this function allows to visit the binary tree r (vector of positive and negative integers)
# in the correct order

# INPUT: the binary tree r (vector of positive and negative integers)
# OUTPUT: either the index of the vector r to be visited next (integer)
#         or -1 if no more nodes in the tree have to be visited

k<-(-1)
l<-length(r);
for (j in 0:(l-1)){
	i<-(l-j)
	p<-r[i]
	if ((p<(-1))&(k==-1)) k<-i;
}
cat(file=fp, "next_op output:\n\tk = ", k, "\n")
return(k)
}