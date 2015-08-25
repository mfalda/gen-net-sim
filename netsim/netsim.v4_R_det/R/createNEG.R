###     createNEG  (2007-04-02)
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

createNEG<-function(M){

# this function could be simplified and put directly in the main code of the function simulatenet.
# It randomply assignes values 1 or -1 to elements (i,j) different from 0 in M, and return the modified
# version of M

#\arguments{
#  \item{M}{(matrix of 0 or 1) : the connectivity matrix containing information on network topology M(i,j)=1 means that gene j regulates
#	     gene i }
#  }

cat(file=fp, "\n***createNEG***\ninput:\n\tM = ", M, "\n")
ind<-which(M==1,arr.ind=TRUE)
L<-dim(ind)[1]
segno<-sample_s(c(-1,1),L,replace=TRUE, chi="createNEG")
M[ind]<-segno
cat(file=fp, "createNEG output:\n\tM = ", M, "\n")
return(M)
}
