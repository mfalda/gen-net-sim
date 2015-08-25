###     MOD1  (2007-12-05)
###
###	Copyright (C) 2007  Di Camillo Barbara
###
###	This program is free software; you can redistribute it and/or
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

MOD1<-function(Ng, Cg=0,max.con)

# MOD1 creates a basic module of type 1    (Figure 1 in the paper)

## INPUT
# Ng: number of nodes in the module    (integer>=1)
# Cg: required average clustering coeff. for the module  (real, non negative, <1)
# max.con: maximum number of regulators a node can have

## OUTPUT
# m: m is the connectivity matrix of the module, with element (i,j)=1 if an arc goes from j to i and 0 otherwise

{
 cat(file=fp, "\n***mod1***\ninput:\n\tNg = ", Ng, "\n\tCg = ", Cg, "\n\tmax_con = ", max.con, "\n")
 m<-matrix(0,Ng,Ng)
 if (Ng>1)
       {x<-seq(1,Ng-1,1)
	ind<-cbind(x,x+1)
	m[ind]<-1
	m[Ng,1]<-1
	s<-sampleB(c(1,-1),1)
	if (s==1) {m[Ng,1]<-0; m[1,Ng]<-1}
       }
else   m[1,1]<-1

if ((Cg>0)&(Ng>2)) m<-triangola(M=m,Cg=Cg,max.con)
cat(file=fp, "mod1 output:\n\tm = ", m, "\n")
return(m)
}
