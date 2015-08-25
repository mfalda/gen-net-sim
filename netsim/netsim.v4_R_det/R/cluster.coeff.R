###     cluster.coeff  (2006-12-12)
###
###	Copyright (C) 2006  Di Camillo Barbara
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


cluster.coeff<-function(W)

## INPUT
# W (real) is the connectivity matrix with element (i,j) different from 0 if an arc goes from j to i and 0 otherwise

## OUTPUT
# a list of   1) the vector of the clustering coefficient of each node (real, non negative)
#             2) the average clustering coefficient (real)


{
 cat(file=fp, "\n***cluster_coeff***\ninput:\n\tW = ", W, "\n")
W<-W+t(W)
ind<-which(W!=0,arr.ind=TRUE)
W[ind]<-1  #matrice simmetrica con solo 0 e 1
diag(W)<-0
N<-dim(W)[1]
Cg<-rep(0,N)
for (i in (1:N))
 {Kg<-sum(W[i,])
  neighbours<-which(W[i,]!=0)
  ng<-0
  L<-length(neighbours)
  if (L>1)
    {for (j in (1:L))
      {ind.j<-neighbours[j]
       for (h in (j:L))
	{ind.h<-neighbours[h]
         if (W[ind.j,ind.h]==1) ng<-ng+1
	}
      }
     Cg[i]<-2*ng/(Kg*(Kg-1))
    }
 }

coeff<-mean(Cg,na.rm=TRUE)
cat(file=fp, "cluster_coeff output:\n\tlist(coeff,Cg) = ", coeff, "\n", Cg, "\n")
return(list(coeff,Cg))
}