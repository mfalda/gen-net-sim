###     connectivitygeometric  (2007-12-05)
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


connectivitygeometric<-function(N=50,k=3, weight.mean=1, weight.sd=0.1)

{
# this function create a directed graph. It gives as output a list of two matrixes M and Mdiscr.
# M is the connectivity matrix with element (i,j)=1 if an arc goes from j to i and 0 otherwise
# (M is the transpose of the standard adjacence matrix)
# Mdiscr is structured as M but contains weights (i,j) belonging to real positive numbers


#\arguments{
#  \item{N}{(integer) Number of genes in the network.}
#  \item{k}{(real) average number of in-degree and out-degree connections.}
#  \item{weight.mean}{(real)  Mean of the Gaussian distribution from which regulatory weigths are sampled. Default=1 (See References).}
#  \item{weight.sd}{(real)  Standard deviation of the Gaussian distribution from which regulatory weigths are sampled. Default=0.1 (See References).}
#  }

 cat(file=fp, "\n***connectivity_geometric***\ninput:\n\tN = ", N,"\n\tk = ", k, "\n\tweight_mean = ", weight.mean, "\n\tweight_sd = ", weight.sd, "\n")
 Mdiscr<-matrix(0,ncol=N,nrow=N)
 r<-sqrt(k/(N*pi))*1.1*sqrt(2)
 x<-runif_s(N, chi="cm")
 y<-runif_s(N, chi="cm")
 for (i in (1:N))
  {d<-sqrt((x[i]-x)^2+(y[i]-y)^2)
   ind<-setdiff(which((d<r)&(d>0)),seq(0,(i-1),1))
   s<-sampleB(c(0,1),length(ind),replace=TRUE)
   ind1<-which(s==1)
   if (length(ind1)>0) Mdiscr[i,ind[ind1]]<-1
   ind0<-which(s==0)
   if (length(ind0)>0)  Mdiscr[ind[ind0],i]<-1
  }

###check that every gene has at least 1 regulator
Sr<-apply(Mdiscr,1,sum)
indL<-which(Sr==0)
L<-length(indL)
while (L>0)
 {ri<-indL[1]
  regulatedind<-which(Mdiscr[,ri]==1)
  if (length(regulatedind)==0)
    {#cat("\n k is low with respect to N. It is suggested to increase k")
     x[ri]<-runif_s(1, chi="cm")
     y[ri]<-runif_s(1, chi="cm")
     d<-sqrt((x[ri]-x)^2+(y[ri]-y)^2)
     ind<-which((d<r)&(d>0))
     s<-sampleB(c(0,1),length(ind),replace=TRUE)
     ind1<-which(s==1)
     if (length(ind1)>0) Mdiscr[ri,ind[ind1]]<-1
     ind0<-which(s==0)
     if (length(ind0)>0)  Mdiscr[ind[ind0],ri]<-1
    }
  else {s<-sampleB(regulatedind,1);  Mdiscr[ri,s]<-1}

  Sr<-apply(Mdiscr,1,sum)
  indL<-which(Sr==0)
  L<-length(indL)
 }



#######assign weights to the connectivity matrix##########
ind<-which(Mdiscr==1,arr.ind=TRUE)
L<-dim(ind)[1]
aus<-abs(rnorm_s(L,weight.mean, weight.sd, chi="cm"))
M<-matrix(0,ncol=N,nrow=N)
M[ind]<-aus
cat(file=fp, "cg output:\n\tlist(M,Mdiscr) = ", M, "\n", Mdiscr, "\n")
return(list(M,Mdiscr))
}