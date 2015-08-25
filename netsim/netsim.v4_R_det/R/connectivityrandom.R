###	(2006-12-12)
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


connectivityrandom<-function(N=50,max.con=12,k=3,weight.mean=1, weight.sd=0.1)

{
 cat(file=fp, "\n***connectivity_random***\ninput:\n\tN = ", N, "\n\tmax_con = ", max.con, "\n\tk = ", k, "\n\tweight_mean = ", weight.mean, "\n\tweight_sd = ", weight.sd, "\n")
# this function create a directed graph. It gives as output a list of two matrixes M and Mdiscr.
# M is the connectivity matrix with element (i,j)=1 if an arc goes from j to i and 0 otherwise
# (M is the transpose of the standard adjacence matrix)
# Mdiscr is structured as M but contains weights (i,j) belonging to real positive numbers

#\arguments{
#  \item{N}{(integer) Number of genes in the network.}
#  \item{max.con}{(integer) Maximum number of regulators that each node (gene) in the network can have. Default=5}
#  \item{k}{(real) average number of in-degree and out-degree connections.}
#  \item{weight.mean}{(real) Mean of the Gaussian distribution from which regulatory weigths are sampled. Default=1 (See References).}
#  \item{weight.sd}{(real) Standard deviation of the Gaussian distribution from which regulatory weigths are sampled. Default=0.1 (See References).}
#  }

 if (max.con>N) max.con<-N
 Mdiscr<-M<-matrix(0,ncol=N,nrow=N)
 p<-k/N
 Pnum<-rep(0,max.con)
 for (j in (1:max.con))
   {a<-N-j+1
    L<-N-a
    F<-1
    if (L>0) {for(i in (0:L)) F<-F*(a+i)}
    Pnum[j]<-(F/factorial(j))*(p^j)*((1-p)^(N-j))
   }

 for (i in (1:N))
  {num<-sampleB(seq(1,max.con,1),1,prob=Pnum)
   ind<-sampleB(seq(1,N,1),num)
   Mdiscr[i,ind]<-1
  }

#######assign weights to the connectivity matrix##########
ind<-which(Mdiscr==1,arr.ind=TRUE)
L<-dim(ind)[1]
aus<-abs(rnorm_s(L,weight.mean, weight.sd, chi="cr"))
M[ind]<-aus
cat(file=fp, "cr output:\n\tlist(M,Mdiscr) = ", M, "\n", Mdiscr, "\n")
return(list(M,Mdiscr))
}
