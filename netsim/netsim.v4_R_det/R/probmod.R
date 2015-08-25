###     probmod  (2007-12-05)
###
###	Copyright (C) 2007  Di Camillo Barbara
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

probmod<-function(M,h,Sin,Sout,STin,STout,Freq.in,Freq.out,toll)
## INPUT
# M: the connectivity matrix of the considered module, with element (i,j)=1 if an arc goes from j to i and 0 otherwise
# h: nodes still to be connected at current iteration  (vector of integers)
# h.new: nodes to be connected at the following iteration  (vector of integers)
# Sin:  vector of integers (sum of the elements on the rows of the connectivity matrix at current iteration )
# Sout: vector of integers (sum of the elements on the coloumns of the connectivity matrix at current iteration)
# STin:  vector of real numbers that, in each position k, contains the expected number of nodes with in-degree k (NA if not constrained)
# STout: vector of real numbers that, in each position k, contains the expected number of nodes with out-degree k (NA if not constrained)
# Freq.in: vector of real numbers that, in each position k, contains the actual number of nodes with in-degree k  (at current iteration)
# Freq.out: vector of real numbers that, in each position k, contains the actual number of nodes with out-degree k  (at current iteration)
# toll:   vector of real numbers that, in each position k, express the absolute tollerance on STin[k] and STout[k]

## OUTPUT
# a list of 4 elements:
# 1)  score proportional to the probability of the module being analized to be selected as a network subgraph at current iteration
# 2)  connectivity matrix of the module, with element (i,j)=1 if an arc goes from j to i and 0 otherwise
# 3)  matrix of scores Sij (real). See the paper for details
# 4)  label(character) indicating if the parameter max.con  (maximum in-degree) has to be increased to allow convergence

{
 cat(file=fp, "\n***probmod***\ninput:\n\tM = ", M, "\n\th = ", h, "\n\tSin = ", Sin, "\n\tSout = ", Sout,"\n\tSTin = ", STin, "\n\tSTout = ", STout, "\n\tFreq_in = ", Freq.in, "\n\tFreq_out = ", Freq.out, "\n\ttoll = ", toll, "\n")
  nh<-length(h)
  ng<-dim(M)[1]
  Sc<-checkIN<-checkOUT<-matrix(0,ng,nh)

  M.in<-apply(M,1,sum)
  M.out<-apply(M,2,sum)
  memory<-matrix(NA,ng,3)


  Pm<-1
  lab<-"NO"
  i<-0

  while (i<ng)
   {i<-i+1
    n.out<-M.out[i]
    n.in<-M.in[i]
    ind1<-intersect(which(memory[,1]==n.in),which(memory[,2]==n.out))
    if (length(ind1)>0) {ind<-memory[ind1,3]; Sc[i,]<-Sc[ind,]; checkIN[i,]<-checkIN[ind,]; checkOUT[i,]<-checkOUT[ind,]}
     else
      {if (n.out!=0)
	{S.out<-Score(S=Sout[h],ST=STout,Freq=Freq.out,n=n.out,toll=toll)
	 indok<-which(S.out!=-Inf)
          checkOUT[i,indok]<-1
          if (length(indok)<1) {Pm<-NA; i<-ng}
          indInf <- setdiff(seq(1, nh, 1), indok)
	  S.out[indInf]<-min(c(0,S.out[indok]))-1
 	}
       else {S.out<-rep(0,nh);checkOUT[i,]<-rep(1,nh)}

       if (n.in!=0)
	 {S.in<-Score(S=Sin[h],ST=STin,Freq=Freq.in,n=n.in,toll=toll)
	  indok<-which(S.in!=-Inf)
          checkIN[i,indok]<-1
          if (length(indok)<1) {Pm<-NA; i<-ng; lab<-"in"}
          indInf <- setdiff(seq(1, nh, 1), indok)
	  S.in[indInf]<-min(c(0,S.in[indok]))-1
      	 }
       else {S.in<-rep(0,nh); checkIN[i,]<-1}
       Sc[i,]<-S.out+S.in

      }

    }#end while (i<ng)


  if (!is.na(Pm))
   {Pm<-sum(Sc)/ng

    #check in
    aus<-apply(checkIN,2,sum)
    L<-length(which(aus!=0))
    if (L<ng) {Pm<-NA; i<-ng; lab<-"in"}
     else
      {rs<-apply(checkIN,1,sum)
       ord.ind<-order(rs)
       I<-vector()
       j<-0
       while (j<ng)
   	{j<-j+1
	 ind<-ord.ind[j]
	 I.add<-which(checkIN[ind,]!=0)
         I<-union(I,I.add)
         if (length(I)<j) {Pm<-NA; j<-ng}
	}
      }


    #check out
    checkOUT[which(checkOUT<0)]<-0
    checkOUT[which(checkOUT>0)]<-1
    aus<-apply(checkOUT,2,sum)
    L<-length(which(aus!=0))
    if (L<ng) {Pm<-NA; i<-ng}
     else
      {rs<-apply(checkOUT,1,sum)
       ord.ind<-order(rs)
       I<-vector()
       j<-0
       while (j<ng)
   	{j<-j+1
	 ind<-ord.ind[j]
	 I.add<-which(checkOUT[ind,]!=0)
         I<-union(I,I.add)
         if (length(I)<j) {Pm<-NA; j<-ng}
	}
      }



   }# end if (!is.na(Pm))
cat(file=fp, "probmod output:\n\tlist(Pm,M,Sc,lab) = ", Pm, "\n", M,"\n", Sc, "\n", lab, "\n")
return(list(Pm,M,Sc,lab))

 }





