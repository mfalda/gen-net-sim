###     Score  (2007-12-05)
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

Score<-function(S,ST,Freq,n,toll)

## INPUT
# S:  integer expressing the in or the out degree of a node i in the module being considered
# ST:  vector of real numbers that, in each position k, contains the expected number of nodes with in-degree or out-degree k (NA if not constrained)
# Freq: vector of real numbers that, in each position k, contains the actual number of nodes with in-degree or out-degree k  (at current iteration)
# n : integer expressing the in or the out degree of a node j in the network
# toll:   vector of real numbers that, in each position k, express the absolute tollerance on ST[k]

## OUTPUT
# Sc: vector of real numbers expressing the scores S_ij (j=1, ..., num) with num = number of nodes in the module being considered
{
 cat(file=fp, "\n***score***\ninput:\n\tS = ", S,"\n\tST = ", ST, "\n\tFreq = ", Freq, "\n\tn = ", n, "\n\ttoll = ", toll, "\n")
 S.new<-S+n
 T1<-ST[S+1]
 T2<-ST[S.new+1]
 OLD1<-Freq[S+1]
 OLD2<-Freq[S.new+1]
 NEW1<-OLD1-1
 NEW2<-OLD2+1
 toll1<-toll[S+1]
 toll2<-toll[S.new+1]
 a<-abs(OLD1-T1)
 b<-abs(NEW1-T1)
 ind<-which(is.na(T1))
 if (length(ind)<length(a))
  {a[ind]<-0
   b[ind]<-0
   m<-apply(cbind(a,b),1,max)
   S1<-(sign(a-b))*m/T1

   ind1<-which((a-toll1)<0)
   ind2<-which((b-toll1)>0)
   indinf<-intersect(ind1,ind2)
   S1[indinf]<-(-Inf)

   S1[ind]<-0
  }
 else S1<-rep(0,length(a))


 a<-abs(OLD2-T2)
 b<-abs(NEW2-T2)
 m<-apply(cbind(a,b),1,max)
 ind<-which(is.na(T2))
 if (length(ind)<length(a))
  {a[ind]<-0
   b[ind]<-0
   m<-apply(cbind(a,b),1,max)
   S2<-(sign(a-b))*m/T2

   ind1<-which((a-toll2)<0)
   ind2<-which((b-toll2)>0)
   indinf<-intersect(ind1,ind2)
   S2[indinf]<-(-Inf)

   S2[ind]<-0
  }
 else S2<-rep(0,length(a))

 ind<-which(T2==0)
 S2[ind]<-(-Inf)

 Sc<-(S1+S2)/2
cat(file=fp, "scoremodular output:\n\tSc = ", Sc, "\n")
 return(Sc)
}
