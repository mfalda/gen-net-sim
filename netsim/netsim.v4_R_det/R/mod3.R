###     MOD3  (2007-12-05)
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


MOD3<-function(Ng.out,Ng.in,STin,STout,max.con,Cg=0)

# MOD3 creates a basic module of type 3    (Figure 1 in the paper)

## INPUT
# Ng.out: number of nodes in the module with role indicated by x in figure 1   (integer>=1)
# Ng.in: number of nodes in the module with role indicated by y in figure 1   (integer>=1)
# STin:  vector of real numbers that, in each position k, contains the expected number of nodes with in-degree k (NA if not constrained)
# STout: vector of real numbers that, in each position k, contains the expected number of nodes with out-degree k (NA if not constrained)
# max.con: maximum number of regulators a node can have
# Cg: required average clustering coeff. for the module  (real, non negative, <1)

## OUTPUT
# m: m is the connectivity matrix of the module, with element (i,j)=1 if an arc goes from j to i and 0 otherwise

{
 cat(file=fp, "\n***mod3***\ninput:\n\tNg_out = ", Ng.out, "\n\tNg_in = ", Ng.in, "\n\tSTin = ", STin, "\n\tSTout = ", STout, "\n\tmax_con = ", max.con, "\n\tCg = ", Cg, "\n")
m<-matrix(0,Ng.out+Ng.in,Ng.out+Ng.in)
N<-Ng.in+Ng.out
S<-sum(STout[2:(Ng.in+1)],na.rm=TRUE)
if (S>0) p<-STout[2:(Ng.in+1)]/S
else p<-rep(1/Ng.in,Ng.in)

M.out<-sampleB(seq(1,Ng.in,1),Ng.out,prob=p,replace=TRUE)

while (sum(M.out)<(N-1))
  {if (S>0)
	{STout<-(STout[1:(Ng.in+1)]/S)*Ng.in
         aus<-max(M.out)+1
   	 Freq.out<-hist(M.out,breaks=seq(0,Ng.in+1,1),right=FALSE,plot=FALSE)$counts/N
         Sc<-Score(S=M.out,ST=STout,Freq=Freq.out,n=1,toll=rep(Inf,(Ng.out+1)))
         indok<-which(Sc!=-Inf)
	 indInf <- setdiff(seq(1, Ng.out, 1), indok)
	 Sc[indInf]<-min(c(0,Sc[indok]))-1
         p<-Sc/sum(Sc)
	}
   else p<-rep(1/Ng.in,Ng.in)

   ind.M<-sampleB(seq(1,Ng.out,1),1,prob=p)

   M.out[ind.M]<-M.out[ind.M]+1
  }

for (j in (1:Ng.out))
  {n.reg<-M.out[j]
   if (j==1)
     {m[1:n.reg,(j+Ng.in)]<-1
      indS<-1:n.reg
      indBS<-setdiff(1:Ng.in,indS)
     }
    else
     {L<-length(indBS)
      Ls<-min(L,n.reg-1)

      ns<-n.reg-Ls
      S<-sum(STin[1:(Ng.out+1)],na.rm=TRUE)
      Sin<-apply(m,1,sum)
      if (S>0)
	{
    STin<-(STin[1:(Ng.out+1)]/S)*Ng.out
         aus<-max(Sin)+1
         Freq.in<-hist(Sin,breaks=seq(0,aus,1),right=FALSE,plot=FALSE)$counts/N
         Sc<-Score(S=Sin[indS],ST=STin,Freq=Freq.in,n=1,toll=rep(Inf,(Ng.out+1)))
         indok<-which(Sc!=-Inf)
	 indInf <- setdiff(seq(1, length(indS), 1), indok)
	 Sc[indInf]<-min(c(0,Sc[indok]))-1
         p.sc<-Sc/sum(Sc)
	}
       else p.sc<-rep(1/length(indS),length(indS))

       ind<-which(Sin[indS]==max.con)
       p.sc[ind]<-0
       p.sc<-p.sc/sum(p.sc)

      ind1<-sampleB(indS,ns,prob=p.sc)
      if (Ls>0)
	{ind<-indBS[1:Ls]
         indS<-c(indS,ind)
         indBS<-setdiff(1:Ng.in,indS)
         ind1<-c(ind1,ind)
	}
      m[ind1,(j+Ng.in)]<-1
     }
  }


if ( (Cg>0)&((Ng.in+Ng.out)>2) ) m<-triangola(M=m,Cg=Cg,max.con)
cat(file=fp, "mod3 output:\n\tm = ", m, "\n")
return(m)
}
