###     probmodund  (2009-12-05)
###
###	Copyright (C) 2006-2009  Di Camillo Barbara
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

probmodund<-function(M,h,Sout,STout,Freq.out,toll)
{
 cat(file=fp, "\n***probmod_und***\ninput:\n\tM = ", M, "\n\th = ", h, "\n\tSout = ", Sout,"\n\tSTout = ", STout, "\n\tFreq_out = ", Freq.out, "\n\ttoll = ", toll, "\n")
  nh<-length(h)
  ng<-dim(M)[1]
  Sc<-checkOUT<-matrix(0,ng,nh)

  M.out<-apply(M,2,sum)
  memory<-matrix(NA,ng,3)


  Pm<-1
  lab<-"NO"
  i<-0

  while (i<ng)
   {i<-i+1
    n.out<-M.out[i]
    ind1<-which(memory[,2]==n.out)
    if (length(ind1)>0) {ind<-memory[ind1,3]; Sc[i,]<-Sc[ind,]; checkOUT[i,]<-checkOUT[ind,]}
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
       Sc[i,]<-S.out
      }
   }#end while (i<ng)

  if (!is.na(Pm))
   {Pm<-sum(Sc)/ng
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
cat(file=fp, "probmod_und output:\n\tlist(Pm,M,Sc,lab) = ", Pm, "\n", M,"\n", Sc, "\n", lab, "\n")
return(list(Pm,M,Sc,lab))

}





