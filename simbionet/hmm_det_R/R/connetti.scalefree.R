###     connetti.scalefree  (2009-04-02)
###
###	Copyright (C) 2006, 2009  Di Camillo Barbara
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


connetti.scalefree<-function(Mdiscr,STout,STin,dist,toll,max.con,und=FALSE)
{
cat(file=fp, "\n***connetti.scalefree***\ninput:\n\tMdiscr = ", Mdiscr,"\n\tSTout = ", STout, "\n\tSTin = ", STin, "\n\tdist = ", dist, "\n\ttoll = ", toll, "\n\tmax.con = ", max.con, "\n\tund = ", und, "\n")
non.connessi<-which(dist==Inf)
  connessi<-which(dist!=Inf)
  N<-length(dist)

  Sout<-apply(Mdiscr,2,sum)
  m<-max(Sout)+1
  Freq.out<-hist(Sout,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts
  Sc.out<-Score(S=Sout,ST=STout,Freq=Freq.out,n=1,toll)
  if (und==TRUE) Sc.out[which(Sout==max.con)]<-(-Inf)

  Sin<-apply(Mdiscr,1,sum)
  m<-max(Sin)+1
  Freq.in<-hist(Sin,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts
  Sc.in<-Score(S=Sin,ST=STin,Freq=Freq.in,n=1,toll)
  Sc.in[which(Sin==max.con)]<-(-Inf)

  co.in<-intersect(connessi,which(Sc.in!=(-Inf)))
  nonco.in<-intersect(non.connessi,which(Sc.in!=(-Inf)))
  co.out<-intersect(connessi,which(Sc.out!=(-Inf)))
  nonco.out<-intersect(non.connessi,which(Sc.out!=(-Inf)))
  while ( ((length(nonco.in)==0)|(length(co.out)==0))&((length(co.in)==0)|(length(nonco.out)==0)) )
    {cat("\n WARNING: tolerance parameters were relaxed to allow connectivity of the graph\n")
     toll<-toll+1
     cat(file=fp, "linea = 50", "\n")
     Sc.out<-Score(S=Sout,ST=STout,Freq=Freq.out,n=1,toll=rep(Inf,N))
    cat(file=fp, "linea = 52", "\n")
     Sc.in<-Score(S=Sin,ST=STin,Freq=Freq.in,n=1,toll=rep(Inf,N))
     co.in<-intersect(connessi,which(Sc.in!=(-Inf)))
     nonco.in<-intersect(non.connessi,which(Sc.in!=(-Inf)))
     co.out<-intersect(connessi,which(Sc.out!=(-Inf)))
     nonco.out<-intersect(non.connessi,which(Sc.out!=(-Inf)))
    }

  choice<-0
  cat(file=fp, "Sc.in = ", Sc.in, "\n")
  cat(file=fp, "Sc.out = ", Sc.out, "\n")
  cat(file=fp, "co.in = ", co.in, "\n")
  cat(file=fp, "nonco.in = ", nonco.in, "\n")
  cat(file=fp, "co.out = ", co.out, "\n")
  cat(file=fp, "nonco.out = ", nonco.out, "\n")
  Sc1.in<-Sc.in[nonco.in]
  Sc1.out<-Sc.out[co.out]
  Sc2.in<-Sc.in[co.in]
  Sc2.out<-Sc.out[nonco.out]
  if ( (length(nonco.in)!=0)&(length(co.out)!=0)&(length(co.in)!=0)&(length(nonco.out)!=0) )
	{p1<-mean(Sc1.in)/2+mean(Sc1.out)/2
         p2<-mean(Sc2.in)/2+mean(Sc2.out)/2
         p<-c(p1,p2)
         if (length(which(p<=0))>0) p<-p-min(p)+1/(N^2)
	 p<-p/sum(p)
         choice<-sampleB(c(1,2),1,prob=p)
	}

  if ( ((length(nonco.in)==0)|(length(co.out)==0)) | (choice==2) )
	{if (length(which(Sc2.in<=0))>0) Sc2.in<-Sc2.in-min(Sc2.in)+1/(N^2)
	 p.in<-Sc2.in/sum(Sc2.in)
	 ind.in<-sampleB(co.in,1,prob=p.in)
         if (length(which(Sc2.out<=0))>0) Sc2.out<-Sc2.out-min(Sc2.out)+1/(N^2)
    cat(file=fp, "Sc2.out = ", Sc2.out, "\n")
    cat(file=fp, "nonco.out = ", nonco.out, "\n")
	 p.out<-Sc2.out/sum(Sc2.out)
	 ind.out<-sampleB(nonco.out,1,prob=p.out)
	}
  if ( ((length(co.in)==0)|(length(nonco.out)==0)) | (choice==1) )
	{if (length(which(Sc1.in<=0))>0) Sc1.in<-Sc1.in-min(Sc1.in)+1/(N^2)
	 p.in<-Sc1.in/sum(Sc1.in)
	 ind.in<-sampleB(nonco.in,1,prob=p.in)
         if (length(which(Sc1.out<=0))>0) Sc1.out<-Sc1.out-min(Sc1.out)+1/(N^2)
	 p.out<-Sc1.out/sum(Sc1.out)
	 ind.out<-sampleB(co.out,1,prob=p.out)
	}

  Mdiscr[ind.in,ind.out]<-1
  if (und==TRUE)  Mdiscr[ind.out,ind.in]<-1
  cat(file=fp, "connetti.scalefree output:\n\tMdiscr = ", Mdiscr, "\n")
 return(Mdiscr)
}