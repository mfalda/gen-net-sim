###     assign.nodes.und  (2009-12-05)
###
###	Copyright (C) 2006-2009  Di Camillo Barbara
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

assign.nodes.und<-function(M,Mdiscr,h,hubs,Sc,Sout,max.con)
{
cat(file=fp, "\n***assign_nodes_und***\ninput:\n\tM = ", M, "\n\tMdiscr = ", Mdiscr, "\n\th = ", h, "\n\thubs = ", hubs, "\n\tSc = ", Sc, "\n\tSout = ", Sout, "\n\tmax_con = ", max.con, "\n")
 nh<-length(h)
 ng<-dim(M)[1]
 or.h<-h
 aus.h<-seq(1,nh,1)
 new.hubs<-vector()
 index<-rep(0,ng)
 M.out<-apply(M,2,sum)
 ORD<-order(M.out,decreasing=TRUE)
 for (j in (1:ng))
  {i<-ORD[j]
   p<-Sc[i,]
   if (length(which(p<=0))>0) p<-p-min(p)+1/(nh^2)
   Sout.h<-Sout[h]
   n<-M.out[i]
   ind<-which(Sout.h>(max.con-n))
   if (length(ind)>0) p[ind]<-0
   p<-p[aus.h]/sum(p[aus.h])
   ind.h<-sampleB(h,size=1,prob=p)
   index[i]<-ind.h
   h<-setdiff(h,ind.h)
   aus.h<-setdiff(aus.h,which(or.h==ind.h))
  }

 for (j in (1:ng))
  {ri<-index[j]
   ind<-which(M[j,]==1)
   co<-index[ind]
   Mdiscr[ri,co]<-Mdiscr[co,ri]<-1
   if (j %in% hubs) new.hubs<-c(new.hubs,ri)
  }

cat(file=fp, "assign_nodes_und output:\n\tlist(Mdiscr,h,new.hubs) = ", Mdiscr, "\n", h, "\n", new.hubs, "\n")
return(list(Mdiscr,h,new.hubs))

}
