###     pathlength  (2009-May-12th)
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

pathlength<-function(W,und=FALSE)
#W matrix with (i,j)!=0 if j regulates i
#und=FALSE if the network is directed
{
 cat(file=fp, "\n***pathlength***\ninput:\n\tW = ", W, "\n\tund = ", und, "\n")
 N<-dim(W)[1]
 if (und==TRUE)  W<-W+t(W)
 W[which(W!=0,arr.ind=TRUE)]<-1
 M<-matrix(NA,N,N)
 dist.tot<-vector()


 for (i in (1:N))
  {distanza<-rep(NA,N)
   #distanza[i]<-0
   visitati<-vector()  #visitati<-i
   conta<-1
   neighbours<-which(W[,i]!=0)#setdiff(which(W[,i]!=0),i)
   L<-length(neighbours)
   while (L>0)
    { distanza[neighbours]<-conta
      #distanza[i]<-0
      visitati<-union(visitati,neighbours)
      conta<-conta+1
      new.neighbours<-vector()
      for (j in (1:L))
         {ind<-neighbours[j]
	        new.neighbours<-c(new.neighbours,which(W[,ind]!=0))
         }
      new.neighbours<-union(new.neighbours,new.neighbours)
      new.neighbours<-setdiff(new.neighbours,visitati)
      neighbours<-new.neighbours
      L<-length(neighbours)
    }

  M[,i]<-distanza

 }

if (und==TRUE) diag(M)<-0

M[which(is.na(M),arr.ind=TRUE)]<-Inf
cat(file=fp, "pathlength output:\n\tM = ", M, "\n")
return(M)

}

