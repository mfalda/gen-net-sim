###     check.conn  (2009-12-05)
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


check.conn<-function(Mdiscr)
{
cat(file=fp, "\n***check_conn***\ninput:\n\tMdiscr = ", Mdiscr, "\n")
N<-dim(Mdiscr)[1]
 Mt<-t(Mdiscr)
 Maus<-Mdiscr+Mt
 ind<-which(Maus!=0,arr.ind=TRUE)
 Maus[ind]<-1
    
 vert<-1
 colore<-c(1,rep(0,N-1))
 dist<-c(0,rep(Inf,N-1))
 grigi<-1
 while (length(grigi)!=0)
  {u<-grigi[1]
   adj<-which(Maus[u,]!=0) #neighbors of u
   L<-length(adj)
   if (L>0) 
     {for (i in (1:L))
	{v<-adj[i]
	 if (colore[v]==0)
	   {colore[v]<-1
	    dist[v]<-dist[u]+1
	    grigi<-c(grigi,v)
	   }	  
	}
      }
   lg<-length(grigi)
   if (lg>=2) grigi<-grigi[2:lg]
    else grigi<-vector()
  }#end while

cat(file=fp, "check_conn output:\n\tdist = ", dist, "\n")
return(dist)
}
