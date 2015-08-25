###     hubs  (2009-May-12th)
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

hubs<-function(mod)
{#hubs defined based on out-degree and ability to reach other nodes
  cat(file=fp, "\n***hubs***\ninput:\n\tmod = ", mod, "\n")
 S<-apply(mod,2,sum)
 ind<-which(S==max(S))
 pl<-pathlength(mod)
 if ( min(diag(pl))!=Inf )  feedback<-TRUE
 else feedback<-FALSE

 if (length(ind)==1) H<-ind
 else
   {diag(pl)<-0
    P<-apply(pl,2,sum)[ind]
    if (min(P)!=Inf) {ind2<-which(P==min(P)); H<-ind[ind2]}
    else H<- ind
   }
 #hubs defined based on global-degree and ability to reach other nodes
 S<-S+apply(mod,1,sum)
 ind<-which(S==max(S))
 Hio<-ind
 cat(file=fp, "hubs output:\n\tlist(H,feedback,Hio) = ", H, "\n", feedback, "\n", Hio, "\n")
 return(list(H,feedback,Hio))
}
