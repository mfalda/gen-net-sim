###     createMOD  (2009-May-12th)
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


createMOD<-function(m=4,auto=FALSE)
#\arguments{
#  \item{m}{maximum size of constitutive modules}
#  \item{auto}{if auto=FALSE modules with autoregulation are not allowed}
#  }

{
##### CHECK PARAMETERS ####
        if (!is.numeric(m))
        stop("'m' must be numeric")
        if ((m<2)|(m>5))  stop("'m' must be greater than 1 and lower than 6")

        if ((auto!=TRUE)&(auto!=FALSE))
        stop("'auto' must be TRUE or FALSE")
##### END CHECK PARAMETERS ####
  cat(file=fp, "\n***createMOD***\ninput:\n\tm = ", m, "\n\tauto = ", auto, "\n")
 MOD<-list()
 b<-1
 path<-"modules/"
 for (i in (2:m))
  {if (auto==FALSE) etich<-paste(.libPaths(), "/HMMb/", path,"M",i,".txt",sep="")
   else  etich<-paste(.libPaths(), "/HMMb/autoreg_",path,"M",i,".txt",sep="")
   M<-as.matrix(read.table(etich))
   codes<-as.vector(M[1,])
   M<-M[-1,]
   L<-dim(M)[2]
   for (j in (1:L))
      { mod<-matrix(M[,j],i,i)
        aus<-length(which((mod==t(mod))==0))
        if (aus==0) simm<-TRUE
        else simm<-FALSE
        H<-hubs(mod)
        MOD[[b]]<-list("code"=codes[j],"net"=mod,"hubs"=H[[1]],"CC"=cluster.coeff(mod)[[1]],"autoreg"=auto,"feedback"=H[[2]],"hubsio"=H[[3]],"SIMM"=simm,"dim.m"=i)
        b<-b+1
      }
  }

return(MOD)

}
