###     triangola  (2006-12-12)
###
###	Copyright (C) 2006  Di Camillo Barbara
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
triangola<-function(M,Cg,max.con)
#this function add random links to modules so to have an average clustering coefficient of the module equal to Cg

## INPUT
# M:  M is the connectivity matrix of the module, with element (i,j)=1 if an arc goes from j to i and 0 otherwise
# Cg:  desired average clustering coeff. of a module (real non negative)
# max.con:  maximum number of regulators (i.e. incoming edges) a node can have  (integer)

## OUTPUT
# M modified so that the module has, in average, a clustering coefficient  to Cg

{
 cat(file=fp, "\n***triangola***\ninput:\n\tM = ", M,"\n\tCg = ", Cg, "\n\tmax_con = ", max.con, "\n")
 Ng<-dim(M)[1]
 Dmem<-diag(M)
 diag(M)<-0
 for (i in (1:Ng))
    {ind<-union(which(M[i,]!=0),which(M[,i]!=0))  #neighbours
     S.i<-length(ind)
     if (S.i>=2)
	{ng<-Cg*S.i*(S.i-1)/2  #average number of links among neighbours

	 ng.s<-round(rnorm_s(1,ng,1, chi="triangola"))

	 if (ng.s>0)
	   {L<-length(ind)
	    conta<-0
	    coord<-matrix(0,1,2)
            k<-0
	    for (j in (1:(L-1)))
	     {for (h in ((j+1):L))
		{aus1<-M[ind[j],ind[h]]
                 aus2<-M[ind[h],ind[j]]
		 if ( (aus1==1)|(aus2==1) ) conta<-conta+1
                  else coord<-rbind(coord,c(ind[j],ind[h]),c(ind[h],ind[j]))
		}
	     }

	    ng.s<-ng.s-conta

            Lr<-ng.s
            while (Lr>0)
             {M.in<-apply(M,1,sum)
	      ind<-which(M.in<max.con)
              if (length(ind)>0)
		{ind.aus<-which(coord[,1]%in%ind)
		 k<-length(ind.aus)
		 if (k>0)
		   {coord<-matrix(coord[ind.aus,],ncol=2)
		    sk<-sampleB(seq(1,k,1),1)
                    coord.ind<-coord[sk,]
		    M[coord.ind[1],coord.ind[2]]<-1
                    ind<-which((coord[, 1] != coord.ind[1])&(coord[, 1] != coord.ind[2]))
		    if (length(ind)>0) {coord<-matrix(coord[ind,],ncol=2); Lr<-Lr-1}
		     else Lr<-0
		   }
                  else Lr<-0
		}
	       else Lr<-0

             }

	   }
	}
    }
diag(M)<-Dmem
cat(file=fp, "triangola output:\n\tM = ", M, "\n")
return(M)
}



