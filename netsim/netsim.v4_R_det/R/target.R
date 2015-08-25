###     target  (2007-04-02)
###
###	Copyright (C) 2006, 2007  Di Camillo Barbara
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

target<-function(n,R,M,N,EXT.IN,EXT.FUN,t)

## INPUT
#n: current estimate of the variables (the gene expression vector) at time t  (real non negative)
# R: list of vectors of positive and negative integers. Each element of Rules is a binary tree codyfing for regulatory interactions
# M: the connectivity matrix (real non negative numbers) with element M(i,j) different from 0 if an arc goes from j to i and 0 otherwise
# N: structured as M but contains 1, or -1 in  position (i,j) depending if j regulates i in a positive or negative fashion
# EXT.IN:  (matrix of real numbers) Connectivity matrix containing parameters of sign and efficiecy of regulation of external inputs on genes. It must have N rows, with N=number of nodes in the network and a number of coloumns equal to the number of external inputs. If equal to NA (default), no external input effect is considered.
# EXT.FUN:  list of functions. (see details in simulateprofiles help in the directory "man"). List of functions describing the external input signal as a function of time. If equal to an empty list (default), no external input effect is considered.
#t: considered time point


## OUTPUT
# target value T of gene expression vector at time t  (real non negative)

{
 cat(file=fp, "\n***target***\ninput:\n\tn = ", n, "\n\tR = ")
 print_list(R)
 cat(file=fp, "\tM = ", M, "\n\tN = ", N, "\n\tt = ", t, "\n")
Ln<-length(n)
 T<-numeric(Ln)  #== rep(0,Ln)
 if (length(EXT.FUN)==0)
  {for(i in 1:Ln)
    {rule<-R[[i]]
	   aus<-boole.result(rule,n,M,N,i)
    if (aus[1]==(-1)) T[i]<-1-aus[2]
     else    T[i]<-aus[2]
    }
   }

 else
  {ind.reg<-which(apply(abs(EXT.IN),1,sum)!=0)
   for(i in 1:Ln)
     {rule<-R[[i]]
      aus<-boole.result(rule,n,M,N,i)
      if (i %in% ind.reg)
        {ind<-which(EXT.IN[i,]!=0)
         L<-length(ind)
         val<-0
         for (j in (1:L))
          {expon<-1/abs(EXT.IN[i,ind[j]])
           valr<-EXT.FUN[[ind[j]]](t)^expon
           signr<-sign(EXT.IN[i,ind[j]])
           val<-val+valr*signr
          }
         val<-val+aus[1]*aus[2]
         if  (val<0) T[i]<-max(0,(1+val))
          else T[i]<-min(1,val)
        }
       else {if (aus[1]==(-1)) T[i]<-1-aus[2]
              else    T[i]<-aus[2]
            }
      }
  }
  cat(file=fp, "target output:\n\tT = ", T, "\n")
return(T)
}