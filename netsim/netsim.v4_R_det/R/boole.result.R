###     boole.result  (2007-04-02)
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


boole.result<-function(r,n,M,N,regulated){
#INPUT
# n: current estimate of the variables (the gene expression vector) at time t  (real non negative)
# r: vector of positive and negative integers. r is a binary tree codyfing for regulatory interactions among the regulators of the regulated gene
# M: the connectivity matrix (real non negative numbers) with element M(i,j) different from 0 if an arc goes from j to i and 0 otherwise
# N: structured as M but contains 1, or -1 in  position (i,j) depending if j regulates i in a positive or negative fashion
# regulated: the index of the regulated gene being considered

#OUTPUT
# a vector of two elements (real) codifying for the sign and the value of the regulation

 cat(file=fp, "\n***boole_result***\ninput:\n\tr = ", r, "\n\tn = ", n, "\n\tM = ", M, "\n\tN = ", N, "\n\tregulated = ", regulated, "\n")
l<-length(r)
signr<-valr<-rep(0,l)
if (is.nan(r[1])) {signr[1]<-1; valr[1]<-0 }
else
 {for(i in 1:l)
   {if(r[i]>0) {expon<-1/abs(M[regulated,r[i]]); valr[i]<-n[r[i]]^expon; signr[i]<-sign(N[regulated,r[i]])}}
  while (next.op(r)!=-1){
	i<-next.op(r)
	if(r[i]==-3){	#OR
		if (signr[2*i]==signr[2*i+1]) {valr[i]<-min(1,sum(valr[2*i],valr[2*i+1])); signr[i]<-signr[2*i]}
		 else  {valr[i]<-min(1,1+signr[2*i]*valr[2*i]+signr[2*i+1]*valr[2*i+1]); signr[i]<-1}
  	r[i]<-0
		r[2*i]<-(-1)
		r[2*i+1]<-(-1)
	}
	if(r[i]==-2){	#AND
		if (signr[2*i]==signr[2*i+1]) {valr[i]<-min(valr[2*i],valr[2*i+1]); signr[i]<-signr[2*i]}
		 else  {valr[i]<-max(0,signr[2*i]*valr[2*i]+signr[2*i+1]*valr[2*i+1]); signr[i]<-1}
		r[i]<-0
		r[2*i]<-(-1)
		r[2*i+1]<-(-1)
	}
   }#end while
  }
cat(file=fp, "boole_result output:\n\tc(signr[1],valr[1]) = ", c(signr[1],valr[1]), "\n")
return(c(signr[1],valr[1]))
}

