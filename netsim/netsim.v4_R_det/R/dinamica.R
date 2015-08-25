###     dinamica  (2007-04-02)
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

dinamica<-function(parms,n0,times){
M<-parms[[1]]
R<-parms[[2]]
N<-parms[[3]]
k<-parms[[4]]
alpha<-parms[[5]]
theta<-parms[[6]]
act.fun<-parms[[7]]
EXT.IN<-parms[[8]]
EXT.FUN<-parms[[9]]
Xmin<-parms[[10]]

res<-signif((log(2)/mean(k))/length(times),1)
S<-max(abs(times/res-round(times/res)))
while (S>1e-8)
   {res<-res/10; S<-max(abs(times/res-round(times/res)))}
passi<-max(times)/res

D<-matrix(n0,length(n0))
for(i in 1:passi){
	n<-D[,i]
	targ<-target(n,R,M,N,EXT.IN,EXT.FUN,res*(i-1))
        if (act.fun=="linear")   targetT<-targ*(1-Xmin)+Xmin
         else  {#norm.f<-(1-theta)^alpha/((1-theta)^alpha+0.5^alpha)
		#aus<-(targ-theta)
		#aus[which(aus<0)]<-0
		#targetT<-(aus^alpha/(aus^alpha+0.5^alpha))/norm.f
		targetT<-(1/(1+exp(-alpha*(targ-theta))))*(1-Xmin)+Xmin
		#aus<-apply(M,1,sum)
		#targetT[which(aus==0)]<-0
	       }

	# Euler Algorithm
	incr<-matrix(res*k*(targetT-n),ncol=1)
	n<-n+incr
	D<-cbind(D,n)
}
L<-length(times)
ind<-rep(0,L)
aus<-seq(0,max(times),res)
for (i in (1:L))
  {ind[i]<-which.min(abs(aus-times[i]))
   if ( (aus[ind[i]]-times[i])>1e-8) stop("error")
  }
D<-D[,ind]

if (length(n0)==1) D<-matrix(D,nrow=1)

return(D)
}
