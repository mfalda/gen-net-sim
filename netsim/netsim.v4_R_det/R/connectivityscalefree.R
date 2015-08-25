###     connectivitymodular  (2006-12-12)
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

connectivityscalefree<-function(N=50,max.con=12,gamma=2.2,r.tol=0.1,a.tol=1,  weight.mean=1, weight.sd=0.1)
{
 cat(file=fp, "\n***connectivity_scalefree***\ninput:\n\tN = ", N, "\n\tmax_con = ", max.con, "\n\tgamma = ", gamma, "\n\tr_tol = ", r.tol, "\n\ta_tol = ", a.tol, "\n\tweight_mean = ", weight.mean, "\n\tweight_sd = ", weight.sd, "\n")
# this function create a directed graph. It gives as output a list of two matrixes M and Mdiscr.
# M is the connectivity matrix with element (i,j)=1 if an arc goes from j to i and 0 otherwise
# (M is the transpose of the standard adjacence matrix)
# Mdiscr is structured as M but contains weights (i,j) belonging to real positive numbers

 #\arguments{
#  \item{N}{(integer)  Number of nodes in the network.}
#  \item{max.con}{(integer)  Maximum number of incoming edges that each node in the network can have. Default=12}
#  \item{gamma}{(real) Parameter of the power law distribution of the degree of connectivity of the nodes in the graph.Default=2.2  equal to the average value observed by (Albert and Barabasi (2000)) in protein networks.}
#  \item{r.tol}{(real) relative tolerance parameter for the power law distribution. Default=0.1 (See Details).}
#  \item{a.tol}{(real) absolute tolerance parameter for the power law distribution. Default=0.1 (See Details).}
#  \item{weight.mean}{(real) Mean of the Gaussian distribution from which regulatory weigths are sampled. Default=1 (See References).}
#  \item{weight.sd}{(real) Standard deviation of the Gaussian distribution from which regulatory weigths are sampled. Default=0.1 (See References).}
#  }



	##### CHECK PARAMETERS ####
        if (!is.numeric(N))
        stop("`N' must be numeric")

        if (!is.numeric(max.con))
        stop("`max.con' must be numeric")

        if (!is.numeric(gamma))
        stop("`gamma' must be numeric")

        if (!is.numeric(r.tol)) stop("`r.tol' must be numeric")
	else
	 {if ((r.tol>1)|(r.tol<0))
          stop("`r.tol' must be greater than 0 and lower than 1")
         }

        if (!is.numeric(a.tol)) stop("`a.tol' must be numeric")
	else
	 {if (a.tol<0)
          stop("`a.tol' must be a positive number")
         }
       ############### END CHECK PARAMETERS ####

 Mdiscr<-matrix(0,ncol=N,nrow=N)
 if (max.con>N) max.con<-N

 Prob<-c(seq(1,N,1)^(-gamma),0)
 Prob<-Prob/(sum(Prob))
 Freq.out<-Freq.in<-c(N,rep(0,N+1))
 STout<-STin<-c(NA,Prob*N)

 p<-Prob[1:max.con]
 STin[(max.con+2):(N+2)]<-0
 STin[2:(max.con+1)]<-N*p/sum(p)
 STin[1]<-0

 aus<-cbind(STout*r.tol,rep(a.tol,length(STout)))
 toll.out<-apply(aus,1,max)
 toll.out[which(STout==0)]<-0		# When STin==0 toll is set to 0
 aus<-cbind(STin*r.tol,rep(a.tol,length(STin)))
 toll.in<-apply(aus,1,max)
 toll.in[which(STin==0)]<-0		# When STin==0 toll is set to 0

 Mdiscr[1,2]<-1

 inthenet<-1
 not.inthenet<-seq(2,N,1)
 NN<-length(not.inthenet)

 while (NN>0)
   { not.regulated<-give.outlink<-vector()
     for (j in (1:NN))
       {i<-not.inthenet[j]
 	Sin<-apply(Mdiscr,1,sum)
   	m<-max(Sin)+1
   	Freq.in[1:m]<-hist(Sin,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts

   	Sout<-apply(Mdiscr,2,sum)
   	m<-max(Sout)+1
  	Freq.out[1:m]<-hist(Sout,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts

   	mx<-length(inthenet)
   	aus.in<-Freq.in[2:(mx+1)]-STin[2:(mx+1)]-toll.in[2:(mx+1)]
   	aus.in[which(is.na(aus.in))]<-1
   	numposs<-which(aus.in<0)
   	Lp<-min(length(numposs),max.con)

	if (Lp>0)
	  {aus.in<-Freq.in[2:(mx+1)]-STin[2:(mx+1)]+toll.in[2:(mx+1)]
   	   aus.in[which(is.na(aus.in))]<-1
	   aus2<-aus.in
	   primi<-which(aus2<0)
      	   if (length(primi)>0)
	     {num<-which.min(aus2)
	      if (num>max.con) num<-0
	     }
	    else
	     {p<-rep(0,Lp)
      	      for (n in (1:Lp)) p[n]<-Score.sf(S=0,ST=STin,Freq=Freq.in,n=n,toll=toll.in)
              indInf<-which(p==(-Inf))
              p[indInf]<-0
              if (min(p)<=0) p<-p-min(p)*(1.01)
              p[indInf]<-0
              if (sum(p)>0) num<-sampleB(seq(1,Lp,1),1,prob=p/sum(p))
	      else num<-0
   	     }
	  }

	if (num==0) not.regulated<-c(not.regulated,i)

	linked<-vector()
   	aus.give.outlink<-union(give.outlink,give.outlink)
   	mem.o<-vector()
   	while (num>0)
     	  {if (length(aus.give.outlink)>0)
	     {o<-aus.give.outlink[1]
	      aus.give.outlink<-setdiff(aus.give.outlink,o)
	      mem.o<-c(mem.o,o)
	      Mdiscr[i,o]<-1
	      num<-num-1
	      linked<-c(linked,o)
	     }
          else
	     {Sout<-apply(Mdiscr,2,sum)
   	      m<-max(Sout)+1
   	      Freq.out[1:m]<-hist(Sout,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts
	      Sc<-Score.sf(S=Sout[inthenet],ST=STout,Freq=Freq.out,n=1,toll=toll.out)
      	      indInf<-which(Sc==(-Inf))
              Sc[indInf]<-0
              if (min(Sc)<=0) Sc<-Sc-min(Sc)*(1.01)
              Sc[indInf]<-0
              ind<-setdiff(which(Sc>0),linked)
              m.c<-length(ind)

              if (m.c>0)
	   	{p.ind<-Sc[ind]/sum(Sc[ind])
            	 s<-sampleB(ind,1,prob=p.ind)
	    	 Mdiscr[i,s]<-1
            	 num<-num-1
	    	 linked<-c(linked,s)
           	} #end if (m.c>0)

              else #(m.c==0)
	        {aus<-Freq.out[2:(mx+1)]-STout[2:(mx+1)]+toll.out[2:(mx+1)]
   	         aus[which(is.na(aus))]<-1
	         available<-which(aus<0)
		 if (length(available)>0) min.out.i<-min(available)
		 else
		   {aus<-Freq.out[2:(mx+1)]-STout[2:(mx+1)]-toll.out[2:(mx+1)]
   	            aus[which(is.na(aus))]<-1
	            available<-which(aus<0)
		    min.out.i<-min(available)
		   }
		 campione<-vector()
	         indici<-inthenet
                 while ((length(campione)==0)&(length(indici)>0))
	     	   {aus.regulated<-sampleB(indici,1)
            	    campione<-setdiff(which(Mdiscr[aus.regulated,]==1),linked)
	    	    if (length(campione)==0) indici<-setdiff(indici,aus.regulated)
	    	   }
		 if (length(campione)>0)
		   {aus.regulator<-sampleB(campione,1)
            	    Mdiscr[i,aus.regulator]<-1
	    	    Mdiscr[aus.regulated,aus.regulator]<-0
	    	    Mdiscr[aus.regulated,i]<-1
		    num<-num-1
	 	    linked<-c(linked,aus.regulator)
		    if (min.out.i>1) give.outlink<-c(give.outlink,rep(i,(min.out.i-1)))
		   }
		  else
		    {num<-0
		     Mdiscr[i,]<-0
		     not.regulated<-c(not.regulated,i)
		    }

	   	} #end else #(m.c==0)
	     }
	  } #end while (num>0)
	if (sum(Mdiscr[i,])>0) inthenet<-c(inthenet,i)

	L.m<-length(mem.o)
        if (L.m>0)
	  {mj<-1
           L.go<-length(give.outlink)
           for (h in (1:L.m))
	     {ind<-which(give.outlink==mem.o[h])[1]
      	      indok<-setdiff(seq(1,L.go,1),ind)
              give.outlink<-give.outlink[indok]
	      L.go<-L.go-1
	     }
          }
       } #end for (j in (1:NN))


    if (length(not.regulated)>0)
     {if (length(intersect(not.inthenet,not.regulated))==NN)
	{#cat("\n The algorithm converged to a topology in which in-degree and out-degree constraint are not compatible. \n tollerance constraints are relaxed.\n")
	 toll1<-which(toll.in!=0)
         toll.in[toll1]<-toll.in[toll1]+1
	 toll1<-which(toll.out!=0)
         toll.out[toll1]<-toll.out[toll1]+1
	}
     }
    not.inthenet<-not.regulated
    NN<-length(not.inthenet)

   } #END while (NN>0)


###check that every gene has at least 1 regulator
Sr<-apply(Mdiscr,1,sum)
ind<-which(Sr==0)
L<-length(ind)
if (L>0)
 {for (i in (1:L))
   {ri<-ind[i]
    num<-1
    Sout<-apply(Mdiscr,2,sum)
    Sc<-Score(S=Sout,ST=STout,Freq=Freq.out,n=1,toll.out)
    ind.Sc<-which(Sc==(-Inf))
    Sc[ind.Sc]<-0
    if (length(which(Sc<=0))>0) Sc<-Sc-min(Sc)+1/(N^2)
    Sc[ind.Sc]<-0
    if (sum(Sc)==0)
	{#cat("\n WARNING: tolerance parameters were relaxed to allow each node to have at least 1 regulator")
	 Sc<-Score(S=Sout,ST=STout,Freq=Freq.out,n=1,toll=rep(Inf,N))
	 ind.Sc<-which(Sc==(-Inf))
  	 Sc[ind.Sc]<-0
  	 if (length(which(Sc<=0))>0) Sc<-Sc-min(Sc)+1/(N^2)
  	 Sc[ind.Sc]<-0
	}
    p.out<-Sc/sum(Sc)
    ind.s<-sampleB(seq(1,N,1),num,prob=p.out)
    Mdiscr[ri,ind.s]<-1
    a1<-Sout[ind.s]
    a2<-a1+1
    Freq.out[a1+1]<-Freq.out[a1+1]-1
    Freq.out[a2+1]<-Freq.out[a2+1]+1
   }
 }


#######assign weights to the connectivity matrix##########
ind<-which(Mdiscr==1,arr.ind=TRUE)
L<-dim(ind)[1]
aus<-abs(rnorm_s(L,weight.mean, weight.sd, chi="csf"))
M<-matrix(0,ncol=N,nrow=N)
M[ind]<-aus

cat(file=fp, "csf output:\n\tlist(M,Mdiscr) = ", M, "\n", Mdiscr, "\n")

return(list(M,Mdiscr))
}

