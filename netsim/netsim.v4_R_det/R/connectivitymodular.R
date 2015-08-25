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


connectivitymodular<-function(N=50,max.con=12,gamma=2.2,INdegree=c("free","out"), Cf.cl=0.4, num.subnet=c(5,5,10),r.tol=0.1,a.tol=1, weight.mean=1, weight.sd=0.1)

{
 cat(file=fp, "\n***connectivity_modular***\ninput:\n\tN = ", N, "\n\tmax_con = ", max.con, "\n\tgamma = ", gamma, "\n\tINdegree = ", INdegree, "\n\tCf_cl = ", Cf.cl, "\n\tnum_subnet = ", num.subnet, "\n\tr_tol = ", r.tol, "\n\ta_tol = ", a.tol, "\n\tweight_mean = ", weight.mean, "\n\tweight_sd = ", weight.sd, "\n")
# this function create a directed graph. It gives as output a list of two matrixes M and Mdiscr.
# M is the connectivity matrix with element (i,j)=1 if an arc goes from j to i and 0 otherwise
# (M is the transpose of the standard adjacence matrix)
# Mdiscr is structured as M but contains weights (i,j) belonging to real positive numbers

#\arguments{
#  \item{N}{(integer) Number of genes in the network.}
#  \item{max.con}{(integer) Maximum number of regulators that each node (gene) in the network can have. Default=5}
#  \item{gamma}{(real) Parameter of the power law distribution of the degree of connectivity of the nodes in the graph.Default=2.2  equal to the average value observed by (Albert and Barabasi (2000)) in protein networks.}
#  \item{INdegree}{(character) The in-degree distribution: either "free", i.e. it is not constrained to follow any distribution or "out", i.e. it is constrained to follow the same distribution of the out-degree (in this latter case the algorithm does not always converge)}
#  \item{Cf.cl}{(real) Average clustering coefficient of each sub-network in the graph (see References for further details).Default=0.4, leading to an average clustering coefficient in the entire network ranging between 0.1 and 0.3 as observed by (Barabasi and Albert (1999)) in protein networks.}
#  \item{num.subnet}{(integer) Vector of 3 elements containing the maximum number of nodes in each sub-network motif of type 1, 2 or 3 respectively (see References for further details). Default=c(5,5,10)}
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

        INdegree<-match.arg(INdegree)

        if (!is.numeric(Cf.cl))
        stop("`max.reg' must be numeric")

        if (!is.numeric(num.subnet))
        stop("`num.subnet' must be numeric")
        if (length(num.subnet)!=3)
        stop("`num.subnet' must have length 3")

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

 aus<-cbind(STout*r.tol,rep(a.tol,length(STout)))
 toll<-apply(aus,1,max)

 STin[(max.con+2):(N+2)]<-0

 if (INdegree=="out") STin[2:(max.con+1)]<-N*Prob[1:max.con]/sum(Prob[1:max.con])
  else STin[2:(max.con+1)]<-NA

 num1<-num.subnet[1]
 num2<-num.subnet[2]
 num3<-num.subnet[3]
 h<-seq(1,N,1)
 h.new<-vector()
 Cf.c<-Cf.cl
 it<-0
 Pm1<-Pm2<-Pm3<-1/3

 while (length(h)>1)
   {Nrim<-length(h)

    Sin<-apply(Mdiscr,1,sum)
    Sout<-apply(Mdiscr,2,sum)
    if (Cf.c>0) CC<-cluster.coeff(Mdiscr)[[1]]    else CC<-0
    if (CC>=Cf.c) CCi<-0 else CCi<-Cf.c

    controllo<-1
    while (controllo==1)
     {if (!is.na(Pm1))  {aus1<-MODULE1(num1,Nrim,Mdiscr,h,h.new,Sin,Sout,STin, STout,Freq.in, Freq.out, max.con, CCi,toll);Pm1<-aus1[[1]]}
      if (!is.na(Pm2))  {aus2<-MODULE2(num2,Nrim,Mdiscr,h,h.new,Sin,Sout,STin, STout,Freq.in, Freq.out, max.con, CCi,toll);Pm2<-aus2[[1]]}
      if (!is.na(Pm3))  {aus3<-MODULE3(num3,Nrim,Mdiscr,h,h.new,Sin,Sout,STin, STout,Freq.in, Freq.out, max.con, CCi,toll);Pm3<-aus3[[1]]}

      prob.mod<-c(Pm1,Pm2,Pm3)
      prob.mod[which(prob.mod==0)]<-NA
      ind<-which(!is.na(prob.mod))
      if (length(ind>0))
        {controllo<-0
	 mn<-min(prob.mod[ind])
         if (mn<=0) prob.mod<-prob.mod-mn+1/9
         prob.mod[which(is.na(prob.mod))]<-0
	 prob.mod<-prob.mod/sum(prob.mod)
         mod.type<-sampleB(c(1,2,3),1,prob=prob.mod)
         if (mod.type==1) {M<-aus1[[2]]; Sc<-aus1[[3]];  hubs<-aus1[[5]]}
	 if (mod.type==2) {M<-aus2[[2]]; Sc<-aus2[[3]];  hubs<-aus2[[5]]}
	 if (mod.type==3) {M<-aus3[[2]]; Sc<-aus3[[3]];  hubs<-aus3[[5]]}
	}

      else
	{Pm1<-Pm2<-Pm3<-1/3
	 if (CCi>0)  {CCi<-max(0,CCi-0.1); Cf.c<-max(0,Cf.c-0.1)}
          else
	    {if ( (aus1[[4]]=="in")&(aus2[[4]]=="in")&(aus3[[4]]=="in")&(max.con<N))
		{max.con<-max.con+1
		 if (INdegree=="out") STin[2:(max.con+1)]<-N*Prob[1:max.con]/sum(Prob[1:max.con])
  		  else STin[2:(max.con+1)]<-NA
		 Cf.c<-Cf.cl
		 if (CC>=Cf.c) CCi<-0 else CCi<-Cf.c
		 #cat("\n WARNING: maximum in.degree set to",max.con,"in iteration",it+1,"connecting the remaining",length(h)+length(h.new),"hubs, because lower in.degree is not compatible with module structure and other parameters setting\n")
		}
	     else
		{if (gamma==0) stop("computation failed! Try with different parameter settings!")
		 gamma<-max(gamma-0.2,0)

                 if (gamma!=0)
		   {Prob<-c(seq(1,N,1)^(-gamma),0)
 		    Prob<-Prob/(sum(Prob))
 		    #p<-Prob[seq(1,max.con,1)]/sum(Prob[seq(1,max.con,1)])
 		    STout<-c(NA,Prob*N)
		    if (INdegree=="out") STin[2:(max.con+1)]<-N*Prob[1:max.con]/sum(Prob[1:max.con])
  			else STin[2:(max.con+1)]<-NA
		    Cf.c<-Cf.cl
		    if (CC>=Cf.c) CCi<-0 else CCi<-Cf.c
		    #cat("\n WARNING: gamma set to",gamma,"in iteration",it+1,"connecting the remaining",length(h)+length(h.new),"hubs, because higher gamma is not compatible with module structure and other parameters setting\n")
		   }
		  else
		   {Prob<-c(rep(1/N,N),0)
 		    #p<-Prob[seq(1,max.con,1)]/sum(Prob[seq(1,max.con,1)])
 		    STout<-c(NA,Prob*N)
		    if (INdegree=="out") STin[2:(max.con+1)]<-N*Prob[1:max.con]/sum(Prob[1:max.con])
  			else STin[2:(max.con+1)]<-NA
		    Cf.c<-Cf.cl
		    if (CC>=Cf.c) CCi<-0 else CCi<-Cf.c
		    #cat("\n WARNING: power law distribution replaced by flat distribution in iteration",it+1,"connecting the remaining",length(h)+length(h.new),"hubs, because power-law distribution is not compatible with module structure and other parameters setting\n")
		   }

		 }

	     }
        }

     }#end while (controllo==1)

    Ng<-dim(M)[1]
    aus<-assign.nodes(M,Mdiscr,h,hubs=hubs,Sc=Sc,Sin,max.con)
    Mdiscr<-aus[[1]]
    h<-aus[[2]]
    h.new<-c(h.new,aus[[3]])
    Sout<-apply(Mdiscr,2,sum)
    m<-max(Sout)+1
    Freq.out[1:m]<-hist(Sout,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts
    Sin<-apply(Mdiscr,1,sum)
    m<-max(Sin)+1
    Freq.in[1:m]<-hist(Sin,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts

    Nrim<-Nrim-Ng
    if (Nrim<=1)
     {it<-it+1
      if ((it==1)&(Nrim==1)) h<-c(h,h.new)
       else h<-h.new
      h.new<-vector()
      if (length(h)>0) h<-h[which(Sin[h]!=max.con)]
      Pm1<-Pm2<-Pm3<-1/3
     }

   }# END WHILE  (length(h)>1)


###check that every gene has at least 1 regulator
Sr<-apply(Mdiscr,1,sum)
ind<-which(Sr==0)
L<-length(ind)
if (L>0)
 {for (i in (1:L))
   {ri<-ind[i]
    num<-1
    Sout<-apply(Mdiscr,2,sum)
    Sc<-Score(S=Sout,ST=STout,Freq=Freq.out,n=1,toll)
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
aus<-abs(rnorm_s(L,weight.mean, weight.sd, chi="cm"))
M<-matrix(0,ncol=N,nrow=N)
M[ind]<-aus
cat(file=fp, "cm output:\n\tlist(M,Mdiscr) = ", M, "\n", Mdiscr, "\n")
return(list(M,Mdiscr))
}



