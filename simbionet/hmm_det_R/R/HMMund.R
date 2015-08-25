###     HMMund  (2009-May-12th)
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


HMMund<-function(N=50,Cf.cl=0.3,gamma=2.2,DEGREE=NULL,MODULES,prior.p.subnet=NULL,max.con=12,sepgraph=TRUE,r.tol=0.1,a.tol=1,iter=1)

{

#\arguments{
#  \item{N}{Number of nodes in the network.}
#  \item{Cf.cl}{Average clustering coefficient of each sub-network in the graph (see References for further details).Default=0.4, leading to an average clustering coefficient in the entire network ranging between 0.1 and 0.3 as observed by (Barabasi and Albert (1999)) in protein networks.}
#  \item{gamma}{Parameter of the power law distribution of the degree of connectivity of the nodes in the graph.Default=2.2  equal to the average value observed by (Albert and Barabasi (2000)) in protein networks.}
#  \item{DEGREE}{Probability distribution of the degree. Default=NULL: power-law distribution will be used. Otherwise a vector of length (N+1) should be provided with, in position k, the average number of nodes with desired degree (k-1)}
#  \item{MODULES}{list of L modules obtained by running function "generatemodules(autoreg=FALSE)" if self-regulation is prohibited; "generatemodules(autoreg=TRUE)" if self-regulation is allowed"}
#  \item{prior.p.subnet}{Vector of L elements containing the a priori probability of having a sub-network motif of type L. Default=equal probability}
#  \item{max.con}{Maximum degree that each node in the network can have. Default=12}
#  \item{sepgraph}{If FALSE networks consisting of separate subgraph are prohibited. Default is TRUE.}
#  \item{r.tol}{relative tolerance parameter for the power law distribution. Default=0.1 (See Details).}
#  \item{a.tol}{absolute tolerance parameter for the power law distribution. Default=0.1 (See Details).}
#  }

   cat(file=fp, "\n***HMMund***\n\n")

   ##### CHECK PARAMETERS ####
        if (!is.numeric(N))
        stop("`N' must be numeric")

        if (!is.numeric(Cf.cl))
        stop("`max.reg' must be numeric")

        if (!is.numeric(gamma))
        stop("`gamma' must be numeric")

        if (!is.null(DEGREE))
         {if (!is.numeric(DEGREE))  stop("DEGREE must be numeric")
          if (length(which(DEGREE<0))>0)  stop("DEGREE must contain positive numbers")
          #~ if (sum(DEGREE,na.rm=TRUE)!=N)  stop("sum(DEGREE) must be equal to N")
          if (length(DEGREE)!=(N+1))  stop("DEGREE must be of length (N+1)")
          if (max(DEGREE,na.rm=TRUE)>N)  stop("DEGREE must have AT MAXIMUM value N")
         }

        if (!is.list(MODULES))
        stop("'MODULES' must be a list")
        L<-length(MODULES)
        CC<-FB<-AU<-MC<-DIM<-SIMM<-rep(0,L)
        for (i in (1:L))
          {CC[i]<-MODULES[[i]]$CC                     #clustering coefficient
           FB[i]<-MODULES[[i]]$feedback               #presence of feedback loop
           AU[i]<-MODULES[[i]]$autoreg                #presence of autoregulation
           MC[i]<-max(apply(MODULES[[i]]$net,1,sum))  #max in.degree
           DIM[i]<-MODULES[[i]]$dim.m
           SIMM[i]<-MODULES[[i]]$SIMM
          }
        if (length(which(AU==TRUE))>0) autoreg<-TRUE
        else autoreg<-FALSE


        if (is.null(prior.p.subnet)) prior.p.subnet<-rep(1/L,L)
        if (!is.numeric(prior.p.subnet))
        stop("`prior.p.subnet' must be numeric")
        if (length(prior.p.subnet)!=L)
        stop("`prior.p.subnet' must have length L=",L)
        if (length(which(is.na(prior.p.subnet)))>0)
        stop("`prior.p.subnet' can't contain NA values")
        if ((length(which(prior.p.subnet>1))>0)|(length(which(prior.p.subnet<0))>0))
        stop("`prior.p.subnet' must contain values greater than 0 and lower than 1")
        if (round((sum(prior.p.subnet))-1,10)!=0)
        stop("`prior.p.subnet' must have sum 1")

        prior.p.subnet[which(SIMM==FALSE)]<-0
        prior.p.subnet<-prior.p.subnet/sum(prior.p.subnet)


        if (!is.numeric(max.con))
        stop("`max.con' must be numeric")
        if (max.con<dim(MODULES[[L]]$net)[1])
          {prior.p.subnet[which(MC>max.con)]<-0
           prior.p.subnet<-prior.p.subnet/sum(prior.p.subnet)
          }


        if ((sepgraph!=TRUE)&(sepgraph!=FALSE))
        stop("'sepgraph' must be TRUE or FALSE")

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

 LG<-rep(-1,N)
 Mdiscr<-matrix(0,ncol=N,nrow=N)
 if (max.con>N) max.con<-N
 Freq.out<-Freq.in<-c(N,rep(0,N+1))

 if (is.null(DEGREE))
  {Prob<-c(NA,seq(1,N,1)^(-gamma),0)
   Prob<-Prob/(sum(Prob,na.rm=TRUE))
  }
 else Prob<-c(DEGREE/(sum(DEGREE,na.rm=TRUE)),0)

 STout<-rep(0,(N+2))
 STout[(max.con+2):(N+2)]<-0
 STout[2:(max.con+1)]<-N*Prob[2:(max.con+1)]/sum(Prob[2:(max.con+1)],na.rm=TRUE)
 aus<-cbind(STout*r.tol,rep(a.tol,length(STout)))
 toll<-apply(aus,1,max)


 h<-seq(1,N,1)
 h.new<-vector()
 Cf.cl.int<-Cf.cl
 it<-0

 CCp<-rep(1,L)
 if (Cf.cl!=0)
   {CCs<-sort(union(CC,CC))
    CCind<-list()
    Lc<-length(CCs)
    p<-rep(1,Lc)
    for (i in (1:Lc)) CCind[[i]]<-which(CC==CCs[i])
   }


 while (length(h)>1)
   {Nrim<-length(h)

    Sout<-apply(Mdiscr,2,sum)
    #################
    Mcc<-cluster.coeff(Mdiscr)[[1]] #mean(CC[ind.possible])
    if (Cf.cl>0)
       {if (Cf.cl.int!=Mcc)
          {segno<-sign(CCs-Mcc)
           indp<-which(segno>0)
           indm<-which(segno<=0)
           if (Cf.cl.int>Mcc)
             {p[indp]<-(1-abs(Cf.cl.int-CCs[indp]))+p[indp]
              p[indm]<-(-abs(Cf.cl.int-CCs[indm]))+p[indm]
             }
           else
             {p[indm]<-(1-abs(Cf.cl.int-CCs[indm]))+p[indm]
              p[indp]<-(-abs(Cf.cl.int-CCs[indp]))+p[indp]
             }
           if (min(p)<=0) p[which(p<0)]<-0.0001 #if (min(p)<=0) p<-p-min(p)+0.001
           p<-p/sum(p)
           for (ccc in (1:Lc)) CCp[(CCind[[ccc]])]<-p[ccc]
           }
        else
          {p<-p+(1-abs(Cf.cl.int-CCs))
           p<-p/sum(p)
           for (ccc in (1:Lc)) CCp[(CCind[[ccc]])]<-p[ccc]
          }
       }
    ppp<-prior.p.subnet*CCp
    if ((max(ppp)==0)|(is.na(max(ppp))))   stop("\n prior.p.subnet values incompatible with other parameter settings! \n")
    ppp<-ppp/sum(ppp,na.rm=TRUE)

    controllo<-1
    controllo2<-1

    if (length(h)<max(DIM)) ppp[which(DIM>length(h))]<-0
    if ((max(ppp)==0)|(is.na(max(ppp)))) {controllo<-0; controllo2<-0}

    while (controllo==1)
     {k<-min(3,length(which(ppp>0)))
      m<-sampleB(1:L,k,prob=ppp)
      mod1<-MODULES[[(m[1])]]$net
      aus1<-probmodund(M=mod1,h=h,Sout=Sout,STout=STout,Freq.out=Freq.out,toll=toll)
      Pm1<-aus1[[1]]

      if (k>1)
        {mod2<-MODULES[[(m[2])]]$net
         aus2<-probmodund(M=mod2,h=h,Sout=Sout,STout=STout,Freq.out=Freq.out,toll=toll)
         Pm2<-aus2[[1]]
        }
      else Pm2<-NA

      if (k>2)
        {mod3<-MODULES[[(m[3])]]$net
         aus3<-probmodund(M=mod3,h=h,Sout=Sout,STout=STout,Freq.out=Freq.out,toll=toll)
         Pm3<-aus3[[1]]
        }
      else Pm3<-NA

      pm<-prob.mod<-c(Pm1,Pm2,Pm3)
      prob.mod[which(is.na(prob.mod))]<-0
      mn<-min(prob.mod)
      if (mn<=0) prob.mod<-prob.mod-mn+1/9
      prob.mod[which(is.na(pm))]<-0
	    ind<-which(prob.mod!=0)
      if (length(ind>0))
        {prob.mod<-prob.mod/sum(prob.mod)
         mod.type<-sampleB(c(1,2,3),1,prob=prob.mod)
         if (mod.type==1) {M<-aus1[[2]]; Sc<-aus1[[3]];  hubs<-MODULES[[(m[1])]]$hubsio}
	       if (mod.type==2) {M<-aus2[[2]]; Sc<-aus2[[3]];  hubs<-MODULES[[(m[2])]]$hubsio}
	       if (mod.type==3) {M<-aus3[[2]]; Sc<-aus3[[3]];  hubs<-MODULES[[(m[3])]]$hubsio}
         controllo<-0
         if (is.na(Pm1))   ppp[[(m[1])]]<-0
         if ((is.na(Pm2))&(k>1))   ppp[[(m[2])]]<-0
         if ((is.na(Pm3))&(k>2))   ppp[[(m[3])]]<-0
         ppp<-ppp/sum(ppp)
        }
      else
        {ppp[m]<-0
         ppp<-ppp/sum(ppp)
         if ((max(ppp)==0)|(is.na(max(ppp))))
          {if (Cf.cl.int>0)
             {Cf.cl.int<-max(0,Cf.cl.int-0.1)
              if (Cf.cl>0)
                 {if (Cf.cl.int!=Mcc)
                    {segno<-sign(CCs-Mcc)
                    indp<-which(segno>0)
                    indm<-which(segno<=0)
                    if (Cf.cl.int>Mcc)
                      {p[indp]<-(1-abs(Cf.cl.int-CCs[indp]))+p[indp]
                       p[indm]<-(-abs(Cf.cl.int-CCs[indm]))+p[indm]
                      }
                    else
                      {p[indm]<-(1-abs(Cf.cl.int-CCs[indm]))+p[indm]
                       p[indp]<-(-abs(Cf.cl.int-CCs[indp]))+p[indp]
                      }
                    if (min(p)<=0) p[which(p<0)]<-0.0001
                    p<-p/sum(p)
                    for (ccc in (1:Lc)) CCp[(CCind[[ccc]])]<-p[ccc]
                    }
                  else
                  {p<-p+(1-abs(Cf.cl.int-CCs))
                  p<-p/sum(p)
                  for (ccc in (1:Lc)) CCp[(CCind[[ccc]])]<-p[ccc]
                  }
                }
              ppp<-prior.p.subnet*CCp
              if ((max(ppp)==0)|(is.na(max(ppp))))   stop("\n prior.p.subnet values incompatible with other parameter settings! \n")
              ppp<-ppp/sum(ppp,na.rm=TRUE)
             }

           else
            {if ((gamma==0)|(!is.null(DEGREE))) stop("computation failed! Try again or consider using different tolerance parameter settings or different topology parameter settings")
		         gamma<-max(gamma-0.2,0)
             if (gamma!=0)
		          {Prob<-c(seq(1,N,1)^(-gamma),0)
 		           Prob<-Prob/(sum(Prob))
 		           # DEGREE==NULL
 		           STout<-rep(0,(N+2))
               STout[(max.con+2):(N+2)]<-0
               STout[2:(max.con+1)]<-N*Prob[2:(max.con+1)]/sum(Prob[2:(max.con+1)],na.rm=TRUE)
               aus<-cbind(STout*r.tol,rep(a.tol,length(STout)))
               toll<-apply(aus,1,max)

               Cf.cl.int<-Cf.cl
		           if (Cf.cl>0)
                 {if (Cf.cl.int!=Mcc)
                    {segno<-sign(CCs-Mcc)
                    indp<-which(segno>0)
                    indm<-which(segno<=0)
                    if (Cf.cl.int>Mcc)
                      {p[indp]<-(1-abs(Cf.cl.int-CCs[indp]))+p[indp]
                       p[indm]<-(-abs(Cf.cl.int-CCs[indm]))+p[indm]
                      }
                    else
                      {p[indm]<-(1-abs(Cf.cl.int-CCs[indm]))+p[indm]
                       p[indp]<-(-abs(Cf.cl.int-CCs[indp]))+p[indp]
                      }
                    if (min(p)<=0) p[which(p<0)]<-0.0001
                    p<-p/sum(p)
                    for (ccc in (1:Lc)) CCp[(CCind[[ccc]])]<-p[ccc]
                    }
                  else
                  {p<-p+(1-abs(Cf.cl.int-CCs))
                  p<-p/sum(p)
                  for (ccc in (1:Lc)) CCp[(CCind[[ccc]])]<-p[ccc]
                  }
                }
                ppp<-prior.p.subnet*CCp
                if ((max(ppp)==0)|(is.na(max(ppp))))   stop("\n prior.p.subnet values incompatible with other parameter settings! \n")
                ppp<-ppp/sum(ppp,na.rm=TRUE)
                cat("\n WARNING: gamma set to",gamma,"in iteration",it+1,"connecting the remaining",length(h)+length(h.new),"hubs, because higher gamma is not compatible with current network structure and other parameters setting\n")
		           }
		         else
		          {Prob<-c(rep(1/N,N),0)
 		           # DEGREE==NULL
		           STout<-rep(0,(N+2))
               STout[(max.con+2):(N+2)]<-0
               STout[2:(max.con+1)]<-N*Prob[2:(max.con+1)]/sum(Prob[2:(max.con+1)],na.rm=TRUE)
               aus<-cbind(STout*r.tol,rep(a.tol,length(STout)))
               toll<-apply(aus,1,max)

               Cf.cl.int<-Cf.cl
		           if (Cf.cl>0)
                 {if (Cf.cl.int!=Mcc)
                    {segno<-sign(CCs-Mcc)
                    indp<-which(segno>0)
                    indm<-which(segno<=0)
                    if (Cf.cl.int>Mcc)
                      {p[indp]<-(1-abs(Cf.cl.int-CCs[indp]))+p[indp]
                       p[indm]<-(-abs(Cf.cl.int-CCs[indm]))+p[indm]
                      }
                    else
                      {p[indm]<-(1-abs(Cf.cl.int-CCs[indm]))+p[indm]
                       p[indp]<-(-abs(Cf.cl.int-CCs[indp]))+p[indp]
                      }
                    if (min(p)<=0) p[which(p<0)]<-0.0001
                    p<-p/sum(p)
                    for (ccc in (1:Lc)) CCp[(CCind[[ccc]])]<-p[ccc]
                    }
                  else
                  {p<-p+(1-abs(Cf.cl.int-CCs))
                  p<-p/sum(p)
                  for (ccc in (1:Lc)) CCp[(CCind[[ccc]])]<-p[ccc]
                  }
                }
              ppp<-prior.p.subnet*CCp
              if ((max(ppp)==0)|(is.na(max(ppp))))   stop("\n prior.p.subnet values incompatible with other parameter settings! \n")
              ppp<-ppp/sum(ppp,na.rm=TRUE)
              cat("\n WARNING: power law distribution replaced by flat distribution in iteration",it+1,"connecting the remaining",length(h)+length(h.new),"hubs, because power-law distribution is not compatible with current network structure and other parameters setting\n")
		          }
             }
	         }#end if (max(ppp,na.rm=TRUE)==0)
        }


     } #end controllo==1

  if (controllo2==1)
   {Ng<-dim(M)[1]
    if (length(hubs)==Ng) hubs<-sampleB(hubs,Ng%/%2)
    aus<-assign.nodes.und(M,Mdiscr,h,hubs=hubs,Sc=Sc,Sout,max.con)
    Mdiscr<-aus[[1]]
    ########livello gerarchico#############
    LG[setdiff(h,aus[[2]])]<-LG[setdiff(h,aus[[2]])]+1
    #######################################
    h<-aus[[2]]
    h.new<-c(h.new,aus[[3]])
    Sout<-apply(Mdiscr,2,sum)
    m<-max(Sout)+1
    Freq.out[1:m]<-hist(Sout,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts

    Nrim<-Nrim-Ng
    if (Nrim<=1)
     {it<-it+1
      if ((it==1)&(Nrim==1)) h<-c(h,h.new)
       else h<-h.new
      h.new<-vector()
      if (length(h)>0) h<-h[which(Sout[h]!=max.con)]
     }
   } #end  if (controllo2==1)
   else  #controllo2==0
    {if (it>0)
      {if (length(h.new)==0) h<-vector()
       else
         {it<-it+1
          h<-h.new
          h.new<-vector()
          if (length(h)>0) h<-h[which(Sout[h]!=max.con)]
         }
       }
      else #it==0
        {if (length(h.new)==0) stop("\n prior.p.subnet values incompatible with other parameter settings! \n")
         else
         {it<-it+1
          h<-c(h,h.new)
          h.new<-vector()
          if (length(h)>0) h<-h[which(Sout[h]!=max.con)]
         }
       }
     }

   }# END WHILE  (length(h)>1)

###check graph connectivity and add links if it is not completely connected


if (sepgraph==FALSE)
{dist<-check.conn(Mdiscr)
 non.connessi<-which(dist==Inf)
 L<-length(non.connessi)
 while (L>0)
  {Mdiscr<-connetti.scalefree(Mdiscr,STout,STout,dist,toll,max.con,und=TRUE)
   dist<-check.conn(Mdiscr)
   non.connessi<-which(dist==Inf)
   L<-length(non.connessi)
  }
}

if (autoreg==FALSE) diag(Mdiscr)<-0
etich<-paste("topology_",iter,".txt",sep="")
write.table(Mdiscr,etich,sep="\t",row.names=FALSE,col.names=FALSE)
etich<-paste("Hierarch_level",iter,".txt",sep="")
write.table(cbind(1:N,LG),etich,sep="\t",row.names=FALSE,col.names=FALSE)

return(Mdiscr)
}


