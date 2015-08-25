###     MODULE3  (2007-12-05)
###
###	Copyright (C) 2007  Di Camillo Barbara
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

MODULE3<-function(num3,Nrim,Mdiscr,h,h.new,Sin,Sout,STin, STout,Freq.in, Freq.out,max.con, Cf.c,toll)
## INPUT
# num3: maximum number of nodes allowed for modules 3 (integer)
# Nrim: number of nodes still to be connected at current iteration  (integer)
# Mdiscr: the connectivity matrix of the network under construction, with element (i,j)=1 if an arc goes from j to i and 0 otherwise
# h: nodes still to be connected at current iteration  (vector of integers)
# h.new: nodes to be connected at the following iteration  (vector of integers)
# Sin:  vector of integers (sum of the elements on the rows of the connectivity matrix at current iteration )
# Sout: vector of integers (sum of the elements on the coloumns of the connectivity matrix at current iteration)
# STin:  vector of real numbers that, in each position k, contains the expected number of nodes with in-degree k (NA if not constrained)
# STout: vector of real numbers that, in each position k, contains the expected number of nodes with out-degree k (NA if not constrained)
# Freq.in: vector of real numbers that, in each position k, contains the actual number of nodes with in-degree k  (at current iteration)
# Freq.out: vector of real numbers that, in each position k, contains the actual number of nodes with out-degree k  (at current iteration)
# max.con: maximum number of regulators that a node can have
# Cf.c:   required average clustering coeff. for the module  (real, non negative, <1)
# toll:   vector of real numbers that, in each position k, express the absolute tollerance on STin[k] and STout[k]

## OUTPUT
# a list of 5 elements:
# 1)  score proportional to the probability of the module being analized to be selected as a network subgraph at current iteration
# 2)  connectivity matrix of the module, with element (i,j)=1 if an arc goes from j to i and 0 otherwise
# 3)  matrix of scores Sij (real). See the paper for details
# 4)  label(character) indicating if the parameter max.con  (maximum in-degree) has to be increased to allow convergence
# 5)  indexes of nodes in the module that serve the role of x nodes ( with role indicated by x in figure 1)

{
 cat(file=fp, "\n***module3***\ninput:\n\tnum3 = ", num3,"\n\tNrim = ", Nrim,"\n\tMdiscr = ", Mdiscr, "\n\th = ", h,"\n\th_new = ", h.new,"\n\tSin = ", Sin, "\n\tSout = ", Sout,"\n\tSTin = ", STin, "\n\tSTout = ", STout, "\n\tFreq_in = ", Freq.in, "\n\tFreq_out = ", Freq.out, "\n\tmax_con = ", max.con, "\n\t Cf_c = ", Cf.c, "\n\ttoll = ", toll, "\n")
 controllo<-0
    while (controllo==0)
        {max.g3<-min(num3,Nrim)
         max.g3.UP<-min(max.con,max.g3-1)

         Ng<-sampleB(seq(2,max.g3,1),1)
         Ng.UP<-sampleB(seq(1,min(max.g3.UP,(Ng-1)),1),1)
         Ng.DOWN<-Ng-Ng.UP

         mm<-MOD3(Ng.UP,Ng.DOWN,STin=STin,STout=STout,max.con,Cf.c)
         aus.p<-probmod(M=mm,h=h,Sin=Sin,Sout=Sout,STin=STin,STout=STout,Freq.in=Freq.in,Freq.out=Freq.out,toll=toll)
	 Pm3<-aus.p[[1]]
	 #Sc<-aus.p[[2]]

         if (!is.na(Pm3))   controllo<-1
	  else {num3<-num3-1; if (num3<2) {controllo<-1; aus.p[[1]]<-Pm3<-NA}}

       }

aus.p[[5]]<-seq(Ng.DOWN+1,Ng.UP+Ng.DOWN,1) #hubs
cat(file=fp, "module3 output:\n\taus.p = ")
print_list(aus.p)
return(aus.p)
}






