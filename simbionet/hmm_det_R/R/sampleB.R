###     sampleB  (2009-12-05)
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

sampleB<-function(x,size,replace=FALSE,prob=NULL)
{
 cat(file=fp, "\n***sampleB***\ninput:\n\tx = ", x, "\n\tsize = ", size, "\n\treplace = ", replace, "\n")
  if (!is.null(prob))
   cat(file=fp, "\tprob = ", prob, "\n")
 if (length(x)==1) res<-rep(x,size)
 else {if (missing(prob)) res<-sample_s(x,size,replace,chi="sampleB")
        else res<-sample_s(x,size,replace,prob, chi="sampleB")
      }
cat(file=fp, "sampleB output:\n\tres = ", res, "\n")
return(res)
}
