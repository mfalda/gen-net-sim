###     sampleB  (2007-12-05)
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

sampleB<-function(x,size,replace=FALSE,prob=NULL)

#I modified the original function "sample_s" since it gave errors in case x had length 1
## INPUT

# x: Either a (numeric, complex, character or logical) vector of more than one element from which to choose, or a positive integer.
# size: non-negative integer giving the number of items to choose.
# replace: TRUE/FALSE Should sampling be with replacement?
#prob: a vector of probability weights for obtaining the elements of the vector being sampled.

## OUTPUT
# a sample_s of the specified size from the elements of x using either with or without replacement


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
