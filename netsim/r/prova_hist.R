dyn.load("hist.dll")

hist_dll <-
    function(x, breaks, right, include_border, naok)
{
    stopifnot(is.vector(x))
    stopifnot(is.vector(breaks))
    vett <- .Call("hist",
                 as.vector(x),
                 as.vector(breaks),
                 as.integer(right),
                 as.integer(include_border),
                 as.integer(naok)
					)
    return(vett)
}


x <- vector()
x1 <- runif(100, 1, 10)
x1
for (xv in x1)
    x <- c(x, round(xv))
x

breaks = seq(1,10)
breaks

l = hist_dll(x, breaks,right=0,include_border=1,naok=0)
l

l1 = hist(x, breaks,plot=FALSE,right=FALSE)$count
l1