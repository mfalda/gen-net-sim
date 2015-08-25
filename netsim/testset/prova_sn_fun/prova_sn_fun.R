library(netsim)

p<-read.table("../parameters.txt")

f <- function(x)
{
	return (1/(1+exp(-10*(x-0.5))))
}

init()

simulatenet(N=5, f.pr.and=f, lambda=p[,1], alpha=p[,2], beta=p[,3], x0=p[,4], Xmin=p[,5], Xmax=p[,6], times=seq(0, 5, 0.1))

cat("ok")