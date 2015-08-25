library("netsim")

N <- 10

net<-simulatenet(N, connectivity="MTM", save=TRUE, itera=5)

Dati <- matrix(c(1.793561,1.793561,1.793561,1.793561,1.793561,4.239327,4.239327,4.239327,4.239327,4.239327,1.031338,1.031338,1.031338,1.031338,1.031338,1.789676,1.789676,1.789676,1.789676,1.789676,7.767749,7.767749,7.767749,7.767749,7.767749), 5, 5)

pesi <- NULL

lambda <-c(1.1056752,0.9747519,1.0007163,1.0061211,1.1354204)

alpha <- c(9.869912,9.955962,10.056057,9.922354,9.950865)

beta <- c(0.4934634,0.4930697,0.4926102,0.4976155,0.4905647)

N<-dim(Dati)[1]
params <- c(4,4,4,4,4)
times<-seq(1,5)
save <- FALSE

ris <- simulateprofiles(weights=pesi, Rules=NULL,lambda=lambda, alpha=alpha, beta=beta, method="rkf45", params=params, itera=3, times=times)
ris
cat("ok")
